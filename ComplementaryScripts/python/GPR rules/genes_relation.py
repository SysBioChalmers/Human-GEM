
import sys
import os
import re
import hashlib
import itertools

'''
Set a functions to parse complexes files and return GPR and/or rules
Note: the job done here have been or will be translated into matlab
'''

def filter_genes(genes_file, keywords):
    '''
        return two dictionnaries:
            - ENSEMBL ID <-> array of UNIPROT ID
            - UNIPROT ID <-> array of ENSEMBL ID

        genes_file: csv export of the sheet GTPRs of the HMR excel file.
        keywords: array of string, list of word use to filter the gene
        select genes if keywords in columns ncbi_desc or uniprot_subunit_struct

        without keywords, 'results' and 'hash_ens_dict' are empty

        GTPR columns
        # Ensembl: Gene ID                              0
        Ensembl: Transcript ID
        Ensembl: Protein ID
        Ensembl: UniProt ID
        Ensembl: HGNC ID
        Ensembl: Entrez ID                              5
        NCBI: EntrezID
        NCBI: chromosome
        NCBI: description
        NCBI: type_of_gene
        NCBI: Other_designations                        10
        NCBI: Gene Summary
        UniProt: UniProt ID
        UniProt: Gene names
        UniProt: Subcellular location
        UniProt: Protein names                          15
        UniProt: Catalytic activity
        UniProt: EC number
        UniProt: Function
        UniProt: Subunit structure
        GENE NAME                                       20
        GENE ID 1
        GENE ID 2
        SHORT NAME
        HMR Reaction Ids based on the old gene Ids
    '''

    GTPR_ENS_ID_INDEX = 0
    GTPR_ENTREZ_ID_INDEX = 6  # use the ncbi col, there are duplicate ids in col index 5

    NCBI_DESCRIPTION_INDEX = 8
    NCBI_SUMMARY_INDEX = 14

    UNIPROT_ID_INDEX = 12
    UNIPROT_GENE_NAME_INDEX = 13
    UNIPROT_PROT_NAME_INDEX = 15
    UNIPROT_SUBSTRUCT_INDEX = 19

    count = 0
    uniprot_ens_dict = {}
    ens_uniprot_dict = {}
    results = {}
    found = False
    hash_ens_dict = {}

    with open(genes_file, 'r') as fh:
        fh.readline()
        fh.readline()
        reaction_dict = {}
        for i, line in enumerate(fh):
            if line[0] == '#':
                continue
            line = line.strip('\n')
            linearr = line.split('\t')

            ens_gene_id = linearr[GTPR_ENS_ID_INDEX].strip()
            uniprot_id = linearr[UNIPROT_ID_INDEX].strip()

            ncbi_desc = linearr[NCBI_DESCRIPTION_INDEX].strip()
            uniprot_prot_name = linearr[UNIPROT_GENE_NAME_INDEX].strip()
            uniprot_subunit_struct= linearr[UNIPROT_SUBSTRUCT_INDEX].strip()

            if ens_gene_id == '-':
                continue

            if ens_gene_id not in ens_uniprot_dict:
                ens_uniprot_dict[ens_gene_id] = []
            if uniprot_id:
                ens_uniprot_dict[ens_gene_id].append(uniprot_id)

            # store entrez id <=> ens id relation
            if uniprot_id and uniprot_id != '-':
                uniprot_ids = [el.strip() for el in uniprot_id.split(';') if el.strip()]
                for el in uniprot_ids:
                    if el not in uniprot_ens_dict:
                        uniprot_ens_dict[el] = []
                    uniprot_ens_dict[el].append(ens_gene_id)

            if keywords:
                # keywords not used anymore
                found = False
                fields = [ncbi_desc, uniprot_subunit_struct]
                for f in fields:
                    for word in keywords:
                        if word.lower() in f:
                            found = True
                            break

                    if found:
                        break

                if found:
                    count += 1
                    # add more information to help manual curration
                    le_hash =  hashlib.sha1(uniprot_subunit_struct).hexdigest()
                    results[ncbi_desc] = [ens_gene_id, uniprot_prot_name, ncbi_desc, le_hash, uniprot_subunit_struct]
                    if le_hash not in hash_ens_dict:
                        hash_ens_dict[le_hash] = []
                    hash_ens_dict[le_hash].append(ens_gene_id)

    print ("Gene match:", count)
    return results, uniprot_ens_dict, ens_uniprot_dict, hash_ens_dict


def parse_CORUM_txt(CORUM_file, uniprot_ens_dict):  # Corum_allComplexes.txt
    ''' 
        CORUM_file: 'Corum_allComplexes.txt' FROM http://mips.helmholtz-muenchen.de/corum/#download

        column indexes:
        ComplexID                                  0
        ComplexName
        Organism
        Synonyms
        Cell line
        subunits(UniProt IDs)                      5
        subunits(Entrez IDs)
        Protein complex purification method
        GO ID
        GO description
        FunCat ID                                  10
        FunCat description
        Complex comment
        PubMed ID
        Subunits comment
        subunits(Protein name)                     15
        subunits(Gene name)
        subunits(Gene name syn)
        SWISSPROT organism
        Disease comment
    '''

    complex_id_ens_dict = {}
    ens_ebi_complex_dict = {}

    entrez_without_ens = 0
    with open(CORUM_file, 'r') as fh:
        fh.readline()
        for i, line in enumerate(fh):
            line = line.strip('\n')
            # print line
            linearr = line.split('\t')
            complex_id = linearr[0]
            complex_name = linearr[1]
            organism = linearr[2]
            subunit_entrez_ids = linearr[6]
            subunit_uniprot_ids = linearr[5]

            if organism != "Human" or subunit_uniprot_ids in ["None", "NULL"]:
                continue

            for uniprot_id in [el.strip() for el in subunit_uniprot_ids.split(';')]:
                if uniprot_id in uniprot_ens_dict:
                    ens_ids = uniprot_ens_dict[uniprot_id] # list
                else:
                    entrez_without_ens += 1
                    continue
                for ens_id in ens_ids:
                    if ens_id not in ens_ebi_complex_dict:
                        ens_ebi_complex_dict[ens_id] = []
                    ens_ebi_complex_dict[ens_id].append([complex_id, complex_name])
                    if complex_id not in complex_id_ens_dict:
                        complex_id_ens_dict[complex_id] = []
                    complex_id_ens_dict[complex_id].append(ens_id)

    print ("gene assign to a CORUM annotated complex:", len(ens_ebi_complex_dict))
    print ("uniprot id w/o ens id:", entrez_without_ens)

    return ens_ebi_complex_dict, complex_id_ens_dict


def parse_EBI_txt(EBI_file, uniprot_ens_dict):  # EBI_complex_homo_sapiens_20171218.tsv
    '''
        EBI_file: EBI_complex_homo_sapiens_20171218.tsv from ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/
        column indexes:
        ComplexID                                  0
        subunits(UniProt IDs)                      4
    '''

    complex_id_ens_dict = {}
    ens_ebi_complex_dict = {}
    entrez_without_ens = 0

    with open(EBI_file, 'r') as fh:
        fh.readline()
        for i, line in enumerate(fh):
            line = line.strip('\n')
            linearr = line.split('\t')
            complex_id = linearr[0]
            subunit_uniprot_ids = linearr[4]

            if not subunit_uniprot_ids:
                continue

            for uniprot_id in [el.strip() for el in subunit_uniprot_ids.split('|')]:
                uniprot_id = uniprot_id.split('(')[0]
                if uniprot_id in uniprot_ens_dict:
                    ens_ids = uniprot_ens_dict[uniprot_id] # list
                else:
                    entrez_without_ens += 1
                    continue

                for ens_id in ens_ids:
                    if ens_id not in ens_ebi_complex_dict:
                        ens_ebi_complex_dict[ens_id] = []
                    ens_ebi_complex_dict[ens_id].append(complex_id)
                    if complex_id not in complex_id_ens_dict:
                        complex_id_ens_dict[complex_id] = []
                    complex_id_ens_dict[complex_id].append(ens_id)

    print ("gene assign to a EBI annotated complex:", len(ens_ebi_complex_dict))
    print ("uniprot id w/o ens id:", entrez_without_ens)

    return ens_ebi_complex_dict, complex_id_ens_dict


def parse_EBI_CORUM_file(both_complex_file, uniprot_ens_dict):  # both_cplx_file.txt from Hao
    '''
        both_complex_file: 'both_cplx_file.txt' from Hao
        The file contains EBI ID, CORUM ID and a unique ID for complexes

        column indexes:
        ComplexID                                  0
        subunits(UniProt IDs)                      3
    '''

    complex_id_ens_dict = {}
    ens_complex_dict = {}
    entrez_without_ens = 0

    with open(both_complex_file, 'r') as fh:
        fh.readline()
        for i, line in enumerate(fh):
            line = line.strip('\n')
            linearr = line.split('\t')
            complex_id = linearr[0]
            subunit_uniprot_ids = linearr[3]

            if not subunit_uniprot_ids:
                continue

            for uniprot_id in [el.strip() for el in subunit_uniprot_ids.split(',')]:
                if uniprot_id in uniprot_ens_dict:
                    ens_ids = uniprot_ens_dict[uniprot_id] # list
                else:
                    entrez_without_ens += 1
                    continue

                for ens_id in ens_ids:
                    if ens_id not in ens_complex_dict:
                        ens_complex_dict[ens_id] = []
                    ens_complex_dict[ens_id].append(complex_id)
                    if complex_id not in complex_id_ens_dict:
                        complex_id_ens_dict[complex_id] = []
                    complex_id_ens_dict[complex_id].append(ens_id)

    print ("gene assign to a MergedComplex ID:", len(ens_complex_dict))
    print ("uniprot id w/o ens id:", entrez_without_ens)

    return ens_complex_dict, complex_id_ens_dict


def build_genes_complex_file(genes_file, CORUM_file):

    '''
        not needed anyore
    '''

    ## filter the genes that may be in a complex
    keywords = ['subunit','complex', 'component', 'hetero']
    gene_filtered, uniprot_ens_dict, ens_uniprot_dict, hash_ens_dict = filter_genes(genes_file, keywords)
    ens_ebi_complex_dict, complex_id_ens_dict = parse_CORUM_txt(CORUM_file, uniprot_ens_dict)


    LINES = []
    # reformat and search for name reference in NSS, SKIPED
    name_id_list = [[v[1], v[0]] for k, v in gene_filtered.items() if v[1]]

    id_name_found_NSS = {}
    for el in gene_filtered:
        ens_id = gene_filtered[el][0]
        if ens_id in ens_ebi_complex_dict:
            list_complex = ens_ebi_complex_dict[ens_id]
            complex_id_str = [complex_id for complex_id, complex_name in list_complex]
            complex_name_str = [complex_name for complex_id, complex_name in list_complex]
            gene_filtered[el] += [complex_id_str, complex_name_str]
        else:
            gene_filtered[el] += ['','']

        for name, ens in name_id_list:
            if ens == ens_id:
                continue

            if re.search('(?:^|[.;(),\s/-])%s(?:$|[.;(),\s/-])' % name, gene_filtered[el][4], re.I):
                if ens not in id_name_found_NSS:
                    id_name_found_NSS[ens] = []
                id_name_found_NSS[ens].append(ens_id)

        LINES.append(gene_filtered[el])

    dict_debug = {}
    LINES = sorted(LINES, key=lambda x : x[2])
    LINES.insert(0, ["ENSEMBL ID", "UNIPROT GENE NAME", "NCBI DESC","HASH NSS","NCBI SUBUNIT STRUCT (NSS)","CORUM COMPLEX ID","CORUM COMPLEX NAME",
                 "GENE IN SAME HASH","GENE IN SAME COMPLEX","NAME REFERENCED IN"])
    for i, row in enumerate(LINES):
        genes_in_same_hash = ''
        genes_in_same_complex = ''
        name_referenced_in = ''
        if i != 0:
            ens_id = row[0]
            le_hash = row[3]
            nss = row[4]
            complex_id_list = row[5]

            if len(hash_ens_dict[le_hash]) > 2:
                genes_in_same_hash = "; ".join({e for e in hash_ens_dict[le_hash] if e != ens_id})
            if complex_id_list:
                r = []
                for cplxid in complex_id_list:
                    if len(complex_id_ens_dict[cplxid]) > 2:
                        for e in {e for e in complex_id_ens_dict[cplxid] if e != ens_id}:
                            r.append(e)
                if r:
                    genes_in_same_complex = "; ".join({e for e in r})
            if ens_id in id_name_found_NSS:
                name_referenced_in = "; ".join(id_name_found_NSS[ens_id])

            if row[-1]:
                row[-1] = "; ".join(row[-1])
                row[-2] = "; ".join(row[-2])
        d = row + [genes_in_same_hash, genes_in_same_complex, name_referenced_in]
        l = "\t".join(d)
        print (l)
        dict_debug[row[0]] = d# store for debug

    '''print "======= debug ============="
    print dict_debug['ENSG00000181090']
    print "complex_ids %s " % dict_debug['ENSG00000181090'][5]
    for cplxid in  dict_debug['ENSG00000181090'][5].split('; '):
        print cplxid, complex_id_ens_dict[cplxid]
    print "same complex_ids %s " % dict_debug['ENSG00000181090'][8]
    same_complex_ids = dict_debug['ENSG00000181090'][8]
    for el in same_complex_ids.split('; '):
        print "checking", el, ens_uniprot_dict[el]
        for el2 in ens_ebi_complex_dict[el]:
            cplxid = el2[0]
            print "%s is in complex %s" % (el, cplxid)
            print cplxid, complex_id_ens_dict[cplxid]
            for el3 in complex_id_ens_dict[cplxid]:
                print el3, ens_uniprot_dict[el3]'''

def get_corum_GR_rules(reactions_file, genes_file, CORUM_file):
    '''
        reactions files: csv export of the sheet RXN of the HMR excel file.
        genes_file: csv export of the sheet GTPRs of the HMR excel file.
        CORUM_file: see parse_CORUM_txt
    '''
    gene_filtered, uniprot_ens_dict, ens_uniprot_dict, hash_ens_dict = filter_genes(genes_file, [])
    ens_ebi_complex_dict, complex_id_ens_dict = parse_CORUM_txt(CORUM_file, uniprot_ens_dict)

    with open(reactions_file, 'r') as fh:
        for line in fh:
            if line[0] == "#":
                continue
            line = line.strip('\n')
            linearr = line.split('\t')
            uniprot_ids = []
            corum_ids = []
            current_ens_ids = linearr[5].split(',')
            same_complex = {}
            test = []

            for ens_id in current_ens_ids:
                ens_id = ens_id.strip()
                if ens_id in ens_uniprot_dict:
                    uni_ids = set(ens_uniprot_dict[ens_id])
                    if len(uni_ids) > 1:
                        print ("error")
                        print (ens_id, ens_uniprot_dict[ens_id])
                        exit()
                    for el in uni_ids:
                        uniprot_ids.append(el)
                if len(current_ens_ids) > 1 and ens_id in ens_ebi_complex_dict:
                    for cpx in ens_ebi_complex_dict[ens_id]:
                        cpx_id = cpx[0]
                        if cpx_id not in same_complex:
                            same_complex[cpx_id] = set()
                        same_complex[cpx_id].add(ens_id)
                    corum_ids.append([e[0] for e in ens_ebi_complex_dict[ens_id]])

            if same_complex:
                for k, v in same_complex.items():
                    if len(v) == 1:
                        del same_complex[k]
            if same_complex:
                list_sets = same_complex.values()
                if len(list_sets) > 1:
                    # merge set if share 1 element
                    while True:
                        break # not merge group where shring one gene id
                        breaked = False
                        for a, b in itertools.combinations(list_sets, 2):
                            if a & b:
                                c = a.union(b)
                                list_sets.remove(a)
                                list_sets.remove(b)
                                list_sets.append(c)
                                breaked = True
                                break
                        if not breaked:
                            break

                same_complex = list_sets

            relation = []
            if same_complex:
                or_relation = []
                processed_id = {}
                for el in same_complex:
                    for e in el:
                        processed_id[e] = None
                    s = "(" + " AND ".join(el) + ")"
                    or_relation.append(s)
                relation = " OR ".join(or_relation)
                for el in current_ens_ids:
                    if el not in processed_id:
                        processed_id[el] = None
                        relation += " OR " + el
            else:
                relation = " OR ".join(current_ens_ids)

            l = [linearr[1], ", ".join(current_ens_ids), ", ".join(uniprot_ids), ", ".join(["/".join(e) for e in corum_ids]), relation]
            print ("\t".join(l))


def get_EBI_GR_rules(reactions_file, genes_file, EBI_file):
    '''
        reactions files:  csv export of the sheet RXN of the HMR excel file.
        genes_file: csv export of the sheet GTPRs of the HMR excel file.
        EBI_file: see parse_EBI_txt
    '''
    gene_filtered, uniprot_ens_dict, ens_uniprot_dict, hash_ens_dict = filter_genes(genes_file, [])
    ens_ebi_complex_dict, complex_id_ens_dict = parse_EBI_txt(EBI_file, uniprot_ens_dict)

    with open(reactions_file, 'r') as fh:
        for line in fh:
            if line[0] == "#":
                continue
            line = line.strip('\n')
            linearr = line.split('\t')
            uniprot_ids = []
            ebi_ids = []
            current_ens_ids = linearr[5].split(',')
            same_complex = {}
            test = []

            for ens_id in current_ens_ids:
                ens_id = ens_id.strip()
                if ens_id in ens_uniprot_dict:
                    uni_ids = set(ens_uniprot_dict[ens_id])
                    if len(uni_ids) > 1:
                        print ("error")
                        print (ens_id, ens_uniprot_dict[ens_id])
                        exit()
                    # store all uniprot ids
                    for el in uni_ids:
                        uniprot_ids.append(el)
                if len(current_ens_ids) > 1 and ens_id in ens_ebi_complex_dict:
                    # group the list of ensembl ids that belong to the same complex id
                    """ same_complex = {
                        "complex0" : [ENSBL000UUUUU, ENSBL000VVVVV],
                        "complex1" : [ENSBL000XXXXX, ENSBL000YYYYY],
                        "complex2" : [ENSBL000ZZZZZ],
                    }
                    """
                    for cpx_id in ens_ebi_complex_dict[ens_id]:
                        if cpx_id not in same_complex:
                            same_complex[cpx_id] = set()
                        same_complex[cpx_id].add(ens_id)
                    # store all ebi ids
                    ebi_ids.append(ens_ebi_complex_dict[ens_id])

            if same_complex:
                """
                    delete comlexId having only one ensembl ID = no relation
                    same_complex = {
                        "complex0" : [ENSBL000UUUUU, ENSBL000VVVVV],
                        "complex1" : [ENSBL000XXXXX, ENSBL000YYYYY],
                    }
                """
                for k, v in same_complex.items():
                    if len(v) == 1:
                        del same_complex[k] # might gives error with python3
            if same_complex:
                list_sets = same_complex.values()
                if len(list_sets) > 1:
                    # merge set if share 1 element
                    while True:
                        break # skip merge group where shring one gene id
                        breaked = False
                        for a, b in itertools.combinations(list_sets, 2):
                            if a & b:
                                c = a.union(b)
                                list_sets.remove(a)
                                list_sets.remove(b)
                                list_sets.append(c)
                                breaked = True
                                break
                        if not breaked:
                            break

                same_complex = list_sets

            relation = []
            if same_complex:
                or_relation = []
                processed_id = {}
                for el in same_complex:
                    for e in el:
                        processed_id[e] = None
                    # all element in the same complex id have a AND relationship
                    """same_complex = {
                        "complex0" : [ENSBL000UUUUU, ENSBL000VVVVV],
                        "complex1" : [ENSBL000XXXXX, ENSBL000YYYYY],
                    } => ["(ENSBL000UUUUU AND ENSBL000VVVVV)", "(ENSBL000XXXXX AND ENSBL000YYYYY)"]
                    """
                    s = "(" + " AND ".join(el) + ")"
                    or_relation.append(s)
                """
                    merge 'AND groups with OR relationship
                    ["(ENSBL000UUUUU AND ENSBL000VVVVV)", "(ENSBL000XXXXX AND ENSBL000YYYYY)"]
                    = >
                    "(ENSBL000UUUUU AND ENSBL000VVVVV) OR (ENSBL000XXXXX AND ENSBL000YYYYY)"
                """
                relation = " OR ".join(or_relation)
                """
                    add other ids, that have not been used yet
                    "(ENSBL000UUUUU AND ENSBL000VVVVV) OR (ENSBL000XXXXX AND ENSBL000YYYYY)"
                    =>
                    "(ENSBL000UUUUU AND ENSBL000VVVVV) OR (ENSBL000XXXXX AND ENSBL000YYYYY) OR ENSBL000ZZZZZ"
                """
                for el in current_ens_ids:
                    if el not in processed_id:
                        processed_id[el] = None
                        relation += " OR " + el
            else:
                relation = " OR ".join(current_ens_ids)

            l = [linearr[1], ", ".join(current_ens_ids), ", ".join(uniprot_ids), ", ".join(["/".join(e) for e in ebi_ids]), relation]
            print ("\t".join(l))


def get_BOTH_GR_rules(reactions_file, genes_file, both_cpl_file):
    '''
        reactions files:  csv export of the sheet RXN of the HMR excel file.
        genes_file: csv export of the sheet GTPRs of the HMR excel file.
        both_cpl_file: file generate by Hao that combines EBI and CORUM complex IDs
    '''

    gene_filtered, uniprot_ens_dict, ens_uniprot_dict, hash_ens_dict = filter_genes(genes_file, [])
    ens_complex_dict, complex_id_ens_dict = parse_EBI_CORUM_file(both_cpl_file, uniprot_ens_dict)

    with open(reactions_file, 'r') as fh:
        for line in fh:
            if line[0] == "#":
                continue
            line = line.strip('\n')
            linearr = line.split('\t')
            uniprot_ids = []
            ebi_ids = []
            current_ens_ids = linearr[5].split(',')
            same_complex = {}
            test = []

            for ens_id in current_ens_ids:
                ens_id = ens_id.strip()
                if ens_id in ens_uniprot_dict:
                    uni_ids = set(ens_uniprot_dict[ens_id])
                    if len(uni_ids) > 1:
                        print ("error")
                        print (ens_id, ens_uniprot_dict[ens_id])
                        exit()
                    for el in uni_ids:
                        uniprot_ids.append(el)
                if len(current_ens_ids) > 1 and ens_id in ens_complex_dict:
                    for cpx_id in ens_complex_dict[ens_id]:
                        if cpx_id not in same_complex:
                            same_complex[cpx_id] = set()
                        same_complex[cpx_id].add(ens_id)
                    ebi_ids.append(ens_complex_dict[ens_id])

            if same_complex:
                # might gives error with python3
                for k, v in same_complex.items():
                    if len(v) == 1:
                        del same_complex[k]
            if same_complex:
                list_sets = same_complex.values()
                if len(list_sets) > 1:
                    # merge set if share 1 element
                    while True:
                        break # skip merge group when sharing one gene id
                        breaked = False
                        for a, b in itertools.combinations(list_sets, 2):
                            if a & b:
                                c = a.union(b)
                                list_sets.remove(a)
                                list_sets.remove(b)
                                list_sets.append(c)
                                breaked = True
                                break
                        if not breaked:
                            break

                same_complex = list_sets

            relation = []
            if same_complex:
                or_relation = []
                processed_id = {}
                for el in same_complex:
                    for e in el:
                        processed_id[e] = None
                    s = "(" + " AND ".join(el) + ")"
                    or_relation.append(s)
                relation = " OR ".join(or_relation)
                for el in current_ens_ids:
                    if el not in processed_id:
                        processed_id[el] = None
                        relation += " OR " + el
            else:
                relation = " OR ".join(current_ens_ids)

            l = [linearr[1], ", ".join(current_ens_ids), ", ".join(uniprot_ids), ", ".join(["/".join(e) for e in ebi_ids]), relation]
            print ("\t".join(l))

if __name__ == "__main__":

    # ens_ebi_complex_dict = parse_CORUM_txt(sys.argv[1], sys.argv[2])
    # add_corum_cols_to_interest_genes_file(sys.argv[3], ens_ebi_complex_dict)

    # build_genes_complex_file(sys.argv[1], sys.argv[2])

    # print ("reactions_file, genes_file, CORUM_file")
    # python genes_relation.py ../RXN-20180122.csv GTPRs_20171212.txt Corum_allComplexes.txt
    # get_corum_GR_rules(sys.argv[1], sys.argv[2], sys.argv[3])

    #print ("reactions_file, genes_file, EBI_file")
    # python genes_relation.py ../RXN-20180122.csv GTPRs_20171212.txt EBI_complex_homo_sapiens_20171218.tsv
    #get_EBI_GR_rules(sys.argv[1], sys.argv[2], sys.argv[3])

    print ("reactions_file, genes_file, both_cpl_file")
    # python genes_relation.py ../RXN-20180122.csv GTPRs_20171212.txt both_cplx_file.txt > both_cpl_genes_rules_results.txt
    get_BOTH_GR_rules(sys.argv[1], sys.argv[2], sys.argv[3])





