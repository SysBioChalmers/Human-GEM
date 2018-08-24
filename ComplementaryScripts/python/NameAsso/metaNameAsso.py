import sys
import parse_HMDB
import parse_ChEBI
import parse_LIPIDMAPS
import pubchem
import parse_MNX
import os

'''
Using the name of the metabolite, search correspoding ID in CHEBI, LIPIDMAPS, HMDB, PUBCHEM, and MNX databases

output 2 files :

- a matrix of all matching IDs, with columns:
Name    HMDB (perfect match)    HMDB (synonyms)    CHEBI (perfect match)    CHEBI (synonyms)    LIPIDMAPS (perfect match)    LIPIDMAPS (synonyms)    PUBCHEM (perfect match)    PUBCHEM (synonyms)    MNX    MNX bigg    MNX chebi    MNX hmdb    MNX kegg    MNX lipidmaps
the MNX external IDs can be compare with the perfect match


- a table of MNX entry with the HMR ID associated with columns:
META ID   MNX_ID  Description Formula Charge  Mass    InChI   SMILES  Source  InChIKey
'''

file_path = os.path.dirname(os.path.abspath(__file__))

# the following input files should be regenerated

# from http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip (v4.0)
# then generated using parse_HMDB.parse_HMDB_XML_file() and write_HMDB_tab()
HMDB_file = os.path.join(file_path, 'input', 'HMDB.tab')

# from ftp://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete.sdf.gz
# then generated using parse_ChEBI.parse_CHEBI_SDF_file() and write_CHEBI_tab()
CHEBI_file = os.path.join(file_path, 'input', 'CHEBI.tab')

# from http://www.lipidmaps.org/resources/downloads/LMSDFDownload12Dec17.tar.gz
# then generated using parse_LIPIDMAPS.parse_LIPIDMAPS_SDF_file() and write_LIPIDMAPS_tab()
LIPIDMAPS_file = os.path.join(file_path, 'input', 'LIPIDMAPS.tab')

# from ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Synonym-filtered.gz
PUBCHEM_file = os.path.join(file_path, 'input', 'CID-Synonym-filtered')
# SQLlite database generated using pubchem.generate_pubchem_db() 
PUBCHEM_db_file = os.path.join(file_path, 'input', 'pubchem.db')
PUBCHEM_db_file = os.path.join('input', 'pubchem.db')

if not os.path.isfile(PUBCHEM_db_file):
    pubchem.generate_pubchem_db(PUBCHEM_file)

# https://www.metanetx.org/mnxdoc/mnxref.html
# download chem_xref.tsv and chem_prop.tsv
# light version contains only identifiers for kegg bigg hmdb chebi pubchem lipidmaps..
MNX_xref_file = os.path.join(file_path, 'input', 'chem_xref_light.tsv')
MNX_prop_file = os.path.join(file_path, 'input', 'chem_prop.tsv')

# from HMA repository, csv export of the metabolite sheet of the excel HMR database
metabolite_file = os.path.join(file_path, 'input', 'human_metaboliteListFromExcel.txt')

MNX_dict = parse_MNX.parse_MNX_xrefs(MNX_xref_file)
MNX_name, MNX_chebi, MNX_hmdb, MNX_bigg, MNX_lipidmaps, MNX_kegg = parse_MNX.convert_MNX_dict(MNX_dict)

CHEBI_dict, CHEBI_secondary, CHEBI_name, CHEBI_IUPAC_name, CHEBI_synonyms, CHEBI_lipidmaps, CHEBI_hmdb, CHEBI_pubchem  \
    = parse_ChEBI.parse_CHEBI_tab(CHEBI_file, synonyms_as_dict=True, lipidmaps_as_dict=True, pubchem_as_dict=True, hmdb_as_dict=True)

HMDB, HMDB_secondary, HMDB_name, HMDB_iupac_name, HMDB_synonyms, HMDB_chebi, HMDB_pubchem = \
    parse_HMDB.parse_HMDB_tab(HMDB_file, synonyms_as_dict=True, chebi_as_dict=True, pubchem_as_dict=True)

LIPIDMAPS, LIPIDMAPS_name, LIPIDMAPS_systematic_name, LIPIDMAPS_synonyms, LIPIDMAPS_chebi, LIPIDMAPS_pubchem, LIPIDMAPS_hmdb = \
    parse_LIPIDMAPS.parse_LIPIDMAPS_tab(LIPIDMAPS_file, synonyms_as_dict=True, chebi_as_dict=True, pubchem_as_dict=True, hmdb_as_dict=True)

conn, cur = pubchem.get_connection_cursor(PUBCHEM_db_file)

id_matrix_output_file = sys.argv[1]
match_id_prop_output_file = sys.argv[2]


# pair hmr meta id - mnx id
match_pair = []

# search only once each name, the metabolite_file contains duplicate names due to metabolite located in multiple compartments
unique_name = set()
c = 0

with open(metabolite_file, 'r') as fh, open(id_matrix_output_file, 'w') as fw: # open metabolite file
    fw.write("Name\tHMDB (perfect match)\tHMDB (synonyms)" + \
          "\tCHEBI (perfect match)\tCHEBI (synonyms)" + \
          "\tLIPIDMAPS (perfect match)\tLIPIDMAPS (synonyms)" + \
          "\tPUBCHEM (perfect match)\tPUBCHEM (synonyms)" + \
          "\tMNX\tMNX bigg\tMNX chebi\tMNX hmdb\tMNX kegg\tMNX lipidmaps\n")
    for line in fh:
        if line[0] == '#':
            continue
        linearr = line.split('\t')
        meta_name = linearr[1][:-3]
        hmr_id = linearr[8]  # e.g. m00002c

        if meta_name in unique_name:
            continue

        unique_name.add(meta_name)
        lower_name = meta_name.lower()
        c += 1

        # HMDB
        HMDB_ids = []
        HMDB_ids_exact = []
        if lower_name in HMDB_name:
            HMDB_ids_exact += HMDB_name[lower_name]
        if lower_name in HMDB_iupac_name:
            HMDB_ids_exact += HMDB_iupac_name[lower_name]
        if lower_name in HMDB_synonyms:
            HMDB_ids += HMDB_synonyms[lower_name]
        HMDB_ids = set(HMDB_ids)
        HMDB_ids_exact = set(HMDB_ids_exact)
        HMDB_ids = HMDB_ids - HMDB_ids_exact

        # CHEBI
        CHEBI_ids = []
        CHEBI_ids_exact = []
        if lower_name in CHEBI_name:
            CHEBI_ids_exact += CHEBI_name[lower_name]
        if lower_name in CHEBI_IUPAC_name:
            CHEBI_ids_exact += CHEBI_IUPAC_name[lower_name]
        if lower_name in CHEBI_synonyms:
            CHEBI_ids += CHEBI_synonyms[lower_name]
        CHEBI_ids = set(CHEBI_ids)
        CHEBI_ids_exact = set(CHEBI_ids_exact)
        CHEBI_ids = CHEBI_ids - CHEBI_ids_exact

        # lipidmaps
        LIPIDMAPS_ids = []
        LIPIDMAPS_ids_exact = []
        if lower_name in LIPIDMAPS_name:
            LIPIDMAPS_ids_exact += LIPIDMAPS_name[lower_name]
        if lower_name in LIPIDMAPS_systematic_name:
            LIPIDMAPS_ids_exact += LIPIDMAPS_systematic_name[lower_name]
        if lower_name in LIPIDMAPS_synonyms:
            LIPIDMAPS_ids += LIPIDMAPS_synonyms[lower_name]
        LIPIDMAPS_ids = set(LIPIDMAPS_ids)
        LIPIDMAPS_ids_exact = set(LIPIDMAPS_ids_exact)
        LIPIDMAPS_ids = LIPIDMAPS_ids - LIPIDMAPS_ids_exact

        # pubchem
        pubchem_result_exact = pubchem.get_full_entry_from_name(lower_name, cur, pos=1)
        pubchem_results = pubchem.get_full_entry_from_name(lower_name, cur)
        if pubchem_results:
            for k in pubchem_result_exact:
				# remove results overlap
                if k in pubchem_results:
                    del pubchem_results[k]

        pubchem_ids_exact = set([str(e) for e in pubchem_result_exact.keys()])
        pubchem_ids = set([str(e) for e in pubchem_results.keys()])

        # MNX
        MNX_ids = []
        MNX_bigg = []
        MNX_chebi = []
        MNX_hmdb = []
        MNX_kegg = []
        MNX_lipidmaps = []
        
        if lower_name in MNX_name:
            MNX_ids += MNX_name[lower_name]
            MNX_ids = set(MNX_ids)
            for MNX_id in MNX_ids:
                MNX_bigg += MNX_dict[MNX_id]['bigg']
                MNX_chebi += MNX_dict[MNX_id]['chebi']
                MNX_hmdb += MNX_dict[MNX_id]['hmdb']
                MNX_kegg += MNX_dict[MNX_id]['kegg']
                MNX_lipidmaps += MNX_dict[MNX_id]['lipidmaps']

                match_pair.append([hmr_id, MNX_id])

        MNX_ids = set(MNX_ids)

        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                meta_name,
                "; ".join(HMDB_ids_exact),
                "; ".join(HMDB_ids),
                "; ".join(CHEBI_ids_exact),
                "; ".join(CHEBI_ids),
                "; ".join(LIPIDMAPS_ids_exact),
                "; ".join(LIPIDMAPS_ids),
                "; ".join(pubchem_ids_exact),
                "; ".join(pubchem_ids),
                "; ".join(MNX_ids),
                "; ".join(MNX_bigg),
                "; ".join(MNX_chebi),
                "; ".join(MNX_hmdb),
                "; ".join(MNX_kegg),
                "; ".join(MNX_lipidmaps),
            ))

print ("Metabolite count", c)

MNX_prop_dict = parse_MNX.parse_MNX_prop(MNX_prop_file)

with open(match_id_prop_output_file, 'w') as fw:
    fw.write("META ID\tMNX_ID\tDescription\tFormula\tCharge\tMass\tInChI\tSMILES\tSource\tInChIKey\n")
    for hmr_id, mnx_id in match_pair:
        fw.write("%s\t%s\t%s\n" % (hmr_id, mnx_id, "\t".join(MNX_prop_dict[mnx_id])))

