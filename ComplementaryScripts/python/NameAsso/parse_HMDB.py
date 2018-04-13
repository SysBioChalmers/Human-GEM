import re

def parse_HMDB_XML_file(file):
    HMDB = {}
    HMDB_secondary = {}
    HMDB_synonyms = {}

    accession = None
    sec_accession = []
    name = None
    description = ''
    synonyms = []
    formula = None
    avg_weight = None
    mono_weight = None
    iupac_name = None
    iupac_trad_name = None
    functions = []
    inchi = None
    kegg = None
    bigg = None
    pubchem = None
    chebi = None

    parse_description = False
    parse_synonyms = False
    parse_sec_accession = False
    parse_function = False

    with open(file, 'r') as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            # print (line)
            if line.startswith("<metabolite>") or line.startswith("</hmdb>"):
                if not accession:
                    continue
                # store metabolite
                accession = accession.strip() if accession else accession
                name = name.strip() if name else name
                iupac_name = iupac_name.strip() if iupac_name else iupac_name

                synonyms = [el.strip() for el in synonyms]

                v = {
                    'accession': accession,
                    'name': name,
                    'description': description.strip() if description else description,
                    'sec_accession': sec_accession,
                    'formula': formula.strip() if formula else formula,
                    'avg_weight': avg_weight.strip() if avg_weight else avg_weight,
                    'mono_weight': mono_weight.strip() if mono_weight else mono_weight,
                    'iupac_name': iupac_name,
                    # 'iupac_trad_name': iupac_trad_name.strip() if iupac_trad_name else iupac_trad_name,
                    'functions': functions,
                    'inchi': inchi.strip() if inchi else inchi,
                    'bigg': bigg.strip() if bigg else bigg,
                    'chebi': chebi.strip() if chebi else chebi,
                    'kegg': kegg.strip() if kegg else kegg,
                    'pubchem': pubchem.strip() if pubchem else pubchem,
                    'synonyms': synonyms
                }

                missing_values = [e for e in ['accession', 'name', 'description'] if not v[e]]
                if len(missing_values):
                    print ("Warning: missing information for HMDB", missing_values)
                    print (v)
                    if accession == "HMDB0013326" and missing_values == ['description']:
                        v['description'] = "Trans-2-dodecenoylcarnitine is classified as a member of the fatty acid esters. Fatty acid esters are carboxylic ester derivatives of a fatty acid. Trans-2-dodecenoylcarnitine is considered to be a practically insoluble (in water) and a weak acidic compound. Trans-2-dodecenoylcarnitine can be found in Cow milk, pasteurized, vitamin A + D added, 0% fats, Cow milk, pasteurized, vitamin A + D added, 1% fats, Cow milk, pasteurized, vitamin A + D added, 2% fats, and Cow milk, pasteurized, vitamin D added, 3.25% fats. Trans-2-dodecenoylcarnitine can be found in blood and urine. Within a cell, Trans-2-dodecenoylcarnitine is primarily located in the extracellular space and near the membrane."
                        print ("fixed")
                    else:
                        continue

                HMDB[accession] = v

                for sec_acc in sec_accession:
                    if sec_acc in HMDB:
                        print ("Error: sec accession '%s' already in HMDB dict (%s)" % (sec_acc, HMDB[sec_acc]["name"]))
                        exit(1)

                    if sec_acc in HMDB_secondary:
                        print ("Error: sec accession '%s' already in HMDB_sec dict (%s)" % (sec_acc, HMDB_secondary[sec_acc]["name"]))
                        exit(1)

                    HMDB_secondary[sec_acc] = HMDB[accession]

                if name in HMDB_synonyms:
                    print ("Warning: name '%s' (%s) already in HMDB_synonyms dict (%s)" % (name, accession, HMDB_synonyms[name]["accession"]))
                    # exit(1)
                HMDB_synonyms[name] = HMDB[accession]

                if iupac_name:
                    if iupac_name in HMDB_synonyms and HMDB_synonyms[iupac_name]["accession"] != accession:
                        print ("warning: iupac_name '%s' (%s) already in HMDB_synonyms dict (%s)" % (iupac_name, accession, HMDB_synonyms[iupac_name]["accession"]))
                        # exit(1)
                    HMDB_synonyms[iupac_name] = HMDB[accession]

                '''if iupac_trad_name:
                    if iupac_trad_name in HMDB_synonyms and HMDB_synonyms[iupac_trad_name]["accession"] != accession:
                        print ("Error: iupac_trad_name '%s' (%s) already in HMDB_synonyms dict (%s)" % (iupac_trad_name, accession, HMDB_synonyms[iupac_trad_name]["accession"]))
                        exit(1)
                    HMDB_synonyms[iupac_trad_name] = HMDB[accession]'''

                for synonym in synonyms:
                    if synonym in HMDB_synonyms and HMDB_synonyms[synonym]["accession"] != accession:
                        print ("Warning: synonym '%s' (%s) already in HMDB_synonyms dict (%s)" % (synonym, accession, HMDB_synonyms[synonym]["accession"]))
                        # exit(1)
                    HMDB_synonyms[synonym] = HMDB[accession]


                #reset
                accession = None
                sec_accession = []
                name = None
                description = ''
                synonyms = []
                formula = None
                avg_weight = None
                mono_weight = None
                iupac_name = None
                iupac_trad_name = None
                functions = []
                inchi = None
                kegg = None
                bigg = None
                pubchem = None
                chebi = None

                parse_description = False
                parse_synonyms = False
                parse_sec_accession = False
                parse_function = False
            else:
                try:
                    if parse_sec_accession:
                        if line.startswith("</secondary_accessions>"):
                            parse_sec_accession = False
                            continue
                        sec_accession.append(re.match('<accession>([^<]+)</accession>', line).group(1))
                    elif parse_synonyms:
                        if line.startswith("</synonyms>"):
                            parse_synonyms = False
                            continue
                        elif line and not line.startswith("</synonym>"):
                            synonyms.append(re.match('<synonym>([^<]*)(?:</synonym>)?', line).group(1))
                    elif parse_description:
                        # print line
                        if line.endswith("</description>"):
                            if line != "</description>":
                                description += re.match('(?:<description>)?([^<>]+)</description>', line).group(1)
                            parse_description = False
                            continue
                        elif line:
                            description += re.match('(?:<description>)?([^<>]+)(</description>)?', line).group(1)
                    elif parse_function:
                        if line.startswith("</biofunctions>"):
                            parse_function = False
                            continue
                        if line != "</biofunction>":
                            functions.append(re.match('<biofunction>([^<]*)(?:</biofunction>)?', line).group(1))

                    elif line.startswith("<accession>"):
                        accession = re.match('<accession>([^<]+)</accession>', line).group(1)
                    elif line.startswith("<secondary_accessions>"):
                        parse_sec_accession = True
                    elif line.startswith("<name>") and not name:
                        name = re.match('<name>([^<]+)</name>', line).group(1)
                    elif line.startswith("<description>") and not description:
                            if line.endswith("</description>"):
                                description = re.match('<description>([^<>]*)</description>', line).group(1)
                            else:
                                description = re.match('<description>([^<>]*)', line).group(1)
                                parse_description = True
                    elif line.startswith("<synonyms>"):
                        parse_synonyms = True
                    elif line.startswith("<chemical_formula>"):
                        formula = re.match('<chemical_formula>([^<]+)</chemical_formula>', line).group(1)
                    elif line.startswith("<average_molecular_weight>"):
                        avg_weight = re.match('<average_molecular_weight>([^<]+)</average_molecular_weight>', line).group(1)
                    elif line.startswith("<monisotopic_molecular_weight>"):
                        mono_weight = re.match('<monisotopic_molecular_weight>([^<]+)</monisotopic_molecular_weight>', line).group(1)
                    elif line.startswith("<iupac_name>"):
                        iupac_name = re.match('<iupac_name>([^<]+)</iupac_name>', line).group(1)
                    # elif line.startswith("<traditional_iupac>"):
                    #     iupac_trad_name = re.match('<traditional_iupac>([^<]+)</traditional_iupac>', line).group(1)
                    elif line.startswith("<biofunctions>"):
                        parse_function = True

                    elif line.startswith("<inchi>"):
                        inchi = re.match('<inchi>([^<]*)</inchi>', line).group(1)
                    elif line.startswith("<kegg_id>"):
                        kegg = re.match('<kegg_id>([^<]*)</kegg_id>', line).group(1)
                    elif line.startswith("<bigg_id>"):
                        bigg = re.match('<bigg_id>([^<]*)</bigg_id>', line).group(1)
                    elif line.startswith("<pubchem_compound_id>"):
                        pubchem = re.match('<pubchem_compound_id>([^<]+)</pubchem_compound_id>', line).group(1)
                    elif line.startswith("<chebi_id>"):
                        chebi = re.match('<chebi_id>([^<]+)</chebi_id>', line).group(1)
                except Exception as e:
                    print (e)
                    print (line)
                    print (i)
                    exit(1)


    return HMDB, HMDB_secondary, HMDB_synonyms


def write_HMDB_tab(output_file, HMDB):
    with open(output_file, 'w') as fw:
        header = ['accession',
                    'name',
                    'description',
                    'sec_accession',
                    'formula',
                    'avg_weight',
                    'mono_weight',
                    'iupac_name',
                    'functions',
                    'inchi',
                    'bigg',
                    'chebi',
                    'kegg',
                    'pubchem',
                    'synonyms']
        fw.write("\t".join(header) + "\n")
        for k, v in HMDB.items():
            for col in header:
                if not v[col]:
                    v[col] = ''
                if col == "sec_accession" or col == "functions":
                    fw.write("; ".join(v[col]) + "\t")
                elif col == "synonyms":
                    fw.write(";; ".join(v[col]) + "\n") # use double ; separator
                else:
                    fw.write(v[col] + "\t")


def parse_HMDB_tab(file, synonyms_as_dict=False, chebi_as_dict=False, pubchem_as_dict=False):
    HMDB = {}
    HMDB_secondary = {}
    HMDB_name = {}
    HMDB_iupac_name = {}
    HMDB_synonyms = {}
    HMDB_chebi = {}
    HMDB_pubchem = {}
    with open(file, 'r') as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            linearr = line.strip("\n").split("\t")
            d = {}
            for i, k in enumerate(header):
                d[k] = linearr[i]
            d['sec_accession'] = [el.strip() for el in d['sec_accession'].split("; ") if el.strip()]
            d['synonyms'] = [el.strip() for el in d['synonyms'].strip().split(";; ") if el.strip()]
            HMDB_id = d["accession"]
            HMDB[HMDB_id] = d

            for el in d['sec_accession']:
                if el not in HMDB_secondary:
                    HMDB_secondary[el] = []
                HMDB_secondary[el].append(HMDB_id)

            if d['name']:
                if d['name'].lower() not in HMDB_name:
                    HMDB_name[d['name'].lower()] = [HMDB_id]
                else:
                    HMDB_name[d['name'].lower()].append(HMDB_id)

            if d['iupac_name']:
                if d['iupac_name'].lower() not in HMDB_iupac_name:
                    HMDB_iupac_name[d['iupac_name'].lower()] = [HMDB_id]
                else:
                    HMDB_iupac_name[d['iupac_name'].lower()].append(HMDB_id)

            if synonyms_as_dict and d['synonyms']:
                for el in d['synonyms']:
                    el = el.lower()
                    if el not in  HMDB_synonyms:
                        HMDB_synonyms[el] = []
                    HMDB_synonyms[el].append(HMDB_id)

            if chebi_as_dict and d["chebi"]:
                if d["chebi"] in HMDB_chebi:
                    HMDB_chebi[d["chebi"]].append(HMDB_id)
                else:
                    HMDB_chebi[d["chebi"]] = [HMDB_id]

            if pubchem_as_dict and d["pubchem"]:
                if d["pubchem"] in HMDB_pubchem:
                    HMDB_pubchem[d["pubchem"]].append(HMDB_id)
                else:
                    HMDB_pubchem[d["pubchem"]] = [HMDB_id]

    return HMDB, HMDB_secondary, HMDB_name, HMDB_iupac_name, HMDB_synonyms, HMDB_chebi, HMDB_pubchem


if __name__ == "__main__":
    import sys
    # HMDB, HMDB_secondary, HMDB_synonyms = parse_HMDB_XML_file(sys.argv[1])
    # write_HMDB_tab(sys.argv[2], HMDB)
    HMDB, HMDB_secondary, HMDB_name, HMDB_iupac_name, HMDB_synonyms, HMDB_chebi, HMDB_pubchem = parse_HMDB_tab(sys.argv[2], pubchem_as_dict=True)
    # parse_HMDB_tab(sys.argv[1])

    if 'acetate' in HMDB_name:
        print (HMDB_name['acetate'])

    if 'acetate' in HMDB_iupac_name:
        print (HMDB_iupac_name['acetate'])

    if 'acetate' in HMDB_synonyms:
        print (HMDB_synonyms['acetate'])

    print (HMDB.values()[0])

    print (HMDB_pubchem['5281780'])