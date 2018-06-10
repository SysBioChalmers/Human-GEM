import re

def parse_LIPIDMAPS_SDF_file(file):
    LIPIDMAPS = {}
    LIPIDMAPS_synonyms = {}

    accession = None
    common_name = None
    systematic_name = None
    synonyms = []
    formula = None
    weight = None
    inchi = None
    pubchem = None
    HMDB = None
    chebi = None
    kegg = None
    status = False

    with open(file, 'r') as fh:
        lines = fh.readlines()
        for i, line in enumerate(lines):
            line = line.strip()

            # print (line)
            if line.startswith("$$$$"):
                if not accession:
                    continue
                # store metabolite
                accession = accession.strip() if accession else accession
                common_name = common_name.strip() if common_name else common_name
                systematic_name = systematic_name.strip() if systematic_name else systematic_name

                v = {
                    'accession': accession,
                    'name': common_name,
                    'systematic_name': systematic_name,
                    'formula': formula.strip() if formula else formula,
                    'weight': weight.strip() if weight else weight,
                    'inchi': inchi.strip() if inchi else inchi,
                    'chebi': chebi.strip() if chebi else chebi,
                    'hmdb': HMDB.strip() if HMDB else HMDB,
                    'kegg': kegg.strip() if kegg else kegg,
                    'pubchem': pubchem.strip() if pubchem else pubchem,
                    'status' : str(status),
                    'synonyms': synonyms
                }

                missing_values = [e for e in ['accession', 'systematic_name'] if not v[e]]
                if len(missing_values):
                    print ("Warning: missing information for LIPIDMAPS", missing_values)
                    print (v)
                    continue

                LIPIDMAPS[accession] = v

                if common_name:
                    if common_name in LIPIDMAPS_synonyms:
                        print ("Warning: name '%s' (%s) already in LIPIDMAPS_synonyms dict (%s)" % \
                            (common_name, accession, LIPIDMAPS_synonyms[common_name]["accession"]))
                        # exit(1)
                    LIPIDMAPS_synonyms[common_name] = LIPIDMAPS[accession]

                if systematic_name:
                    if systematic_name in LIPIDMAPS_synonyms and LIPIDMAPS_synonyms[systematic_name]["accession"] != accession:
                        print ("warning: systematic_name '%s' (%s) already in LIPIDMAPS_synonyms dict (%s)" % \
                            (systematic_name, accession, LIPIDMAPS_synonyms[systematic_name]["accession"]))
                        # exit(1)
                    LIPIDMAPS_synonyms[systematic_name] = LIPIDMAPS[accession]

                for synonym in synonyms:
                    if synonym in LIPIDMAPS_synonyms and LIPIDMAPS_synonyms[synonym]["accession"] != accession:
                        print ("Warning: synonym '%s' (%s) already in LIPIDMAPS_synonyms dict (%s)" % (synonym, accession, LIPIDMAPS_synonyms[synonym]["accession"]))
                        # exit(1)
                    LIPIDMAPS_synonyms[synonym] = LIPIDMAPS[accession]

                #reset
                accession = None
                common_name = None
                systematic_name = None
                synonyms = []
                formula = None
                weight = None
                inchi = None
                pubchem = None
                HMDB = None
                chebi = None
                kegg = None
                status = False
            else:
                try:
                    if line.startswith("> <LM_ID>"):
                        accession = lines[i+1]
                    elif line.startswith("> <COMMON_NAME>"):
                        common_name = lines[i+1]
                    elif line.startswith("> <SYSTEMATIC_NAME>"):
                        systematic_name = lines[i+1]
                    elif line.startswith("> <SYNONYMS>"):
                        synonyms = [e.strip() for e in lines[i+1].strip().split(';')]
                    elif line.startswith("> <EXACT_MASS>"):
                        weight = lines[i+1]
                    elif line.startswith("> <FORMULA>"):
                        formula = lines[i+1]
                    elif line.startswith("> <PUBCHEM_CID>"):
                        pubchem = lines[i+1]
                    elif line.startswith("> <KEGG_ID>"):
                        kegg = lines[i+1]
                    elif line.startswith("> <CHEBI_ID>"):
                        chebi = lines[i+1]
                    elif line.startswith("> <HMDBID>"):
                        HMDB = lines[i+1]
                    elif line.startswith("> <INCHI>"):
                        inchi = lines[i+1]
                    elif line.startswith("> <STATUS>"):
                        status = True if lines[i+1].strip().startswith("Active") else False
                except Exception as e:
                    print (e)
                    print (line)
                    print (i)
                    exit(1)

    return LIPIDMAPS


def write_LIPIDMAPS_tab(output_file, LIPIDMAPS):
    with open(output_file, 'w') as fw:
        header = ['accession',
                    'name',
                    'systematic_name',
                    'formula',
                    'weight',
                    'inchi',
                    'chebi',
                    'kegg',
                    'pubchem',
                    'hmdb',
                    'status',
                    'synonyms']
        fw.write("\t".join(header) + "\n")
        for k, v in LIPIDMAPS.items():
            for col in header:
                if not v[col]:
                    v[col] = ''
                if col == "synonyms":
                    fw.write(";; ".join(v[col]) + "\n") # use double ; separator
                else:
                    fw.write(v[col] + "\t")


def parse_LIPIDMAPS_tab(file, synonyms_as_dict=False, chebi_as_dict=False, pubchem_as_dict=False, hmdb_as_dict=False):
    LIPIDMAPS = {}
    LIPIDMAPS_name = {}
    LIPIDMAPS_systematic_name = {}
    LIPIDMAPS_synonyms = {}
    LIPIDMAPS_chebi = {}
    LIPIDMAPS_pubchem = {}
    LIPIDMAPS_hmdb = {}

    with open(file, 'r') as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            linearr = line.strip("\n").split("\t")
            d = {}
            for i, k in enumerate(header):
                d[k] = linearr[i]
            d['synonyms'] = [el.strip() for el in d['synonyms'].strip().split(";; ") if el.strip()]
            LIPIDMAPS_id = d["accession"]
            LIPIDMAPS[LIPIDMAPS_id] = d

            if synonyms_as_dict:
                for el in d['synonyms']:
                    el = el.lower()
                    if el not in LIPIDMAPS_synonyms:
                        LIPIDMAPS_synonyms[el] = []
                    LIPIDMAPS_synonyms[el].append(LIPIDMAPS_id)

            if d['name']:
                if d['name'].lower() not in LIPIDMAPS_name:
                    LIPIDMAPS_name[d['name'].lower()] = [LIPIDMAPS_id]
                else:
                    LIPIDMAPS_name[d['name'].lower()].append(LIPIDMAPS_id)

            if d['systematic_name']:
                if d['systematic_name'].lower() not in LIPIDMAPS_systematic_name:
                    LIPIDMAPS_systematic_name[d['systematic_name'].lower()] = [LIPIDMAPS_id]
                else:
                    LIPIDMAPS_systematic_name[d['systematic_name'].lower()].append(LIPIDMAPS_id)

            if chebi_as_dict and d["chebi"]:
                if d["chebi"] in LIPIDMAPS_chebi:
                    LIPIDMAPS_chebi[d["chebi"]].append(LIPIDMAPS_id)
                else:
                    LIPIDMAPS_chebi[d["chebi"]] = [LIPIDMAPS_id]

            if pubchem_as_dict and d["pubchem"]:
                if d["pubchem"] in LIPIDMAPS_pubchem:
                    LIPIDMAPS_pubchem[d["pubchem"]].append(LIPIDMAPS_id)
                else:
                    LIPIDMAPS_pubchem[d["pubchem"]] = [LIPIDMAPS_id]

            if hmdb_as_dict and d["hmdb"]:
                if d["hmdb"] in LIPIDMAPS_hmdb:
                    LIPIDMAPS_hmdb[d["hmdb"]].append(LIPIDMAPS_id)
                else:
                    LIPIDMAPS_hmdb[d["hmdb"]] = [LIPIDMAPS_id]

    return LIPIDMAPS, LIPIDMAPS_name, LIPIDMAPS_systematic_name, LIPIDMAPS_synonyms, LIPIDMAPS_chebi, LIPIDMAPS_pubchem, LIPIDMAPS_hmdb


if __name__ == "__main__":
    import sys
    # LIPIDMAPS = parse_LIPIDMAPS_SDF_file(sys.argv[1])
    # write_LIPIDMAPS_tab(sys.argv[2], LIPIDMAPS)
    LIPIDMAPS, LIPIDMAPS_name, LIPIDMAPS_systematic_name, LIPIDMAPS_synonyms, LIPIDMAPS_chebi, LIPIDMAPS_pubchem, LIPIDMAPS_hmdb = parse_LIPIDMAPS_tab(sys.argv[2], hmdb_as_dict=True)
    # print (LIPIDMAPS_hmdb)

    print (LIPIDMAPS_hmdb['HMDB0030705'])