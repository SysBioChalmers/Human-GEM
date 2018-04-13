import re

def parse_CHEBI_SDF_file(file):
    CHEBI = {}
    CHEBI_secondary = {}
    CHEBI_synonyms = {}

    accession = None
    star = None
    sec_accession = []
    name = None
    IUPAC_name = []
    synonyms = []
    formula = None
    weight = None
    mono_weight = None
    inchi = None
    pubchem = []
    HMDB = []
    lipidmaps = []
    kegg = []

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
                name = name.strip() if name else name

                v = {
                    'accession': accession,
                    'sec_accession': sec_accession,
                    'star': star.strip() if star else star,
                    'name': name,
                    'iupac_name': IUPAC_name,
                    'formula': formula.strip() if formula else formula,
                    'weight': weight.strip() if weight else weight,
                    'mono_weight': mono_weight.strip() if mono_weight else mono_weight,
                    'inchi': inchi.strip() if inchi else inchi,
                    'hmdb': HMDB,
                    'kegg': kegg,
                    'pubchem': pubchem,
                    'lipidmaps': lipidmaps,
                    'synonyms': synonyms
                }

                missing_values = [e for e in ['accession', 'name'] if not v[e]]
                if len(missing_values):
                    print ("Warning: missing information for HMDB", missing_values)
                    print (v)
                    continue

                CHEBI[accession] = v

                for sec_acc in sec_accession:
                    if sec_acc in CHEBI:
                        print ("Error: sec accession '%s' already in CHEBI dict (%s)" % (sec_acc, CHEBI[sec_acc]["name"]))
                        exit(1)

                    if sec_acc in CHEBI_secondary:
                        print ("Error: sec accession '%s' already in CHEBI_sec dict (%s)" % (sec_acc, CHEBI_secondary[sec_acc]["name"]))
                        exit(1)

                    CHEBI_secondary[sec_acc] = CHEBI[accession]

                if name:
                    if name in CHEBI_synonyms:
                        print ("Warning: name '%s' (%s) already in CHEBI_synonyms dict (%s)" % \
                            (name, accession, CHEBI_synonyms[name]["accession"]))
                        # exit(1)
                    CHEBI_synonyms[name] = CHEBI[accession]

                if IUPAC_name:
                    for el in IUPAC_name:
                        if el in CHEBI_synonyms and CHEBI_synonyms[el]["accession"] != accession:
                            print ("warning: IUPAC_name '%s' (%s) already in CHEBI_synonyms dict (%s)" % \
                                (el, accession, CHEBI_synonyms[el]["accession"]))
                            # exit(1)
                        CHEBI_synonyms[el] = CHEBI[accession]

                for synonym in synonyms:
                    if synonym in CHEBI_synonyms and CHEBI_synonyms[synonym]["accession"] != accession:
                        print ("Warning: synonym '%s' (%s) already in CHEBI_synonyms dict (%s)" % (synonym, accession, CHEBI_synonyms[synonym]["accession"]))
                        # exit(1)
                    CHEBI_synonyms[synonym] = CHEBI[accession]

                #reset
                accession = None
                star = None
                sec_accession = []
                name = None
                IUPAC_name = []
                synonyms = []
                formula = None
                weight = None
                mono_weight = None
                inchi = None
                pubchem = []
                HMDB = []
                lipidmaps = []
                kegg = []
            else:
                try:
                    if line.startswith("> <ChEBI ID>"):
                        accession = lines[i+1].lstrip("CHEBI:")
                    if line.startswith("> <Star>"):
                        star = lines[i+1]
                    elif line.startswith("> <ChEBI Name>"):
                        name = lines[i+1]
                    elif line.startswith("> <Secondary ChEBI ID>"):
                        l = i
                        while lines[l+1].startswith("CHEBI"):
                            sec_accession.append(lines[l+1].lstrip("CHEBI:").strip())
                            l += 1
                    elif line.startswith("> <IUPAC Names>"):
                        l = i
                        while not lines[l+1].startswith(">") and lines[l+1].strip():
                            IUPAC_name.append(lines[l+1].strip())
                            l += 1
                    elif line.startswith("> <Synonyms>"):
                        l = i
                        while not lines[l+1].startswith(">") and lines[l+1].strip():
                            synonyms.append(lines[l+1].strip())
                            l += 1
                    elif line.startswith("> <Mass>"):
                        weight = lines[i+1]
                    elif line.startswith("> <Monoisotopic Mass>"):
                        mono_weight = lines[i+1]
                    elif line.startswith("> <Formulae>"):
                        formula = lines[i+1]
                    elif line.startswith("> <PubChem Database Links>"):
                        l = i
                        while not lines[l+1].startswith(">") and lines[l+1].strip():
                            if lines[l+1].startswith("CID"):
                                pubchem.append(lines[l+1].lstrip("CID:").strip())
                            l += 1
                    elif line.startswith("> <KEGG COMPOUND Database Links>"):
                        l = i
                        while not lines[l+1].startswith(">") and lines[l+1].strip():
                            kegg.append(lines[l+1].strip())
                            l += 1
                    elif line.startswith("> <HMDB Database Links>"):
                        l = i
                        while lines[l+1].startswith("HMDB") and lines[l+1].strip():
                            HMDB.append(lines[l+1].strip())
                            l += 1
                    elif line.startswith("> <LIPID MAPS instance Database Links>"):
                        l = i
                        while lines[l+1].startswith("LM") and lines[l+1].strip():
                            lipidmaps.append(lines[l+1].strip())
                            l += 1
                    elif line.startswith("> <InChI>"):
                        inchi = lines[i+1]
                except Exception as e:
                    print (e)
                    print (line)
                    print (i)
                    exit(1)


    print (CHEBI['17336']) # multi lipidsmaps, secid and kegg
    print (CHEBI['17374'])  # multi pubchem
    print (CHEBI['5445']) # 17 synonyms

    return CHEBI


def write_CHEBI_tab(output_file, CHEBI):
    with open(output_file, 'w') as fw:
        header = ['accession',
                    'sec_accession', # list
                    'star',
                    'name',
                    'iupac_name', # list
                    'formula',
                    'weight',
                    'mono_weight',
                    'inchi',
                    'hmdb', # list
                    'kegg', # list
                    'pubchem', # list
                    'lipidmaps', # list
                    'synonyms'] # list
        fw.write("\t".join(header) + "\n")
        for k, v in CHEBI.items():
            for col in header:
                if not v[col]:
                    v[col] = ''
                if col in ["sec_accession", "iupac_name", 'hmdb', 'kegg', 'pubchem', 'lipidmaps']:
                    fw.write("; ".join(v[col]) + "\t")
                elif col == "synonyms":
                    fw.write(";; ".join(v[col]) + "\n") # use double ; separator
                else:
                    fw.write(v[col] + "\t")


def parse_CHEBI_tab(file, synonyms_as_dict=False, lipidmaps_as_dict=False, pubchem_as_dict=False, hmdb_as_dict=False):
    CHEBI = {}
    CHEBI_secondary = {}
    CHEBI_name = {}
    CHEBI_IUPAC_name = {}
    CHEBI_synonyms = {}
    CHEBI_lipidmaps = {}
    CHEBI_pubchem = {}
    CHEBI_hmdb = {}

    with open(file, 'r') as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            linearr = line.strip("\n").split("\t")
            d = {}
            for i, k in enumerate(header):
                d[k] = linearr[i]
            d['sec_accession'] = [el.strip() for el in d['sec_accession'].split("; ") if el.strip()]
            d['iupac_name'] = [el.strip() for el in d['iupac_name'].split("; ") if el.strip()]
            d['lipidmaps'] = [el.strip() for el in d['lipidmaps'].split("; ") if el.strip()]
            d['kegg'] = [el.strip() for el in d['kegg'].split("; ") if el.strip()]
            d['pubchem'] = [el.strip() for el in d['pubchem'].split("; ") if el.strip()]
            d['hmdb'] = [el.strip() for el in d['hmdb'].split("; ") if el.strip()]
            d['synonyms'] = [el.strip() for el in d['synonyms'].strip().split(";; ") if el.strip()]
            CHEBI_id = d["accession"]
            CHEBI[CHEBI_id] = d

            if d['sec_accession']:
                for el in d['sec_accession']:
                    if el not in CHEBI_secondary:
                        CHEBI_secondary[el] = []
                    CHEBI_secondary[el].append(CHEBI_id)

            if d['name']:
                if d['name'].lower() not in CHEBI_name:
                    CHEBI_name[d['name'].lower()] = [CHEBI_id]
                else:
                    CHEBI_name[d['name'].lower()].append(CHEBI_id)

            if d['iupac_name']:
                for el in d['iupac_name']:
                    el = el.lower()
                    if el not in CHEBI_IUPAC_name:
                        CHEBI_IUPAC_name[el] = [CHEBI_id]
                    else:
                        CHEBI_IUPAC_name[el].append(CHEBI_id)

            if synonyms_as_dict:
                for el in d['synonyms']:
                    el = el.lower()
                    if el not in CHEBI_synonyms:
                        CHEBI_synonyms[el] = []
                    CHEBI_synonyms[el].append(CHEBI_id)

            if lipidmaps_as_dict and d["lipidmaps"]:
                for el in d["lipidmaps"]:
                    if el in CHEBI_lipidmaps:
                        CHEBI_lipidmaps[el].append(CHEBI_id)
                    else:
                        CHEBI_lipidmaps[el] = [CHEBI_id]

            if pubchem_as_dict and d["pubchem"]:
                for el in d["pubchem"]:
                    if el in CHEBI_pubchem:
                        CHEBI_pubchem[el].append(CHEBI_id)
                    else:
                        CHEBI_pubchem[el] = [CHEBI_id]

            if hmdb_as_dict and d["hmdb"]:
                for el in d["hmdb"]:
                    if el in CHEBI_hmdb:
                        CHEBI_hmdb[el].append(CHEBI_id)
                    else:
                        CHEBI_hmdb[el] = [CHEBI_id]

    return CHEBI, CHEBI_secondary, CHEBI_name, CHEBI_IUPAC_name, CHEBI_synonyms, CHEBI_lipidmaps, CHEBI_hmdb, CHEBI_pubchem


if __name__ == "__main__":
    import sys
    CHEBI = parse_CHEBI_SDF_file(sys.argv[1])
    write_CHEBI_tab(sys.argv[2], CHEBI)
    CHEBI, CHEBI_secondary, CHEBI_name, CHEBI_IUPAC_name, CHEBI_synonyms, CHEBI_lipidmaps, CHEBI_hmdb, CHEBI_pubchem = parse_CHEBI_tab(sys.argv[2], hmdb_as_dict=True)
    print (CHEBI_hmdb)