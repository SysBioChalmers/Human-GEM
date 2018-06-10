
def parse_MNX_xrefs(file):
    res = {}
    with open(file, 'r') as fh:
        #XREF   MNX_ID  Evidence    Description
        for i, line in enumerate(fh):
            if line[0] == "#":
                continue
            linearr = line.split('\t')
            linearr[3] = linearr[3].strip()
            if line[:4] == "MNXM":
                if linearr[0] != linearr[1]:
                    print (line, i)
                    exit(1)

            elif linearr[1][:4] == "MNXM":
                name = [e.strip() for e in linearr[3].split('|') if e.strip()]
                MNX_id = linearr[1]
                # if MNX_id in res:
                #     print "%s in dict" % (MNX_id)
                #    exit(1)

                if MNX_id not in res:
                    res[MNX_id] = {
                        'bigg': [],
                        'chebi': [],
                        'hmdb': [],
                        'kegg': [],
                        'lipidmaps': [],
                        'names' : set()
                    }
                for el in name:
                     res[MNX_id]['names'].add(el)
                if line[:4] == "bigg":
                    res[MNX_id]['bigg'].append(linearr[0][5:].strip())
                elif line[:5] == "chebi":
                    res[MNX_id]['chebi'].append(linearr[0][6:].strip())
                elif line[:4] == "hmdb":
                    res[MNX_id]['hmdb'].append(linearr[0][5:].strip())
                elif line[:4] == "kegg":
                    res[MNX_id]['kegg'].append(linearr[0][5:].strip())
                elif line[:9] == "lipidmaps":
                    res[MNX_id]['lipidmaps'].append(linearr[0][10:].strip())
                else:
                    continue

    '''print (res['MNXM11971'])
    print (res['MNXM10232'])
    print (res['MNXM480330'])
    print (res['MNXM6'])
    print (res['MNXM587669'])'''

    return res

def parse_MNX_prop(file):
    prop = {}
    with open(file, 'r') as fh:
        #MNX_ID Description Formula Charge  Mass    InChI   SMILES  Source  InChIKey
        for i, line in enumerate(fh):
            if line[0] == "#":
                continue
            linearr = line.split('\t')
            linearr[-1] = linearr[-1].strip()
            prop[linearr[0]] = linearr[1:]

    # for k, v in prop.items()[:10]:
    #    print k, v

    return prop



def convert_MNX_dict(d):

    MNX_name = {}
    MNX_chebi = {}
    MNX_hmdb = {}
    MNX_bigg = {}
    MNX_lipidmaps = {}
    MNX_kegg = {}
    for key, value in d.items():
        if value['names']:
            for name in value['names']:
                namel = name.lower()
                if namel not in MNX_name:
                    MNX_name[namel] = [key]
                else:
                    MNX_name[namel].append(key)

        if value['chebi']:
            for chebi_id in value['chebi']:
                if chebi_id not in MNX_chebi:
                    MNX_chebi[chebi_id] = [key]
                else:
                    MNX_chebi[chebi_id].append(key)

        if value['bigg']:
            for bigg_id in value['bigg']:
                if bigg_id not in MNX_bigg:
                    MNX_bigg[bigg_id] = [key]
                else:
                    MNX_bigg[bigg_id].append(key)

        if value['hmdb']:
            for hmdb_id in value['hmdb']:
                if hmdb_id not in MNX_hmdb:
                    MNX_hmdb[hmdb_id] = [key]
                else:
                    MNX_hmdb[hmdb_id].append(key)

        if value['lipidmaps']:
            for lipidmaps_id in value['lipidmaps']:
                if lipidmaps_id not in MNX_lipidmaps:
                    MNX_lipidmaps[lipidmaps_id] = [key]
                else:
                    MNX_lipidmaps[lipidmaps_id].append(key)

        if value['kegg']:
            for kegg_id in value['kegg']:
                if kegg_id not in MNX_kegg:
                    MNX_kegg[kegg_id] = [key]
                else:
                    MNX_kegg[kegg_id].append(key)

    return MNX_name, MNX_chebi, MNX_hmdb, MNX_bigg, MNX_lipidmaps, MNX_kegg


if __name__ == "__main__":
    # https://www.metanetx.org/mnxdoc/mnxref.html
    # download chem_xref.tsv and chem_prop.tsv
    import sys
    MNX_dict = parse_MNX_xrefs(sys.argv[1])
    MNX_name, MNX_chebi, MNX_hmdb, MNX_bigg, MNX_lipidmaps, MNX_kegg = convert_MNX_dict(MNX_dict)

    print ("==============================")
    v = MNX_name.keys()[0]
    print (MNX_dict[MNX_name[v][0]])
    print ("==============================")
    v = MNX_chebi.keys()[0]
    print (MNX_dict[MNX_chebi[v][0]])
    print ("==============================")
    v = MNX_hmdb.keys()[0]
    print (MNX_dict[MNX_hmdb[v][0]])
    print ("==============================")
    v = MNX_bigg.keys()[0]
    print (MNX_dict[MNX_bigg[v][0]])
    print ("==============================")
    v = MNX_lipidmaps.keys()[0]
    print (MNX_dict[MNX_lipidmaps[v][0]])
    print ("==============================")
    v = MNX_kegg.keys()[0]
    print (MNX_dict[MNX_kegg[v][0]])

    MNX_prop_dict = parse_MNX_prop(sys.argv[2])

    print (MNX_prop_dict[MNX_dict.keys()[0]])