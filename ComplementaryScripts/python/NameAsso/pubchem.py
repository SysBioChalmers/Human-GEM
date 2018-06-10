import re
import os
import sqlite3

def generate_pubchem_db(file):
    print ("Generating pubchem database")
    print ("Estimated time: 10-20minutes, output 'pubchem.db' >13GB")
    dbfile = os.path.join(os.path.dirname(__file__), 'pubchem.db')

    if os.path.isfile(dbfile):
        os.unlink(dbfile)

    con = sqlite3.connect(dbfile)
    cur = con.cursor()
    cur.execute("CREATE TABLE alias ("
                "alias TEXT collate nocase NOT NULL, "
                "CID INTEGER NOT NULL, "
                "pos INTEGER NOT NULL)")

    pCID = 0
    pos = 1
    with open(file, 'r') as fh:
        for i, line in enumerate(fh):
            if i != 0 and i % 1000000 == 0:
                print ("Processed %s lines" % i)
            CID, term = line.strip().split('\t')
            if CID != pCID:
                pos = 1

            cur.execute("INSERT INTO alias (alias, CID, pos) VALUES (?, ?, ?)", (term, CID, pos))
            pos += 1
            pCID = CID

    con.commit()
    cur.execute("CREATE INDEX alias_index on alias (alias collate nocase)")
    cur.execute("CREATE INDEX cid_index on alias (cid)")
    con.commit()
    # test 
    '''cur.execute("SELECT * FROM alias WHERE alias = '(2E)-octadecenoyl-CoA'")
    cur.execute("SELECT * FROM alias WHERE alias = 'DL-2-amino-3-mercaptopropionic acid'")
    results = cur.fetchall()
    for r in results:
        print (r)
    cur.execute("SELECT * FROM alias")
    results = cur.fetchall()
    for r in results:
        print (r)'''
    cur.close()
    con.close()


def get_connection_cursor(file):

    def match(expr, item):
        return re.match(expr, item) is not None

    con = sqlite3.connect(file)
    con.create_function("REGEXP", 2, match)
    cur = con.cursor()

    return con, cur


def get_results_as_dict(res):
    result_dict = {}
    for alias, cid, pos in res:
        if cid not in result_dict:
            result_dict[cid] = {
                'kegg': {},
                'chebi': {},
                'lipidmaps': {},
                'hmdb': {},
                'synonyms': {},
                'name': []
            }
        d = result_dict[cid]
        if re.match("C[0-9]{5}", alias):
            d['kegg'][alias] = pos
        elif re.match("CHEBI:[0-9]+", alias):
            d['chebi'][alias[6:]] = pos
        elif re.match("LM[A-Z]{2}[A-Z0-9]+[0-9]{2}", alias):
            d['lipidmaps'][alias] = pos
        elif re.match("HMDB[0-9]{5,}", alias):
            d['hmdb'][alias] = pos
        else:
            d['synonyms'][alias.lower()] = pos

        if pos == 1:
            d['name'].append(alias)

    return result_dict

def get_entry_from_CID(CID, cur, pos=None):
    try:
        int(CID)
    except:
        print ("Error: invalid CID '%s'" % CID)
        exit(1)
    if pos:
        cur.execute("SELECT * FROM alias WHERE CID=%s and pos=%s" % (CID, int(pos)))
    else:
        cur.execute("SELECT * FROM alias WHERE CID=%s" % CID)
    res = cur.fetchall()
    return get_results_as_dict(res)


def get_entry_from_name(name, cur, pos=None):
    if pos:
        cur.execute("SELECT * FROM alias WHERE alias='%s' and pos=%s" % (name.replace('\'', '\'\''), int(pos)))
    else:
        cur.execute("SELECT * FROM alias WHERE alias='%s'" % name.replace('\'', '\'\''))
    res = cur.fetchall()

    return get_results_as_dict(res)


def get_full_entry_from_name(name, cur, pos=None):
    if pos:
        cur.execute("SELECT * FROM alias WHERE alias='%s' and pos=%s" % (name.replace('\'', '\'\''), int(pos)))
    else:
        cur.execute("SELECT * FROM alias WHERE alias='%s'" % name.replace('\'', '\'\''))
    res = cur.fetchall()
    full_res = {}
    for alias, cid, pos in res:
        r = get_entry_from_CID(cid, cur)
        full_res[cid] = r[cid]

    return full_res


def get_entry_from_kegg(kegg_id, cur):
    return get_entry_from_name(kegg_id, cur)

def get_entry_from_lipidmaps(lipidmaps_id, cur):
    return get_entry_from_name(lipidmaps_id, cur)

def get_entry_from_hmdb(hmdb_id, cur):
    return get_entry_from_name(hmdb_id, cur)

def get_entry_from_chebi(chebi_id, cur):
    chebi_id_arr = chebi_id.split(":")
    if len(chebi_id_arr) == 2:
        chebi_id = "CHEBI:%s" % str(int(chebi_id_arr[1]))
    elif len(chebi_id_arr) == 1:
        chebi_id = "CHEBI:%s" % str(int(chebi_id_arr[0]))
    else:
        print ("invalid chebi id '%s'" % chebi_id)
        exit(1)
    return get_entry_from_name(chebi_id, cur)


if __name__ == "__main__":
    import sys
    # generate_pubchem_db(sys.argv[1])
    conn, cur = get_connection_cursor(sys.argv[1])
    print (get_entry_from_name('Methyl 2-[5-(2-methoxyphenyl)tetrazol-2-yl]acetate', cur, pos=None))
    print (get_entry_from_CID(166486, cur, pos=None))

