import mysql.connector
import os
import sys
import argparse
import requests
import re
from collections import OrderedDict

"""
Fetch and write in OUTPUT FILE gene annotations of the input YAML file. The file must contains Ensembl IDs as gene IDs.
Input:
 - the YAML file format of the GEM
 - version of Ensembl, e.g. 103
 - version of Ensembl, e.g. 103
 - version of human assembly, e.g. 38
 - name of the output file

Steps:
- Fetch all the IDs (gene_stable_id, transcript_stable_id, translation_stable_id, uniprot_id, gene_name, ncbi_id)
    and the description and the synonyms from Ensembl Human database
- Fetch the list of Primary assembly gene IDs from the same database
- Fetch the gene ids in the input YAML FILE
- Write the OUTPUT FILE with the list of IDs only existing in YAML FILE and excluding the no-primary assembly IDs

Output:
 - data written in the TSV output file

Require mysql-connector-python package.
"""

def get_latest_ensembl_info():
    """
    Retrieves - using the REST api - the latest version of the Ensembl (e.g. 103)
    and the ending number version of the latest Human genome version supported (e.g. 38 of GRCh38).
    The two numbers are used to connect to the MySQL Ensembl database.
    """
    def fetch_data_json(url):
        r = requests.get(url, headers={"Content-Type": "application/json"})
        if r.status_code != 200:
            r.raise_for_status()
            exit(1)
        return r.json()

    resp = fetch_data_json("https://rest.ensembl.org/info/data/?")
    release_version = resp["releases"][0]
    print('Ensembl release version found: %d' % release_version)

    resp = fetch_data_json("https://rest.ensembl.org/info/genomes/taxonomy/Homo sapiens?")
    assembly_default = resp[0]["assembly_default"]
    genome_version = int(re.search('[0-9]+', assembly_default).group(0))
    print('Ensembl genome version found: %d' % genome_version)

    return release_version, genome_version

def get_human_database_name(ensembl_version, genome_version):
    return 'homo_sapiens_core_%s_%s' % (ensembl_version, genome_version)


def get_ensembl_db_connection(database_name):
    conn = mysql.connector.connect(host="ensembldb.ensembl.org",
                                   user="anonymous",
                                   database=database_name)
    return conn


def get_yaml_gene_ids(yaml_file):
    """
    Get the list of gene IDs in the model using quick parsing
    The cobra package could be used as well, but it might be overkill
    """

    model_gene_ids_dict = OrderedDict()  # used it as OrderedSet
    try:
        with open(yaml_file, "r") as fh:
            pg = False
            for line in fh:
                line = line.strip()
                if line == '- genes:':
                    pg = True
                elif pg and line.startswith('- id:'):
                    model_gene_ids_dict[line.split(": ", 1)[1].strip('"')] = None
    except Exception as e:
        print(e)
        print("Error: cannot parse yaml file")
        exit(1)
    return model_gene_ids_dict


def retrieve_ensembl_gene_annotations(connection):
    """
    The SQL query was provided by Ensembl helpdesk (ticket #323315).
    The annotations (the IDs but not the last 2 fields 'description' and 'synonym') can be also obtained \
      using the biomart service on the Ensembl website using this url:
    #  http://www.ensembl.org/biomart/martview/825e18e7c897cd07721eddc59c127c4d?
        VIRTUALSCHEMANAME=default&ATTRIBUTES=
        hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|
        hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|
        hsapiens_gene_ensembl.default.feature_page.ensembl_peptide_id|
        hsapiens_gene_ensembl.default.feature_page.uniprotswissprot|
        hsapiens_gene_ensembl.default.feature_page.external_gene_name|
        hsapiens_gene_ensembl.default.feature_page.entrezgene&FILTERS=&VISIBLEPANEL=
        resultspanel
    """

    # Fetch all the Human (not human-GEM!) genes
    cur = connection.cursor()
    cur.execute("set session group_concat_max_len = 4096")
    cur.execute("""SELECT DISTINCT
        gene.stable_id AS gene_stable_id,
        group_concat(DISTINCT transcript.stable_id SEPARATOR ';') AS transcript_stable_id,
        group_concat(DISTINCT translation.stable_id SEPARATOR ';') AS translation_stable_id,
        group_concat(DISTINCT ups.display_label SEPARATOR ';') AS uniprot_id,
        group_concat(DISTINCT xref.display_label SEPARATOR ';') AS gene_name,
        group_concat(DISTINCT eg.dbprimary_acc SEPARATOR ';') AS ncbi_id,
        gene.description as description,
        (select group_concat(external_synonym.synonym SEPARATOR ';')
         from external_synonym where external_synonym.xref_id = gene.display_xref_id
        ) as synonyms
    FROM
        gene
        -- gene_name is a special case of xrefs, as it is directly linked
        -- instead of being linked through the object_xref table
        JOIN xref ON gene.display_xref_id=xref.xref_id
        JOIN transcript ON gene.gene_id=transcript.gene_id
        JOIN translation ON transcript.transcript_id=translation.transcript_id
        -- join Uniprot IDs to the translations, as above
        LEFT JOIN (
            SELECT DISTINCT
                translation.translation_id,
                xref.display_label
            FROM
                     translation
                JOIN object_xref ON (
                    translation.translation_id=object_xref.ensembl_id AND
                    object_xref.ensembl_object_type='Translation'
                )
                JOIN xref ON object_xref.xref_id=xref.xref_id
                JOIN external_db ON xref.external_db_id=external_db.external_db_id
            WHERE
                external_db.db_name='Uniprot/SWISSPROT') AS ups ON translation.translation_id=ups.translation_id
        -- join NCBI gene IDs to the genes, as above
        LEFT JOIN (
            SELECT DISTINCT
                gene.gene_id,
                xref.dbprimary_acc
            FROM
                     gene
                JOIN object_xref ON (
                    gene.gene_id=object_xref.ensembl_id AND
                    object_xref.ensembl_object_type='Gene'
                )
                JOIN xref ON object_xref.xref_id=xref.xref_id
                JOIN external_db ON xref.external_db_id=external_db.external_db_id
            WHERE
                external_db.db_name='EntrezGene') AS eg ON gene.gene_id=eg.gene_id
        GROUP BY 1""")

    return cur.fetchall()


def get_primary_assembly_ids(connection):
    """
    Get primary assembly gene IDs
    the query was validated by Ensembl helpdesk (ticket #320660)
    """
    # get primary assembly gene ids
    # the query was confirmed by Ensembl helpdesk (ticket #320660) to be a correct way
    # to get the list primary assembly genes ids
    cur = connection.cursor()
    cur.execute("select gene.stable_id from gene, seq_region where \
                gene.seq_region_id = seq_region.seq_region_id and \
                seq_region. name regexp '^([[:digit:]]+|MT|X|Y)$'")
    return {row[0] for row in cur.fetchall()}


def create_gene_annotation_file(yaml_file, output_file, ensembl_version=None, genome_version=None, gene_ids_list=None):
    if not ensembl_version or not genome_version:
        ensembl_version, genome_version = get_latest_ensembl_info()
    database_name = get_human_database_name(ensembl_version, genome_version)
    return create_annotation_file(yaml_file, database_name, output_file, gene_ids_list=gene_ids_list)


def create_annotation_file(yaml_file, database_name, output_file, gene_ids_list=None):
    """
    Returns the list of gene IDs written in the output file
    """
    if gene_ids_list:
        model_gene_ids_dict = OrderedDict({gid: None for gid in gene_ids_list})
    else:
        model_gene_ids_dict = get_yaml_gene_ids(yaml_file)

    connection = get_ensembl_db_connection(database_name)
    data = retrieve_ensembl_gene_annotations(connection)
    genes_data_dict = {}
    for row in data:
        row = ["" if e is None else e for e in row]  # replace None (null) by empty string
        genes_data_dict[row[0]] = row

    primary_assembly = get_primary_assembly_ids(connection)
    connection.close()
    new_list = []
    try:
        with open(output_file, 'w') as fw:
            fw.write("genes\tgeneENSTID\tgeneENSPID\tgeneUniProtID\tgeneSymbols\t"
                     "geneEntrezID\tgeneNames\tgeneAliases\n")
            for gid in model_gene_ids_dict:
                if gid not in genes_data_dict:
                    print("Error: could not retrieve Ensembl annotations for gene ID '%s'" % gid)
                    continue
                ensembl_gene_data = genes_data_dict[gid]
                data = ['' if e is None else e for e in ensembl_gene_data]  # replace None (null) by empty string
                gene_stable_id, transcript_stable_ids, translation_stable_ids, \
                    uniprot_ids, gene_symbol, ncbi_ids, gene_desc, gene_synonyms = data
                if gene_stable_id not in primary_assembly:
                    print("Warning: gene ID '%s' is not on the primary assembly" % gid)
                    continue

                new_list.append(gid)

                # sorted the IDs using python. It seems to significantly increase the query time when performed
                # at the database level, in the group_concat().
                # transcript_stable_ids = ";".join(sorted(transcript_stable_ids.split(";")))
                # translation_stable_ids = ";".join(sorted(translation_stable_ids.split(";")))
                gene_desc = gene_desc.split("[Source:")[0].strip()  # remove the source part from the description
                # remove the version history. UniProt IDs without versioning do not contain any '.'.
                uniprot_ids = ";".join([e.split(".")[0] for e in uniprot_ids.split(";")])

                fw.write('"%s"\t"%s"\t"%s"\t"%s"\t"%s"\t"%s"\t"%s"\t"%s"\n' %
                         (gene_stable_id, transcript_stable_ids, translation_stable_ids,
                             uniprot_ids, gene_symbol, ncbi_ids, gene_desc, gene_synonyms))
    except Exception as e:
        print(e)
        exit(1)

    return new_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fetch Human genes information through the Ensembl\'s MySQL database')
    parser.add_argument('YAML file', action="store", help="Human-GEM YAML file")
    parser.add_argument('output file', action="store", help="path of the output file")
    parser.add_argument('--ensembl-version', action="store", type=int,
                        help="version of Ensembl, is used to connect to an Ensembl's database. e.g. 103")
    parser.add_argument('--genome-version', action="store", type=int,
                        help="version the assembly, is used to connect to an Ensembl's database. e.g. 38")
    parser.add_argument('--database', action="store", dest="db_name",
                        help="Ensembl's database name to connect with, if provided 'Ensembl version' and"
                             " 'genome version' are ignored")

    args = vars(parser.parse_args())
    yaml_file = args['YAML file']
    ensembl_version = args['ensembl_version']
    genome_version = args['genome_version']
    if (not ensembl_version and genome_version) or (not genome_version and ensembl_version):
        print("Warning: both ENSEMBL VERSION and GENOME VERSION must be provided in order to connect to the database.")
        print("The provided value '%d' is ignored, the lastest available Ensembl realease and"
              " genome version will be used." % (ensembl_version if ensembl_version else genome_version))
    output_file = args['output file']
    db_name = args['db_name']

    if not db_name:
        if not ensembl_version or not genome_version:
            ensembl_version, genome_version = get_latest_ensembl_info()
        db_name = get_human_database_name(ensembl_version, genome_version)

    create_annotation_file(yaml_file, db_name, output_file)
