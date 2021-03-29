import MySQLdb
import os
import sys
import argparse

"""
- Fetch all the IDs (gene_stable_id, transcript_stable_id, translation_stable_id, uniprot_id, gene_name, ncbi_id)
  from ensembl Human database
- Fetch the list of Primary assembly gene IDs from the same database
- Write the OUTPUT File with the list of IDS, excluding the no-primary assembly IDs

Input:
 - version of ensembl, e.g. 95
 - version of human assembly, e.g. 38
 - name of the output file

Output:
 - the IDs mapping written in the output file
"""

parser = argparse.ArgumentParser(description='Ensembl IDs fetcher')
parser.add_argument('ensembl version', action="store", type=int,
 help="version of Ensembl, is used to connect to an Ensembl's database. e.g. 95")
parser.add_argument('genome version', action="store", type=int,
 help="version the assembly, is used to connect to an Ensembl's database. e.g. 38")
parser.add_argument('output file', action="store", help="path of the output file")
parser.add_argument('--database', action="store", dest="dbname", help="Ensembl's database name to connect with, if provided 'ensembl version' and 'genome version' are ignored")

args = vars(parser.parse_args())
ensembl_version = args['ensembl version']
genome_version = args['genome version']
output_file = args['output file']
dbname = args['dbname']

database_name = 'homo_sapiens_core_%s_%s' % (ensembl_version, genome_version)
if dbname:
    database_name = dbname

connN = MySQLdb.connect(host="ensembldb.ensembl.org",
                        user="anonymous",
                        db=database_name)
curN = connN.cursor()

"""
# retrieve all the ids
# the SQL query was provided by ensembl helpdesk (ticket #323315)
# the content can be also obtain using the biomart service on the ensembl website by this url:
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

curN.execute("""SELECT DISTINCT
    gene.stable_id AS gene_stable_id,
    transcript.stable_id AS transcript_stable_id,
    translation.stable_id AS translation_stable_id,
    ups.display_label AS uniprot_id,
    xref.display_label AS gene_name,
    eg.dbprimary_acc AS ncbi_id
FROM
    gene
    -- gene_name is a special case of xrefs, as it is directly linked
    -- instead of being linked through the object_xref table
    JOIN xref ON gene.display_xref_id=xref.xref_id
    JOIN transcript ON gene.gene_id=transcript.gene_id
    JOIN translation ON transcript.transcript_id=translation.transcript_id
    -- join Uniprot IDs to the translations, as above
    LEFT JOIN (
        SELECT
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
        SELECT
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
            external_db.db_name='EntrezGene') AS eg ON gene.gene_id=eg.gene_id""")

data = curN.fetchall()
ids_data = []
for row in data:
    row = ["" if e == None else e for e in row] # replace None (null) by empty string
    ids_data.append(row)

# sort rows
ids_data.sort(key=lambda x: "%s%s%s%s%s%s" % (x[0], x[1], x[2], x[3], x[4], x[5]))

# get primary assembly gene ids
# the query was confirmed by Ensembl helpdesk (ticket #320660) to be a correct way to get the list primary assembly genes ids
curN.execute("select gene.stable_id from gene, seq_region where gene.seq_region_id = seq_region.seq_region_id and seq_region. name regexp '^([[:digit:]]+|MT|X|Y)$'")
data = curN.fetchall()

primary_assembly = set()
for row in data:
    primary_assembly.add(row[0])

if False:
    # build the primary assembly - alternate map
    # can convert no-PA id to PA ID
    # not used, all the no-primary assembl gene are removed
    curN.execute("select alt_allele_group_id, group_concat(gene.stable_id) from alt_allele, gene where alt_allele.gene_id = gene.gene_id group by alt_allele_group_id")
    data = curN.fetchall()

    primary_assembly_map = {}
    for row in data:
        #print row[0]
        pa = None
        not_pa = []
        for gene_id in row[1].split(','):  # all gene_id in the current alt_allele_group_id, one of them may be the primary assembly id
            if gene_id in primary_assembly:
                if pa:
                    # check if there is only on PA from each alt_allele_group_id
                    print "Error: 2 primary assembly genes found: " + pa + " / " + gene_id
                    exit(1)
                pa = gene_id
            else:
                not_pa.append(gene_id)
        if pa:
            for gene_id in not_pa:
                primary_assembly_map[gene_id] = pa

try:
    with open(output_file, 'w') as fw:
        fw.write("Gene_stable_ID\tTranscript_stable_ID\tProtein_stable_ID\tUniProtKB_Swiss_Prot_ID\tGene_name\tNCBI_gene_ID\n")
        for ids in ids_data:
            ids = ["" if e == None else e for e in ids] # replace None (null) by empty string
            gene_stable_id, transcript_stable_id, translation_stable_id, uniprot_id, gene_name, ncbi_id = ids
            if gene_stable_id not in primary_assembly:
                continue
            fw.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_stable_id, transcript_stable_id, translation_stable_id, uniprot_id, gene_name, ncbi_id))
except Exception as e:
    print (e)

