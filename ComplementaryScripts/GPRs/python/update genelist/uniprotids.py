# Code to generate Uniprot ID lists from each GPR of each HMR reaction, started: 2018-03-14

# 1.1 GENERATE DICTIONARY FOR GENES VERSUS UNIPROT PROTEIN IDs
# Let's generate a dictionary with gene list versus UNIPROT ID columns in the GTPRs page; i.e., (Ensembl: Gene ID) and (UniProt: UniProt ID) columns of GTPRs page
 

geneuniprotid = open(r'C:\Users\kocabas\Desktop\geneuniprot.txt', 'r')

trans = {}
for line in geneuniprotid:
	parts=line.rstrip()
	gene=parts[0:15]
	id=parts[16:26]
	trans[gene]=id
	

b='ENSG00000004779' in trans
#print(b)
#print(trans)

geneuniprotid.close()

# 1.2 CHANGE GPRs with UNIPROT IDs
# change list elements using dictionary 


with open(r"C:\Users\kocabas\Desktop\ngpr.csv", "r") as ins: # ngpr is abbreviation for new gprs
	uniprot = []
	for line in ins:
		ens = line.rstrip().split(',')
		for i in range(len(ens)):
			current_ens_id = ens[i]
			ens[i] = trans[current_ens_id]
			#if current_ens_id != trans[current_ens_id]:
				#print ("change %s %s" % (current_ens_id, trans[current_ens_id]))
			
		uniprot.append(ens)

print(uniprot)

output_file = r'C:\Users\kocabas\Desktop\uniprot_ids.txt'
with open(output_file, 'w') as fw:
	for el in uniprot:
		fw.write(','.join(el)+'\n')