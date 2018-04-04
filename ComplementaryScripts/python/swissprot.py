# Cofactors of HMR proteins were retrieved from Swissprot database on 2018-04-03.
# So swissprot.csv file was generated using the data retrieved from Swissprot database and it contains human proteins' Uniprot IDs and their cofactors
# Therefore as a first step, a dictionary was generated using  uniprot_ids and cofactors taken from swissprot.csv file.

d = {}
f=open(r"C:\Users\kocabas\Desktop\swissprot.csv", "r")

for line in f:
	line = line.strip('\n')# remove line break at the end of the line
	array = line.split('\t') # create columns
	uniprot_id = array[0]
	cofactor_col = array[6]
	d[uniprot_id]=cofactor_col.strip()
#print(d['Q9UKJ8'])


# Uniprot ids of each HMR reaction were exchanged with cofactors using the dictionary.
# STEP 2: Change gpr_uniprotids with cofactors and generate gpr_cofactors

with open(r"C:\Users\kocabas\Desktop\gpruniprotids.csv", "r") as ins:
	cofactor = []
	for line in ins:
		uni = line.rstrip().split(',')
		for i in range(len(uni)):
			id = uni[i].strip()
			if id in d:
				if d[id]:
					cofactor_str = d[id]
					line = line.replace(id, '[' + cofactor_str + ']')
				else:
					line = line.replace(id, '')
		cofactor.append(line)


#STEP 3: Write the output to a file

output_file = r'C:\Users\kocabas\Desktop\cofactors.txt'
with open(output_file, 'w') as fw:
	for el in cofactor:
		fw.write(el)

