# Cofactors of HMR proteins were retrieved from Swissprot database on 2018-04-03 in the form of excell file.
# The excell file was converted to a csv file, which is swissprot.csv file.
# Columns 1 and 7 of Swissprot.csv file have Uniprot IDs and cofactors for HMR proteins, respectively.
# Therefore, columns 1 and 7 of Swissprot.csv file were respectively used as keys and values in the process of generating dictionary d .
#array[0] is column 1 of Swissprot.csv file; array[1] is column 7 of Swissprot.csv file.

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
# gpruniprotids.csv file is mathematically a matrix; a list of uniprot ids for each HMR reaction
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

