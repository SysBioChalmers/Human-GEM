# Enzyme regulation data of HMR proteins were retrieved from Uniprot database on 2018-04-23 in the form of excell file.
# The excell file was converted to a txt file, which is C:\Users\kocabas\Desktop\regulation.txt file.
# Columns 1 and 2 of C:\Users\kocabas\Desktop\regulation.txt file have "Uniprot IDs" and "Enzyme regulation data" for human proteins, respectively.
# Therefore, columns 1 and 2 of C:\Users\kocabas\Desktop\regulation.txt file were respectively used as keys and values in the process of generating dictionary d.
# array[0] and array[1] are columns 1 and 2 of C:\Users\kocabas\Desktop\regulation.txt file respectively.

d = {}
f=open(r"C:\Users\kocabas\Desktop\regulation.txt", "r")

for line in f:
	line = line.strip('\n')# remove line break at the end of the line
	array = line.split('\t') # create columns
	uniprot_id = array[0]
	regulation_col = array[1]
	d[uniprot_id]=regulation_col.rstrip()
#print(d['Q9UKM7'])
d['']=''


# Uniprot ids of each HMR reaction were exchanged with enzyme regulation using the dictionary.
# gpruniprotids.csv file is mathematically a matrix; a list of uniprot ids for each HMR reaction
# STEP 2: Change gpr_uniprotids with enzyme regulation and generate gpr_enzymeregulation

with open(r"C:\Users\kocabas\Desktop\gpruniprotids.csv", "r") as ins:
	enzymeregulation = []
	for line in ins:
		uni = line.rstrip().split(',')
		for id in uni:
			id = id.strip()
			if id in d:
			#if d[id]:
				enzymeregulation_str = d[id]
				line = line.replace(id, '[' + enzymeregulation_str + ']')
			#else:
			#	line = line.replace(id, '')
		enzymeregulation.append(line)


#STEP 3: Write the output to a file

output_file = r'C:\Users\kocabas\Desktop\enzymeregulations.txt'
with open(output_file, 'w') as fw:
	for el in enzymeregulation:
		fw.write(el)

