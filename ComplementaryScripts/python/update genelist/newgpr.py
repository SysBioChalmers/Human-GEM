# UPDATE GPRs of HMR based on new gene list
# 2018 January 4-11
# https://github.com/numpy/numpy

import sys
import re
import csv
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))


'''
Take the old gene list as input: 
Convert the old gene list into new gene list using the curated list in the excel file
Write the new gene list into a csv file.
'''


# TASK 1
# The first file is  old genes versus new genes, i.e., new_old_genes.txt
# Open the first file and convert it to a dictionary

with open(os.path.join(dir_path, r'new_old_genes.txt'), "r") as f:
    d = dict(x.rstrip().split(None, 1) for x in f if x.rstrip())

# print(d)
# print (d['LOC642502'])

# TASK 2
# The second file is the column of old GPRs, i.e., old_gpr.csv
# Convert the CSV file into a MATRIX 
9
with open(os.path.join(dir_path, r"old_gene_list.csv"), "r") as ins:
    old_gene_list = []
    for line in ins:
        ens = line.rstrip().split(',')
        for i in range(len(ens)):
            current_old_ens_id = ens[i]
            if not current_old_ens_id:
                continue
            ens[i] = d[current_old_ens_id]
            # if current_old_ens_id != d[current_old_ens_id]:
                # print ("change %s %s" % (current_old_ens_id, d[current_old_ens_id]))
            
        old_gene_list.append(ens)

# print(old_gene_list)

output_file = os.path.join(dir_path, r'new_gene_list.csv')
with open(output_file, 'w') as fw:
    for el in old_gene_list:
        fw.write(','.join(el)+'\n')

