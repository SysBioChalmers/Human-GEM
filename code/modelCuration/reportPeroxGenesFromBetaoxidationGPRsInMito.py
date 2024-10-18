# 1. read gene compartment from genes.tsv and load model
# 2. get the reactions in fatty acid oxidation subsystem and mitochondria compartment
# 3. go through the reactions and check the compartment of each genes in the GPR
# 4. report the genes that are only expressed in peroxisome
import csv
import cobra
import pandas as pd

# load genes.tsv and read info into a dict for genes and their compartments
with open("../../model/genes.tsv", 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    geneCompDict = {row[0]: row[8] for row in reader}

# load model
model = cobra.io.load_yaml_model('../../model/Human-GEM.yml')

# collect reactions from "Fatty acid oxidation" subsystem and in [m] compartment
subsys = 'Fatty acid oxidation'
targetComp = {'m'}
for r in model.reactions:
    if subsys in r.subsystem and targetComp == r.compartments:
        # go through collected reactions find genes that exist only in peroxisome
        for g in r.genes:
            if geneCompDict[g.id] == 'Peroxisome':
                print(f'{r.id} | {g.id} | {r.build_reaction_string(True)} | {r.gene_reaction_rule} | {r.gene_name_reaction_rule}')

