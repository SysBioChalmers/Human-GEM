"""Fetch Human-GEM reaction names from KEGG
Original file is located at
    https://colab.research.google.com/drive/17X0Qx0H4pwjZjLLWHnpp5ac2daH9hOxs
"""

import requests
import re
import yaml
import pandas

"""Get all the KEGG reactions via their API, and save the result to a file."""

KEGG_REACTIONS = 'kegg_reactions.txt'
HG_YAML = '../model/Human-GEM.yml'
F_YAML = '../model/curated-Human-GEM.yml'

with open(KEGG_REACTIONS,'w') as f:
  r = requests.get('http://rest.kegg.jp/list/reaction/')
  f.write(r.text)

"""Extract the KEGG reactions as key-value pairs."""

raw_reactions = open(KEGG_REACTIONS, 'r')
raw_reaction_lines = raw_reactions.readlines()

reaction_id = re.compile('(?:^rn\:)(R\d+)')
reaction_name = re.compile('(?:\t)([^;]+)(?:;)')
kegg_reactions = {}
for line in raw_reaction_lines:
  try:
    kegg_reactions[reaction_id.search(line).group(1)] = reaction_name.search(line).group(1)
  except:
    kegg_reactions[reaction_id.search(line).group(1)] = ''
# print(kegg_reactions[])

"""Fetch Human-GEM reactions from the TSV annotation."""

hg_annotation = pandas.read_csv('../model/reactions.tsv', sep='\t', index_col=0)

""" Traverse the YAML, and for each line that looks like a reaction definition, extract the reaction identifier, and get the matching KEGG id. Then, change the next line that contains the reaction name to the name provided by KEGG."""

with open(HG_YAML, 'r') as inputf:
  with open(F_YAML, 'w') as outputf:
    count = 0
    count_blank = 0
    while True:
      reaction_id = re.compile('(?:^      - id: ")(MAR\d+)')
      reaction_name = re.compile('(?:^      - name: ")()("$)')
      try:
        line = inputf.readline()
        r_id = reaction_id.search(line).group(1)
        outputf.write(line)
        line = inputf.readline()
        r_name = reaction_name.search(line).group(1)
        kegg_id = hg_annotation.loc[r_id]['rxnKEGGID']
        if kegg_id and r_name == "":
            if kegg_reactions[kegg_id] == "":
                count_blank = count_blank + 1
            else:
                line = '      - name: "' + kegg_reactions[kegg_id] + '"\n'
                count = count + 1
      except:
        None
      outputf.write(line)
      if not line:
        break
    print('Reaction names adopted from KEGG: ' + str(count))
    print('Blank names also blank in KEGG: ' + str(count_blank))
