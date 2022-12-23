# collect subcellular localization for existing enzymes from Uniprot/Swissprot and Cell Atlas
# compartment info from both sources were combined and integrated into genes.tsv

import pandas as pd
import csv

# import subcellular location from SwissProt annotation
SwissProt_tsv = pd.read_table("~/Downloads/SwissProt_20221115.tsv")
SwissProt_proteins = SwissProt_tsv['Entry'].to_list()
SwissProt_subCellularlocation = SwissProt_tsv['Subcellular location [CC]'].to_list()


# build a dict for key word mapping of compartment assignment
swissprot_keywords = {
    #'Extracellular': ['',''],   # e
    'Peroxisome': ['peroxisome'],   # x
    'Mitochondria': ['mitochondrion {','mitochondrion matrix'],   # m
    'Cytosol': ['cytoplasm','cytosol','cytoskeleton'],   # c
    'Lysosome': ['lysosome','lysosome lumen','lysosome membrane'],   # l
    'Endoplasmic reticulum': ['endoplasmic reticulum','Endoplasmic reticulum lumen','endoplasmic reticulum membrane'],   # r
    'Golgi apparatus': ['golgi apparatus','golgi apparatus lumen','golgi apparatus membrane'],   # g
    'Nucleus': ['nucleus','nucleolus'],   # n
    'Inner mitochondria': ['mitochondrion inner membrane','mitochondrion intermembrane space']    # i
}


# store swissport compartment info as a dict with a key of uniprot id
comps_from_swissprot = {}
for i, text in enumerate(SwissProt_subCellularlocation):
    out_list = []  # output compartmetns as a list
    comps_from_swissprot[SwissProt_proteins[i]] = out_list
    if not pd.isna(text):
        for key in swissprot_keywords:
            # check if cellular location text contains any keyword from a compartment
            if any(word in text.lower() for word in swissprot_keywords[key]):
                out_list.append(key)
        comps_from_swissprot[SwissProt_proteins[i]] = out_list


# load Human-GEM gene annotation file
HumanGenes_tsv = pd.read_table("../../model/genes.tsv")
Human_genes = HumanGenes_tsv['genes'].to_list()
Human_proteins = HumanGenes_tsv['geneUniProtID'].to_list()


# save SwissProt compartment information to tsv file
swissprot_compList = ['']*len(Human_genes)
for i in range(len(Human_genes)):
    if not pd.isna(Human_proteins[i]):
        if Human_proteins[i] in SwissProt_proteins:
            swissprot_compList[i] = ';'.join(comps_from_swissprot[Human_proteins[i]])
swissprot_compartments = pd.DataFrame()
swissprot_compartments['genes'] = Human_genes
swissprot_compartments['geneUniProtID'] = Human_proteins
swissprot_compartments['compartments'] = swissprot_compList
swissprot_compartments.to_csv('../../data/modelCuration/Swissprot_compartments.tsv', sep="\t", index=False)


# import Cell Atlas compartment info
cell_atlas_compartments = pd.read_table("../../data/modelCuration/CellAtlasCompartments_science_2017.tsv")
cellAtlas_comps = cell_atlas_compartments['Subcellular location'].to_list()
cellAtlas_ensembl_id = cell_atlas_compartments['Ensembl'].to_list()
cellAtlas_uniprot_id = cell_atlas_compartments['Uniprot'].to_list()

# investigate key words
compartment_names = []
for comps in cellAtlas_comps:
    compartment_names.extend(comps.split(','))
unique_compartments = set(compartment_names)
unique_compartments
# Cell Atlas compartments were sorted and the following ones were discarded:
# 'Microtubules', 'Midbody ring', 'Midbody', 'Cytokinetic bridge', 'Microtubule ends', 'Mitotic spindle'
# 'Intermediate filaments', 'Actin filaments', 'Focal adhesion sites', 'Cleavage furrow', 
# 'Lipid droplets', 'Vesicles', 'Endosomes', 'Plasma membrane', 'Cell Junctions'
# 'Centrosome', 'Centriolar satellite'

# construct a dict for key word mapping to Cell Atlas compartments
cellatlas_keywords = {
    #'Extracellular': [''],   # e
    'Peroxisome': ['Peroxisomes'],   # x
    'Mitochondria': ['Mitochondria'],   # m
    'Cytosol': ['Cytosol', 'Cytoplasmic bodies', 'Rods & Rings', 'Aggresome'],   # c
    'Lysosome': ['Lysosomes'],   # l
    'Endoplasmic reticulum': ['Endoplasmic reticulum'],   # r
    'Golgi apparatus': ['Golgi apparatus'],   # g
    'Nucleus': ['Nuclear membrane', 'Nucleoli rim', 'Nucleoli', 'Nuclear bodies', 'Nucleoli fibrillar center', 'Nucleoplasm', 'Kinetochore', 'Mitotic chromosome', 'Nuclear speckles'],   # n
    #'Inner mitochondria': ['']    # i, Cell Atlas does provide such location info
}


# store CellAtlas compartment info to a dict with keys as ensembl ids
geneComps_from_cell_atlas = {}
for gene in Human_genes:
    out_list = []  # store compartments as a list
    if gene not in cellAtlas_ensembl_id:
        geneComps_from_cell_atlas[gene] = out_list
    else:
        gene_ind = cellAtlas_ensembl_id.index(gene)
        for key in cellatlas_keywords:
            # check if cellular location text contains any keyword from a compartment
            if any(word in cellAtlas_comps[gene_ind] for word in cellatlas_keywords[key]):
                out_list.append(key)
        geneComps_from_cell_atlas[gene] = out_list


# integrate compartment info from cellAtlas and SwissProt with following rules:
# 1. output two columns: "compartments" and "compDataSource"
# 2. use one source if another has no assigned compartments
# 3. union compartments if both source are provided
swissprot_comps = pd.read_table("../../data/modelCuration/Swissprot_compartments.tsv")
geneComps_from_swissprot = swissprot_comps['compartments'].to_list()

# combine compartment info from two data sources
geneComps_combined = []
source = ['']*len(Human_genes)
for i, gene in enumerate(Human_genes):
    out_list = []  # save compartmetns as a list
    if pd.isna(geneComps_from_swissprot[i]) and geneComps_from_cell_atlas[gene] != []:
        out_list = geneComps_from_cell_atlas[gene]
        source[i] = 'CellAtlas'
    elif not pd.isna(geneComps_from_swissprot[i]) and geneComps_from_cell_atlas[gene] == []:
        out_list = geneComps_from_swissprot[i].split(';')
        source[i] = 'SwissProt'
    elif not pd.isna(geneComps_from_swissprot[i]) and geneComps_from_cell_atlas[gene] != []:
        union = set(geneComps_from_swissprot[i].split(';')+ geneComps_from_cell_atlas[gene])
        out_list = list(union)
        source[i] = 'SwissProt;CellAtlas'
        if 'Mitochondria' not in geneComps_from_swissprot[i].split(';') and 'Inner mitochondria' in geneComps_from_swissprot[i].split(';'):
            if 'Mitochondria' in geneComps_from_cell_atlas[gene]:
                out_list.remove('Mitochondria')
        # manual inspection of integration
        print(gene+': '+geneComps_from_swissprot[i]+'\t'+';'.join(geneComps_from_cell_atlas[gene])+'\t'+';'.join(out_list))
    geneComps_combined.append(';'.join(out_list))


# append columns to genes.tsv
HumanGenes_tsv['compartments'] = geneComps_combined
HumanGenes_tsv['compDataSource'] = source
HumanGenes_tsv[:0].to_csv('../../model/genes.tsv', sep="\t", index=False)
HumanGenes_tsv.to_csv('../../model/genes.tsv', sep="\t", index=False, quoting=csv.QUOTE_NONNUMERIC, header=False, mode="a")

