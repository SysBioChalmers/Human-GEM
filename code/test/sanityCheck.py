# -*- coding: utf-8 -*-

"""Test functions"""

import pandas as pd
import cobra


def get_column_from_tsv(tsv_file, column_id, to_list=True):
    """
    read a column from a tsv file and convert the content into a list by default
    """
    tsv_content = pd.read_table(tsv_file)
    if to_list:
        return tsv_content[column_id].to_list()
    else:
        return tsv_content[column_id]


def load_yml(yml_file):
    """
    import yml model file, and output rxn and met lists
    """
    model = cobra.io.load_yaml_model(yml_file)

    # get lists of rxns, mets and genes
    modelRxns  = list(map(lambda element : element.id, model.reactions))
    modelMets  = list(map(lambda element : element.id, model.metabolites))
    modelGenes = list(map(lambda element : element.id, model.genes))

    return modelRxns, modelMets, modelGenes


def checkRxnAnnotation(rxns):
    """
    check consistency of rxn lists between model and annotation file
    """
    rxnList = get_column_from_tsv("model/reactions.tsv", "rxns")
    spontaneous = get_column_from_tsv("model/reactions.tsv", "spontaneous", False)
    assert rxnList == rxns, "Reaction annotation mismatch!"
    assert pd.api.types.is_numeric_dtype(spontaneous), "Spontaneous column should be in numeric!"


def checkMetAnnotation(mets):
    """
    check consistency of met lists between model and annotation file
    """
    metList = get_column_from_tsv("model/metabolites.tsv", "mets")
    assert metList == mets, "Metabolite annotation mismatch!"


def checkGeneAnnotation(genes):
    """
    check consistency of gene lists between model and annotation file
    """
    geneList = get_column_from_tsv("model/genes.tsv", "genes")
    assert geneList == genes, "Gene annotation mismatch!"


if __name__ == "__main__":
    rxns, mets, genes = load_yml("model/Human-GEM.yml")
    checkRxnAnnotation(rxns)
    checkMetAnnotation(mets)
    checkGeneAnnotation(genes)
    print("Everything passed")

