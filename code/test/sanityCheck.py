# -*- coding: utf-8 -*-

"""Test functions"""

import pandas as pd
import cobra


def get_column_from_tsv(tsv_file, column_id, to_list=True):
    """
    read a column from a tsv file and convert the content into a list by default
    """
    tsv_content = pd.read_table(tsv_file)
    desired_column = tsv_content[column_id]
    if to_list:
        return desired_column.to_list()
    return desired_column


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
    rxnDeprecated = get_column_from_tsv("data/deprecatedIdentifiers/deprecatedReactions.tsv", "rxns")
    assert not bool(set(rxnList) & set(rxnDeprecated)), "Deprecated reaction reused!"
    assert rxnList == rxns, "Reaction annotation mismatch!"
    assert pd.api.types.is_numeric_dtype(spontaneous), "Spontaneous column should be in numeric!"


def checkMetAnnotation(mets):
    """
    check consistency of met lists between model and annotation file
    """
    metList = get_column_from_tsv("model/metabolites.tsv", "mets")
    metDeprecated = get_column_from_tsv("data/deprecatedIdentifiers/deprecatedMetabolites.tsv", "mets")
    assert not bool(set(metList) & set(metDeprecated)), "Deprecated metabolite reused!"
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
    print("All checks have passed.")
