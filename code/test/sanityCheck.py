# -*- coding: utf-8 -*-

"""Test functions"""

import pandas as pd
import cobra


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
    rxnAssoc = pd.read_table("model/reactions.tsv")
    rxnList = rxnAssoc['rxns'].to_list()
    assert rxnList == rxns, "Reaction annotation mismatch!"
    assert pd.api.types.is_numeric_dtype(
        rxnAssoc['spontaneous']), "Spontaneous column should be in numeric!"


def checkMetAnnotation(mets):
    """
    check consistency of met lists between model and annotation file
    """
    metAssoc = pd.read_table("model/metabolites.tsv")
    metList = metAssoc['mets'].to_list()
    assert metList == mets, "Metabolite annotation mismatch!"


def checkGeneAnnotation(genes):
    """
    check consistency of gene lists between model and annotation file
    """
    geneAssoc = pd.read_table("model/genes.tsv")
    geneList = geneAssoc['genes'].to_list()
    assert geneList == genes, "Gene annotation mismatch!"


if __name__ == "__main__":
    rxns, mets, genes = load_yml("model/Human-GEM.yml")
    checkRxnAnnotation(rxns)
    checkMetAnnotation(mets)
    checkGeneAnnotation(genes)
    print("Everything passed")

