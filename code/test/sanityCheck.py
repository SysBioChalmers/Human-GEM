# -*- coding: utf-8 -*-

"""Test functions"""

import pandas as pd
import cobra
from collections import Counter


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

    return modelRxns, modelMets, modelGenes, model


def checkRxnAnnotation(rxns):
    """
    check consistency of rxn lists between model and annotation file
    """
    rxnList = get_column_from_tsv("model/reactions.tsv", "rxns")
    spontaneous = get_column_from_tsv("model/reactions.tsv", "spontaneous", False)
    rxnDeprecated = get_column_from_tsv("data/deprecatedIdentifiers/deprecatedReactions.tsv", "rxns")
    rxnMisused = set(rxns).intersection(set(rxnDeprecated))
    assert len(rxnMisused) == 0, "Deprecated reaction(s) used: {}".format(rxnMisused)
    rxnMisused = set(rxns).difference(set(rxnList))
    assert len(rxnMisused) == 0, "Reaction(s) not annotated in reactions.tsv: {}".format(rxnMisused)
    rxnMisused = set(rxnList).difference(set(rxns))
    assert len(rxnMisused) == 0, "Annotated reactions are missing from the model: {}".format(rxnMisused)  
    assert pd.api.types.is_numeric_dtype(spontaneous), "Spontaneous column should be in numeric!"


def checkMetAnnotation(mets):
    """
    check consistency of met lists between model and annotation file
    """
    metList = get_column_from_tsv("model/metabolites.tsv", "mets")
    metDeprecated = get_column_from_tsv("data/deprecatedIdentifiers/deprecatedMetabolites.tsv", "mets")
    metMisused = set(mets).intersection(set(metDeprecated))
    assert len(metMisused) == 0, "Deprecated metabolite(s) used: {}".format(metMisused)
    metMisused = set(mets).difference(set(metList))
    assert len(metMisused) == 0, "Metabolite(s) not annotated in metabolites.tsv: {}".format(metMisused)
    metMisused = set(metList).difference(set(mets))
    assert len(metMisused) == 0, "Annotated metabolite(s) are missing from the model: {}".format(metMisused)  

def checkGeneAnnotation(genes):
    """
    check consistency of gene lists between model and annotation file
    """
    geneList = get_column_from_tsv("model/genes.tsv", "genes")
    geneMisused = set(genes).difference(set(geneList))
    assert len(geneMisused) == 0, "Gene(s) not annotated in genes.tsv: {}".format(geneMisused)
    geneMisused = set(geneList).difference(set(genes))
    assert len(geneMisused) == 0, "Annotated gene(s) are missing from the model: {}".format(geneMisused)  


def find_unused_entities(model, entity_type):
    """
    collect unused metabolites or genes if exist in the model
    """

    # collect all metabolites actually used by reactions
    entities_used = []

    for reaction in model.reactions:
        entities = reaction.genes if entity_type == "genes" else reaction.metabolites
        # Loop through each entity and add if not already there
        for entity in entities:
            if entity not in entities_used:
                entities_used.append(entity)

    # go through entities in the model and collect unused ones
    unused_entities = []
    all_entities = model.genes if entity_type == "genes" else model.metabolites
    for entity in all_entities:
        if entity not in entities_used:
            unused_entities.append(entity)

    return unused_entities


def checkUnusedEntities(model, entity_type):
    """
    check if unused genes or metabolites exist in the model
    """

    # collect unused entites
    unused_entities = find_unused_entities(model, entity_type)
    assert len(unused_entities) == 0, f"Found unused {entity_type}: {unused_entities}"


def checkDupRxn(model):
    """
    Check for duplicate reactions in the model
    """

    reaction_equations = [rxn.build_reaction_string(use_metabolite_names=False) for rxn in model.reactions]
    duplicate_reactions = [reaction for reaction, count in Counter(reaction_equations).items() if count > 1]
    dup_list = [model.reactions[idx].id for idx, val in enumerate(reaction_equations) if val in duplicate_reactions]

    if duplicate_reactions:
        output = f"The following {len(dup_list)} reactions are duplicates, please check: " + ';'.join(dup_list)
        print(output)
        
    assert len(duplicate_reactions) == 0, "Found duplicated reactions!"


if __name__ == "__main__":
    rxns, mets, genes, model = load_yml("model/Human-GEM.yml")
    checkRxnAnnotation(rxns)
    checkMetAnnotation(mets)
    checkGeneAnnotation(genes)
    checkUnusedEntities(model, "metabolites")
    checkUnusedEntities(model, "genes")
    checkDupRxn(model)
    print("All checks have passed.")
