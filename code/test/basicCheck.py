# -*- coding: utf-8 -*-

"""Test functions"""

import json
import cobra


def load_yml(yml_file):
    """
    import yml model file, and output rxn and met lists
    """
    model = cobra.io.load_yaml_model(yml_file)

    # get list of rxns
    modelRxns = []
    for r in model.reactions:
        modelRxns.append(r.id)

    # get list of mets
    modelMets = []
    for m in model.metabolites:
        modelMets.append(m.id)

    return modelRxns, modelMets


def checkRxnAnnotation(rxns):
    """
    check consistency of rxn lists between model and annoation file
    """
    rxnJsonFile = "data/annotation/humanGEMRxnAssoc.JSON"
    with open(rxnJsonFile) as rxnJson:
        rxnAssoc = json.load(rxnJson)
    rxnList = rxnAssoc.get('rxns', None)
    assert rxnList == rxns, "Reaction annoation mismatch!"


def checkMetAnnotation(mets):
    """
    check consistency of met lists between model and annoation file
    """
    metJsonFile = "data/annotation/humanGEMMetAssoc.JSON"
    with open(metJsonFile) as metJson:
        metAssoc = json.load(metJson)
    metList = metAssoc.get('mets', None)
    assert metList == mets, "Metabolite annoation mismatch!"


if __name__ == "__main__":
    rxns, mets = load_yml("model/Human-GEM.yml")
    checkRxnAnnotation(rxns)
    checkMetAnnotation(mets)
    print("Everything passed")

