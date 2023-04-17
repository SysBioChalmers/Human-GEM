% This function is for adding new rxns annotated for glycolysis genes of GPRs identified by GPT.
rxnStruct = importTsvFile('../../data/modelCuration/addRxnGly_20230414.tsv');
MetStruct = importTsvFile('../../data/modelCuration/addMetGly_20230414.tsv');

model = importYaml('Human-GEM.yml');

rxnsToAdd.rxns = rxnStruct.rxnID;
rxnsToAdd.eccodes = rxnStruct.ECNumber;
rxnsToAdd.equations = rxnStruct.Equation;
rxnsToAdd.rxnNames = rxnStruct.rxnID;
rxnsToAdd.lb = zeros(length(rxnStruct.rxnID),1);
rxnsToAdd.ub = repmat(1000,length(rxnStruct.rxnID),1);
idx = find(cell2mat(rxnStruct.Reversibility) == 1);
rxnsToAdd.lb(idx) = -1000;
rxnsToAdd.grRules = rxnStruct.grRules;
rxnsToAdd.subSystems = rxnStruct.Subsystems;
rxnsToAdd.rxnConfidenceScores = zeros(length(rxnStruct.rxnID),1);
[~,idx] = ismember(MetStruct.MetID,model.mets);
metsToAdd.mets = MetStruct.MetID(~idx,1);
metsToAdd.metNames = MetStruct.Name(~idx,1);
metsToAdd.compartments = MetStruct.Comp(~idx,1);
metsToAdd.metFormulas = MetStruct.Fromula(~idx,1);
metsToAdd.metCharges = MetStruct.Charge(~idx,1);
tmp = [];
for i = 1:length(metsToAdd.metCharges)
    if ~isempty(metsToAdd.metCharges{i})
        tmp(i,1) = str2double(metsToAdd.metCharges{i});
    else
        tmp(i,1) = NaN;
    end
end
metsToAdd.metCharges = tmp;

[newModel, modelChanges] = addMetabolicNetwork(model,rxnsToAdd,metsToAdd);

% export the new mets, new reactions to the table
structure = importTsvFile('../../model/metabolites.tsv');
[~,idx] = ismember(modelChanges.mets.mets,MetStruct.MetID);
MetStruct.emptyID(1:length(MetStruct.MetID),1) = {''};
% metID metNames	metFormulas	metCharges	compartments	metKEGGID	metPubChemID	metChEBIID	metMetaCycID	metMetaNetXID
metNocomp = cellfun(@(s) s(1:8), MetStruct.MetID(idx,1), 'UniformOutput', false);
structure.mets = [structure.mets;MetStruct.MetID(idx,1)];
structure.metBiGGID = [structure.metBiGGID;MetStruct.emptyID(idx,1)];
structure.metKEGGID = [structure.metKEGGID;MetStruct.KEGG(idx,1)];
structure.metsNoComp = [structure.metsNoComp;metNocomp];
structure.metHMDBID = [structure.metHMDBID;MetStruct.emptyID(idx,1)];
structure.metChEBIID = [structure.metChEBIID;MetStruct.CHEBI(idx,1)];
structure.metPubChemID = [structure.metPubChemID;MetStruct.emptyID(idx,1)];
structure.metLipidMapsID = [structure.metLipidMapsID;MetStruct.emptyID(idx,1)];
structure.metEHMNID = [structure.metEHMNID;MetStruct.emptyID(idx,1)];
structure.metHepatoNET1ID = [structure.metHepatoNET1ID;MetStruct.emptyID(idx,1)];
structure.metRecon3DID = [structure.metRecon3DID;MetStruct.emptyID(idx,1)];
structure.metMetaNetXID = [structure.metMetaNetXID;MetStruct.MetaNetx(idx,1)];
structure.metHMR2ID = [structure.metHMR2ID;MetStruct.emptyID(idx,1)];
structure.metRetired = [structure.metRetired;MetStruct.emptyID(idx,1)];
% sort the order based on the model.mets
[~,idx] = ismember(model.mets,structure.mets);
structure.mets = structure.mets(idx);
structure.metsNoComp = structure.metsNoComp(idx);
structure.metBiGGID = structure.metBiGGID(idx,1);
structure.metKEGGID = structure.metKEGGID(idx,1);
structure.metHMDBID = structure.metHMDBID(idx,1);
structure.metChEBIID = structure.metChEBIID(idx,1);
structure.metPubChemID = structure.metPubChemID(idx,1);
structure.metLipidMapsID = structure.metLipidMapsID(idx,1);
structure.metEHMNID = structure.metEHMNID(idx,1);
structure.metHepatoNET1ID = structure.metHepatoNET1ID(idx,1);
structure.metRecon3DID = structure.metRecon3DID(idx,1);
structure.metMetaNetXID = structure.metMetaNetXID(idx,1);
structure.metHMR2ID = structure.metHMR2ID(idx,1);
structure.metRetired = structure.metRetired(idx,1);
exportTsvFile(structure, '../../model/metabolites.tsv')

structure = importTsvFile('../../model/reactions.tsv');
[~,idx] = ismember(modelChanges.rxns.rxns,rxnStruct.rxnID);
rxnStruct.emptyID(1:length(rxnStruct.rxnID),1) = {''};
structure.rxns = [structure.rxns;rxnStruct.rxnID(idx,1) ];
structure.rxnKEGGID = [structure.rxnKEGGID;rxnStruct.KEGG(idx,1) ];
structure.rxnBiGGID= [structure.rxnBiGGID;rxnStruct.BiGG(idx,1) ];
structure.rxnEHMNID = [structure.rxnEHMNID;rxnStruct.emptyID(idx,1) ];
structure.rxnHepatoNET1ID = [structure.rxnHepatoNET1ID;rxnStruct.emptyID(idx,1) ];
structure.rxnREACTOMEID = [structure.rxnREACTOMEID;rxnStruct.Reactome(idx,1) ];
structure.rxnRecon3DID = [structure.rxnRecon3DID;rxnStruct.emptyID(idx,1) ];
structure.rxnMetaNetXID = [structure.rxnMetaNetXID;rxnStruct.MetaNetx(idx,1) ];
structure.rxnHMR2ID = [structure.rxnHMR2ID;rxnStruct.emptyID(idx,1) ];
structure.rxnRatconID = [structure.rxnRatconID;rxnStruct.emptyID(idx,1) ];
structure.rxnTCDBID = [structure.rxnTCDBID;rxnStruct.emptyID(idx,1) ];
structure.spontaneous = [structure.spontaneous;zeros(length(rxnStruct.rxnID(idx,1)),1) ];
structure.rxnRheaID = [structure.rxnRheaID;rxnStruct.RHEA(idx,1) ];
structure.rxnRheaMasterID = [structure.rxnRheaMasterID;rxnStruct.emptyID(idx,1) ];
structure.rxnRetired = [structure.rxnRetired;rxnStruct.emptyID(idx,1) ];

exportTsvFile(structure, '../../model/reactions.tsv')

