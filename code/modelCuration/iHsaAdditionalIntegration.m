%
% FILE NAME:    iHsaAdditionalIntegration.m
% 
% PURPOSE: This script integrates additional components from the iHsa model
%          (EM Blais, et al. 2017 Nat Commun) into HumanGEM that were not
%          integrated previously. These changes include:
%               - 67 new reactions
%               - 60 new metabolites
%               - 688 new reaction KEGG associations
%               - 1 new gene
%


%% Load model and annotation files

% load Human-GEM model
load('HumanGEM.mat');
ihuman_orig = ihuman;  % keep copy of original version
ihuman = removeConflictingInchiStrings(ihuman);  % remove conflicting inchis

% load metabolite and reaction annotation data
metAssoc = jsondecode(fileread('../../ComplementaryData/annotation/humanGEMMetAssoc.JSON'));
rxnAssoc = jsondecode(fileread('../../ComplementaryData/annotation/humanGEMRxnAssoc.JSON'));

% verify that HumanGEM and annotation structures are aligned
if ~isequal(ihuman.mets, metAssoc.mets) || ~isequal(ihuman.rxns, rxnAssoc.rxns)
    error('HumanGEM is not synced with the metAssoc and/or rxnAssoc structure!');
end

% remove some deprecated fields if they still exist
removeFields = intersect(fieldnames(ihuman),{'proteins';'prRules';'rxnProtMat';'priorCombiningGrRules'});
if ~isempty(removeFields)
    ihuman = rmfield(ihuman, removeFields);
end


%% Add new rxnRatconID field to reaction annotation file

% load iHsa model
load('../../ComplementaryData/iHsa/iHsa.mat');  % loads "iHsa" structure

% add reaction Ratcon IDs to rxn association structure
[inModel,rxnInd] = ismember(ihuman.rxns, iHsa.rxnHMRID);
ids = repmat({''}, numel(ihuman.rxns), 1);
ids(inModel) = iHsa.rxns(rxnInd(inModel));
rxnAssoc.rxnRatconID = ids;


%% Add new metabolites from iHsa

% load new metabolite information from file
fid = fopen('../../ComplementaryData/iHsa/iHsaMetsToAdd.tsv');
metData = textscan(fid,'%s%s%s%d%s%s%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% verify that none of the metabolites exist in the current model
if any(ismember(regexprep(metData{1},'.$',''),ihuman.mets))
    error('One or more metabolite IDs to be added already exist in the model.');
end

% add metabolites to the model
metsToAdd = {};
metsToAdd.mets = metData{1};
metsToAdd.metNames = metData{2};
metsToAdd.metFormulas = metData{3};
metsToAdd.metCharges = metData{4};
metsToAdd.compartments = metData{5};
ihuman = addMets(ihuman,metsToAdd);

% add metabolites to the annotation structure
numOrigMets = numel(metAssoc.mets);
numNewMets = numel(metsToAdd.mets);
newMetInd = (numOrigMets+1:numOrigMets+numNewMets)';
f = fieldnames(metAssoc);
for i = 1:numel(f)
    metAssoc.(f{i})(newMetInd) = {''};  % initialize fields with empty entry
end
metAssoc.mets(newMetInd) = metsToAdd.mets;
metAssoc.metsNoComp(newMetInd) = regexprep(metsToAdd.mets,'.$','');
metAssoc.metPubChemID(newMetInd) = metData{6};
metAssoc.metKEGGID(newMetInd) = metData{7};
metAssoc.metChEBIID(newMetInd) = metData{8};


%% Add new genes from iHsa

newGenes = {'ENSG00000183463'};
if any(ismember(newGenes,ihuman.genes))
    error('One or more genes to be added already exist in the model.');
end

% append new genes to list of model genes
ihuman.genes = [ihuman.genes; newGenes];

% add new columns to rxnGeneMat will be updated after the new reactions are added below.
ihuman.rxnGeneMat(:, end+1:end+numel(newGenes)) = 0;


%% Add new reactions from iHsa

% load new reaction information from file
fid = fopen('../../ComplementaryData/iHsa/iHsaRxnsToAdd.tsv');
rxnData = textscan(fid,'%s%s%s%s%s%s%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% verify that none of the reactions exist in the current model
if any(ismember(rxnData{1},ihuman.rxns))
    error('One or more reactions to be added already exist in the model.');
end

% add reactions to the model
rxnsToAdd = {};
rxnsToAdd.rxns = rxnData{1};
rxnsToAdd.equations = rxnData{2};
rxnsToAdd.eccodes = rxnData{3};
rxnsToAdd.subSystems = cellfun(@(s) {{s}}, rxnData{4});
rxnsToAdd.grRules = rxnData{5};
rxnsToAdd.rxnReferences = rxnData{6};
ihuman = addRxns(ihuman, rxnsToAdd, 3);

% ensure that genes and rxnGeneMat are synced with new grRules
[ihuman.genes, ihuman.rxnGeneMat] = getGenesFromGrRules(ihuman.grRules);

% add reactions to the annotation structure
numOrigRxns = numel(rxnAssoc.rxns);
numNewRxns = numel(rxnData{1});
newRxnInd = (numOrigRxns+1:numOrigRxns+numNewRxns)';
f = fieldnames(rxnAssoc);
for i = 1:numel(f)
    rxnAssoc.(f{i})(newRxnInd) = {''};  % initialize fields with empty entry
end
rxnAssoc.rxns(newRxnInd) = rxnData{1};
rxnAssoc.rxnKEGGID(newRxnInd) = rxnData{7};
rxnAssoc.rxnRatconID(newRxnInd) = rxnData{8};


%% Add new reaction KEGG IDs from iHsa

% import data from file
fid = fopen('../../ComplementaryData/iHsa/iHsaRxnAssocToAdd.tsv');
rxnData = textscan(fid,'%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% update reaction KEGG IDs with iHsa information
[~,rxnInd] = ismember(rxnData{1},ihuman.rxns);
rxnAssoc.rxnKEGGID(rxnInd) = rxnData{2};


%% Finalize and document changes, and export files

% update unconstrained field
[~,boundaryComp] = ismember('Boundary',ihuman.compNames);
ihuman.unconstrained(ihuman.metComps == boundaryComp) = 1;

% verify that HumanGEM and annotation structures are aligned
if ~isequal(ihuman.mets, metAssoc.mets) || ~isequal(ihuman.rxns, rxnAssoc.rxns)
    error('HumanGEM is not synced with the metAssoc and/or rxnAssoc structure!');
end

% export reaction and metabolite annotation structures to JSON
jsonStr = jsonencode(rxnAssoc);
fid = fopen('../../ComplementaryData/annotation/humanGEMRxnAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);
jsonStr = jsonencode(metAssoc);
fid = fopen('../../ComplementaryData/annotation/humanGEMMetAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);

% identify and document model changes
modelChanges = docModelChanges(ihuman_orig,ihuman);
writeModelChanges(modelChanges,'../../ComplementaryData/iHsa/iHsaAdditionalIntegration_modelChanges.tsv');

% export HumanGEM
exportHumanGEM(ihuman,'HumanGEM','../../',{'mat','yml'},false,false);

% clear unneeded variables
clearvars -except ihuman_orig ihuman modelChanges




