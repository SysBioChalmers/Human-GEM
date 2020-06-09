%
% FILE NAME:    createNewHumanBiomassRxn.m
% 
% PURPOSE: This script adds a new biomass reaction to Human-GEM, which is
%          based on literature, databases, and previous models. The biomass
%          reaction is formulated such that its flux corresponds to 1 gram
%          per gram dry weight of biomass (1 g/gDW).
%


%% Load model and annotation files

% load Human-GEM model
load('HumanGEM.mat');
ihuman_orig = ihuman;  % keep copy of original version

% load metabolite and reaction annotation data
metAssoc = jsondecode(fileread('../../ComplementaryData/annotation/humanGEMMetAssoc.JSON'));
rxnAssoc = jsondecode(fileread('../../ComplementaryData/annotation/humanGEMRxnAssoc.JSON'));
changeNotes = {};

% verify that HumanGEM and annotation structures are aligned
if ~isequal(ihuman.mets, metAssoc.mets) || ~isequal(ihuman.rxns, rxnAssoc.rxns)
    error('HumanGEM is not synced with the metAssoc and/or rxnAssoc structure!');
end


%% Remove temporary HepG2 biomass reactions and associated metabolites
% This reaction was implemented temporarily, and is quite specific to HepG2
% cells. Therefore, it should not remain as a permanent reaction in
% HumanGEM, and will be removed here.

% list of reactions associated with deprecated HepG2 biomass reaction
remRxns = {'biomass_HepG2';'HMR_10016';'HMR_10017';'HMR_10018';'HMR_10019';
           'HMR_10020';'HMR_10021';'HMR_10022'};
if ~all(ismember(remRxns,ihuman.rxns))
    error('One or more reactions to be removed are not present in the model!');
end

% list of metabolites associated with deprecated HepG2 biomass reaction
remMetNames = {'fatty acid biomass pool'; 'heparan sulfate'; 'metabolite pool';
               'phosphatidyl pool'; 'phosphatidate'; 'protein pool'};
if ~all(ismember(remMetNames,ihuman.metNames))
    error('One or more metabolites to be removed are not present in the model!');
end
remMetIDs = ihuman.mets(ismember(ihuman.metNames,remMetNames));

% remove reactions from the model and annotation structure
remRxnInd = getIndexes(ihuman,remRxns,'rxns');
ihuman = removeReactionsFull(ihuman,remRxns);
f = fieldnames(rxnAssoc);
for i = 1:numel(f)
    rxnAssoc.(f{i})(remRxnInd) = [];
end

% remove metabolites from the model and annotation structure
remMetInd = getIndexes(ihuman,remMetIDs,'mets');
ihuman = removeMets(ihuman,remMetInd);
f = fieldnames(metAssoc);
for i = 1:numel(f)
    metAssoc.(f{i})(remMetInd) = [];
end

% note changes
changeNotes = [changeNotes; [remRxns, repmat({'reaction was associated with deprecated HepG2 biomass reaction and should therefore be deleted'},numel(remRxns),1)]];
changeNotes = [changeNotes; [remMetIDs, repmat({'metabolite was associated with deprecated HepG2 biomass reaction and should therefore be deleted'},numel(remMetIDs),1)]];


%% Add new metabolites for new biomass reaction

% add metabolites to the model
metsToAdd = {};
metsToAdd.mets = {'m10012c'; 'm10013c'; 'm10014c'; 'm10015c'};
metsToAdd.metNames = {'cofactor_pool_biomass';'protein_pool_biomass';'lipid_pool_biomass';'metabolite_pool_biomass'};
metsToAdd.compartments = 'c';
ihuman = addMets(ihuman,metsToAdd);

% add metabolites to the annotation structure
numOrigMets = numel(metAssoc.mets);
numNewMets = numel(metsToAdd.mets);
newMetInd = (numOrigMets+1:numOrigMets+numNewMets)';
f = fieldnames(metAssoc);
for i = 1:numel(f)
    metAssoc.(f{i})(newMetInd) = {''};
end
metAssoc.mets(newMetInd) = metsToAdd.mets;
metAssoc.metsNoComp(newMetInd) = regexprep(metsToAdd.mets,'.$','');


%% Add new reactions

% load new reaction information from file
fid = fopen('../../ComplementaryData/modelCuration/newHumanBiomassRxns.tsv');
rxnData = textscan(fid,'%s%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% verify that none of the reactions exist in the current model
if any(ismember(rxnData{1},ihuman.rxns))
    error('One or more reactions to be added already exist in the model.');
end

% add reactions to the model
rxnsToAdd = {};
rxnsToAdd.rxns = rxnData{1};
rxnsToAdd.rxnNames = rxnData{2};
rxnsToAdd.equations = rxnData{3};
rxnsToAdd.subSystems = repmat({{'Artificial reactions'}},numel(rxnsToAdd.rxns),1);
ihuman = addRxns(ihuman, rxnsToAdd, 3);
ihuman.priorCombiningGrRules(end+1:end+numel(rxnsToAdd.rxns)) = {''};

% add reactions to the annotation structure
numOrigRxns = numel(rxnAssoc.rxns);
numNewRxns = numel(rxnData{1});
newRxnInd = (numOrigRxns+1:numOrigRxns+numNewRxns)';
f = fieldnames(rxnAssoc);
for i = 1:numel(f)
    rxnAssoc.(f{i})(newRxnInd) = {''};
end
rxnAssoc.rxns(newRxnInd) = rxnData{1};


%% Finalize and document changes, and export files

% make the new biomass reaction the default objective, and close all other
% biomass reactions (lb = ub = 0)
ihuman = setParam(ihuman,'obj','biomass_human',1);
otherBiomassRxns = startsWith(ihuman.rxns,'biomass_') & ~ismember(ihuman.rxns,'biomass_human');
ihuman.ub(otherBiomassRxns) = 0;

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
modelChanges = docModelChanges(ihuman_orig,ihuman,changeNotes);
writeModelChanges(modelChanges,'../../ComplementaryData/modelCuration/newHumanBiomassRxn_modelChanges.tsv');

% export HumanGEM
exportHumanGEM(ihuman,'HumanGEM','../../',{'mat','yml'},false,false);

% clear unneeded variables
clearvars -except ihuman_orig ihuman modelChanges






