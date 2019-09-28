%
% FILE NAME:    curateLipidPools.m
% 
% DATE CREATED: 2019-09-07
%     MODIFIED: 2019-09-27
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: This script revises the formulation of reactions involving lipid
%          and cholesterol-ester pool formation and breakdown in HumanGEM
%          to achieve a more generic and mass-balanced model. In the
%          process, the following changes were made:
%
%             > 6 new pool reactions
%             > 2 new transport reaction
%             > 16 new pool metabolites (3 unique, ignoring compartment)
%             > 23 pool metabolites were updated to the new metabolites
%
%          Full details on model changes can be found in the generated
%          model change log files (lipidPool_modelChanges_mets.tsv and 
%          lipidPool_modelChanges_rxns.tsv)
%


%% Load Human-GEM and annotation files

% load Human-GEM model
load('HumanGEM.mat');
ihuman_orig = ihuman;  % keep copy of original version

% load metabolite and reaction annotation data
metAssoc = jsondecode(fileread('humanGEMMetAssoc.JSON'));
rxnAssoc = jsondecode(fileread('humanGEMRxnAssoc.JSON'));

% verify that HumanGEM and annotation structures are aligned
if ~isequal(ihuman.mets, metAssoc.mets) || ~isequal(ihuman.rxns, rxnAssoc.rxns)
    error('HumanGEM is not synced with the metAssoc and/or rxnAssoc structure!');
end


%% Add new metabolites

% first need to rename an existing "fatty acid pool" metabolite to avoid confusion
ihuman.metNames(ismember(ihuman.metNames,'fatty acid pool')) = {'fatty acid biomass pool'};

% load new metabolite information from file
fid = fopen('../../ComplementaryData/modelCuration/lipidPools/newLipidPoolMets.tsv');
metData = textscan(fid,'%s%s%s%s%d%s%s%s%s%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% verify that none of the metabolite IDs (ignoring compartment) exist in the current model
if any(startsWith(ihuman.mets, regexprep(metData{1},'.$','')))
    error('One or more metabolite IDs to be added already exist in the model.');
end

% add metabolites to the model
metsToAdd = {};
metsToAdd.mets = metData{1};
metsToAdd.metNames = metData{2};
metsToAdd.compartments = metData{3};
metsToAdd.metFormulas = metData{4};
metsToAdd.metCharges = metData{5};
ihuman = addMets(ihuman, metsToAdd);

% add metabolites to the metAssoc structure
numOrigMets = numel(metAssoc.mets);
numNewMets = numel(metData{1});
newMetInd = (numOrigMets+1:numOrigMets+numNewMets)';
f = fieldnames(metAssoc);
for i = 1:numel(f)
    % initialize structure with empty entries
    metAssoc.(f{i})(newMetInd) = {''};
end
metAssoc.mets(newMetInd) = metData{1};
metAssoc.metsNoComp(newMetInd) = regexprep(metData{1},'.$','');
metAssoc.metKEGGID(newMetInd) = metData{6};
metAssoc.metBiGGID(newMetInd) = metData{7};
metAssoc.metChEBIID(newMetInd) = metData{8};
metAssoc.metMetaNetXID(newMetInd) = metData{9};
metAssoc.metLipidMapsID(newMetInd) = metData{10};



%% Add new pool reactions

% load new reaction information from file
fid = fopen('../../ComplementaryData/modelCuration/lipidPools/newLipidPoolRxns.tsv');
rxnData = textscan(fid,'%s%s%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% verify that none of the reactions exist in the current model
if any(ismember(rxnData{1},ihuman.rxns))
    error('One or more reactions to be added already exist in the model.');
end

% add reactions to the model
rxnsToAdd = {};
rxnsToAdd.rxns = rxnData{1};
rxnsToAdd.rxnNames = rxnData{2};
rxnsToAdd.subSystems = cellfun(@(s) {{s}},rxnData{3});
rxnsToAdd.equations = rxnData{4};
ihuman = addRxns(ihuman, rxnsToAdd, 3);

% add reactions to the rxnAssoc structure
numOrigRxns = numel(rxnAssoc.rxns);
numNewRxns = numel(rxnData{1});
newRxnInd = (numOrigRxns+1:numOrigRxns+numNewRxns)';
f = fieldnames(rxnAssoc);
for i = 1:numel(f)
    % initialize structure with empty entries
    rxnAssoc.(f{i})(newRxnInd) = {''};
end
rxnAssoc.rxns(newRxnInd) = rxnData{1};


%% Update cholesterol-ester pool metabolites
% The "cholesterol-ester pool" metabolites will be updated to the new 
% "cholesterol-ester plasma pool" metabolites ONLY in the old cholesterol-
% ester pool formation/degradation reactions (HMR_3537 and HMR_3622). All
% other instances of the "cholesterol-ester pool" metabolite will remain
% unchanged.

% HMR_3537 (lysosome)
oldMetInd = getIndexes(ihuman, 'cholesterol-ester pool[l]', 'metscomps');
newMetInd = getIndexes(ihuman, 'cholesterol-ester plasma pool[l]', 'metscomps');
rxnInd = getIndexes(ihuman,'HMR_3537','rxns');
ihuman.S(oldMetInd, rxnInd) = 0;
ihuman.S(newMetInd, rxnInd) = -1;

% HMR_3622 (endoplasmic reticulum)
oldMetInd = getIndexes(ihuman, 'cholesterol-ester pool[r]', 'metscomps');
newMetInd = getIndexes(ihuman, 'cholesterol-ester plasma pool[r]', 'metscomps');
rxnInd = getIndexes(ihuman,'HMR_3622','rxns');
ihuman.S(oldMetInd, rxnInd) = 0;
ihuman.S(newMetInd, rxnInd) = -1;


%% Update fatty acid, 1-acylglycerol-3P, and acyl-CoA pool metabolites
% The old version of these pool metabolites will be updated to their new
% versions in ALL reactions EXCEPT those that involve their formation/
% degradation. The formation/degradation reactions involving the old pool
% metabolites will remain unchanged.

% load list of existing pool reactions related to lipid pools
fid = fopen('../../ComplementaryData/modelCuration/lipidPools/oldLipidPoolRxns.tsv');
filedata = textscan(fid,'%s','Headerlines',1);
fclose(fid);
oldPoolRxns = filedata{1};
oldPoolRxnInd = getIndexes(ihuman, oldPoolRxns, 'rxns');

% load array specifying the conversion from old lipid pool mets to new mets
fid = fopen('../../ComplementaryData/modelCuration/lipidPools/lipidPoolMetReplacements.tsv');
filedata = textscan(fid,'%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);
poolMets = [filedata{1}, filedata{2}];

% replace old pool metabolites with new pool metabolites
for i = 1:size(poolMets,1)
    metComps = unique(ihuman.metComps(ismember(ihuman.metNames, poolMets(i,1))));
    for j = 1:numel(metComps)
        oldMetInd = find(ismember(ihuman.metNames,poolMets(i,1)) & ihuman.metComps == metComps(j));
        newMetInd = find(ismember(ihuman.metNames,poolMets(i,2)) & ihuman.metComps == metComps(j));
        rxnInd = setdiff(find(ihuman.S(oldMetInd,:) ~= 0), oldPoolRxnInd);
        if ~isempty(rxnInd)
            stoichCoeffs = ihuman.S(oldMetInd, rxnInd);
            ihuman.S(oldMetInd, rxnInd) = 0;
            ihuman.S(newMetInd, rxnInd) = stoichCoeffs;
        end
    end
end

% Rename the old pool metabolites, specifying that their composition is
% specific to liver tissue (see Mardinoglu et al. Nat Commun 2014, PMID:24419221)
for i = 1:size(poolMets,1)
    metInd = getIndexes(ihuman, poolMets(i,1), 'metnames');
    ihuman.metNames(metInd) = regexprep(ihuman.metNames(metInd), 'pool$', 'pool (liver tissue)');
end


%% Identify and document changes to the model

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

% determine and document model changes
modelChanges = docModelChanges(ihuman_orig,ihuman);
writeModelChanges(modelChanges,'../../ComplementaryData/modelCuration/lipidPools/lipidPool_modelChanges.tsv');

% export HumanGEM
exportHumanGEM(ihuman,'HumanGEM','../../',{'mat','yml'},false,false);

% clear unneeded variables
clearvars -except ihuman_orig ihuman modelChanges


