%
% FILE NAME:    miscModelCurationScript_20181004.m
% 
% DATE CREATED: 2018-10-02
%     MODIFIED: 2018-10-04
% 
% PROGRAMMER:   Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: This script detects humanGEM v0.4.0 for addtional duplicate 
%          pairs by ignoring the reactions direction. These new pairs
%          are subjected to manual check and then added to the rxnAssoc
%      



%% Load important items

% load latest version of HumanGEM
load('humanGEM.mat');  % v0.4.0
ihuman.rxnEqns = constructEquations(ihuman);   % construct rxn equations

% load rxnAssoc.mat
load('ComplementaryScripts/modelIntegration/rxnAssoc.mat');
rxnAssoc_orig = rxnAssoc;  % make backup


%% Identify duplicate reaction pairs by ignoring reaction direction

% get all sets of duplicated reactions
[~,~,~,dupRxnSets] = checkDuplicateRxn_mod(ihuman,'FR');

% ignore sets that contain more than two duplicates
if ~all(cellfun(@numel,dupRxnSets) == 2)
    fprintf('\t* Some reactions were duplicated more than once; these will be ignored for now.\n');
    dupRxnSets(cellfun(@numel,dupRxnSets) ~= 2) = [];  % remove these sets
end
fprintf('\t* A total of %u duplicate reaction pairs were found.\n',length(dupRxnSets));

% convert to matrix format, and get associated rxn names
dupRxnInds = vertcat(dupRxnSets{:});
dupRxnNames = ihuman.rxns(dupRxnInds);

% check that duplicate sets come from different models (i.e., one from
% HMR, and one from Recon3D), otherwise remove
rem_sets = sum(startsWith(dupRxnNames,'HMR_'),2) ~= 1;
dupRxnInds(rem_sets,:) = [];
dupRxnNames(rem_sets,:) = [];
fprintf('\t* %u duplicate pairs remain after removing those from the same source model (HMR or Recon3D).\n',size(dupRxnInds,1));

% organize pairs so that first entry is always the HMR rxn
% can be deleted
flipRow = ~startsWith(dupRxnNames(:,1),'HMR_');
dupRxnNames(flipRow,:) = fliplr(dupRxnNames(flipRow,:));

% check if any of the associations already exist in rxnAssoc 
% (join columns by arbitrary string "***" to enable use of ismember)
rem_sets = ismember(join(dupRxnNames,'***'),join([rxnAssoc.rxnHMRID,rxnAssoc.rxnRecon3DID],'***'));
dupRxnInds(rem_sets,:) = [];
dupRxnNames(rem_sets,:) = [];
fprintf('\t* %u duplicate pairs remain after removing those already present in rxnAssoc.mat.\n',size(dupRxnInds,1));

% print out for manual curation
for i=1:length(dupRxnInds)
    ind1 = find(strcmp(ihuman.rxns,dupRxnNames{i,1}));
    ind2 = find(strcmp(ihuman.rxns,dupRxnNames{i,2}));
    fprintf('%s\t%s\n%s\t%s\n\n',dupRxnNames{i,1},ihuman.rxnEqns{ind1},dupRxnNames{i,2},ihuman.rxnEqns{ind2});
end

% get bounds for associated rxns
dupRxnLB = ihuman.lb(dupRxnInds);
dupRxnUB = ihuman.ub(dupRxnInds);

% add information to rxnAssoc structure
rxnAssoc.rxnHMRID = [rxnAssoc.rxnHMRID; dupRxnNames(:,1)];
rxnAssoc.lbHMR = [rxnAssoc.lbHMR; dupRxnLB(:,1)];
rxnAssoc.ubHMR = [rxnAssoc.ubHMR; dupRxnUB(:,1)];
rxnAssoc.rxnRecon3DID = [rxnAssoc.rxnRecon3DID; dupRxnNames(:,2)];
rxnAssoc.lbRecon3D = [rxnAssoc.lbRecon3D; dupRxnLB(:,2)];
rxnAssoc.ubRecon3D = [rxnAssoc.ubRecon3D; dupRxnUB(:,2)];

%% Save rxnAssoc.mat structure
save('rxnAssoc.mat','rxnAssoc');
% print out summary
fprintf('\t* %u duplicate reaction pairs were added to rxnAssoc.mat.\n\n',size(dupRxnInds,1));


%% Incoporate curated grRules based on the CORUM enzyme complexes database

newModel = addCuratedComplexRulesToModel(ihuman, 'curated_CORUM_grRules_20180924');

% Get the index of reactions with modified grRules
indModifiedRxn = find(~strcmp(ihuman.grRules, newModel.grRules));
fprintf('A total of %d grRules are changed with curation from CORUM database.\n',length(indModifiedRxn));



%% Save model
ihuman = newModel;
save('ihuman.mat','ihuman');



