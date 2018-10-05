%
% FILE NAME:    miscModelCurationScript_20181005.m
% 
% DATE CREATED: 2018-10-02
%     MODIFIED: 2018-10-05
% 
% PROGRAMMER:   Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: This script undertakes two tasks: 
%          
%          1. Detect humanGEM v0.4.0 for addtional duplicate reaction pairs
%          by ignoring the reactions direction. These new pairs are
%          subjected to manual check, subsequently added to the rxnAssoc.mat
%          structure and the ones from Recon3D are then removed;
%          
%          2. Incoporate enzyme complex information from CORUM database
%          to humanGEM grRules that include only "OR" relations, followed
%          by manual verification from UniProt and NCBI databases and a
%          redundant removal step (deleting subunit encoding genes from
%          isoenzymes)
%

%% Load model and rxnAssoc

% load latest version of HumanGEM
load('humanGEM.mat');  % v0.4.0
rxnEqns = constructEquations(ihuman);   % construct reaction equations

% load rxnAssoc.mat
load('rxnAssoc.mat');



%% Identify duplicate reaction pairs by ignoring reaction direction

% get the sets of duplicated reactions
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

% check that duplicate sets and make sure one is from
% HMR, and the other from Recon3D
rem_sets = sum(startsWith(dupRxnNames,'HMR_'),2) ~= 1;
dupRxnInds(rem_sets,:) = [];
dupRxnNames(rem_sets,:) = [];
fprintf('\t* %u duplicate pairs remain after removing those from the same source model (HMR or Recon3D).\n',size(dupRxnInds,1));

% organize pairs so that first entry is always the HMR rxn, not needed
% flipRow = ~startsWith(dupRxnNames(:,1),'HMR_');
% dupRxnNames(flipRow,:) = fliplr(dupRxnNames(flipRow,:));

% check if any of the associations already exist in rxnAssoc 
rem_sets = ismember(join(dupRxnNames,'***'),join([rxnAssoc.rxnHMRID,rxnAssoc.rxnRecon3DID],'***'));
dupRxnInds(rem_sets,:) = [];
dupRxnNames(rem_sets,:) = [];
fprintf('\t* %u duplicate pairs remain after removing those already present in rxnAssoc.mat.\n',size(dupRxnInds,1));

% print out for manual curation
% for i=1:length(dupRxnInds)
%    ind1 = find(strcmp(ihuman.rxns,dupRxnNames{i,1}));
%    ind2 = find(strcmp(ihuman.rxns,dupRxnNames{i,2}));
%    fprintf('%s\t%s\n%s\t%s\n\n',dupRxnNames{i,1},rxnEqns{ind1},dupRxnNames{i,2},rxnEqns{ind2});
% end

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
save('../modelIntegration/rxnAssoc.mat','rxnAssoc');
% print out summary
fprintf('\t* %u duplicate reaction pairs were added to rxnAssoc.mat.\n\n',size(dupRxnInds,1));

% Update the rxnRecon3DID field in the model
for i = 1:length(dupRxnInds)
    [~,ind] = ismember(dupRxnNames(i,1),ihuman.rxns);
    if ~isempty(ihuman.rxnRecon3DID{ind})
        % ensure entries are unique, and separted by semicolon if multiple
        updatedEntry = join(unique([strsplit(ihuman.rxnRecon3DID{ind},';'),dupRxnNames(i,2)]),';');
    else
        updatedEntry = dupRxnNames(i,2);
    end
    ihuman.rxnRecon3DID(ind) = updatedEntry;
end

% update .metFrom field by filling the blank entries
% derived from adding boundary metabolite step
indEmpty = find(cellfun('isempty', ihuman.metFrom));
ihuman.metFrom(indEmpty) = {'Recon3D'};
fprintf('A total of %d blank entries were filled in the .metFrom field.\n\n',length(indEmpty));

% delete duplicated reactions
model = removeReactionsFull(ihuman,dupRxnNames(:,2));
fprintf('A total of %d reactions were deleted from the model.\n\n',length(dupRxnNames));


%% Incoporate curated grRules based on the CORUM enzyme complexes database

% add curated grRules from automatic incoporation of enzyme
% complexes info in CORUM database and manual checking, as
% well as additional redundant gene removal step
newModel = addCuratedComplexRulesToModel(model, 'curated_CORUM_grRules_20180924');

% Get the index of reactions with modified grRules
indModifiedRxn = find(~strcmp(model.grRules, newModel.grRules));
fprintf('A total of %d grRules are changed with curation based on CORUM database.\n',length(indModifiedRxn));



%% Save model and clear variables
ihuman = newModel;
save('../../ModelFiles/mat/humanGEM.mat','ihuman');
clear;

