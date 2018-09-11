%
% FILE NAME:    curateDuplicatedRxns.m
% 
% DATE CREATED: 2018-09-07
%     MODIFIED: 2018-09-07
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: A script for identifying, investigating, and merging duplicated
%          reactions in the HumanGEM model. 
%          *** Note that this script is not meant to be run in its
%              entirety, but instead serves as a collection of steps
%              involved in the process of identifying and addressing
%              reaction duplications in the model.


%% Load model and identify duplicated reactions (ignoring rxn direction)

% load latest version of HumanGEM
load('humanGEM.mat');

% get all sets of duplicated reactions
[~,~,~,dupRxnSets] = checkDuplicateRxn_mod(ihuman,'FR');

% construct reaction equations
eqns = constructEquations(ihuman);

% retrieve equations for duplicated rxns
dupRxnEqns = arrayfun(@(i) [eqns(dupRxnSets{i});{''}],(1:length(dupRxnSets))','UniformOutput',false);

% organize rxn indices in a similar format (convert to strings)
dupRxnInds = arrayfun(@(i) [arrayfun(@num2str,dupRxnSets{i},'UniformOutput',false)';{''}],(1:length(dupRxnSets))','UniformOutput',false);

% obtain rxn names
dupRxnNames = arrayfun(@(i) [ihuman.rxns(dupRxnSets{i});{''}],(1:length(dupRxnSets))','UniformOutput',false);

% obtain and organize bounds for each rxn (convert to strings)
dupRxnLBs = arrayfun(@(i) [arrayfun(@num2str,ihuman.lb(dupRxnSets{i}),'UniformOutput',false);{''}],(1:length(dupRxnSets))','UniformOutput',false);
dupRxnUBs = arrayfun(@(i) [arrayfun(@num2str,ihuman.ub(dupRxnSets{i}),'UniformOutput',false);{''}],(1:length(dupRxnSets))','UniformOutput',false);

% obtain and organize grRules for each rxn
dupRxnRules = arrayfun(@(i) [ihuman.grRules(dupRxnSets{i});{''}],(1:length(dupRxnSets))','UniformOutput',false);

% assemble all information into a single cell array
dupRxnData = {'index','rxn','eqn','LB','UB','grRule'};
for i = 1:length(dupRxnSets)
    dupRxnData = [dupRxnData; [dupRxnInds{i}, dupRxnNames{i}, dupRxnEqns{i}, dupRxnLBs{i}, dupRxnUBs{i}, dupRxnRules{i}]];
end

% write data to file
writecell(dupRxnData,'dupRxnData.txt',true,'\t');



% this next part depends on all the duplicated reactions simply being the
% reverse form of each other
isrev = [];
for i = 1:length(dupRxnSets)
    % check if they actually are the reverse of each other
    isrev(i,1) = isequal(ihuman.S(:,dupRxnSets{i}(1)),-ihuman.S(:,dupRxnSets{i}(2)));
end
if all(isrev) && all(cellfun(@numel,dupRxnSets) == 2)
    fprintf('All duplicated reaction sets are just the reverse of one another.\n');
    
    % assess which duplicated reactions have equal bounds or grRules
    equal_bounds = false(length(dupRxnSets),1);  % initialize variable
    equal_rules = false(length(dupRxnSets),1);  % initialize variable
    for i = 1:length(dupRxnSets)
        equal_bounds(i) = (ihuman.lb(dupRxnSets{i}(1)) == -ihuman.ub(dupRxnSets{i}(2))) ...
                       && (ihuman.ub(dupRxnSets{i}(1)) == -ihuman.lb(dupRxnSets{i}(2)));
        equal_rules(i) = strcmp(ihuman.grRules(dupRxnSets{i}(1)),ihuman.grRules(dupRxnSets{i}(2)));
    end
end








%% OLD SCRIPT
% 
% % Below is an old version of a previous script used for analyzing duplicate
% % reactions ("findDuplicateRxns.m"). It will be updated/deleted in the near
% % future.
% 
% 
% %% load latest Recon4 model
% x = load('ComplementaryScripts/ModelIntegration/HMR3_02.mat');
% f = fields(x);
% model = x.(f{1});  % ugly way to rename model without knowing variable name in advance
% 
% 
% 
% %% Add new field to model: rxnRecon3DID
% % This is to keep track of which reactions from Recon3D were "merged" with
% % each of the HMR reactions when generating the combined model (Recon4).
% 
% % load Recon3D-HMR rxn mapping information (loads "Recon3D" variable)
% load('ComplementaryScripts/RxnAssociation/Recon3Rxns2HMR.mat');
% model.rxnRecon3DID = repmat({''},size(model.rxns));
% 
% % get unique list of all HMR rxn IDs mapped to the Recon3D rxns
% HMRrxns = unique(flattenCell(cellfun(@(r) strsplit(r,';'),Recon3D.rxnHMRID,'UniformOutput',false),true));
% HMRrxns(ismember(HMRrxns,{''})) = [];  % remove empty cell if present
% 
% % associate Recon3D rxnIDs to each of the reactions in Recon4
% for i = 1:length(HMRrxns)
%     r3ind = contains(Recon3D.rxnHMRID,HMRrxns{i});
%     r4ind = ismember(model.rxns,HMRrxns{i});
%     model.rxnRecon3DID{r4ind} = strjoin(Recon3D.rxns(r3ind),'; ');
% end
% 
% 
% %% Identify duplicated reactions in Recon4
% 
% % specify names of metabolites to ignore (e.g., protons)
% ignoreMets = {'H+'};
% ignoreMetInds = ismember(model.metNames,ignoreMets);
% if any(~ismember(ignoreMets,model.metNames))
%     error('At least one of the specified ignoreMets was not found in the model.');
% end
% 
% % Set the stoich matrix entries for ignored mets to zero, except in
% % reactions where the ignored met appears in multiple compartments. This is
% % to avoid flagging reactions involving, for example, proton transport, as
% % the same as reactions without proton transport.
% model_temp = model;
% ignoreRxnInds = true(size(model.rxns));
% for i = 1:length(ignoreMets)
%     ind = ismember(model.metNames,ignoreMets(i));
%     ignoreRxnInds(sum(model.S(ind,:) ~= 0,1) > 1) = false;
% end
% model_temp.S(ignoreMetInds,ignoreRxnInds) = 0;
% 
% % Check for duplicated reactions based on S-matrix.
% % Note that reactions with stoichiometry that differs only by scalar
% % multiplication will be treated as identical, and reactions that are
% % identical but in reverse direction from one another will NOT be treated
% % as identical.
% [~,~,~,duplicateRxnSets] = checkDuplicateRxn_mod(model_temp,'S',false);
% 
% % generate reaction equations for manual inspection of reactions
% rxnEqns = constructEquations(model);
% 
% % organize duplicateRxnSets into easily comparable rxn eqn format to
% % facilitate manual curation
% duplicateRxnInfo = {};
% for i = 1:length(duplicateRxnSets)
%     duplicateRxnInfo = [duplicateRxnInfo;[{''},{''}];[model.rxns(duplicateRxnSets{i}),rxnEqns(duplicateRxnSets{i})]];
% end
% 




