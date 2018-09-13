%
% FILE NAME:    curateDuplicatedRxns.m
% 
% DATE CREATED: 2018-09-07
%     MODIFIED: 2018-09-13
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: A script for identifying, investigating, and exporting
%          information on duplicated reactions in the model.
%          *** Note that this script is not meant to be run in its
%              entirety, but instead serves as a collection of steps
%              involved in the process of identifying and addressing
%              reaction duplications in the model.




%% Load model and identify duplicated reactions (ACCOUNT FOR rxn direction)

% load latest version of HumanGEM
load('humanGEM.mat');

% get all sets of duplicated reactions
[~,~,~,dupRxnSets] = checkDuplicateRxn_mod(ihuman,'S');

% skip next step if no duplicate sets, or they involve more than two rxns
if ~isempty(dupRxnSets) && all(cellfun(@numel,dupRxnSets) == 2)
    
    % convert to matrix format, and get associated rxn names
    dupRxnInds = vertcat(dupRxnSets{:});
    dupRxnNames = ihuman.rxns(dupRxnInds);
    
    % check that duplicate sets come from different models (i.e., one from
    % HMR, and one from Recon3D), otherwise remove
    rem_sets = sum(startsWith(dupRxnNames,'HMR_'),2) ~= 1;
    dupRxnInds(rem_sets,:) = [];
    dupRxnNames(rem_sets,:) = [];
    
    % organize pairs so that first entry is always HMR rxn
    flipRow = ~startsWith(dupRxnNames(:,1),'HMR_');
    dupRxnNames(flipRow,:) = fliplr(dupRxnNames(flipRow,:));
    
    % load rxnAssoc.mat
    load('ComplementaryScripts/modelIntegration/rxnAssoc.mat');
    
    % check if any of the associations already exist in rxnAssoc 
    % (join columns by arbitrary string "***" to enable use of ismember)
    rem_sets = ismember(join(dupRxnNames,'***'),join([rxnAssoc.rxnHMRID,rxnAssoc.rxnRecon3DID],'***'));
    dupRxnInds(rem_sets,:) = [];
    dupRxnNames(rem_sets,:) = [];
    
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
    
    % save new rxnAssoc.m structure
    % NOTE: THIS SCRIPT/SAVE WAS LAST RUN ON 2018-09-13
    save('rxnAssoc','rxnAssoc');
    
end




%% Load model and identify duplicated reactions (IGNORING rxn direction)

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



% This next part depends on all the duplicated reactions simply being the
% reverse form of each other. First, remove sets with more than two rxns
rem_sets = (cellfun(@numel,dupRxnSets) ~= 2);
dupRxnSets(rem_sets) = [];

% next check that they are all reverses of each other
isrev = [];
for i = 1:length(dupRxnSets)
    isrev(i,1) = isequal(ihuman.S(:,dupRxnSets{i}(1)),-ihuman.S(:,dupRxnSets{i}(2)));
end

% remove those that are not reverses of each other.
dupRxnSets(~isrev) = [];


% now process the remaining sets

% evaluate some properties about each set of duplicated reactions
equal_bounds = false(length(dupRxnSets),1);  % initialize variable
equal_rules = false(length(dupRxnSets),1);  % initialize variable
from_diff_model = false(length(dupRxnSets),1);  %initialize variable
for i = 1:length(dupRxnSets)
    % check if the rxns have equal bounds
    equal_bounds(i) = (ihuman.lb(dupRxnSets{i}(1)) == -ihuman.ub(dupRxnSets{i}(2))) ...
        && (ihuman.ub(dupRxnSets{i}(1)) == -ihuman.lb(dupRxnSets{i}(2)));
    
    % check if the rxns have the same grRules
    equal_rules(i) = strcmp(ihuman.grRules(dupRxnSets{i}(1)),ihuman.grRules(dupRxnSets{i}(2)));
    
    % check if the rxns come from different sources (i.e., one from HMR, and one from Recon3D)
    from_diff_model(i) = sum(startsWith(ihuman.rxns(dupRxnSets{i}),'HMR_')) == 1;
end


% update rxnAssoc.mat ONLY with the reaction pairs that came from two
% different models (one from HMR, one from Recon3D).
load('ComplementaryScripts/modelIntegration/rxnAssoc.mat');  % load rxnAssoc structure
addRxnAssocInds = dupRxnSets(from_diff_model);  % select duplicate rxn sets to add

% ensure that the HMR rxn is listed first, and the Recon3D rxn second
for i = 1:length(addRxnAssocInds)
    if ~startsWith(ihuman.rxns(addRxnAssocInds{i}(1)),'HMR_')
        addRxnAssocInds{i} = fliplr(addRxnAssocInds{i});
    end
end

% convert index variable to a matrix
addRxnAssocInds = vertcat(addRxnAssocInds{:});

% retrieve info about rxn names
addRxnAssocNames = ihuman.rxns(addRxnAssocInds);

% check if any of the associations already exist in rxnAssoc
% (join columns by arbitrary string "***" to enable use of ismember)
exclude_ind = ismember(join(addRxnAssocNames,'***'),join([rxnAssoc.rxnHMRID,rxnAssoc.rxnRecon3DID],'***'));
addRxnAssocInds(exclude_ind,:) = [];
addRxnAssocNames(exclude_ind,:) = [];

% retrieve associated bounds ('lb' and 'ub') for each reaction
% NOTE: need to switch directions (LB = -UB, and UB = -LB) for the rev rxn
addRxnAssocLB = [ihuman.lb(addRxnAssocInds(:,1)), -ihuman.ub(addRxnAssocInds(:,2))];
addRxnAssocUB = [ihuman.ub(addRxnAssocInds(:,1)), -ihuman.lb(addRxnAssocInds(:,2))];

% add information to rxnAssoc structure
rxnAssoc.rxnHMRID = [rxnAssoc.rxnHMRID; addRxnAssocNames(:,1)];
rxnAssoc.lbHMR = [rxnAssoc.lbHMR; addRxnAssocLB(:,1)];
rxnAssoc.ubHMR = [rxnAssoc.ubHMR; addRxnAssocUB(:,1)];
rxnAssoc.rxnRecon3DID = [rxnAssoc.rxnRecon3DID; addRxnAssocNames(:,2)];
rxnAssoc.lbRecon3D = [rxnAssoc.lbRecon3D; addRxnAssocLB(:,2)];
rxnAssoc.ubRecon3D = [rxnAssoc.ubRecon3D; addRxnAssocUB(:,2)];

% save new rxnAssoc.m structure
% NOTE: THIS SCRIPT/SAVE WAS LAST RUN ON 2018-09-13
save('rxnAssoc','rxnAssoc');








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




