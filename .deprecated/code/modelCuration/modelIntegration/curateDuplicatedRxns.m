%
% FILE NAME:    curateDuplicatedRxns.m
% 
% PURPOSE: This script updates the humanGEM model and the rxnAssoc.mat
%          reaction association structure, with the following steps:
%           1) Add boundary metabolites to humanGEM
%               - Boundary metabolites are added to balance reactions that
%                 contain only one extracellular metabolite, to be
%                 consistent with RAVEN model formulations.
%               - In both "metFrom" and "rxnFrom" fields, entries of
%                  "reducedRecon3D" are changed to "Recon3D".
%           2) Duplicate reactions (considering direction) are identified
%               - New duplicated reactions are identified, focusing only on
%                 those that are written in the same direction, and
%                 originate from different models (i.e., one from HMR, and
%                 one from Recon3D)
%               - New reaction duplicate pairs are added to the
%                 rxnAssoc.mat structure.
%           3) Duplicate reactions (ignoring direction) are identified
%               - New duplicated reactions are identified, focusing only on
%                 those that are the reverse of one another, and originate
%                 from different models (i.e., one from HMR, and the other
%                 from Recon3D).
%               - New reaction duplicate pairs are added to the
%                 rxnAssoc.mat structure.
%           4) Duplicate reactions are removed from the humanGEM model
%               - All pairs of duplicate reactions that were identified in
%                 steps 2 and 3 are removed from the model structure.
%               - New entries are added to the rxnRecon3DID field of the
%                 model, based on these new duplicate reaction pairs.
%           5) Update GPR rules and associated model fields in humanGEM
%               - The combineModelGPRs function is run, taking the newly
%                 modified humanGEM model as input. This updates the GPRs,
%                 by combining information from grRules from iHsa and
%                 Recon3D into the current model.
%           6) Save new humanGEM model and rxnAssoc.mat structure
%               - The newly updated humanGEM model, as well as the updated
%                 rxnAssoc.mat structure, are saved. These files are saved
%                 as humanGEM_new.mat and rxnAssoc_new.mat, respectively,
%                 to avoid unintended file overwrites. These should be used
%                 to manually overwrite the files when it is determined
%                 that they are ready.
%


%% Load important items

% load latest version of HumanGEM
load('humanGEM.mat');
ihuman_orig = ihuman;  % save original so we can compare later if needed

% load rxnAssoc.mat
load('ComplementaryScripts/modelIntegration/rxnAssoc.mat');
rxnAssoc_orig = rxnAssoc;  % save original so we can compare later if needed



%% Add boundary metabolites to HumanGEM

prevMetNum = length(ihuman.mets);  % record number of mets before addition

% add boundary metabolites
ihuman = addBoundaryMets(ihuman);
% Boundary metabolites were added to 1639 reactions.
% New (boundary) versions of 1198 metabolites were added to the model.

newMetNum = length(ihuman.mets);

% update the .metFrom field, since it has blank entries for the new mets
if prevMetNum ~= newMetNum
    ihuman.metFrom(prevMetNum+1:end) = {'Recon3D'};
end
% also go ahead and change any "reducedRecon3D" entries to just "Recon3D"
ihuman.metFrom(ismember(ihuman.metFrom,'reducedRecon3D')) = {'Recon3D'};

% update the .rxnFrom field by changing any "reducedRecon3D" entries to "Recon3D"
ihuman.rxnFrom(ismember(ihuman.rxnFrom,'reducedRecon3D')) = {'Recon3D'};

% clear intermediate variables
clear('prevMetNum','newMetNum')


%% Identify duplicated reactions, ACCOUNTING FOR rxn direction

% get all sets of duplicated reactions
fprintf('Searching for duplicate reactions, ACCOUNTING FOR reaction direction...\n');
[~,~,~,dupRxnSets] = checkDuplicateRxn_mod(ihuman,'S');

% check if any duplicate sets were actually found
if isempty(dupRxnSets)
    fprintf('\t* No duplicate rxns were found.\n\n');
else
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
    flipRow = ~startsWith(dupRxnNames(:,1),'HMR_');
    dupRxnNames(flipRow,:) = fliplr(dupRxnNames(flipRow,:));
    
    % check if any of the associations already exist in rxnAssoc 
    % (join columns by arbitrary string "***" to enable use of ismember)
    rem_sets = ismember(join(dupRxnNames,'***'),join([rxnAssoc.rxnHMRID,rxnAssoc.rxnRecon3DID],'***'));
    dupRxnInds(rem_sets,:) = [];
    dupRxnNames(rem_sets,:) = [];
    fprintf('\t* %u duplicate pairs remain after removing those already present in rxnAssoc.mat.\n',size(dupRxnInds,1));
    
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
    
    % print some info
    if isempty(dupRxnInds)
        fprintf('\t* rxnAssoc.mat was not modified.\n\n');
    else
        fprintf('\t* %u duplicate reaction pairs were added to rxnAssoc.mat.\n\n',size(dupRxnInds,1));
    end
    
end

% clear intermediate variables
clear('dupRxnInds','dupRxnLB','dupRxnNames','dupRxnSets','dupRxnUB','flipRow','rem_sets');

% there are 9 pairs of duplicate reactions identified in this section

%% Identify duplicated reactions, IGNORING rxn direction

% get all sets of duplicated reactions
fprintf('Searching for duplicate reactions, IGNORING reaction direction...\n');
[~,~,~,dupRxnSets] = checkDuplicateRxn_mod(ihuman,'FR');

% The optional code below can be run to export a report of the identified
% reaction duplicates, but is skipped by default.
if ( false )
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
    writecell2file(dupRxnData,'dupRxnData.txt',true,'\t');
end


% check if any duplicate sets were actually found
if isempty(dupRxnSets)
    fprintf('\t* No duplicate rxns were found.\n\n');
else
    
    % ignore sets that contain more than two duplicates
    if ~all(cellfun(@numel,dupRxnSets) == 2)
        fprintf('\t* Some reactions were duplicated more than once; these will be ignored for now.\n');
        dupRxnSets(cellfun(@numel,dupRxnSets) ~= 2) = [];  % remove these sets
    end
    fprintf('\t* A total of %u duplicate reaction pairs were found.\n',length(dupRxnSets));
    
    % Only consider duplicated reactions that are the reverse of one another.
    isrev = [];
    for i = 1:length(dupRxnSets)
        isrev(i,1) = isequal(ihuman.S(:,dupRxnSets{i}(1)),-ihuman.S(:,dupRxnSets{i}(2)));
    end
    
    % remove those that are not reverses of each other.
    dupRxnSets(~isrev) = [];
    fprintf('\t* %u duplicate pairs remain after removing those that are NOT the reverse of one another.\n',length(dupRxnSets));
    
    % Evaluate some properties about each set of duplicated reactions. Some
    % of the information will not be used, but may be nice to have.
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
    
    % Only condsider reaction pairs that came from two different source 
    % models (i.e., one from HMR, one from Recon3D).
    addRxnAssocInds = dupRxnSets(from_diff_model);  % select duplicate rxn sets to add
    fprintf('\t* %u duplicate pairs remain after removing those from the same source model (HMR or Recon3D).\n',length(addRxnAssocInds));
    
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
    fprintf('\t* %u duplicate pairs remain after removing those already present in rxnAssoc.mat.\n',size(addRxnAssocInds,1));
    
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
    
    % print some info
    if isempty(addRxnAssocInds)
        fprintf('\t* rxnAssoc.mat was not modified.\n\n');
    else
        fprintf('\t* %u duplicate reaction pairs were added to rxnAssoc.mat.\n\n',size(addRxnAssocInds,1));
    end
    
end

% clear intermediate variables
clear('addRxnAssocInds','addRxnAssocLB','addRxnAssocUB','addRxnAssocNames','i',...
    'dupRxnSets','equal_bounds','equal_rules','exclude_ind','from_diff_model','isrev');

% there are 976 pairs of duplicate reactions identified in this section

%% Remove duplicate reactions from the model

% identify rxns to be deleted (only new Recon3D IDs added to rxnAssoc)
delRxnNames = rxnAssoc.rxnRecon3DID(length(rxnAssoc_orig.rxnRecon3DID)+1:end);
delRxnInds = ismember(ihuman.rxns,delRxnNames);

% also get a list of the corresponding HMR rxn IDs that will be kept
keepRxnNames = rxnAssoc.rxnHMRID(length(rxnAssoc_orig.rxnHMRID)+1:end);

% before deleting rxns, update the rxnRecon3DID field in the model
for i = 1:length(keepRxnNames)
    [~,ind] = ismember(keepRxnNames(i),ihuman.rxns);
    if ~isempty(ihuman.rxnRecon3DID{ind})
        % ensure entries are unique, and separted by semicolon if multiple
        updatedEntry = join(unique([strsplit(ihuman.rxnRecon3DID{ind},';'),delRxnNames(i)]),';');
    else
        updatedEntry = delRxnNames(i);
    end
    ihuman.rxnRecon3DID(ind) = updatedEntry;
end

% delete model reactions
ihuman = removeReactions(ihuman, delRxnNames, true, true, true);
fprintf('A total of %u reactions were deleted from the model.\n\n',length(delRxnNames));

% manually update fields that are not recognized by the removeReactions function
ihuman.rxnRecon3DID(delRxnInds) = [];
ihuman.priorCombiningGrRules(delRxnInds) = [];

% don't actually need to update these fields, because they will be
% overwritten by the "combineModelGPRs" function anyway.
% ihuman.prRules(delRxnInds) = [];
% ihuman.rxnProtMat(delRxnInds,:) = [];



%% Update GPR rules based on newly identified reaction associations

% use existing function, which will automatically load and use the iHsa and
% Recon3D models
ihuman = combineModelGPRs(ihuman);


%% Save files

% save new model and rxnAssoc.mat structure
save('../../model/Human-GEM.mat','ihuman');
save('rxnAssoc.mat','rxnAssoc');








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




