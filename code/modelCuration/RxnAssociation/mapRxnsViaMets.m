function results = mapRxnsViaMets(model,mnx,mapRxns)
%mapRxnsViaMets  Map model reactions to MNX database via metabolite MNXIDs.
%
% mapRxnsViaMets searches an MNX database structure for reactions that
% share the same set of metabolites as reactions in model, and returns the
% corresponding MNX reaction IDs and reaction equations. The mapping is
% performed by looking for any MNX reactions that contain the same set of
% metabolite MNX IDs as a model reaction, meaning that STOICHIOMETRY
% COEFFICIENTS ARE IGNORED.
%
% This function can also handle model metabolites associated with multiple
% MNX IDs. In this case, only one of the multiple met IDs for a particular
% metabolite need to match a metabolite in the corresponding MNX reaction.
%
% USAGE:
%
%   results = mapRxnsViaMets(model,mnx,mapRxns)
%
% INPUT:
%
%   model    Model structure.
%
%   mnx      (Optional, will be generated if not provided) MNX database
%            structure, generated using the following command: 
%            mnx = buildMNXmodel('both');
%
%   mapRxns  (Optional, Default = all reactions) A logical vector of the
%            same size as model.rxns, indicating which reactions should be
%            mapped. Reactions corresponding to mapRxns entries that are
%            FALSE will be skipped.
%
% OUTPUT:
%
%   results  A results structure with the following fields:
%
%       rxnMNXID: A column cell array containing the MNX IDs mapped to each
%                 reaction. For reactions to which multiple MNX IDs were
%                 mapped, the corresponding rxnMNXID entry will be a nested
%                 cell of MNX IDs.
%
%      rxnMNXeqn: A column cell array containing the reaction equations
%                 from the MNX database that correspond to the MNX ID(s) 
%                 mapped to each model reaction.
%
%          notes: A column of strings providing information on if and how
%                 each model reaction was mapped to any MXN ID(s).
%
%
% Jonathan Robinson, 2018-05-24


% handle input arguments
if nargin < 3
    mapRxns = true(size(model.rxns));
elseif ~isequal(size(mapRxns),size(model.rxns))
    error('mapRxns input must have same dimensions as model.rxns');
end
if nargin < 2
    mnx = [];
end

% initialize results structure
results.rxnMNXID = {};
results.rxnMNXeqn = {};
results.notes = repmat({''},size(model.rxns));
results.notes(~mapRxns) = {'IGNORED'};


%% Prepare model and MNX structures for processing

fprintf('Pre-processing model and MNX structures... ');

% if the model.metMNXID field contains multiple columns, convert it to a 
% single column of nested cells
if size(model.metMNXID,2) > 1
    model.metMNXID = nestCell(model.metMNXID,true);
elseif all(cellfun(@ischar,model.metMNXID)) && any(contains(model.metMNXID,';'))
    empty_inds = cellfun(@isempty,model.metMNXID);
    model.metMNXID = cellfun(@(id) strsplit(id,';\s*','DelimiterType','RegularExpression'),model.metMNXID,'UniformOutput',false);
    model.metMNXID(empty_inds) = {''};  % to deal with cells that should be empty, but are recognized as non-empty
end

% generate binary (logical) version of S-matrix indicating which mets are
% involved in each rxn
S = (full(model.S) ~= 0);

% do not try to map model rxns that contain a metabolite without an MNXID
met_noID = cellfun(@isempty,model.metMNXID);
rxn_noID = any(S(met_noID,:),1)';
results.notes(mapRxns & rxn_noID) = {'SKIPPED: one or more mets lacking MNXID'};
mapRxns(rxn_noID) = false;

% Create an additional field in model that lists all possible metMNXIDs
% associated with each reaction.
for i = 1:length(model.rxns)
    model.rxnMetsMNX{i,1} = unique(horzcat(model.metMNXID{S(:,i)}));
end


% For faster processing later on, remove all rxns in the MNX structure that
% contain at least one metabolite that is not present in the model.

% get list of all metMNXIDs in the model (for reactions to be mapped)
met_inds = any(S(:,mapRxns),2);
allMetMNXIDs = unique(horzcat(model.metMNXID{met_inds}))';
allMetMNXIDs(cellfun(@isempty,allMetMNXIDs)) = [];

% load MNX database structure if not provided as input
if isempty(mnx)
    fprintf('MNX database structure not provided. It will be loaded.\n\n');
    mnx = buildMNXmodel('both');
end

% flatten mnx.rxnMets (single column -> multi column cell array)
mnx.rxnMets = flattenCell(mnx.rxnMets,true);

% find MNX rxns with mets that aren't present in model
mismatches = ~ismember(mnx.rxnMets,allMetMNXIDs);
mismatches(cellfun(@isempty,mnx.rxnMets)) = false;
del_rxns = any(mismatches,2);

% restore mnx.rxnMets to its original format of a column of nested cells
mnx.rxnMets = nestCell(mnx.rxnMets,true);

% remove reactions from MNX structure
nrxns = length(mnx.rxns);
f = fields(mnx);
for i = 1:length(f)
    if iscell(mnx.(f{i})) && size(mnx.(f{i}),1) == nrxns
        mnx.(f{i})(del_rxns,:) = [];
    end
end

% initialize rxn MNX IDs
rxnMNXID = repmat({''},size(model.rxns));

fprintf('Done.\n');


%% STAGE 1: Find exact matches to MNX reactions that are balanced

fprintf('Searching for exact, balanced rxn matches... ');

% subset MNX structure to include only the balanced reactions
mnx_bal = mnx;
unbal = ~ismember(mnx_bal.rxnBalanced,'true');
mnx_bal.rxns(unbal) = [];
mnx_bal.rxnMets(unbal) = [];

% map reactions
rxnMNXID(mapRxns) = findMNXrxns(model,mnx_bal,mapRxns);

% find newly mapped rxns
newMapped = mapRxns & ~cellfun(@isempty,rxnMNXID);
results.notes(newMapped) = {'MATCH: exact'};
mapRxns(newMapped) = false;  % update to ignore rxns that are now mapped

% clear MNX structure to free up some memory
clear mnx_bal

fprintf('Done.\n');


%% STAGE 2: Ignore protons (H+) and water (H2O) in reaction equations

fprintf('Searching for matches, ignoring balance status, protons (H+), and water (H2O)... ');

% remove H+ and H2O from all model reaction equations
rem_ind = ismember(model.metFormulas,{'H2O','H'});
model.S(rem_ind,:) = 0;

% also remove from MNX reaction equations
rem_ind = ismember(mnx.metFormulas,{'H2O','H'});
rem_ind(contains(mnx.metNames,{'(.)','deuterium','tritium','hydride'},'IgnoreCase',true)) = false;  % exclude other variants
rem_IDs = mnx.mets(rem_ind);
mnx.rxnMets = flattenCell(mnx.rxnMets,true);  % flatten cell
mnx.rxnMets(ismember(mnx.rxnMets,rem_IDs)) = {''};  % remove IDs
mnx.rxnMets = nestCell(mnx.rxnMets,true);  % re-nest cell

% map reactions
rxnMNXID(mapRxns) = findMNXrxns(model,mnx,mapRxns);

% find newly mapped rxns
newMapped = mapRxns & ~cellfun(@isempty,rxnMNXID);
results.notes(newMapped) = {'MATCH: ignored H+ and H2O'};

fprintf('Done.\n');


%% Finalize output

fprintf('Finalizing output... ');

results.notes(cellfun(@isempty,rxnMNXID) & mapRxns) = {'NO MATCHES FOUND'};

% retrieve MNX reaction equations for all mapped MNXIDs
rxnMNXeqn = repmat({''},size(model.rxns));
for i = 1:length(rxnMNXID)
    if ~isempty(rxnMNXID{i})
        rxnMNXeqn{i} = mnx.rxnEqnNames(ismember(mnx.rxns,rxnMNXID{i}));
    end
end

% add data to results structure
results.rxnMNXID = rxnMNXID;
results.rxnMNXeqn = rxnMNXeqn;

fprintf('Done.\n');

end




%% Additional functions

function rxnIDs = findMNXrxns(model,mnx,mapRxns)
%Finds MNX rxns with same set of mets as model rxns.

% convert stoich matrix to logicals, where all non-zero coeffs are TRUE
S = (full(model.S) ~= 0);

% exclude reactions that are not to be mapped
S(:,~mapRxns) = [];
model.rxns(~mapRxns) = [];
model.rxnMetsMNX(~mapRxns) = [];

% get number of mets in each reaction, for model and MNX
Nmets_per_rxn_model = sum(S,1)';
Nmets_per_rxn_mnx = cellfun(@numel,mnx.rxnMets);

% initialize output
rxnIDs = repmat({''},size(model.rxns));

% iterate through model reactions, searching for matching MNX reactions
h = waitbar(0,'Scanning reactions...');
for i = 1:length(model.rxns)
    
    waitbar(i/length(model.rxns),h);
    
    % find MNX rxns that have the same number of mets as the model rxn
    match = (Nmets_per_rxn_mnx == Nmets_per_rxn_model(i));
    if ~any(match)
        continue
    end
    
    % find MNX rxns whose mets match at least one of those associated with the model rxn
    match(match) = cellfun(@(m) all(ismember(m,model.rxnMetsMNX{i})),mnx.rxnMets(match));
    
    % now check if each of the mets in the model rxn are contained in the MNX rxns
    met_inds = find(S(:,i));
    for j = 1:length(met_inds)

        % stop checking metabolites if no further matches remain
        if ~any(match)
            break
        end
        
        % remove matching MNX rxns that don't contain the current model met
        match(match) = cellfun(@(m) any(ismember(model.metMNXID{met_inds(j)},m)),mnx.rxnMets(match));
    end
    
    % obtain IDs of all MNX rxns with a matching set of mets as the model rxn
    rxnIDs{i} = mnx.rxns(match);
    
end
close(h);

end

