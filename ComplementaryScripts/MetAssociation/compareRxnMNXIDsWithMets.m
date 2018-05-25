function results = compareRxnMNXIDsWithMets(model,mnx)
%compareRxnMNXIDsWithMets  Check consistency between rxn and met MNXIDs.
%
% compareRxnMNXIDsWithMets determines if the reaction MNXID associations
% have distributed their information to the metabolite MNXID associations.
% For each metabolite in the model, the set of reactions in which it
% participates is collected. The rxnMNXIDs assigned to this set of
% reactions is then obtained, and the corresponding reaction information is
% retreieved from the MNX database. If there is no overlap between the set
% of metabolite MNX IDs participating in the reaction, and the set of 
% MNXIDs assigned to the model metabolite participating in the reaction,
% the reaction and metabolite will be flagged.
%
% USAGE:
%
%   results = compareRxnMNXIDsWithMets(model,mnx);
%
% INPUT:
%
%   model    A model structure containing reaction and metabolite MNX ID
%            association fields ("rxnMNXID" and "metMNXID", respectively).
%
%   mnx      (Optional) An MNX database structure, containing reaction-
%            related information retrieved from the MNX database, generated
%            using the following command: mnx = buildMNXmodel('rxn');
%            By default, the function will automatically run the above
%            command to regenerate the MNX database structure (slower).
%
% OUTPUT:
%
%   results  A results structure containing information on the flagged
%            reactions and the associated mets and MNXIDs, organized as a
%            cell array.
%
%
% Jonathan Robinson 2018-05-25


% handle input arguments
if nargin < 2
    mnx = buildMNXmodel('rxn');
end

% if metMNXID and/or rxnMNXID field contains multiple columns, consolidate 
% into a single column with nested cell entries
if size(model.metMNXID,2) > 1
    model.metMNXID = nestCell(model.metMNXID,true);
end
if size(model.rxnMNXID,2) > 1
    model.rxnMNXID = nestCell(model.rxnMNXID,true);
end

% initialize results structure
results = {'model met','MNXIDs assigned to model met','flagged model rxn','flagged MNXIDs mapped to rxn'};

% iterate through each of the model metabolites
h = waitbar(0,'Processing metabolites...');
for i = 1:length(model.mets)
    
    % identify all reactions in which the metabolite is involved
    rxn_ind = find(model.S(i,:) ~= 0);
    
    % iterate through each of the reactions
    for r = rxn_ind
        if isempty(model.rxnMNXID{r})
            % skip reactions that have not been mapped to any MNXIDs
            continue
        end
        
        % obtain MNXIDs of mets participating in each of those MNX rxns
        [~,MNX_rxnInds] = ismember(model.rxnMNXID{r},mnx.rxns);
        MNX_metIDs = mnx.rxnMets(MNX_rxnInds);
        
        % identify rxns where none of the mets are associated with the
        % current model metabolite
        flag_MNX_rxn = cellfun(@(m) ~any(ismember(m,model.metMNXID{i})),MNX_metIDs);
        
        % add flagged MNX rxn and associated model met and rxn info to
        % results structure
        if any(flag_MNX_rxn)
            if isempty(model.metMNXID{i})
                results = [results; [model.mets(i), {''}, ...
                                     model.rxns(r), strjoin(mnx.rxns(MNX_rxnInds(flag_MNX_rxn)),'; ')]];
            else
                results = [results; [model.mets(i), strjoin(model.metMNXID{i},'; '), ...
                                     model.rxns(r), strjoin(mnx.rxns(MNX_rxnInds(flag_MNX_rxn)),'; ')]];
            end
        end
    end
    waitbar(i/length(model.mets),h);
end
close(h);

if size(results,1) == 1
    fprintf('No reactions were flagged!\n');
end










