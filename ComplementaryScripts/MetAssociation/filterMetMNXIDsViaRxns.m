function [fmodel,removed] = filterMetMNXIDsViaRxns(model,mnx,ignoreComp)
%filterMetMNXIDsViaRxns  Remove met MNXID associations based on their rxns.
%
% filterMetMNXIDsViaRxns determines which met MNXIDs (if any) should be 
% removed for mets associated with multiple MNXIDs. This is accomplished by 
% finding all model rxns that involve the met, obtaining their associated 
% MNX rxn IDs, and retrieving those rxns from the MNX database. Any MNXIDs 
% associated with the met that are not included in that set of rxns from 
% the MNX database will be removed.
% 
%
% USAGE:
%
%   [fmodel,removed] = filterMetMNXIDsViaRxns(model);
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
%   ignoreComp   (Optional, Default FALSE) If TRUE, metabolite compartments
%                will be ignored. In this case, identical mets of different
%                compartments will be lumped together, and when searching
%                for rxns involving the metabolite, the compartment will be
%                ignored.
%
% OUTPUT:
%
%   fmodel   A filtered model, which has removed metabolite MNXIDs that
%            were not found in any of the MNX reactions associated with the
%            model and that metabolite.
%
%   removed  A structure containing more detailed information on the met
%            MNXID associations that were removed.
%
%
% Jonathan Robinson 2018-05-26


% handle input arguments
if nargin < 2 || isempty(mnx)
    mnx = buildMNXmodel('rxn');
end
if nargin < 3
    ignoreComp = false;
end

% initialize "removed" structure and "fmodel" outputs
fmodel = model;
removed.mets = {};
removed.metMNXID = {};

% if metMNXID and/or rxnMNXID field contains multiple columns, consolidate 
% into a single column with nested cell entries
if size(model.metMNXID,2) > 1
    model.metMNXID = nestCell(model.metMNXID,true);
end
if size(model.rxnMNXID,2) > 1
    model.rxnMNXID = nestCell(model.rxnMNXID,true);
end

if ( ignoreComp )
    
    % check if model has already merged met compartments, or if "mets"
    % contains non-unique elements
    if length(unique(model.mets)) ~= length(model.mets)
        error('Input model "mets" field must contain unique (non-repeated) elements.');
    elseif any(regexp(model.mets{1},'\d$'))
        % Note: this test is specific to models whose met IDs end in
        % numbers (e.g., "m00001"), with compartments appended to the end
        % of the ID (e.g., "m00001c" or "m00001[c]").
        fprintf('\nIt appears that the model compartments have already been merged.\n');
        fprintf('The "ignoreComp" flag will be ignored, since it will have no effect.\n');
        S = model.S;
        ignoreComp = false;
    else
        
        % strip compartment label from model met ID
        if endsWith(model.mets{1},']')
            % compartment name is formatted as "m00001[c]"
            model.mets = regexprep(model.mets,'\[\w\]$','');
        else
            % compartment name is formatted as "m00001c"
            model.mets = regexprep(model.mets,'.$','');
        end
        
        % check if ignoring compartments will actually do anything
        if length(unique(model.mets)) == length(model.mets)
            fprintf('\nIt appears that the model compartments have already been merged.\n');
            fprintf('The "ignoreComp" flag will be ignored, since it will have no effect.\n');
            S = model.S;
            ignoreComp = false;
        else
            
            fprintf('\nMetabolite compartments will be ignored. NOTE: This process will merge metMNXIDs\n');
            fprintf('of mets that are identical except for their compartment. If any mets are associated\n');
            fprintf('with compartment-specific metMNXIDs, this is not recommended.\n');
            
            % If ignoring compartments, convert the stoich matrix into a 
            % binary met-rxn association matrix (i.e., set all nonzero 
            % entries = 1), and combine all associations for metabolites 
            % that are identical except for their compartment. Also combine
            % their metMNXIDs.
            S = (model.S ~= 0);
            [~,uniq_ind,met_groups] = unique(model.mets);
            h = waitbar(0,'Merging model metabolites across compartments...');
            for i = 1:max(met_groups)
                ind = (met_groups == i);
                S(ind,:) = repmat(any(S(ind,:),1),sum(ind),1);
                model.metMNXID(ind) = repmat({unique(horzcat(model.metMNXID{ind}))},sum(ind),1);
                waitbar(i/max(met_groups),h);
            end
            close(h);
            
            % now merge mets of different compartments into single met
            model.mets_nocomp = model.mets;  % first save these for indexing later
            model.mets = model.mets(uniq_ind);
            S = S(uniq_ind,:);
            model.metMNXID = model.metMNXID(uniq_ind);
            
        end 
        
    end

else
    S = model.S;
end

% Identify all mets with two or more MNXID associations. Mets with one or
% zero MNXIDs will be ignored.
multi_ind = find(cellfun(@numel,model.metMNXID) > 1);

% Iterate through mets with multiple MNXIDs, and determine which MNXIDs
% should be removed.
h = waitbar(0,'Processing metabolites...');
for i = 1:length(multi_ind)
    
    m = multi_ind(i);
    
    % get indices of all rxns in which the met participates
    rxn_ind = S(m,:) ~= 0;
    
    % obtain a unique list of all the rxn MNXIDs associated with those
    % rxns, as well as their indices in the MNX database structure
    MNX_rxnIDs = unique(horzcat(model.rxnMNXID{rxn_ind}));
    MNX_rxn_ind = ismember(mnx.rxns,MNX_rxnIDs);
    
    if isempty(MNX_rxnIDs)
        % If none of the rxns are associated to any rxn MNXIDs, just skip
        % this metabolite. Otherwise, all met MNXIDs associated with the
        % metabolite will be removed, which is not so helpful.
        continue
    end
    
    % obtain unique set of all mets participating in the MNX rxns
    MNX_metIDs = unique([mnx.rxnMets{MNX_rxn_ind}]);
    
    % determine which of the MNXIDs currently associated to the model met
    % are not found in any of the MNX rxns
    rem_mets = ~ismember(model.metMNXID{m},MNX_metIDs);
    
    if ~any(rem_mets)
        % if no met MNXID associations should be removed, continue to next met
        continue
    elseif all(rem_mets)
        % This shouldn't happen, especially if "compareRxnMNXIDsWithMets" 
        % is run first, and all the proper metMNXIDs have been added to the
        % model. However, it is still possible, and if it does occur, we
        % do not want to remove all MNXIDs associated to the met, so NONE
        % will be removed. However, the information will be reported in the
        % "removed" structure, so that these cases can be manually
        % inspected by the user.
        removed.mets = [removed.mets; model.mets(m)];
        removed.metMNXID = [removed.metMNXID; {'ALL SHOULD BE REMOVED!'}];
    else
        % remove one or more (but not all) MNXID associations from the met
        removed.mets = [removed.mets; model.mets(m)];
        removed.metMNXID = [removed.metMNXID; {model.metMNXID{m}(rem_mets)}];
        model.metMNXID{m}(rem_mets) = [];
    end
    
    waitbar(i/length(multi_ind),h);
end
close(h);

if any(cellfun(@(x) ismember({'ALL SHOULD BE REMOVED!'},x),removed.metMNXID))
    fprintf('\n*** NOTE: There was at least one case where ALL the MNXIDs associated\n');
    fprintf('          with a metabolite were not found in any of the involved rxns,\n');
    fprintf('          and therefore NONE were removed. See "removed" structure for\n');
    fprintf('          further information.\n\n');
end


% update filtered model (fmodel) metMNXID field
if ( ignoreComp )
    % if compartment was ignored (and thus mets were merged), distribute
    % the updated metMNXID assignments to all compartment-versions of each
    % metabolite.
    for i = 1:length(model.mets)
        ind = ismember(model.mets_nocomp,model.mets(i));
        fmodel.metMNXID(ind) = model.metMNXID(i);
    end
else
    fmodel.metMNXID = model.metMNXID;
end





