function results = compareRxnMNXIDsWithMets(model,mnx,ignoreComp)
%compareRxnMNXIDsWithMets  Check consistency between rxn and met MNXIDs.
%
% compareRxnMNXIDsWithMets determines if the model rxn MNXID associations
% have distributed their information to the model met MNXID associations.
% For each metabolite in the model, the set of reactions in which it
% participates is identified. The rxnMNXIDs assigned to this set of
% reactions is then obtained, and the corresponding reaction information is
% retreieved from the MNX database. If there is no overlap between the set
% of metabolite MNX IDs participating in the reaction, and the set of 
% MNXIDs assigned to the metabolite in the model, the reaction and 
% metabolite will be flagged and reported in the results structure.
%
% USAGE:
%
%   results = compareRxnMNXIDsWithMets(model,mnx,ignoreComp);
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
%   results  A results structure containing information on the flagged
%            reactions and the associated mets and MNXIDs, organized as a
%            cell array, containing the following column headers:
%               'model met'
%               'metMNXIDs assigned to met'
%               'flagged model rxn'
%               'flagged rxnMNXIDs mapped to rxn'
%


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
        else
            
            fprintf('\nMetabolite compartments will be ignored.\n');
            
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

% initialize results structure
results = {'model met','metMNXIDs assigned to met','flagged model rxn','flagged rxnMNXIDs mapped to rxn'};

% iterate through each of the model metabolites
h = waitbar(0,'Processing metabolites...');
for i = 1:length(model.mets)
    
    % identify all reactions in which the metabolite is involved
    rxn_ind = find(S(i,:) ~= 0);
    
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
        
        % add flagged MNX rxn and associated model met and rxn info to the
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


