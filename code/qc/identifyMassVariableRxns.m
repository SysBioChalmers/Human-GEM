function [massVariableSets,massConsistentSets] = identifyMassVariableRxns(model,minMets,maxMets)
%identifyMassVariableRxns  Find rxn sets with potential mass inconsistency.
%
%   This function searches a model for sets of reactions that involve the
%   same set of metabolites except for one, which indicates that one or
%   more of the reactions in the set is mass-imbalanced. The metabolites
%   that differ among the reactions in a set will be compared in their
%   metFormulas; if they are different, the set will be classified as
%   "mass variable reactions", otherwise the set will be classified as
%   "mass consistent reactions".
%
% USAGE:
%
%   [massVariableSets,massConsistentSets] = identifyMassVariableRxns(model,minMets,maxMets);
%
% INPUT:
%
%   model    Model structure.
%
%   minMets  (Optional, default = 4) the minimum number of metabolites a
%            reaction must contain to be considered in the analysis.
%
%   maxMets  (Optional, default = 100) the maximum number of metabolites a
%            reaction can contain to be considered in the analysis.
%
% OUTPUT:
%
%   massVariableSets    Sets of reactions that differ by one metabolite,
%                       where the differing metabolites do NOT all have the
%                       same metFormula. The massVariableSets output is
%                       formatted as an R x N logical matrix, where R is
%                       the number of reactions in the model, and N is the
%                       number of sets idenfied. Each column represents a
%                       reaction set, where the members of that set are
%                       indicated by a 1 ("true").
%
%   massCosistentSets   Sets of reactions that differ by one metabolite,
%                       where the differing metabolites all have identical
%                       metFormulas. The massCosistentSets output is
%                       formatted as an R x N logical matrix, where R is
%                       the number of reactions in the model, and N is the
%                       number of sets idenfied. Each column represents a
%                       reaction set, where the members of that set are
%                       indicated by a 1 ("true").
%


% handle input arguments
if nargin < 3
    maxMets = 100;
end
if nargin < 2 || isempty(minMets)
    minMets = 4;
end

% initialize variables
massVariable_rxn_sets = [];
massConsistent_rxn_sets = [];

% iterate through each number of metabolites that the reactions can contain
for n = minMets:maxMets
    
    % find all reactions involving "n" metabolites
    rxn_ind = sum(model.S ~= 0,1)' == n;
    
    % find all mets involved in these reactions
    met_ind = any(model.S(:,rxn_ind) ~= 0,2);
    
    % extract subset of stoich matrix for these reactions and metabolites
    s = model.S(met_ind,rxn_ind);
    
    % keep track of original reaction indices
    orig_rxn_ind = (1:length(model.rxns))';
    orig_rxn_ind = orig_rxn_ind(rxn_ind);
    
    % calculated the hamming distance between these reactions
    rxn_dist = squareform(pdist(s','hamming'));
    
    % we are interested in reactions that differ by only one metabolite, which
    % would correspond to a hamming distance of 2/Nmets
    d = 2/sum(met_ind);
    
    % find all cases where at least 2 other rxns differ by this distance
    check_rxns = sum(rxn_dist == d) >= 2;
    
    % obtain unique list of similar reaction sets (do this by retrieving rxns
    % with the specified hamming distance, as well as a distance of zero, so
    % that the reaction set includes the checked reaction itself).
    rxn_sets = unique(ismember(rxn_dist(check_rxns,:),[0,d]),'rows');
    
    
    % remove any sets that are just a subset of another
    i = 0;
    while ~isempty(rxn_sets)
        i = i+1;
        r = rxn_sets - rxn_sets(i,:);
        if sum(all(r >= 0,2)) > 1
            rxn_sets(i,:) = [];
            i = 0;
        elseif i == size(rxn_sets,1)
            break
        end
    end
        
    % convert to logical that corresponds to original rxn indexing
    rxn_sets = rxn_sets';  % transpose matrix
    rxn_sets_orig = false(length(model.rxns),size(rxn_sets,2));
    for i = 1:size(rxn_sets,2)
        rxn_sets_orig(orig_rxn_ind(rxn_sets(:,i)),i) = true;
    end
    
    % for each reaction set, determine which differ by metabolites that vary in
    % their formula
    ignore_sets = [];
    for i = 1:size(rxn_sets_orig,2)
        diff_mets = sum(model.S(:,rxn_sets_orig(:,i)) ~= 0, 2) == 1;
        if length(unique(model.metFormulas(diff_mets))) == 1
            ignore_sets = [ignore_sets; i];
        end
    end
    
    % append these reaction sets to the existing sets
    if ~isempty(ignore_sets)
        massConsistent_rxn_sets = [massConsistent_rxn_sets, rxn_sets_orig(:,ignore_sets)];
        rxn_sets_orig(:,ignore_sets) = [];
    end
    massVariable_rxn_sets = [massVariable_rxn_sets, rxn_sets_orig];
    
end

% assign output
massVariableSets = logical(massVariable_rxn_sets);
massConsistentSets = logical(massConsistent_rxn_sets);



