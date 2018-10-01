%
% FILE NAME:    constrainVariableMassReactions.m
% 
% DATE CREATED: 2018-09-28
%     MODIFIED: 2018-09-28
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script to identify and constrain all reactions that involve the
%          same reactants producing different compounds, such that the
%          different compounds vary in mass. For example:
%
%           LCAT39e: cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosahexenoylglycerophosphocholine (Delta 4, 7, 10, 13, 16, 19), Sn1-Lpc (22:6)[s]
%            LCAT5e: cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosadienoylglycerophosphocholine (Delta 11,14)[s]
%           LCAT31e: cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Octadeca-Trienoylglycerophosphocholine, Sn1-Lpc (18:3, Delta 6, 9, 12)[s]
%
%          These reactions have the same reactants (cholesterol and PC-LD
%          pool), but differ in only one product, which results in a PC-LD
%          pool of variable mass (since some of these reactions are
%          reversible, and some of the mets involved can be broken down
%          into their components).


massVariable_rxn_sets = [];
massConsistent_rxn_sets = [];

for n = 4:100
    
    % find all reactions with X metabolites
    rxn_ind = sum(ihuman.S ~= 0,1)' == n;
    
    % find all mets involved in these reactions
    met_ind = any(ihuman.S(:,rxn_ind) ~= 0,2);
    
    % extract subset of stoich matrix for these reactions and metabolites
    s = ihuman.S(met_ind,rxn_ind);
    
    % keep track of original reaction indices
    orig_rxn_ind = (1:length(ihuman.rxns))';
    orig_rxn_ind = orig_rxn_ind(rxn_ind);
    
    % calculated the hamming distance between these reactions
    rxn_dist = squareform(pdist(s','hamming'));
    
    % we are interested in reactions that differ by only one metabolite, which
    % would correspond to a hamming distance of 2/Nmets
    d = 2/sum(met_ind);
    
    % find all cases where at least 2 other rxns differ by this distance
    check_rxns = find(sum(rxn_dist == d) >= 2);
    
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
    rxn_sets_orig = false(length(ihuman.rxns),size(rxn_sets,2));
    for i = 1:size(rxn_sets,2)
        rxn_sets_orig(orig_rxn_ind(rxn_sets(:,i)),i) = true;
    end
    
    % for each reaction set, determine which differ by metabolites that vary in
    % their formula
    ignore_sets = [];
    for i = 1:size(rxn_sets_orig,2)
        diff_mets = sum(ihuman.S(:,rxn_sets_orig(:,i)) ~= 0, 2) == 1;
        if length(unique(ihuman.metFormulas(diff_mets))) == 1
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

massVariable_rxn_sets = logical(massVariable_rxn_sets);
massConsistent_rxn_sets = logical(massConsistent_rxn_sets);










