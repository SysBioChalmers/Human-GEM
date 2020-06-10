function new_model = addBoundaryMets(model, exch_only)
%addBoundaryMets  Add boundary metabolites to a model structure.
%
%   Boundary metabolites are pseudometabolites that exist at the system
%   boundary, and are used for the import/export of mass into/out of the
%   system. 
%   For example: "glc[extracellular] <==> glc[boundary]"
%
%   Some models do not use boundary metabolites, but instead formulate
%   these exchange reactions (or "demand/DM" or "sink" reactions) without
%   an explicit product (or reactant).
%   For example: "glc[extracellular] <==> "
%
%   This function identifies all reactions that involve the conversion of 
%   a single metabolite into nothing (or vice versa), and balances the
%   reaction with the same metabolite in a different (boundary) comparment.
%   The function will also add any new boundary metabolites to the model if
%   they did not yet exist.
%
% USAGE:
%
%   new_model = addBoundaryMets(model, exch_only);
%
%
% INPUTS:
%
%   model       Model structure.
%
%   exch_only   (Optional, default = FALSE) 
%               If TRUE, only exchange rxns (i.e., into/out of the
%               extracellular compartment) will be considered.
%               If FALSE, reactions involving mets in other compartments
%               will also be considered. For example:
%
%               "met[nucleus] -->" becomes "met[nucleus] --> met[boundary]"
%
%
% OUTPUTS:
%
%   new_model   New model structure, with added boundary metabolites.
%

% handle input arguments
if nargin < 2
    exch_only = false;
end

% add a boundary compartment to the model if it does not yet exist
if ~ismember('boundary',lower(model.compNames))
    model.compNames(end+1) = {'Boundary'};
    if ismember('x',model.comps)
        error('Consider renaming the "x" compartment, or modify this function to use another letter.');
    else
        model.comps(end+1) = {'x'};
    end
end

% get index of boundary and extracellular compartments
[~,bound_comp_ind] = ismember('boundary',lower(model.compNames));
[~,xcell_comp_ind] = ismember('extracellular',lower(model.compNames));
if (xcell_comp_ind == 0) && (exch_only)
    error('No compartments named "Extracellular" were found.');
end

% add an "unconstrained" field to the model if it does not yet exist
if ~isfield(model,'unconstrained')
    model.unconstrained = double(model.metComps == bound_comp_ind);
elseif ~all(model.unconstrained == double(model.metComps == bound_comp_ind))
    warning('Ensure that the "unconstrained" field is properly updated.');
end

% find all reactions that involve only a single metabolite
if (exch_only)
    % if only considering exchange (extracellular) rxns, ensure that the
    % single metabolite is in the extracellular compartment
    xcell_met_ind = (model.metComps == xcell_comp_ind);
    unbal_rxn_inds = find((sum(model.S ~= 0) == 1) & (sum(model.S(xcell_met_ind,:) ~= 0) == 1))';
else
    unbal_rxn_inds = find(sum(model.S ~= 0) == 1)';
end

% obtain the names of the mets participating in these rxns
[unbal_met_inds,~] = find(model.S(:,unbal_rxn_inds) ~= 0);
unbal_mets = unique(model.metNames(unbal_met_inds));

% get list of existing boundary metabolite names
bound_mets = unique(model.metNames(model.metComps == bound_comp_ind));

% identify new boundary metabolites that must be added to the model
add_bound_mets = unbal_mets(~ismember(unbal_mets,bound_mets));

% Need to construct the associated metID for each of these new mets, which
% will simply be the same as the ID for the non-boundary version of the
% metabolite, with the compartment replaced.
[~,ref_met_ind] = ismember(add_bound_mets,model.metNames);
if all(endsWith(model.mets,']'))
    add_bound_met_IDs = regexprep(model.mets(ref_met_ind),'\[.\]$','[x]');
elseif any(endsWith(model.mets,']'))
    error('All metabolite IDs must use the same format for describing the compartment.');
else
    add_bound_met_IDs = regexprep(model.mets(ref_met_ind),'.$','x');
end

% add new boundary mets to the model
metsToAdd.mets = add_bound_met_IDs;
metsToAdd.metNames = add_bound_mets;
metsToAdd.compartments = 'x';
metsToAdd.unconstrained = ones(size(add_bound_mets));
new_model = addMets(model,metsToAdd);

% now add the boundary mets to the model S-matrix
S = new_model.S;
for i = 1:length(unbal_rxn_inds)
    % find the boundary met corresponding to the current reaction
    met_ind = ismember(new_model.metNames,new_model.metNames(unbal_met_inds(i))) & (new_model.metComps == bound_comp_ind);
    
    % balance the reaction by giving the newly added boundary met a stoich
    % coeff that is the negative of the original sole metabolite
    new_model.S(met_ind,unbal_rxn_inds(i)) = -S(unbal_met_inds(i),unbal_rxn_inds(i));
end

% print some results
fprintf('\nBoundary metabolites were added to %u reactions.\n',length(unbal_rxn_inds));
fprintf('New (boundary) versions of %u metabolites were added to the model.\n\n',length(metsToAdd.mets));







