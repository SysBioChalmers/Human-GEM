function new_model = addBoundaryMets(model)
%addBoundaryMets  Add boundary metabolites to a model structure.
%
%   Boundary metabolites are pseudometabolites that exist at the system
%   boundary, and are used for the import/export of mass into/out of the
%   system. 
%   For example: "glc(boundary) <==> glc(extracellular)"
%
%   Some model formats (e.g., Cobra) do not use boundary metabolites, but
%   instead formulate these exchange reactions (or "demand/DM" or "sink"
%   reactions) without an explicit reactant. 
%   For example: "glc(extracellular) <==> "
%
%   This function identifies all reactions that are formulated with an 
%   exchange of an extracellular metabolite into nothing as a reaction with
%   a boundary metabolite, and adds any new boundary metabolites to the 
%   model, if they did not yet exist.
%
% USAGE:
%
%   new_model = addBoundaryMets(model);
%
% INPUTS:
%
%   model    Model structure.
%
% OUTPUTS:
%
%   new_model   New model structure, with added boundary metabolites.
%
%
% Jonathan Robinson, 2018-09-12


% add a boundary compartment to the model if it does not yet exist
if ~ismember('boundary',lower(model.compNames))
    model.compNames(end+1) = {'Boundary'};
    if ismember('x',model.comps)
        error('Consider renaming the "x" compartment, or modify this function to use another letter.');
    else
        model.comps(end+1) = {'x'};
    end
end

% get index of boundary compartment
[~,bound_comp_ind] = ismember('boundary',lower(model.compNames));

% add an "unconstrained" field to the model if it does not yet exist
if ~isfield(model,'unconstrained')
    model.unconstrained = double(model.metComps == bound_comp_ind);
elseif ~all(model.unconstrained == double(model.metComps == bound_comp_ind))
    warning('Ensure that the "unconstrained" field is properly updated.');
end

% find all reactions that involve only a single metabolite
unbal_rxn = find(sum(model.S ~= 0) == 1)';

% obtain the names of the mets participating in these rxns
[unbal_met_inds,~] = find(model.S(:,unbal_rxn) ~= 0);
unbal_mets = unique(model.metNames(unbal_met_inds));

% get list of existing boundary metabolite names
bound_mets = unique(model.metNames(model.metComps == bound_comp_ind));

% identify new boundary metabolites that must be added to the model
add_bound_mets = unbal_mets(~ismember(unbal_mets,bound_mets));


% add new boundary mets to the model
metsToAdd.metNames = add_bound_mets;
metsToAdd.compartments = 'x';
metsToAdd.unconstrained = ones(size(add_bound_mets));
new_model = addMets(model,metsToAdd);

% now add the boundary mets to the model S-matrix
for i = 1:length(unbal_rxn)
    % find the boundary met corresponding to the current reaction
    met_ind = ismember(new_model.metNames,new_model.metNames(unbal_met_inds(i))) & (new_model.metComps == bound_comp_ind);
    
    % balance the reaction by giving the newly added boundary met a stoich
    % coeff that is the negative of the original sole metabolite
    new_model.S(met_ind,unbal_rxn(i)) = -new_model.S(unbal_met_inds(i),unbal_rxn(i));
end








