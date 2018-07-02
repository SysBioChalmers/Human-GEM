function model_new = addMetCompsField(model)
%addMetCompsField  Add "metComps" field to a model based on "mets" field.
%
% USAGE:
%
%   model_new = addMetCompsField(model);
%
% INPUTS:
%
%   model   A genome-scale metabolic model. Must contain a "mets" field,
%           which should specify metabolite compartment information at the
%           end of each entry (e.g., h2o_c or h2o[c]).
%
% OUTPUTS:
%
%   model_new   The input model with a new "metComps" field added, which is
%               numeric vector indicating the compartment index of each
%               metabolite. If the input model does not contain a "comps"
%               field, this field will also be generated and added to the
%               output model_new.
%
%
% Jonathan L. Robinson, 2018-07-02


% determine format of met compartment abbreviations
if endsWith(model.mets{1},']')
    % compartment abbrev is contained within brackets at end of met ID
    metCompAbbrevs = regexp(model.mets,'\[(\w)\]$','tokens');
    metCompAbbrevs = cellfun(@(a) a{1}{1},metCompAbbrevs,'UniformOutput',false);
else
    % assume compartment abbrev is last character of met ID
    metCompAbbrevs = regexp(model.mets,'\w$','match');
    metCompAbbrevs = cellfun(@(a) a{1},metCompAbbrevs,'UniformOutput',false);
end

% check "comps" field
if ~isfield(model,'comps')
    % if field does not exist, create one
    model.comps = unique(metCompAbbrevs);
elseif any(~ismember(metCompAbbrevs,model.comps))
    error('One or more compartments taken from met ID do not occur in "comps".');
end

% generate metComps field
[~,model.metComps] = ismember(metCompAbbrevs,model.comps);

% assign output
model_new = model;


