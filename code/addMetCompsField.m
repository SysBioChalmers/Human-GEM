function model_new = addMetCompsField(model)
%addMetCompsField
%   Add "metComps" field to a model based on "mets" field.
%
%   model       A genome-scale metabolic model. Must contain a "mets" field,
%               which should specify metabolite compartment information at
%               the end of each element (e.g., h2o_c or h2o[c])
%
%   model_new   The output model with a new "metComps" field, which is
%               numeric vector indicating the compartment index of each
%               metabolite. If the input model does not contain a "comps"
%               field, this field will also be generated and added to the
%               output model_new
%
%   Usage: model_new = addMetCompsField(model);
%


metCompAbbrevs=cell(numel(model.mets),1);
metCompAbbrevs(:)={''};

% deal with the compartment abbrev one-by-one
for i=1:numel(model.mets)
    % determine format of met compartment abbreviations
    if endsWith(model.mets{i},']')
        % compartment abbrev is contained within brackets at end of met ID
        tmp = regexp(model.mets{i},'\[(\w+)\]$','tokens');
        metCompAbbrevs{i} = tmp{1};
    elseif contains(model.mets{1},'_')
        % compartment abbrev is at the end of metID and separated by underscore
        metCompAbbrevs{i} = regexp(model.mets{i},'\_(\w+)$','tokens');
    else
        % compartment abbrev is last character(s) of met ID without separator
        metCompAbbrevs{i} = regexp(model.mets{i},'[a-z]+$','match');
    end
end
metCompAbbrevs=cellfun(@char,metCompAbbrevs,'UniformOutput',false);

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

end
