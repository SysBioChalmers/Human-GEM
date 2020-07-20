%
% FILE NAME:    removeBoundaryCompartment.m
%
% PURPOSE: This script removes the boundary compartment [x] and all
%          metabolites within this compartment from the model, as discussed
%          in #172

% load model
ihuman = importHumanYaml('../../model/Human-GEM.yml');

% delete unconstrained (boundary) metabolites
model_del = simplifyModel(ihuman);
del_mets = setdiff(ihuman.mets, model_del.mets);
fprintf('\nRemoved %u boundary metabolites from the model.\n\n', numel(del_mets));

% remove the boundary compartment itself and update the metComps field
metCompSymbols = model_del.comps(model_del.metComps);
[~,comp_ind] = ismember('boundary', lower(model_del.compNames));
model_del.comps(comp_ind) = [];
model_del.compNames(comp_ind) = [];
[~,model_del.metComps] = ismember(metCompSymbols, model_del.comps);


% write new model to yml
writeHumanYaml(model_del, '../../model/Human-GEM.yml');


% import metabolite annotations
metAssoc = jsondecode(fileread('../../data/annotation/humanGEMMetAssoc.JSON'));

% remove boundary metabolites
del_ind = ismember(metAssoc.mets, del_mets);
f = fieldnames(metAssoc);
for i = 1:numel(f)
    metAssoc.(f{i})(del_ind) = [];
end

% verify that annotation structure is aligned with model
if ~isequal(model_del.mets, metAssoc.mets)
    error('Model and metabolite annotation structures not synced!');
end

% export metabolite annotations
jsonStr = jsonencode(metAssoc);
fid = fopen('../../data/annotation/humanGEMMetAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);





