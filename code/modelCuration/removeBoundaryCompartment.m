%
% FILE NAME:    removeBoundaryCompartment.m
% 
% DATE CREATED: 2020-05-21
%      UPDATED: 2020-05-21
% 
% PROGRAMMER:   Jonathan Robinson
%               National Bioinformatics Infrastructure Sweden
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: This script removes the boundary compartment [x] and all
%          metabolites within this compartment from the model.
%

% load model
ihuman = importHumanYaml('../../modelFiles/yml/HumanGEM.yml');

% delete unconstrained (boundary) metabolites
model_del = simplifyModel(ihuman);
fprintf('\nRemoved %u boundary metabolites from the model.\n\n', numel(ihuman.mets) - numel(model_del.mets));

% remove the boundary compartment itself and update the metComps field
metCompSymbols = model_del.comps(model_del.metComps);
[~,comp_ind] = ismember('boundary', lower(model_del.compNames));
model_del.comps(comp_ind) = [];
model_del.compNames(comp_ind) = [];
[~,model_del.metComps] = ismember(metCompSymbols, model_del.comps);


% write new model to yml
writeHumanYaml(model_del, '../../modelFiles/yml/HumanGEM.yml');










