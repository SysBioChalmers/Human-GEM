function [outModel, modelChanges]=gapfill4Biomass(model,refModel)
% gapfill4Biomass
%   Fill gaps for the input model using "fitTasks" and based on the reference
%   model so that the resulting model can grow on Ham's media
%
%   Input:
%   model           model structure
%   refModel        reference model from which to retrieve reactions
%
%
%   Output:
%   outModel        model structure with additional reactions necessary for
%                   biomass formation
%   modelChanges    A cell structure containing the original and new reaction
%                   and/or metabolite properties. Reaction-related changes,
%                   if found, will be reported in the modelChanges.rxns
%                   field, whereas metabolite-related changes are reported
%                   in the modelChanges.mets field
%                   
%
%   Usage: [outModel, modelChanges]=gapfill4Biomass(model,refModel)
%
%

% handle input arguments
if nargin < 2
    error('No enough input.');
end

% confirm the relations of input model and reference model - share 50% rxns
if length(intersect(model.rxns, refModel.rxns)) < 0.5*length(model.rxns)
   error('Please check if a correct reference model is used, its relation to the input model is uncertain.');
end

% save the initial model
model_orig = model;

% add boundary mets
model = addBoundaryMets(model);
refModel = addBoundaryMets(refModel);


%% gap-filling

% load metabolic task for growth under Ham's media
[ST, I] = dbstack('-completenames');
path = fileparts(ST(I).file);
essentialTasks = fullfile(path,'../data/metabolicTasks','metabolicTasks_Essential.xlsx');
taskStruct = parseTaskList(essentialTasks);
%taskStruct = taskStruct(end);

% gap-filling with fitTask
fprintf('Start gap-filling for essential tasks...\n');
[outModel, addedRxns]=fitTasks(model,refModel,[],true,[],taskStruct);

% confirm the growth
taskReport = checkTasks(outModel,[],false,false,false,taskStruct);
if ~all(taskReport.ok)
    error('Gap-filling process failed in enabling the biomass object function to carry flux.');
end

% clean/remove fields introducted by above steps (e.g. fitTasks)
fieldsToRemove = {'rxnFrom', 'metFrom', 'geneFrom'};
outModel = rmfield(outModel, intersect(fieldsToRemove,fieldnames(outModel)));
if ismember('id',fieldnames(outModel))
    outModel.id = '';
end


%% empty newly introduced grRules

% clean the grRules associated with gap-filling reactions
gapfilledRxns = setdiff(outModel.rxns, model.rxns);
[~, ind] = ismember(gapfilledRxns, outModel.rxns);
outModel.grRules(ind) = {''};

% re-generate gene and rxnGeneMat fields
[outModel.genes, outModel.rxnGeneMat] = getGenesFromGrRules(outModel.grRules);


%% prepare output

% remove boundary mets
outModel = simplifyModel(outModel);

% record model changes
modelChanges = docModelChanges(model_orig, outModel);


