function eGenes = getTaskEssentialGenes(INIT_output, refModel, taskStruct)
% Identify genes essential for different tasks in different tINIT models.
%
% Input:
%
%   INIT_output     structure containing tINIT models (or information
%                   necessary to regenerate tINIT models) for which gene
%                   essentiality will be evaluated.
%
%   refModel        reference model from which the tINIT models were
%                   generated.
%
%   taskStruct      metabolic task structure
%
%
% Output:
%
%   eGenes          results structure with the following fields:
%       taskList    list of metabolic tasks that were tested
%       tissues     list of tissues (model IDs) corresponding to each model
%       geneList    cell array of the list of genes from each model
%       essentialGenes   cell array with one entry per model, where each
%                        entry is a logical matrix with rows corresponding
%                        to genes (in geneList) and columns to tasks (in
%                        taskList). Entries in the matrix are true when a
%                        gene is essential for a task, and false otherwise.
%
%
% Usage:
%
%   eGenes = getTaskEssentialGenes(INIT_output, refModel, taskStruct);
%
%
% Jonathan Robinson, 2019-10-26



if ~isfield(INIT_output, 'model')
    % if INIT_output does not contain complete models, they first need to
    % be regenerated
    models = regen_tINIT_model(refModel, INIT_output);
else
    models = INIT_output.model(:);
end

% ensure that model ID is set to the tissue names
for i = 1:numel(models)
    if ~isfield(INIT_output,'tissues')
        models{i}.id = INIT_output.id{i};
    else
        models{i}.id = INIT_output.tissues{i};
    end
end

% replace Ensembl IDs with gene symbols
for i = 1:length(models)
    if any(startsWith(models{i}.genes,'ENSG000'))
        idMapping = [models{i}.genes, models{i}.geneShortNames];
        [grRules,genes,rxnGeneMat] = replaceGrRules(models{i}.grRules,idMapping);
        models{i}.grRules = grRules;
        models{i}.genes = genes;
        models{i}.rxnGeneMat = rxnGeneMat;
    end
end

% iterate through the different models
parfor i = 1:length(models)
    
    m = models{i};
    
    % determine essential genes for each task
    [~,essentialGenes] = checkTasksGenes(m,[],false,false,true,taskStruct);
    
    % collect results
    tissues{i,1} = m.id;
    geneList{i,1} = m.genes;
    allEssentials{i,1} = essentialGenes;
    
end

% gather results into eGenes structure
eGenes = {};
eGenes.taskList = {taskStruct(:).description}';
eGenes.tissues = tissues;
eGenes.geneList = geneList;
eGenes.essentialGenes = allEssentials;



