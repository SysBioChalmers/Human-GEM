function eGenes = getTaskEssentialGenes(models_file_start, chunk)
% Identify genes essential for different tasks in different models.
% This code is largely copied from the Human1 paper
%
% Input:
%
%   models_file_start  Start of path of .mat file containing a cell array of tINIT models.
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

models_file = [models_file_start num2str(chunk) '.mat' ];
out_filename = [models_file_start 'eGenes-' num2str(chunk) '.mat'];

% verify that the cluster node has sufficient CPUs available
x = parcluster('local');
if x.NumWorkers < 8
    error('Insufficient CPUs available (%u) for parallel processing!', x.NumWorkers); %for some reason, Vera often gives 8 cores, which is enough.
end

%for debugging:
%models_file = '../data/init_models_depmap15.mat';

% load models and metabolic tasks
if isnumeric(models_file)
    models_file = strcat('init_models-', num2str(models_file), '.mat');
end
x = load(models_file);
fns = fieldnames(x);
models = getfield(x,fns{1});

% specify paths for cluster use
addpath(genpath('../components/RAVEN'));
addpath(genpath('../components/COBRA'));
addpath(genpath('../components/Human-GEM'));

%setRavenSolver('gurobi');

taskStruct = parseTaskList('metabolicTasks_Essential.txt');


% initialize variables
n = numel(models);
tissues = cell(n,1);
geneList = cell(n,1);
allEssentials = cell(n,1);

% iterate through models in parallel for-loop
%parpool(8);  % assume 8 CPUs
parfor i = 1:length(models)
    % determine essential genes for each task
    [~,essentialGenes] = checkTasksGenes(closeModel(models{i}),[],false,false,true,taskStruct);
    
    % collect results
    tissues{i,1} = models{i}.id;
    geneList{i,1} = models{i}.genes;
    allEssentials{i,1} = essentialGenes;
    
end

% gather results into eGenes structure
eGenes = {};
eGenes.taskList = {taskStruct(:).description}';
eGenes.tissues = tissues;
eGenes.geneList = geneList;
eGenes.essentialGenes = allEssentials;

% save eGenes structure
save(out_filename, 'eGenes');

