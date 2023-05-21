function eGenes = getTaskEssentialGenesCluster(models_file_start, chunk)
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


% specify paths for cluster use
addpath(genpath('../components/RAVEN'));
addpath(genpath('../components/COBRA'));
addpath(genpath('../components/Human-GEM'));

cd '../components/Human-GEM/code/DepMapGeneEss'

eGenes = getTaskEssentialGenes(models_file_start, chunk)

