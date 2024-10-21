function [eGenes, INIT_output] = estimateEssentialGenes(model, dataFile, taskStruct, useGeneSymbol)
% generate tINIT models and estimate essential genes
%
% Input:
%
%   model          reference human or animal model
%
%   dataFile       (opt, default Hart2015_RNAseq.txt)
%
%   taskStruct     metabolic task structure (opt, default is Essential tasks)
%
%   useGeneSymbol  use gene symbols as ids and in grRules (opt, default TRUE)
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
%   INIT_output     structure containing tINIT models (or information
%                   necessary to regenerate tINIT models) for which gene
%                   essentiality will be evaluated.
%
% Note: This function may take long computation time.
%


if nargin < 2
    % load Hart et al. RNA-Seq cell line data
    dataFile = 'Hart2015_RNAseq.txt';
end

if nargin < 3
    taskStruct = parseTaskList('metabolicTasks_Essential.txt');
end

if nargin < 4
    useGeneSymbol = true;
end

% replace gene IDs with gene symbols
if useGeneSymbol
    idMapping = [model.genes, model.geneShortNames];
    [grRules,genes,rxnGeneMat] = replaceGrRules(model.grRules,idMapping);
    model.grRules = grRules;
    model.genes = genes;
    model.rxnGeneMat = rxnGeneMat;
end

% pre-process RNA-Seq data
disp('Step 1: preprocess and preliminary step')
tmp = readtable(dataFile);
arrayData.genes = tmp.gene;
arrayData.tissues = tmp.Properties.VariableNames(2:end)';
arrayData.levels = table2array(tmp(:,2:end));


% Run some preliminary steps
[~,deletedDeadEndRxns] = simplifyModel(model,true,false,true,true,true);
cModel = removeReactions(model,deletedDeadEndRxns,false,true);
[taskReport, essentialRxnMat] = checkTasks(cModel,[],true,false,true,taskStruct);


% add pre-processing results to arrayData structure
arrayData.deletedDeadEndRxns = deletedDeadEndRxns;
arrayData.taskReport = taskReport;
arrayData.essentialRxnMat = essentialRxnMat;
arrayData.threshold = 1;
        

% run tINIT 
disp('Step 2: get tissue models')
model = addBoundaryMets(model);
params = {};
INIT_output = {};
    
for i = 1:length(arrayData.tissues)
    disp(['Tissue ', num2str(i), ' out of ',  num2str(length(arrayData.tissues)),': ', arrayData.tissues{i}])
        
    % First try to run tINIT with shorter time limit. If it fails, then
    % try again with a longer time limit.
    try
        params.TimeLimit = 1000;
        init_model = getINITModel2(model,arrayData.tissues{i},[],[],arrayData,[],true,[],true,true,taskStruct,params);
    catch
        params.TimeLimit = 5000;
        init_model = getINITModel2(model,arrayData.tissues{i},[],[],arrayData,[],true,[],true,true,taskStruct,params);
    end
        
    init_model.id = arrayData.tissues{i};
    INIT_output.id{i,1} = init_model.id;
    INIT_output.model{i,1} = init_model;
            
end
    
disp('Step 3: get essential genes')
% get essential genes for each model and task
eGenes = getTaskEssentialGenes(INIT_output, model, taskStruct);
eGenes.refModel = model;
    
end
