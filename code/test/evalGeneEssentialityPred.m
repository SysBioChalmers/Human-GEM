function [metrics, modelEssential] = evalGeneEssentialityPred(model, expData, metabTasks, modelEssential)
% Compare GEM-based gene essentiality predictions with experimental data
%
% Inputs:
%
%   model       Genome-scale metabolic model structure.
%
%               NOTE: The model should include boundary metabolites! These
%               can be added using the addBoundaryMets function.
%
%   expData     A two-column cell array, where the first column contains
%               gene names or IDs (of the same type as those used in the
%               model), and the second column indicates whether the gene
%               was found as essential or non-essential.
%
%               The second column can be text (e.g., 'essential', 'non') or
%               numeric (0 = non-essential, 1 = essential). For example:
%
%               {'gene1'   'essential'  }       {'gene1'   [1]}
%               {'gene2'   'non'        }       {'gene2'   [0]}
%               {'gene3'   'conditional'}       {'gene3'   [1]}
%               {...       ...          }       {...       ...}
%
%               NOTE: If essentiality is provided as strings, only genes
%               starting with "non" will be treated as non-essential. This
%               means that genes annotated as "essential" or "conditional"
%               will both be treated as essential.
%
%               NOTE: ALL genes that were tested for essentiality should be
%               included in expData. If expData contains only essential
%               genes, then it will be assumed that all genes in the genome
%               (and thus all genes in the model) were tested.
%
%   metabTasks  Either the filename or structure containing the metabolic
%               tasks that will be tested when using the model to predict
%               gene essentiality. The task structure can be loaded using
%               the RAVEN "parseTaskList" function.
%
%               Genes will be predicted as essential if their deletion
%               impairs ANY of the metabolic tasks in metabTasks. If the
%               model can perform ALL metabolic tasks upon the deletion of
%               a gene, that gene will be predicted as non-essential.
%
%   modelEssential  (Optional) A list of model-predicted essential genes.
%                   If provided, the model gene essentiality prediction
%                   will be skipped and the provided list used instead.
%
% Outputs:
%
%   metrics     A result structure with the following fields:
%                   sensitivity
%                   specificity
%                   accuracy
%                   F1 statistic
%                   Matthew's Correlation Coefficient (MCC)
%                   p-value associated with a hypergeometric test of the
%                           true and false positives and negatives
%                           (2x2 contingency table).
%
%   modelEssential   List of model-predicted essential genes.
%

if nargin < 4
    modelEssential = [];
end


%% Pre-process expData

% extract information from expData and convert essentiality to 0, 1
ex.genes = expData(:,1);
if length(ex.genes) > length(unique(ex.genes))
    error('Gene names or IDs in expData are not unique. Duplicated entries must be removed.');
end
if isnumeric(expData{1,2})
    ex.essentiality = double(cell2mat(expData(:,2)));
else
    ex.essentiality = double(~startsWith(lower(expData(:,2)), 'non'));
end

% check that gene IDs/names in the model are the same type as in expData
if sum(ismember(ex.genes, model.genes)) < 3
    error('The gene name or ID type in expData seem to differ from those used in the model.');
end

% define essential and non-essential gene lists
ex.essential = ex.genes(ex.essentiality == 1);
if all(ex.essentiality == 1)
    fprintf('NOTE: All genes in expData are marked as essential; it will therefore be assumed that all genes in the genome were tested.\n\n');
    ex.nonessential = setdiff(model.genes, ex.essential);
    ex.genes = [ex.genes; ex.nonessential];
else
    ex.nonessential = ex.genes(ex.essentiality == 0);
end


%% Run gene essentiality predictions with the model

if isempty(modelEssential)
    
    if ischar(metabTasks)
        taskStruct = parseTaskList(metabTasks);
    else
        taskStruct = metabTasks;
    end
    
    % first confirm that the model can perform all the tasks
    taskReport = checkTasks(model,[],false,false,false,taskStruct);
    if ~all(taskReport.ok)
        fprintf('\nThe provided model could not perform the following tasks:\n');
        fprintf('\t> %s\n', taskReport.description{~taskReport.ok});
        fprintf('\n');
        error('Model failed task(s) before gene deletion!');
    end
    
    % get essential genes for each task in taskStruct
    [~,essentialGeneMat] = checkTasksGenes(model,[],false,false,true,taskStruct);
    
    % essential genes are counted as those that are essential for ANY task
    essentialGeneVect = any(essentialGeneMat > 0, 2);
    pred.essential = model.genes(essentialGeneVect);
    pred.nonessential = setdiff(model.genes, pred.essential);
    
else
    pred.essential = modelEssential;
    pred.nonessential = setdiff(model.genes, modelEssential);
end


%% Evaluate prediction performance

TP = sum(ismember(pred.essential, ex.essential));        % true positives
TN = sum(ismember(pred.nonessential, ex.nonessential));  % true negatives
FP = sum(ismember(pred.essential, ex.nonessential));     % false positives
FN = sum(ismember(pred.nonessential, ex.essential));     % false negatives
  
% calculate some metrics
sensitivity = TP./(TP + FN);
specificity = TN./(TN + FP);
accuracy = (TP + TN)./(TP + TN + FP + FN);
F1 = 2*TP./(2*TP + FP + FN);
MCC = ((TP.*TN) - (FP.*FN))./sqrt((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN));
[~, p_hyper] = fishertest([TP, FP; FN, TN], 'tail', 'right');

% combine metrics into output structure
metrics.sensitivity = sensitivity;
metrics.specificity = specificity;
metrics.accuracy = accuracy;
metrics.F1 = F1;
metrics.MCC = MCC;
metrics.p_hyper = p_hyper;

% assign output
modelEssential = pred.essential;


end


