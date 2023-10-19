function results = evaluateHart2015Essentiality(eGenes)
% Evaluate and compare experimental fitness genes with essentiality results
% predicted from 5 cell-line specific GEMs from Hart 2015 datasets
%
% Input:
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
% Output:
%
%
%       results    evaluation with values of TP, TN, FP, FN, accuracy,
%                  sensitivity, specificity, F1, MCC
%
% Example usage:
%
%   results = evaluateEssentialityHart2015(eGenes)
%


%% load Hart2015 essentiality data

% Table S2 from the Hart2015 datasets:
% In this study, essential genes are, by definition, a subset of fitness genes
x = readtable('Hart2015_TableS2.xlsx');
    
% remove duplicated rows (MARCH1 and MARCH2 genes)
ind1 = find(ismember(x.Gene,'MARCH1'),1,'last');
ind2 = find(ismember(x.Gene,'MARCH2'),1,'last');
x([ind1;ind2],:) = [];
    
% extract gene list and cell types
genes = x.Gene;
celltypes = regexprep(upper(x.Properties.VariableNames(3:7)'),'BF_','');
    
% for each cell type, determine the "fitness" genes, defined as those with
% a 5% FDRs were chosen as thresholds (the values are extracted from the
% supporting information of the Hart 2015 paper. Genes observed in 3 or	more
% of the 5 TKO screens (n=1,580) were considered as core fitness genes
bf_data = table2array(x(:,3:7));
bf_thresh = {'HCT116',  1.57
    'HELA',   15.47
    'GBM',     3.20
    'RPE1',    6.84
    'DLD1',    3.57};

% construct fitness matrix according to threshold values
fitness_mat = zeros(numel(genes),numel(celltypes));
for i = 1:numel(celltypes)
    [~,ind] = ismember(celltypes{i},bf_thresh(:,1));
    fitness_mat(:,i) = bf_data(:,i) > bf_thresh{ind,2};
end
fitness_mat(isnan(bf_data)) = NaN;  % mark missing values differently than non-hits

% also get the genes that were essential in all 5 cell lines ("all")
celltypes(end+1) = {'all'};
fitness_mat(:,end+1) = all(fitness_mat == 1,2);
    
% put Hart2015 essentiality data into data structure for comparison
expdata = {};
expdata.genes = genes;
expdata.tissues = celltypes;
expdata.essential = fitness_mat;


%% Load and organize model-predicted gene essentiality results

% re-organize essentialGenes logical matrices into a uniform 3D matrix,
% where rows represent genes, columns are tasks, and the 3rd dimension is
% different cell lines.
eGenes.essentialMat = false(numel(eGenes.refModel.genes), numel(eGenes.taskList), numel(eGenes.tissues));
for i = 1:numel(eGenes.tissues)
    [hasMatch,ind] = ismember(eGenes.refModel.genes,eGenes.geneList{i});
    eGenes.essentialMat(hasMatch,:,i) = eGenes.essentialGenes{i}(ind(hasMatch),:);
end


%% compare model predictions to experimental essentiality data

% specify model-determined essential genes
modelPred = squeeze(any(eGenes.essentialMat,2));  % essential for any task

tissues = eGenes.tissues;
tissues(end+1) = {'all'};  % also include "all" category for Hart2015 dataset

% calculate true and false positives and negatives
[TP,TN,FP,FN,Penr] = deal(nan(numel(tissues),1));  % initialize variables
for i = 1:numel(tissues)
    
    [~,tissue_ind] = ismember(tissues(i), expdata.tissues);
    if tissue_ind == 0
        continue  % a few tissues are missing from the DepMap dataset
    end
    
    modelGenes = eGenes.refModel.genes;
    if strcmpi(tissues{i},'all')
        modelEssential = modelGenes(sum(modelPred,2) == size(modelPred,2));
    else
        modelEssential = modelGenes(modelPred(:,i));
    end
    modelNonEssential = setdiff(modelGenes,modelEssential);
    
    expGenes = expdata.genes(~isnan(expdata.essential(:,tissue_ind)));
    expEssential = expdata.genes(expdata.essential(:,tissue_ind) == 1);
    expNonEssential = setdiff(expGenes,expEssential);
    
    TP(i) = sum(ismember(modelEssential, expEssential));        % true positives
    TN(i) = sum(ismember(modelNonEssential, expNonEssential));  % true negatives
    FP(i) = sum(ismember(modelEssential, expNonEssential));     % false positives
    FN(i) = sum(ismember(modelNonEssential, expEssential));     % false negatives    
    
    Penr(i) = EnrichmentTest(intersect(modelGenes,expGenes), intersect(modelEssential,expGenes), intersect(expEssential,modelGenes));
end

% calculate some metrics
sensitivity = TP./(TP + FN);
specificity = TN./(TN + FP);
accuracy = (TP + TN)./(TP + TN + FP + FN);
F1 = 2*TP./(2*TP + FP + FN);
MCC = ((TP.*TN) - (FP.*FN))./sqrt((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN));  % Matthews correlation coefficient
PenrAdj = adjust_pvalues(Penr,'Benjamini');

% get results for cell types
results = [{'cellLine','TP','TN','FP','FN','accuracy','sensitivity','specificity','F1','MCC','Penr','logPenr','PenrAdj','logPenrAdj'};
          [tissues, num2cell([TP, TN, FP, FN, accuracy, sensitivity, specificity, F1, MCC, Penr, -log10(Penr), PenrAdj, -log10(PenrAdj)])]];


end

function [penr,pdep] = EnrichmentTest(pop,sample,successes)
% Caculates p-value of enrichment of successes in a sample.
%
% Evaluates the significance of enrichment (and depletion) of SUCCESSES in
% a SAMPLE drawn from population POP using the hypergeometric test.

x = numel(intersect(successes,sample));  % calc # of successes in sample
m = numel(pop);  % calc size of population
k = numel(intersect(successes,pop));
n = numel(sample);

penr = hygecdf(x-1,m,k,n,'upper');
pdep = hygecdf(x,m,k,n);

end

