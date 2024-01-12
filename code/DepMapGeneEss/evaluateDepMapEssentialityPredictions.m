%This code is largely copied from the single cell modeling paper, and there in turn copied from the Human1 paper ("An atlas of human metabolism").
%The code compares the essential genes predicted by the model with experimental data from CRISPR screens.

%cd C:/Work/MatlabCode/components/human-GEM/Human-GEMDepMapEval/Human-GEM/code/DepMapGeneEss %to be used if not running the whole file, may need to change this
cd(fileparts(which(mfilename))); %to be used if running the whole file

nModelComps = 1; %change to 2 if you want to compare 2 models (i.e., 2 versions of Human-GEM), 3 for three etc.
nChunksPerModel = 2;%Change to 40 for a full evaluation


%% load CRISPR screen data
    
x = readtable('data/Achilles_gene_effect.csv');

genes = x.Properties.VariableDescriptions(2:end)';
genes = regexprep(genes,"Original column heading: '",'');
%genes = regexprep(genes,"\s\(\d+\)$",'');
genes = regexprep(genes,"\s\(\d+\)'",'');

celltypes = x.DepMap_ID;

% convert to matrix and transpose
x = table2array(x(:,2:end))';

% collect results
expdata = {};
expdata.genes = genes;
expdata.tissues = celltypes;
expdata.gene_effect = x;

%clear x genes celltypes col_types cell_col i


%% Load and organize model-predicted gene essentiality results

% load model analysis results
%   Results are organized in the "eGenes" structure, with the fields:
%       .taskList   list of metabolic tasks tested
%       .tissues    list of tissue (cell line) names corresponding to each
%                   of the tINIT models
%       .geneList   cell array of cell arrays, where each entry contains a
%                   list of gene names contained in each of the tINIT
%                   models
%       .essentialGenes   cell array of logical matrices, with genes as
%                         rows (corresponding to genes in ".geneList") and
%                         metabolic tasks as columns (corresponding to
%                         tasks in ".taskList"). If a gene is essential for
%                         a task, its corresponding entry will be 1,
%                         otherwise it will be zero.
%                         

% load reference model (Human-GEM)
x = load('data/prepDataGeneSymbols.mat');
refModel = x.prepData.refModel;

allFiles = cell(nModelComps,1);
allFiles{1} = cell(nChunksPerModel,1);
%allFiles{2} = cell(nChunksPerModel,1);
%allFiles{3} = cell(nChunksPerModel,1);
allOutFiles = {'data/geneEss_model1.txt'};

%allOutFiles = {'data/geneEss_model1.txt';
%            'data/geneEss_model2.txt';
%            'data/geneEss_model3.txt'
%            };
for i = 1:nChunksPerModel
    allFiles{1}{i} = ['data/depmap_models-eGenes-' num2str(i) '.mat'];
%    allFiles{2}{i} = ['../data/DepMap/ftINIT/depmap_models_newalg-eGenes-' num2str(i) '.mat'];
%    allFiles{3}{i} = ['../data/DepMap/ftINIT2/depmap_models_newalg-eGenes-' num2str(i) '.mat'];
end


for j = 1:length(allFiles)
        
    % load and merge gene essentiality results
    files = allFiles{j};
    eGenesAll = {};
    for i = 1:numel(files)
        disp(i)
        load(files{i})
        if i == 1
            eGenesAll = eGenes;
        else
            eGenesAll.tissues = [eGenesAll.tissues; eGenes.tissues];
            eGenesAll.geneList = [eGenesAll.geneList; eGenes.geneList];
            eGenesAll.essentialGenes = [eGenesAll.essentialGenes; eGenes.essentialGenes];
        end
    end
    eGenes = eGenesAll;

    % re-organize essentialGenes logical matrices into a uniform 3D matrix,
    % where rows represent genes, columns are tasks, and the 3rd dimension is
    % different cell lines.
    for i = 1:numel(eGenes.tissues)
        [hasMatch,ind] = ismember(refModel.genes, eGenes.geneList{i});
        eGenes.essentialMat(hasMatch,:,i) = eGenes.essentialGenes{i}(ind(hasMatch),:);
    end

    % clear intermediate variables
    clear hasMatch i ind refModelInd resultsFileName eGenesAll files f x


    %% compare model predictions to experimental essentiality data

    % % To scan several thresholds for the depmap dataset
    % scan.threshold = (linspace(-1,0,21))';
    % scan.mcc = zeros(numel(eGenes.tissues),numel(scan.threshold));
    % for k = 1:numel(scan.threshold)
    %     DepMapThresh = scan.threshold(k);


    % specify model-determined essential genes
    modelPred = squeeze(any(eGenes.essentialMat,2));  % essential for any task
    % modelPred = squeeze(any(eGenes.essentialMat(:,end,:),2));  % essential for only biomass production

    % specify essentiality threshold (for DepMap data only)
    threshold = -0.6;

    tissues = eGenes.tissues;
    expdata.essential = double(expdata.gene_effect < threshold);

    % calculate true and false positives and negatives
    [TP,TN,FP,FN,Penr] = deal(nan(numel(tissues),1));  % initialize variables
    for i = 1:numel(tissues)

        [~,tissue_ind] = ismember(tissues(i), expdata.tissues);
        if tissue_ind == 0
            error('bad');
        end

        modelGenes = refModel.genes;
        if strcmpi(tissues{i},'all')
            modelEssential = modelGenes(sum(modelPred,2) == size(modelPred,2));
        else
            modelEssential = modelGenes(modelPred(:,i));
        end
        modelNonEssential = setdiff(modelGenes, modelEssential);

        expGenes = expdata.genes(~isnan(expdata.essential(:, tissue_ind)));
        expEssential = expdata.genes(expdata.essential(:, tissue_ind) == 1);
        expNonEssential = setdiff(expGenes, expEssential);

        TP(i) = sum(ismember(modelEssential, expEssential));        % true positives
        TN(i) = sum(ismember(modelNonEssential, expNonEssential));  % true negatives
        FP(i) = sum(ismember(modelEssential, expNonEssential));     % false positives
        FN(i) = sum(ismember(modelNonEssential, expEssential));     % false negatives    

        Penr(i) = enrichment_test(intersect(modelGenes,expGenes), intersect(modelEssential,expGenes), intersect(expEssential,modelGenes));
    end

    % calculate some metrics
    sensitivity = TP./(TP + FN);
    specificity = TN./(TN + FP);
    accuracy = (TP + TN)./(TP + TN + FP + FN);
    F1 = 2*TP./(2*TP + FP + FN);
    MCC = ((TP.*TN) - (FP.*FN))./sqrt((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN));  % Matthews correlation coefficient
    PenrAdj = adjust_pvalues(Penr, 'Benjamini');

    %% Write results to file

    results = [{'cellLine','TP','TN','FP','FN','accuracy','sensitivity','specificity','F1','MCC','Penr','logPenr','PenrAdj','logPenrAdj'};
              [tissues, num2cell([TP, TN, FP, FN, accuracy, sensitivity, specificity, F1, MCC, Penr, -log10(Penr), PenrAdj, -log10(PenrAdj)])]];

    % include threshold without sign or decimal (e.g., -0.6 = '06')
    %outfile = strcat('essentiality_prediction_results_', regexprep(num2str(abs(threshold)),'\.',''), '.txt');      
    %outfile = fullfile(proj_dir, 'results', 'gene_essentiality', outfile);
    outfile = allOutFiles{j};
    writecell(results, outfile, 'Delimiter', '\t');
end


%% Now get model statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load the models


allFilesTempl = {'data/depmap_models-'};
%allFilesTempl = {'data/depmap_models-';
%            '../data/DepMap/ftINIT/depmap_models_newalg-';
%            '../data/DepMap/ftINIT2/depmap_models_newalg-'
%            };
        
stats = cell(length(allFilesTempl),1);

for j = 1:length(allFilesTempl)
    %get the stats
    stat = struct();
    stat.numRxns = []; %will be extended in the loop
    stat.numRxnsWithGpr = []; %will be extended in the loop
    for i = 1:nChunksPerModel
        disp(num2str(i))
        if j == 1
            x = load([allFilesTempl{j} num2str(i) '.mat']).depmap_models;
        else
            x = load([allFilesTempl{j} num2str(i) '.mat']).depmap_models;
        end
        for k = 1:length(x)
            m = x{k};
            stat.numRxns = [stat.numRxns;length(m.rxns)];
            stat.numRxnsWithGpr = [stat.numRxnsWithGpr;sum(~strcmp(m.grRules,''))];
        end
    end

    stats{j} = stat;
end

save('data/ModelStatData.mat', 'stats')


