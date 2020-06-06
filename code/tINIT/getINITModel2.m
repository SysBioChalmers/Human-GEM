function [model, metProduction, essentialRxnsForTasks, addedRxnsForTasks, deletedDeadEndRxns, deletedRxnsInINIT, taskReport] = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT)
% getINITModel2
%   Generates a model using the INIT algorithm, based on proteomics and/or
%   transcriptomics and/or metabolomics and/or metabolic tasks.
%   This is the newer version, which has updated handling of gene rules to
%   differentiate between isozymes and enzyme complexes.
%
%   refModel            a model structure. The model should be in the
%                       closed form (no exchange reactions open). Import
%                       using import(filename,false). If the model is not
%                       loaded using importModel, it might be that there
%                       is no "unconstrained" field. In that case,
%                       manually add the field like:
%                       model.unconstrained = false(numel(model.mets),1);
%   tissue              tissue to score for. Should exist in either
%                       hpaData.tissues or arrayData.tissues
%   celltype            cell type to score for. Should exist in either
%                       hpaData.celltypes or arrayData.celltypes for this
%                       tissue (opt, default is to use the max values
%                       among all the cell types for the tissue.
%   hpaData             HPA data structure from parseHPA (opt if arrayData
%                       is supplied, default [])
%   arrayData           gene expression data structure (opt if hpaData is
%                       supplied, default [])
%       genes           cell array with the unique gene names
%       tissues         cell array with the tissue names. The list may not
%                       be unique, as there can be multiple cell types per
%                       tissue
%       celltypes       cell array with the cell type names for each tissue
%       levels          GENESxTISSUES array with the expression level for
%                       each gene in each tissue/celltype. NaN should be
%                       used when no measurement was performed
%       threshold       a single value or a vector of gene expression 
%                       thresholds, above which genes are considered to be
%                       "expressed". (opt, by default, the mean expression
%                       levels of each gene across all tissues in arrayData
%                       will be used as the threshold values)
%       singleCells     binary value selecting whether to use the
%                       single-cell algorithm to identify expressed genes.
%                       If used, specify cell subpopulations in CELLTYPES
%                       (opt, default [])
%       plotResults     true if single cell probability distributions
%                       should be plotted (opt, default = false)
%   metabolomicsData    cell array with metabolite names that the model
%                       should produce (opt, default [])
%   removeGenes         if true, low-abundance genes will be removed from
%                       grRules, unless they are the only gene associated 
%                       with a reaction, or a subunit of an enzyme complex
%                       (see "removeLowScorGenes" function for details).
%                       If false, grRules will not be modified; however,
%                       genes that were associated only with removed 
%                       reactions will not be present in the final model.
%                       (opt, default true).
%   taskFile            a task list in Excel format. See parseTaskList for
%                       details (opt, default [])
%   useScoresForTasks   true if the calculated reaction scored should be 
%                       used as weights when fitting to tasks (opt, default
%                       true)
%   printReport         true if a report should be printed to the screen
%                       (opt, default true)
%   taskStructure       task structure generated from parseTaskList. Can be
%                       used as an alternative way to define tasks when Excel
%                       sheets are not suitable. Overrides taskFile (opt,
%                       default [])
%   params              parameter structure as used by getMILPParams. This
%                       is for the INIT algorithm. For the the MILP 
%                       problems solved to fit tasks, see paramsFT (opt,
%                       default [])
%   paramsFT            parameter structure as used by getMILPParams. This
%                       is for the fitTasks step. For the INIT algorithm,
%                       see params (opt, default [])
%
%   model                   the resulting model structure
%   metProduction           array that indicates which of the
%                           metabolites in metabolomicsData that could be
%                           produced. Note that this is before the
%                           gap-filling process to enable defined tasks. To
%                           see which metabolites that can be produced in
%                           the final model, use canProduce.
%                           -2: metabolite name not found in model
%                           -1: metabolite found, but it could not be produced
%                           1: metabolite could be produced
%   essentialRxnsForTasks   cell array of the reactions which were
%                           essential to perform the tasks
%   addedRxnsForTasks       cell array of the reactions which were added in
%                           order to perform the tasks
%   deletedDeadEndRxns      cell array of reactions deleted because they
%                           could not carry flux (INIT requires a
%                           functional input model)
%   deletedRxnsInINIT       cell array of the reactions which were deleted by
%                           the INIT algorithm
%   taskReport              structure with the results for each task
%   	id                  cell array with the id of the task
%       description         cell array with the description of the task
%       ok                  boolean array with true if the task was successful
%       essential           cell array with cell arrays of essential
%                           reactions for the task
%       gapfill             cell array of cell arrays of reactions included
%                           in the gap-filling for the task
%
%   This is the main function for automatic reconstruction of models based
%   on the INIT algorithm (PLoS Comput Biol. 2012;8(5):e1002518). Not all
%   settings are possible using this function, and you may want to call the
%   functions scoreComplexModel, runINIT and fitTasks individually instead.
%
%   NOTE: Exchange metabolites should normally not be removed from the model
%   when using this approach, since checkTasks/fitTasks rely on putting specific
%   constraints for each task. The INIT algorithm will remove exchange metabolites
%   if any are present. Use importModel(file,false) to import a model with
%   exchange metabolites remaining.
%
%   Usage: [model, metProduction, essentialRxnsForTasks, addedRxnsForTasks,...
%               deletedDeadEndRxns, deletedRxnsInINIT, taskReport] = ...
%               getINITModel2(refModel, tissue, celltype, hpaData, arrayData,...
%               metabolomicsData, removeGenes, taskFile, useScoresForTasks, ...
%               printReport, taskStructure, params, paramsFT);
%


if nargin < 5
    arrayData = [];
end
if nargin < 6
    metabolomicsData = [];
end
if nargin < 7 || isempty(removeGenes)
    removeGenes = true;
end
if nargin < 8
    taskFile = [];
end
if nargin < 9 || isempty(useScoresForTasks)
    useScoresForTasks = true;
end
if nargin < 10 || isempty(printReport)
    printReport = true;
end
if nargin < 11
    taskStructure = [];
end
if nargin < 12
    params = [];
end
if nargin < 13
    paramsFT = [];
end

%Check that the model is in the closed form
if ~isfield(refModel,'unconstrained')
    EM = 'Exchange metabolites should normally not be removed from the model when using getINITModel. Use importModel(file,false) to import a model with exchange metabolites remaining (see the documentation for details)';
    dispEM(EM);
end

%Create the task structure if not supplied
if any(taskFile) && isempty(taskStructure)
    taskStructure = parseTaskList(taskFile);
end

% sc-tINIT to identify confidence levels of gene expression
if isfield(arrayData,'singleCells')
    if arrayData.singleCells == 1
        % Check to ensure cell type is defined
        if ~isfield(arrayData,'celltypes')
            dispEM('arrayData must contain cell type information if sc-tINIT is to be used','false');   
        end
        if ~ismember(upper(celltype),upper(arrayData.celltypes))
            dispEM('The cell type name does not match');   
        end
        
        % Analyze only cell type of interest
        J = strcmpi(arrayData.celltypes,celltype);
        
        % Analyze only genes included in the reference model
        I = ismember(arrayData.genes,refModel.genes);
        
        % Convert expression data to population fractions
        binary_levels = arrayData.levels(I,J) ~= 0; % Classify each gene as detected (1) or not (0) in each cell
        cell_count_levels = sum(binary_levels,2); % Number of cells expressing each transcript
        cell_frac_levels = cell_count_levels/size(binary_levels,2); % Number of cells expressing each transcript

        % Bin cell_frac_counts manually
        x = 0:.01:1;
        for i = 1:length(x)
            cell_frac_count(i) = sum(round(cell_frac_levels,2) == x(i));
        end
        
        % Fit four beta distributions
        cell_frac_count(cell_frac_count == 0) = NaN; % Remove zeros from optimization
        cell_frac_count(1) = NaN; % Remove non-expressed genes from optimization
        x_lim = 1; % Somewhat arbitrary parameter to fit left tail of distr.
        myfun = @(par) nansum((cell_frac_count(1:find(x >= x_lim,1)) - ...
            abs(par(1))*betapdf(x(1:find(x >= x_lim,1)),abs(par(2)),abs(par(3))) - ...
            abs(par(4))*betapdf(x(1:find(x >= x_lim,1)),abs(par(5)),abs(par(6))) - ...
            abs(par(7))*betapdf(x(1:find(x >= x_lim,1)),abs(par(8)),abs(par(9))) - ...
            abs(par(10))*betapdf(x(1:find(x >= x_lim,1)),abs(par(11)),abs(par(12)))).^2);
        
        par0 = [4,2,100,7,2,30,7,5,20,5,15,20];
        opts = optimset('Display','off');
        [par,f_val] = fminsearch(myfun,par0,opts);
        par = abs(par);
        
        % Plot results
        if (isfield(arrayData,'plotResults'))
            if arrayData.plotResults == true
                figure(); hold on; plot(x,cell_frac_count,'ko','MarkerSize',5);
                plot(x,abs(par(1))*betapdf(x,abs(par(2)),abs(par(3))),'b-','LineWidth',1)
                plot(x,abs(par(4))*betapdf(x,abs(par(5)),abs(par(6))),'b-','LineWidth',1)
                plot(x,abs(par(7))*betapdf(x,abs(par(8)),abs(par(9))),'b-','LineWidth',1)
                plot(x,abs(par(10))*betapdf(x,abs(par(11)),abs(par(12))),'b-','LineWidth',1)
                plot(x,abs(par(1))*betapdf(x,abs(par(2)),abs(par(3))) + ...
                    abs(par(4))*betapdf(x,abs(par(5)),abs(par(6))) + ...
                    abs(par(7))*betapdf(x,abs(par(8)),abs(par(9))) + ...
                    abs(par(10))*betapdf(x,abs(par(11)),abs(par(12))),'-','Color',[.5 .5 .5],'LineWidth',2)
                xlabel('Expression Probability');ylabel('# of genes');set(gca,'FontSize',14,'LineWidth',1.25);
                title('Expression prediction','FontSize',18,'FontWeight','bold')
            end
        end
        
        % Score genes based on population expression (p = .05)
        exprs_cutoff_1 = find(cdf('beta',x,par(2),par(3)) >.95,1)-1; % Find index of no confidence genes
        exprs_cutoff_2 = find(cdf('beta',x,par(5),par(6)) >.95,1)-1; % Find index of low confidence genes
        exprs_cutoff_3 = find(cdf('beta',x,par(8),par(9)) >.95,1)-1; % Find index of low confidence genes
        exprs_cutoffs = sort([exprs_cutoff_1,exprs_cutoff_2,exprs_cutoff_3]);
        gene_scores = cell_frac_levels*0;
        gene_scores(cell_frac_levels <= x(exprs_cutoffs(1))) = 4; % Not detected
        gene_scores(logical((cell_frac_levels >= x(exprs_cutoffs(1))).*(cell_frac_levels < x(exprs_cutoffs(2))))) = 3; % Low detection
        gene_scores(logical((cell_frac_levels >= x(exprs_cutoffs(2))).*(cell_frac_levels < x(exprs_cutoffs(3))))) = 2; % Medium detection
        gene_scores(cell_frac_levels > x(exprs_cutoffs(3))) = 1; % High detection

        % Replace hpaData with singleCellData
        if printReport == true
            dispEM('Single cell data is not currently compatible with HPA data. \n         Replacing hpaData with single cell-based scoring.',false);
        end
        hpaData.genes = arrayData.genes;
        hpaData.tissues = arrayData.tissues;
        hpaData.celltypes = arrayData.celltypes;
        hpaData.levels = [{'High'},{'Medium'},{'Low'},{'None'}];
        hpaData.gene2Level = zeros(length(arrayData.genes),length(arrayData.celltypes));
        for i = 1:length(find(J))
            find_var = find(J,i);
            hpaData.gene2Level(I,find_var(end)) = gene_scores;
        end
        
        % Remove arrayData from the analysis (Might be a bad idea)
        clear arrayData
        arrayData = [];
    end
end


if printReport == true
    if any(celltype)
        fprintf(['***Generating model for: ' tissue ' - ' celltype '\n']);
    else
        fprintf(['***Generating model for: ' tissue '\n']);
    end
    if ~isempty(hpaData)
        fprintf('-Using HPA data\n');
    end
    if ~isempty(arrayData)
        fprintf('-Using array data\n');
    end
    if ~isempty(metabolomicsData)
        fprintf('-Using metabolomics data\n');
    end
    if ~isempty(taskFile) || ~isempty(taskStructure)
        fprintf('-Using metabolic tasks\n');
    end
    fprintf('\n');
    
    printScores(refModel,'Reference model statistics',hpaData,arrayData,tissue,celltype);
end

%Remove dead-end reactions to speed up the optimization and to
%differentiate between reactions removed by INIT and those that are
%dead-end
if isfield(arrayData,'deletedDeadEndRxns')                  % temporary shortcut, should be revised.
    deletedDeadEndRxns = arrayData.deletedDeadEndRxns;      % temporary shortcut, should be revised.
elseif isfield(hpaData,'deletedDeadEndRxns')                % temporary shortcut, should be revised.
    deletedDeadEndRxns = hpaData.deletedDeadEndRxns;        % temporary shortcut, should be revised.
else                                                        % temporary shortcut, should be revised.
    [~, deletedDeadEndRxns]=simplifyModel(refModel,true,false,true,true,true);
end                                                         % temporary shortcut, should be revised.
cModel = removeReactions(refModel,deletedDeadEndRxns,false,true);

%Store the connected model like this to keep track of stuff
if printReport == true
    printScores(cModel,'Pruned model statistics',hpaData,arrayData,tissue,celltype);
end

%If tasks have been defined, then go through them and get essential
%reactions
if ~isempty(taskStructure)

    if isfield(arrayData,'taskReport') && isfield(arrayData,'essentialRxnMat')  % temporary shortcut, should be revised.
        taskReport = arrayData.taskReport;                                      % temporary shortcut, should be revised.
        essentialRxnMat = arrayData.essentialRxnMat;                            % temporary shortcut, should be revised.
    elseif isfield(hpaData,'taskReport') && isfield(hpaData,'essentialRxnMat')  % temporary shortcut, should be revised.
        taskReport = hpaData.taskReport;                                        % temporary shortcut, should be revised.
        essentialRxnMat = hpaData.essentialRxnMat;                              % temporary shortcut, should be revised.
    else                                                                        % temporary shortcut, should be revised.
        [taskReport, essentialRxnMat]=checkTasks(cModel,[],printReport,true,true,taskStructure);
    end                                                                         % temporary shortcut, should be revised.
    
    essentialRxnsForTasks = cModel.rxns(any(essentialRxnMat,2));
    
    %Find metabolites present in taskStruct. We want to avoid removing
    %these metabolites from the final model (even though they may no longer
    %participate in any reacitons) so that the final model is still able to
    %complete all of the tasks without any errors.
    taskMets = union(vertcat(taskStructure.inputs),vertcat(taskStructure.outputs));
    taskMets = union(taskMets, parseRxnEqu(vertcat(taskStructure.equations)));
    modelMets = strcat(cModel.metNames,'[',cModel.comps(cModel.metComps),']');
    [inModel,metInd] = ismember(taskMets,modelMets);
    essentialMetsForTasks = cModel.mets(metInd(inModel));
    
    %Remove tasks that cannot be performed
    taskStructure(taskReport.ok == false) = [];
    if printReport == true
        printScores(removeReactions(cModel,setdiff(cModel.rxns,essentialRxnsForTasks),true,true),'Reactions essential for tasks',hpaData,arrayData,tissue,celltype);
    end
else
    essentialRxnsForTasks = {};
    essentialMetsForTasks = {};
end

% Score the connected model
rxnScores = scoreComplexModel(cModel,hpaData,arrayData,tissue,celltype);

%Run the INIT algorithm. The exchange reactions that are used in the final
%model will be open, which doesn't fit with the last step. Therefore
%delete reactions from the original model instead of taking the output. The
%default implementation does not constrain reversible reactions to only
%carry flux in one direction. Runs without the constraints on reversibility
%and with all output allowed. This is to reduce the complexity of the
%problem.
[~, deletedRxnsInINIT, metProduction] = runINIT(simplifyModel(cModel),rxnScores,metabolomicsData,essentialRxnsForTasks,0,true,false,params);
initModel = removeReactions(cModel,deletedRxnsInINIT,false,true);

% remove metabolites separately to avoid removing those needed for tasks
unusedMets = initModel.mets(all(initModel.S == 0,2));
initModel = removeMets(initModel, setdiff(unusedMets,essentialMetsForTasks));

if printReport == true
    printScores(initModel,'INIT model statistics',hpaData,arrayData,tissue,celltype);
    printScores(removeReactions(cModel,setdiff(cModel.rxns,deletedRxnsInINIT),true,true),'Reactions deleted by INIT',hpaData,arrayData,tissue,celltype);
end

%The full model has exchange reactions in it. fitTasks calls on fillGaps,
%which automatically removes exchange metabolites (because it assumes that
%the reactions are constrained when appropriate). In this case the
%uptakes/outputs are retrieved from the task sheet instead. To prevent
%exchange reactions being used to fill gaps, they are delete from the
%reference model here.
initModel.id = 'INITModel';

%If gaps in the model should be filled using a task list
if ~isempty(taskStructure)
    %Remove exchange reactions and reactions already included in the INIT
    %model
    refModelNoExc = removeReactions(refModel,union(initModel.rxns,getExchangeRxns(refModel)),true,true);
    
    %At this stage the model is fully connected and most of the genes with
    %good scores should have been included. The final gap-filling should
    %take the scores of the genes into account, so that "rather bad"
    %reactions are preferred to "very bad" reactions. However, reactions
    %with positive scores will be included even if they are not connected
    %in the current formulation. Therefore, such reactions will have to be
    %assigned a small negative score instead.
    if useScoresForTasks == true
        refRxnScores = scoreComplexModel(refModelNoExc,hpaData,arrayData,tissue,celltype);
        [outModel,addedRxnMat] = fitTasks(initModel,refModelNoExc,[],true,min(refRxnScores,-0.1),taskStructure,paramsFT);
    else
        [outModel,addedRxnMat] = fitTasks(initModel,refModelNoExc,[],true,[],taskStructure,paramsFT);
    end
    if printReport == true
        printScores(outModel,'Functional model statistics',hpaData,arrayData,tissue,celltype);
        printScores(removeReactions(outModel,intersect(outModel.rxns,initModel.rxns),true,true),'Reactions added to perform the tasks',hpaData,arrayData,tissue,celltype);
    end
    
    addedRxnsForTasks = refModelNoExc.rxns(any(addedRxnMat,2));
else
    outModel = initModel;
    addedRxnMat = [];
    addedRxnsForTasks = {};
end

% The model can now perform all the tasks defined in the task list.
model = outModel;

% If requested, attempt to remove negative-score genes from the model, 
% depending on their role (isozyme or complex subunit) in each grRule.
% See the "removeLowScoreGenes" function more more details, and to adjust
% any default parameters therein.
if ( removeGenes )
    [~, geneScores] = scoreComplexModel(model,hpaData,arrayData,tissue,celltype);
    model = removeLowScoreGenes(model,geneScores);
end


% At this stage the model will contain some exchange reactions but probably
% not all (and maybe zero). This can be inconvenient, so all exchange
% reactions from the reference model are added, except for those which
% involve metabolites that are not in the model.

% First delete any included exchange reactions in order to prevent the order
% from changing
model = removeReactions(model,getExchangeRxns(model));

% Create a model with only the exchange reactions in refModel
excModel = removeReactions(refModel,setdiff(refModel.rxns,getExchangeRxns(refModel)),true,true);

% Find the metabolites there which are not exchange metabolites and which do
% not exist in the output model
I = ~ismember(excModel.mets,model.mets) & excModel.unconstrained == 0;

% Then find those reactions and delete them
[~, J] = find(excModel.S(I,:));
excModel = removeReactions(excModel,J,true,true);

% Merge with the output model
model = mergeModels({model;excModel},'metNames');
model.id = 'INITModel';
model.description = ['Automatically generated model for ' tissue];
if any(celltype)
    model.description = [model.description ' - ' celltype];
end

if printReport == true
    printScores(model,'Final model statistics',hpaData,arrayData,tissue,celltype);
end

% Add information about essential reactions and reactions included for
% gap-filling and return a taskReport
if ~isempty(taskStructure)
    I = find(taskReport.ok); %Ignore failed tasks
    for i = 1:numel(I)
        taskReport.essential{I(i),1} = cModel.rxns(essentialRxnMat(:,I(i)));
        taskReport.gapfill{I(i),1} = refModelNoExc.rxns(addedRxnMat(:,i));
    end
else
    taskReport = [];
end

%Fix grRules and reconstruct rxnGeneMat
[grRules,rxnGeneMat] = standardizeGrRules(model,true);
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;
end

%This is for printing a summary of a model
function [rxnS, geneS] = printScores(model,name,hpaData,arrayData,tissue,celltype)
[a, b] = scoreComplexModel(model,hpaData,arrayData,tissue,celltype);
rxnS = mean(a);
geneS = mean(b,'omitnan');
fprintf([name ':\n']);
fprintf(['\t' num2str(numel(model.rxns)) ' reactions, ' num2str(numel(model.genes)) ' genes\n']);
fprintf(['\tMean reaction score: ' num2str(rxnS) '\n']);
fprintf(['\tMean gene score: ' num2str(geneS) '\n']);
fprintf(['\tReactions with positive scores: ' num2str(100*sum(a>0)/numel(a)) '%%\n\n']);
end
