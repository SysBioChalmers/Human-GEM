function [taskReport, essentialGenes, taskStructure]=checkTasksGenes(model,inputFile,printOutput,printOnlyFailed,getEssential,taskStructure)
% checkTasksGenes
%   Performs a set of simulations as defined in a task file. This function
%   is identical to "checkTasks", except it allows for the determination of
%   essential genes, rather than essential reactions.
%
%   model           a model structure
%   inputFile       a task list in Excel format. See the function
%                   parseTaskList for details (opt if taskStructure is
%                   supplied)
%   printOutput     true if the results of the test should be displayed
%                   (opt, default true)
%   printOnlyFailed true if only tasks that failed should be displayed
%                   (opt, default false)
%   getEssential    true if the essential genes should be calculated for
%                   all the tasks. (opt, default false)
%   taskStructure   structure with the tasks, as from parseTaskList. If
%                   this is supplied then inputFile is ignored (opt)
%
%   taskReport          structure with the results
%       id              cell array with the id of the task
%       description     cell array with the description of the task
%       ok              boolean array with true if the task was successful
%   essentialGenes      MxN matrix with the essential genes (M) for each
%                       task (N). An element is true if the corresponding
%                       gene is essential in the corresponding task.
%                       Failed tasks and SHOULD FAIL tasks are ignored.
%                       If getEssential is false, then essentialRxns will 
%                       be returned as: false(nGenes,nTasks).
%   taskStructure       structure with the tasks, as from parseTaskList
%
%   This function is used for defining a set of tasks for a model to
%   perform. The tasks are defined by defining constraints on the model,
%   and if the problem is feasible, then the task is considered successful.
%   In general, each row can contain one constraint on uptakes, one
%   constraint on outputs, one new equation, and one change of reaction
%   bounds. If more bounds are needed to define the task, then several rows
%   can be used for each task.
%
%   Usage: [taskReport, essentialGenes, taskStructure]=checkTasksGenes(model,inputFile,...
%           printOutput,printOnlyFailed,getEssential,taskStructure)
%


if nargin<3 || isempty(printOutput)
    printOutput=true;
end
if nargin<4 || isempty(printOnlyFailed)
    printOnlyFailed=false;
end
if nargin<5 || isempty(getEssential)
    getEssential=false;
end

%Prepare the input model
model.b=zeros(numel(model.mets),2);

modelMets=upper(strcat(model.metNames,'[',model.comps(model.metComps),']'));
if ~isfield(model,'unconstrained')
    EM='Exchange metabolites should normally not be removed from the model when using checkTasks. Inputs and outputs are defined in the task file instead. Use importModel(file,false) to import a model with exchange metabolites remaining';
    dispEM(EM,false);
end

%Parse the task file
if nargin<6
    taskStructure=parseTaskList(inputFile);
end

essentialGenes=false(numel(model.genes),numel(taskStructure));

tModel=model;
taskReport=[];
for i=1:numel(taskStructure)
    taskReport.id{i,1}=taskStructure(i).id;
    taskReport.description{i,1}=taskStructure(i).description;
    %Set the inputs
    if ~isempty(taskStructure(i).inputs)
        [I, J]=ismember(upper(taskStructure(i).inputs),modelMets);
        J=J(I); %Only keep the ones with matches
        K=ismember(upper(taskStructure(i).inputs),'ALLMETS');
        L=~cellfun('isempty',strfind(upper(taskStructure(i).inputs),'ALLMETSIN'));
        %Check that all metabolites are either real metabolites or
        %ALLMETS/ALLMETSIN
        if ~all(I|K|L)
            fprintf(['ERROR: Could not find all inputs in "[' taskStructure(i).id '] ' taskStructure(i).description '"\n']);
            taskReport.ok(i,1)=false;
            tModel=model;
            continue;
        end
        if numel(J)~=numel(unique(J))
            EM=['The constraints on some input(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time'];
            dispEM(EM);
        end
        %If all metabolites should be added
        if any(K)
            %Check if ALLMETS is the first metabolite. Otherwise print a
            %warning since it will write over any other constraints that
            %are set
            if K(1)==0
                EM=['ALLMETS is used as an input in "[' taskStructure(i).id '] ' taskStructure(i).description '" but it it not the first metabolite in the list. Constraints defined for the metabolites before it will be over-written'];
                dispEM(EM,false);
            end
            %Use the first match of ALLMETS. There should only be one, but
            %still..
            tModel.b(:,1)=taskStructure(i).UBin(find(K,1))*-1;
        end
        %If metabolites in a specific compartment should be used
        if any(L)
            L=find(L);
            for j=1:numel(L)
                %The compartment defined
                compartment=upper(taskStructure(i).inputs{L(j)}(11:end-1));
                %Check if it exists in the model
                C=find(ismember(upper(model.comps),compartment));
                if any(C)
                    %Match to metabolites
                    tModel.b(model.metComps==C,1)=taskStructure(i).UBin(L(j))*-1;
                else
                    EM=['The compartment defined for ALLMETSIN in "[' taskStructure(i).id '] ' taskStructure(i).description '" does not exist'];
                    dispEM(EM);
                end
            end
        end
        %Then add the normal constraints
        if any(J)
            tModel.b(J,1)=taskStructure(i).UBin(I)*-1;
            tModel.b(J,2)=taskStructure(i).LBin(I)*-1;
        end
    end
    %Set the outputs
    if ~isempty(taskStructure(i).outputs)
        [I, J]=ismember(upper(taskStructure(i).outputs),modelMets);
        J=J(I); %Only keep the ones with matches
        K=ismember(upper(taskStructure(i).outputs),'ALLMETS');
        L=~cellfun('isempty',strfind(upper(taskStructure(i).outputs),'ALLMETSIN'));
        %Check that all metabolites are either real metabolites or
        %ALLMETS/ALLMETSIN
        if ~all(I|K|L)
            fprintf(['ERROR: Could not find all outputs in "[' taskStructure(i).id '] ' taskStructure(i).description '"\n']);
            taskReport.ok(i,1)=false;
            tModel=model;
            continue;
        end
        if numel(J)~=numel(unique(J))
            EM=['The constraints on some output(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time'];
            dispEM(EM);
        end
        %If all metabolites should be added
        if any(K)
            %Check if ALLMETS is the first metabolite. Otherwise print a
            %warning since it will write over any other constraints that
            %are set
            if K(1)==0
                EM=['ALLMETS is used as an output in "[' taskStructure(i).id '] ' taskStructure(i).description '" but it it not the first metabolite in the list. Constraints defined for the metabolites before it will be over-written'];
                dispEM(EM,false);
            end
            %Use the first match of ALLMETS. There should only be one, but
            %still..
            tModel.b(:,2)=taskStructure(i).UBout(find(K,1));
        end
        %If metabolites in a specific compartment should be used
        if any(L)
            L=find(L);
            for j=1:numel(L)
                %The compartment defined
                compartment=upper(taskStructure(i).outputs{L(j)}(11:end-1));
                %Check if it exists in the model
                C=find(ismember(upper(model.comps),compartment));
                if any(C)
                    %Match to metabolites
                    tModel.b(model.metComps==C,2)=taskStructure(i).UBout(L(j));
                else
                    EM=['The compartment defined for ALLMETSIN in "[' taskStructure(i).id '] ' taskStructure(i).description '" does not exist'];
                    dispEM(EM);
                end
            end
        end
        %Then add the normal constraints
        if any(J)
            tModel.b(J,1)=taskStructure(i).LBout(I);
            tModel.b(J,2)=taskStructure(i).UBout(I);
        end
    end
    %Add new rxns
    if ~isempty(taskStructure(i).equations)
        rxn.equations=taskStructure(i).equations;
        rxn.lb=taskStructure(i).LBequ;
        rxn.ub=taskStructure(i).UBequ;
        rxn.rxns=strcat({'TEMPORARY_'},num2str((1:numel(taskStructure(i).equations))'));
        %Allow for new metabolites to be added. This is because it should
        %be possible to add, say, a whole new pathway
        tModel=addRxns(tModel,rxn,3,[],true);
    end
    %Add changed bounds
    if ~isempty(taskStructure(i).changed)
        tModel=setParam(tModel,'lb',taskStructure(i).changed,taskStructure(i).LBrxn);
        tModel=setParam(tModel,'ub',taskStructure(i).changed,taskStructure(i).UBrxn);
    end
    
    %Solve and print
    sol=solveLP(tModel);
    if ~isempty(sol.x)
        if ~taskStructure(i).shouldFail
            taskReport.ok(i,1)=true;
            if printOnlyFailed==false && printOutput==true
                fprintf(['PASS: [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
            end
            %Calculate the essential reactions
            if getEssential==true
                if ~isfield(tModel,'rules')
                    % the getEssentialGenes function requires the model to
                    % have a "rules" field
                    tModel = generateRules(tModel);
                end
                [~, taskEssential]=getEssentialGenes(tModel);
                essentialGenes(taskEssential,i)=true;
            end
        else
            taskReport.ok(i,1)=false;
            if printOutput==true
                fprintf(['PASS (should fail): [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
            end
        end
    else
        if ~taskStructure(i).shouldFail
            taskReport.ok(i,1)=false;
            if printOutput==true
                fprintf(['FAIL: [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
            end
        else
            taskReport.ok(i,1)=true;
            if printOnlyFailed==false && printOutput==true
                fprintf(['FAIL (should fail): [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
            end
        end
    end
    if taskStructure(i).printFluxes && ~isempty(sol.x)
        sol=solveLP(tModel,1);
        if ~isempty(sol.x)
            printFluxes(tModel,sol.x,false,10^-6,[],'%rxnID (%eqn):%flux\n');
            fprintf('\n');
        end
    end
    tModel=model;
end
end
