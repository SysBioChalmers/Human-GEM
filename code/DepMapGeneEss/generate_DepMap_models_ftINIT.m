function [] = generate_DepMap_models_ftINIT(chunk, use1Plus1)
% generate tINIT models from RNA-Seq profiles
%
% Note: This function was designed specifically for use on the cluster
%
% Input:
%
%   chunk     integer ranging from 1 to 10, or 'test' to run a test that
%             generates only one model.
%
%   use1Plus1 if TRUE, run ftINIT with 1+1
%             (Default = FALSE)
%

if nargin < 2
    use1Plus1 = false;
end

%setRavenSolver('gurobi');

pwd()
%% run ftINIT

%we replace the arraydata with the one from depmap
load('./data/arrayDataDepMap.mat')
arrayData = arrayDataDepMap; %no need to add anything


% load inputs generated using "prepare_tINIT_inputs" function
load('./data/prepDataGeneSymbols.mat');

nChunks = 40;
indx = round(linspace(0, numel(arrayData.tissues), nChunks+1))';
model_chunks = arrayfun(@(i) (indx(i)+1:indx(i+1))', (1:nChunks)', 'UniformOutput', false);
if strcmpi(chunk, 'test')
    model_indx = 1;
else
    model_indx = model_chunks{chunk};
end

% extract subset of array data for building models
arrayData.tissues = arrayData.tissues(model_indx);
arrayData.levels = arrayData.levels(:, model_indx);

nModels = length(model_indx);

depmap_models = cell(nModels, 1);

if use1Plus1
   INITSteps = '1+1';
else
   INITSteps = '1+0';
end   

for i = 1:nModels
     if (isempty(depmap_models{i}))
         disp(['running model: ' num2str(i)])
         tic 
         mres = ftINIT(prepData,arrayData.tissues{i},[],[],arrayData,{},getHumanGEMINITSteps(INITSteps),false,true,[]);
         toc
         mres.id = arrayData.tissues{i};
         depmap_models{i,1} = mres;
     end
end

% save results
if isnumeric(chunk)
	if use1Plus1
		filename = strcat('./data/depmap_models_1p1-',num2str(chunk));
	else
		filename = strcat('./data/depmap_models-',num2str(chunk));
	end
else
    filename = 'depmap_models-TEST';
end
save(filename,'depmap_models');

end
