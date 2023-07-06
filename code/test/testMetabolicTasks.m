function status = testMetabolicTasks(taskType)
% Test for metabolic tasks
%
% Input   
%
%   taskType      "essential" for checking viability tasks, or
%                 "verification" for detecting mass and energy balance
%
% Usage: status = testMetabolicTasks(taskType)
%

if nargin < 1
    disp('::error::Task type is missing!');
    exit(1)
end

% Get model path
[ST, I]=dbstack('-completenames');
modelPath=fileparts(fileparts(fileparts(ST(I).file)));

% Import yaml model
ymlFile=fullfile(modelPath,'model','Human-GEM.yml');
ihuman = importYaml(ymlFile, true);

% parse metabolic tasks
model = addBoundaryMets(ihuman);
if taskType == "essential"
    taskFile=fullfile(modelPath,'data','metabolicTasks','metabolicTasks_Essential.txt');
elseif taskType == "verification"
    taskFile=fullfile(modelPath,'data','metabolicTasks','metabolicTasks_VerifyModel.txt');
else
    warning('::error::Unknown task type is provided.');
    exit(1)
end
verificationTasks = parseTaskList(taskFile);

% evaluate task performance
verificationTaskReport=checkTasks(model,[],false,true,false,verificationTasks);
if all(verificationTaskReport.ok)
    fprintf('Suceeded with %s tasks.\n', taskType)
    status = 1;
else
    warning('::error::Failed in %s tasks.', taskType);
    exit(1)
end
