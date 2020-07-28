function status = testYamlConversion
% test the functions for yaml import/export see if the conversion process
% changes the model content
%
% Usage: status = testYamlConversion
%


% Get model path
[ST, I]=dbstack('-completenames');
modelPath=fileparts(fileparts(fileparts(ST(I).file)));

% Import yaml model
ymlFile=fullfile(modelPath,'model','Human-GEM.yml');
model = importHumanYaml(ymlFile, true);

% make sure there is no intermediate Yaml file under the current folder
warning('off', 'MATLAB:DELETE:FileNotFound')
if exist('testYamlConversion.yml','file')
    delete testYamlConversion.yml;
end


% export to yml and then import back
writeHumanYaml(model,'testYamlConversion.yml');
importedHumanGEM = importHumanYaml('testYamlConversion.yml', true);

% remove intermediate Yaml file
delete testYamlConversion.yml;

% compare the imported model from yaml with the original one
if isequal(model, importedHumanGEM)
    % model conversion between Matlab and Yaml files is successful
    status = 1;
else
    error('There are problems during the conversion between Matlab and Yaml files');
end

