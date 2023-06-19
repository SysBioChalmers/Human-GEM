function status = testYamlConversion
% test the functions for yaml import/export see if the conversion process
% changes the model content
%

% Get model path
[ST, I]=dbstack('-completenames');
modelPath=fileparts(fileparts(fileparts(ST(I).file)));

% Import yaml model
ymlFile=fullfile(modelPath,'model','Human-GEM.yml');
model = importYaml(ymlFile, true);

% make sure there is no intermediate Yaml file under the current folder
warning('off', 'MATLAB:DELETE:FileNotFound')
if exist('testYamlConversion.yml','file')
    delete testYamlConversion.yml;
end

% export to yml and then import back
try
    exportYaml(model,'testYamlConversion.yml');
    importedHumanGEM = importYaml('testYamlConversion.yml', true);

    % remove intermediate Yaml file
    delete testYamlConversion.yml;

    % compare the imported model from yaml with the original one
    if ~isequal(model, importedHumanGEM)
        disp('::error::Re-imported model is diffrent from export');
        exit(1);
    end
catch
    disp('::error::There are problems during the conversion import and export');
    exit(1)
end

% model conversion between Matlab and Yaml files is successful
disp('The conversion was successful.')
exit(0)