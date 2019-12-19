function status = testYamlConversion()
%test the functions for yaml import/export if the conversion process changes
%the model content


% load HumanGEM
load('HumanGEM.mat');

warning('off', 'MATLAB:DELETE:FileNotFound')
if exist('testYamlConversion.yml','file')
    delete testYamlConversion.yml;
end


% export to yml and then import back
writeHumanYaml(ihuman,'testYamlConversion.yml');
importedHumanGEM = importHumanYaml('testYamlConversion.yml');

% remove intermediate Yaml file
delete testYamlConversion.yml;


% compare the imported model from yaml with the original one
if isequal(ihuman, importedHumanGEM)
    %model conversion between Matlab and Yaml files is successful
    status = 1;
else
    %There is problems during the conversion between Matlab and Yaml files
    status = 0;
end

