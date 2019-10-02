function status = testYamlConversion()
%test the functions for yaml import/export if the conversion process changes
%the model content


% load HumanGEM
load('humanGEM.mat');

if exist('testYamlConversion.yml','file')==2
    delete testYamlConversion.yml
end

% export to yml and then import back
writeHumanYaml(ihuman,'testYamlConversion.yml');
importedHumanGEM = importHumanYaml('testYamlConversion.yml');

% remove intermediate Yaml file
delete testYamlConversion.yml

%===this section will be removed after resolving the temporary fields
% extract shared fields
ihumanFields = fieldnames(ihuman);
importedModelFields = fieldnames(importedHumanGEM);
sharedFields = intersect(ihumanFields, importedModelFields);

% trim off unique fields
HumanGEM = rmfield(ihuman, setdiff(ihumanFields, sharedFields));
importedHumanGEM = rmfield(importedHumanGEM, setdiff(importedModelFields, sharedFields)); 
%========

% compare the imported model from yaml with the original one
if isequal(HumanGEM, importedHumanGEM)
    %model conversion between Matlab and Yaml files is successful
    status = 1;
else
    %There is problems during the conversion between Matlab and Yaml files
    status = 0;
end

