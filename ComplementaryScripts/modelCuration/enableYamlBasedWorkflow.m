%
% FILE NAME:    enableYamlBasedWorkflow.m
% 
% DATE CREATED: 2019-12-17
%     MODIFIED:
% 
% PROGRAMMER:   Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: This script is to adjust the Matlab structure of HumanGEM model
%          so that a Yaml-based workflow could be implemented, as proposed
%          in #27
%


%% Load model and annotation files

% load Human-GEM model
load('HumanGEM.mat');


%% The following fields are to be deprecated

% the information in fields "compOutside", "metMiriams", "rxnMiriams" is 
% either duplicated or unnecessary, so they will be deleted

% remove the three fields
ihuman = rmfield(ihuman, {'compOutside','rxnMiriams','metMiriams'});


%% The following fields need to be adjusted

% In "rxnReferences" field, there are elements containing quotation marks
% that need to be removed

% get the index of elements that contain quotation marks
indRef = find(contains(ihuman.rxnReferences, '"'));
fprintf('A total of %d elements are found with quotation marks.\n', length(indRef));
ihuman.rxnReferences(indRef) = regexprep(ihuman.rxnReferences(indRef),'"','');
%ihuman.rxnReferences(indRef)

% Some elements in "rxnConfidenceScores" field have values in 'NaN', which
% may lead to inconsistency during mat-yml format converstion. Here these
% elements are set to zero

indNaN = isnan(ihuman.rxnConfidenceScores);
ihuman.rxnConfidenceScores(indNaN) = 0;


%% test the conversion

% export to yml and then import back
writeHumanYaml(ihuman,'testYamlConversion.yml');
importedHumanGEM = importHumanYaml('testYamlConversion.yml');

% remove intermediate Yaml file
delete testYamlConversion.yml

% compare the imported model from yaml with the original one
if isequal(ihuman, importedHumanGEM)
    fprintf('The model conversion between Mat and Yaml is successful.\n');
else
    fprintf('There is problems in the conversion between Mat and Yaml files.\n');
end


% export HumanGEM
exportHumanGEM(ihuman,'HumanGEM','../../',{'mat','yml'},false,false);


