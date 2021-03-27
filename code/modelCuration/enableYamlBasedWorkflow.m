%
% FILE NAME:    enableYamlBasedWorkflow.m
% 
% PURPOSE: This script is to adjust the Matlab structure of HumanGEM model
%          so that a Yaml-based workflow could be implemented, as proposed
%          in #27
%


%% Load model and annotation files

% load Human-GEM model
load('HumanGEM.mat');


%% The following fields are to be deprecated

% the "compOutside" field stores static information that can be deduced
% from common sense, no need to keep in the model


% check "rxnMiriams" field
% get index of non-empty elements
index = find(~cellfun(@isempty, ihuman.rxnMiriams));
rxnMiriams = ihuman.rxnMiriams(index);

% check the rxnMiriams types
allNames_rxnMiriams  = (cellfun(@(x)x.name, rxnMiriams,'UniformOutput', 0));
rxnMiriamsNames = vertcat(allNames_rxnMiriams{:});
unique(rxnMiriamsNames)
% only one type "pmid" was found in "rxnMiriams" field, its informaiton is
% duplicated to the "rxnReferences" field


% check "metMiriams" field
% get index of non-empty elements
ind = find(~cellfun(@isempty, ihuman.metMiriams));
metMiriams = ihuman.metMiriams(ind);

% check the metMiriams types
allNames_metMiriams  = (cellfun(@(x)x.name, metMiriams,'UniformOutput', 0));
metMiriamsNames = vertcat(allNames_metMiriams{:});
unique(metMiriamsNames)
% hmdb
% kegg.compound
% kegg.glycan
% lipidmaps
% obo.chebi:CHEBI
% obo.chebi:CHEBi
% obo.chebi:ChEBI
% there are four types (HMDB, KEGG, LipidMaps, ChEBI) of identifiers are
% found in "metMiriams" field, this information has been continuesly
% curated and stored in annotation files: "humanGEMRxnAssoc.JSON" and
% "humanGEMMetAssoc.JSON"


% since the information in fields "compOutside", "metMiriams", "rxnMiriams"
% is either duplicated or unnecessary, so they will be deleted
ihuman = rmfield(ihuman, {'compOutside','rxnMiriams','metMiriams'});


%% The following fields need to be adjusted

% In "rxnReferences" field, there are elements containing quotation marks
% that need to be removed

% get the index of elements that contain quotation marks
indRef = find(contains(ihuman.rxnReferences, '"'));
fprintf('A total of %d elements are found with quotation marks.\n', length(indRef));
ihuman.rxnReferences(indRef)
ihuman.rxnReferences(indRef) = regexprep(ihuman.rxnReferences(indRef),'"','');


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


