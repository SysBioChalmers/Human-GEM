%
% FILE NAME:    miscModelCurationScript_20190323.m
%
% PURPOSES: This script conducts a number of tasks: 
%           1. Add Metadata to the annotation field of Human1
%           2. Reformat EC-number in eccodes field
%           3. Remove rxnComps field
%           4. Turn `version` into a blank field
%           5. Initialize rxnConfidenceScores field with zero
%


%% Load the model

if ~exist('ihuman','var')
    load('humanGEM.mat');  % version 1.0.0
end

%% Add Metadata to the annotation field of Human1

% Generate the Metadata as a structure
annotation.defaultLB=-1000;
annotation.defaultUB=1000;
annotation.taxonomy='9606';
annotation.note='Human genome-scale metabolic models are important tools for the study of human health and diseases, by providing a scaffold upon which different types of data can be analyzed. This is the latest version of human-GEM, which is a genome-scale model of the generic human cell. The objective of human-GEM is to serve as a community model for enabling integrative and mechanistic studies of human metabolism.';
annotation.sourceUrl='https://github.com/SysBioChalmers/human-GEM';
annotation.authorList='Jonathan Robinson, Hao Wang, Pierre-Etienne Cholley, Pınar Kocabaş';
annotation.email='nielsenj@chalmers.se';
annotation.organization='Chalmers University of Technology';

ihuman.annotation = annotation;
ihuman.description = 'Generic genome-scale metabolic model of Homo sapiens';


%% Reformat EC-number

% Multiple EC numbers should be separated by a semicolon and a blank space
% "; " according to RAVEN issue #184. Some incorrectly formatted eccodes
% elements, which are separated by " or " as reported in #93, are fixed here.
eccodes = regexprep(ihuman.eccodes,'\s*or\s*',';');

% Consistengly add a blank space after each semicolon
eccodes = regexprep(eccodes,';','; ');

% Update to the model
ihuman.eccodes = eccodes;


%% Emtpy version field

% Turn `version` into a blank field to retain a simple and clear work flow
ihuman.version = '';


%% Remove rxnComps field

% Take away the rxnComps field, according to RAVEN #184
ihuman = rmfield(ihuman, 'rxnComps');


%% Initialize rxnConfidenceScores field

% Assign elements in rxnConfidenceScores field with 0, based on #48
ihuman.rxnConfidenceScores(:) = 0;


%% Save the model files

writeHumanYaml(ihuman, 'humanGEM');
movefile('Human-GEM.yml','../../model/');
save('../../model/Human-GEM.mat','ihuman');
