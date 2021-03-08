function [animalGEM, speciesSpecNetwork, gapfillNetwork]=updateAnimalGEM(orthologPairs,rxnsToAdd,metsToAdd,modelId)
% updateAnimalGEM
%   Generate a model by using the Human-GEM as a template and taking into
%   account species-specific pathways/reactions
%
%   Input:
%   orthologPairs        an Nx2 cell array of the ortholog pairs, where the
%                        first column contains gene IDs from the reference
%                        organism, and the second includes gene IDs of the
%                        query organism
%   rxnsToAdd            the structure of species-specific reactions
%   metsToAdd            the structure of species-specific metabolites
%   modelId              model id
%
%
%   Output:
%   animalGEM            an updated animal GEM
%   speciesSpecNetwork   a structure containing species specific reactions
%                        and metabolites
%   gapfillNetwork        a structure containing gap-filled reactions and
%                        metabolites for achieving essential tasks
%
%   NOTE: The detail specification for rxnsToAdd and metsToAdd structures
%         are described in addMetabolicNetwork function.
%
%   Usage: [animalGEM, speciesSpecNetwork, gapfillNetwork]=updateAnimalGEM(orthologPairs,rxnsToAdd,metsToAdd,modelId)
%


% handle input arguments
if nargin < 1
    error('Missing ortholog pairs!');
end

% check orthologPairs and ensure the structure is correct
if ~iscell(orthologPairs) || size(orthologPairs, 2) ~= 2
    error('The data structure of ortholog pairs is incorrect!');
end

if nargin < 4
    modelId = '';
end


%% Load Human-GEM as template

% get model path
[ST, I]=dbstack('-completenames');
modelPath=fileparts(fileparts(ST(I).file));

% Human-GEM has to be included into Matlab path
matFile=fullfile(modelPath,'model','Human-GEM.mat');
ymlFile=fullfile(modelPath,'model','Human-GEM.yml');
if isfile(matFile)
    fprintf('Load Human-GEM from Matlab file.\n');
    load(matFile);
elseif isfile(ymlFile)
    fprintf('Load Human-GEM from Yaml file.\n');
    ihuman = importYaml(ymlFile);
else
    error('ERROR: No model file is found!');
end

% convert gene identifiers from Ensembl ids to gene symbols
[grRules,genes,rxnGeneMat] = translateGrRules(ihuman.grRules,'Name','ENSG');
ihuman.grRules    = grRules;
ihuman.genes      = genes;
ihuman.rxnGeneMat = rxnGeneMat;


%% get animal GEM with updated ortholog pairs and species-specific network

% get ortholog-GEM based on provide ortholog pairs
orthologGEM = getModelFromOrthology(ihuman, orthologPairs);

% integrate species-specific metabolic network
if ~iscell(rxnsToAdd.subSystems{1})
    rxnsToAdd.subSystems = cellfun(@(s) {{s}}, rxnsToAdd.subSystems);
end
[animalGEM, speciesSpecNetwork] = addMetabolicNetwork(orthologGEM, rxnsToAdd, metsToAdd);


%% Gap-filling
[animalGEM, gapfillNetwork]=gapfill4EssentialTasks(animalGEM,ihuman);



%% post-gapfilling procedures

% Use MA reactions identifiers 
rxnAssoc = importTsvFile('reactions.tsv');

% check reaction annotation structure
if ~isequal(rxnAssoc.rxns, ihuman.rxns)
    error('There is inconsistency found in the reaction annotaiton file.');
else
    [hits, ind] = ismember(animalGEM.rxns, rxnAssoc.rxns);
    animalGEM.rxns(hits) = rxnAssoc.rxnMAID(ind(hits));
end
% this replacement would be unnecessary after the implementation of MA ids
% as Human-GEM rxn ids

% add ids
if ~isempty(modelId)
    animalGEM.id = modelId;
end




