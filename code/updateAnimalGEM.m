function [animalGEM, speciesSpecNetwork, gapfillNetwork]=updateAnimalGEM(orthologPairs,rxnsToAdd,metsToAdd,modelId,resetBiomass)
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
%   resetBiomass         reset biomass objective function to "biomass_components"
%                        which is constituted by generic componenets that
%                        suppose to occur in a eukaryotic cell (opt, default TRUE)
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
    % load Human-GEM Matlab file
    load(matFile);
elseif isfile(ymlFile)
    % Load Human-GEM Yaml file
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
[animalGEM, gapfillNetwork]=gapfill4EssentialTasks(animalGEM,ihuman,resetBiomass);
animalGEM.b = animalGEM.b(:,1);   % ensure b field in single column


%% post-gapfilling procedures

% add ids
if ~isempty(modelId)
    animalGEM.id = modelId;
end




