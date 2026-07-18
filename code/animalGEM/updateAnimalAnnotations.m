function [rxnAssoc, metAssoc] = updateAnimalAnnotations(animalGEM, speciesSpecRxns, speciesSpecMets);
% updateAnimalAnnotations
%   Update animal GEM annotation by integrating annotiation of Human-GEM
%   and species-specific rxns/mets
%
%   Input:
%   animalGEM            a given animal GEM
%
%   speciesSpecRxns      a structure of species-specific reactions
%   speciesSpecMets      a structure of species-specific metabolites
%
%
%   Output:
%   rxnAssoc             a structure containing reaction annotations
%   metAssoc             a structure containing metabolite annotations
%                       
%   NOTE: The detail specification for speciesSpecRxns and speciesSpecMets
%          structures are described in addMetabolicNetwork function.
%         
%   Usage: [rxnAssoc, metAssoc] = updateAnimalAnnotations(animalGEM, speciesSpecRxns, speciesSpecMets)
%


% handle input arguments
if nargin < 3 || nargin < 2
    error('Missing species-specific annotations!');
end

% checks


%% Load Human-GEM annotation files

% get model path
[ST, I]=dbstack('-completenames');
modelPath=fileparts(fileparts(ST(I).file));

% Human-GEM has to be included into Matlab path
rxnAnnotation=fullfile(modelPath,'model','reactions.tsv');
metAnnotation=fullfile(modelPath,'model','metabolites.tsv');
if isfile(rxnAnnotation) && isfile(metAnnotation)
    fprintf('Load Human-GEM annotation files.\n');
    rxnAssocHuman = importTsvFile(rxnAnnotation);
    metAssocHuman = importTsvFile(metAnnotation);
else
    error('ERROR: No annotation file(s) are found!');
end


%% integrate human and species-specific annotations

% append species-specific annotiaons
rxnAssoc = integrateAnnotation(rxnAssocHuman, speciesSpecRxns, 'rxns');
metAssoc = integrateAnnotation(metAssocHuman, speciesSpecMets, 'mets');

% remove non-existing elements and order by the fields in model
rxnAssoc = sortCleanArray(rxnAssoc, animalGEM.rxns, 'rxns');
metAssoc = sortCleanArray(metAssoc, animalGEM.mets, 'mets');


end


%% sub-functions
% integrate human and species-specific annotations
function joinedArray= integrateAnnotation(baseArray, uniqueArray, fieldName);

    joinedArray = baseArray;  % human annotation data
    nNewElements = length(uniqueArray.(fieldName));
    baseArrayFields = fieldnames(baseArray);
    uniqueArrayFields = fieldnames(uniqueArray);

    % combine annotation information
    for i = 1:length(baseArrayFields)
        n = baseArrayFields{i};
        if ismember(n,intersect(baseArrayFields, uniqueArrayFields))
            joinedArray.(n) = [baseArray.(n); uniqueArray.(n)];
        else
            if isnumeric(joinedArray.(n)(1))
                joinedArray.(n) = [baseArray.(n); zeros(nNewElements,1)];
            else
                joinedArray.(n) = [baseArray.(n); repmat({''},nNewElements,1)];
            end
        end
    end

end

% sort and clean annotation cell array according to model rxns/mets
function sortedArray= sortCleanArray(inputArray, inputCell, fieldName);

% model components (rxns or mets) must be included in the inputArray
if all(ismember(inputCell, inputArray.(fieldName)))
    idToRemove = setdiff(inputArray.(fieldName), inputCell);
    
    % identifiy human rxns/mets that are not in animal GEM
    if ~isempty(idToRemove)
        [~, indexToRemove] = ismember(idToRemove, inputArray.(fieldName));
        inputArray.(fieldName)(indexToRemove) = [];  %remove non-existing elements
    end
    
    % get the index for aligning annoation array with model rxns/mets
    [~, sortedIndex] = ismember(inputCell, inputArray.(fieldName));
    sortedArray = {};
    sortedArray.(fieldName) = inputArray.(fieldName)(sortedIndex);
    
    arrayFields = setdiff(fieldnames(inputArray), fieldName);
    for i = 1:length(arrayFields)
        n = arrayFields{i};
        if ~isempty(idToRemove)
	        inputArray.(n)(indexToRemove) = [];  %remove non-existing elements
        end
        sortedArray.(n) = inputArray.(n)(sortedIndex);
    end
else
    error('There are model components that are from unknown source.'); 
end

end
