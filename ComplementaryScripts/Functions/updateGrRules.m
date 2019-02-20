function [model] = updateGrRules(fileName, nHeaderLines, colNewGrRules)
% updateGrRules
%   Update specific grRules curation results into humanGEM model. Other
%   modified fields include genes and rxnGeneMat
%
% Input:
%   fileName       the file with curated grRules curation
%   nHeaderLines   the number of heading lines including column names row
%   colNewGrRules  the column number of curated grRules
%
% Output:
%   model          an updated model structure
%
% NOTE: this input file with curated grRules should follow a defined format:
% i) a plaintext file in tab delimitted format; ii) has five columns
% for storing the curation information; iii) the first column must be a
% unique list of rxn ids whose grRules are to be changed; iv) the third column
% has the newly curated grRules perferably (otherwise its column number need
% to be specified in the argument colNewGrRules)
%
% Usage: [model] = updateGrRules(fileName, nHeaderLines, targetColumn)
%
% Hao Wang, 2019-02-20
%

% handel input
if nargin < 3
    colNewGrRules = 3;
end
if nargin < 2
    nHeaderLines = 1;
end

% Load the grRules curation results from input file
inputFile=fullfile('../../ComplementaryData/modelCuration/',fileName);
if ~exist(inputFile,'file')
    error('The file with curated grRules cannot be located. Please specify the correct filename.\n');
else
    fid = fopen(inputFile,'r');
    tmp = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',nHeaderLines);
    fclose(fid);
end

% Get the rxn ids and new grRules
rxnIDs = tmp{1};
newGrRules = tmp{colNewGrRules};


% Load humanGEM model
load('humanGEM.mat');


% Update curated grRules
[~,rxn_ind] = ismember(rxnIDs,ihuman.rxns);
ihuman.grRules(rxn_ind) = newGrRules;


% Update other fields and save the model
[genes,rxnGeneMat] = getGenesFromGrRules(ihuman.grRules);
ihuman.genes = genes;
ihuman.rxnGeneMat = rxnGeneMat;
save('../../ModelFiles/mat/humanGEM.mat','ihuman');

model = ihuman;

end

