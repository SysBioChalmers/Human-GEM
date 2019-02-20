function [model] = updateGrRules(fileName,nHeaderLines,colNewGrRules,autoSave)
% updateGrRules
%   Update specific grRules curation results into humanGEM model. Other
%   modified fields include genes and rxnGeneMat
%
% Input:
%   fileName       name of the file that has curated grRules information
%
%   nHeaderLines   the number of heading lines including column names row
%
%   colNewGrRules  the column number of curated grRules
%
%   autoSave       if TRUE, save changes to .mat model file (opt, default FALSE)
%                  
%
% Output:
%   model          an updated model structure
%
% NOTE: this input file with curated grRules should follow defined rules:
% i) it has to a tab delimitted plaintext file placed under subfolder
% "~/ComplementaryData/modelCuration/"; ii) has five columns for storing
% the curation information; iii) the first column must be a unique list of
% rxn ids whose grRules are to be changed; iv) the third column includes
% the newly curated grRules, perferably (otherwise its column number need
% to be specified in the argument colNewGrRules)
%
% Usage: [model] = updateGrRules(fileName,nHeaderLines,colNewGrRules,autoSave)
%
% Hao Wang, 2019-02-20
%

% handel input
if nargin < 4
    autoSave = false;
end
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
model = ihuman;

% Save changes to .mat model file
if autoSave
    save('../../ModelFiles/mat/humanGEM.mat','ihuman');
end

end
