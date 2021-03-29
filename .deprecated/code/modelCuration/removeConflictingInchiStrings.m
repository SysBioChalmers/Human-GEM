function corrected_model = removeConflictingInchiStrings(model)
% Removes inchi strings that conflict with their corresponding metFormula
%
% Input:
%
%   model   model structure
%
%
% Output:
%
%   corrected_model   model structure with "inchis" field updated such that
%                     all entries conflicting with the corresponding entry
%                     in the "metFormulas" field have been removed.
%
%
% Usage:
%
%   corrected_model = removeConflictingInchiStrings(model);
%


if ~isfield(model,'inchis')
    error('Could not find field named "inchis" in the model (case-sensitive).');
end
    
% identify non-empty inchis entries
ind = find(~cellfun(@isempty, model.inchis));

% extract formulas from inchis strings, which are composed of segments
% separated by delimiter (/) and formula is from the 2nd segment
inchiSplit = cellfun(@(i) strsplit(i,'/'), model.inchis(ind), 'UniformOutput', false);
inchiFormulas = cellfun(@(i) i{2}, inchiSplit, 'UniformOutput', false);
inchiFormulas = regexprep(inchiFormulas,'[^a-zA-Z0-9]','');  % remove special characters

% extract elemental composition of formulas
metFormulas = model.metFormulas(ind);
combinedFormulas = [inchiFormulas; metFormulas];
[~,elMat] = parseFormulas(combinedFormulas);

% compare formula composition and remove inchi entries that disagree
inchiMat = elMat(1:numel(ind),:);
formulaMat = elMat(numel(ind)+1:end,:);
mismatch = any((inchiMat - formulaMat) ~= 0, 2);
model.inchis(ind(mismatch)) = {''};

if any(mismatch)
    fprintf('Removed %u conflicting InChI strings.\n', sum(mismatch));
else
    fprintf('No InChI conflicts were identified.\n');
end

corrected_model = model;


