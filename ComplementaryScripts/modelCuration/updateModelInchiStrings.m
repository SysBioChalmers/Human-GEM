function updated_model = updateModelInchiStrings(model)
% Removes inchi strings that conflict with their corresponding metFormula
%
% Input:
%
%   model           model structure
%
%
% Output:
%
%   updated_model   model structure with updated "inchis" field
%
%
% Usage:
%
%   updated_model = updateModelInchiStrings(model);
%
%
% Jonathan Robinson, 2019-10-28


if ~isfield(model,'inchis')
    error('Could not find field named "inchis" in the model (case-sensitive).');
end
    
% identify non-empty inchis entries and extract formulas
ind = find(~cellfun(@isempty, model.inchis));
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

updated_model = model;


