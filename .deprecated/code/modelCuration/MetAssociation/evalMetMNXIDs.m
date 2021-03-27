function results = evalMetMNXIDs(model,mnx)
%evalMetMNXIDs  Evaluate MNX ID assignments based on met formula and name.
%
% USAGE:
%
%   results = resolveMetIDconflicts(model,mnx)
%
% INPUTS:
%
%   model    Model structure containing ID conflicts to be resolved.
%
%   mnx      MetaNetX Database structure to which model will be compared.
%
% OUTPUTS:
%
%   results  Results structure with the fields listed below:
%                  mets   List of model met IDs. Only mets with an
%                         associated MNX ID are included, and mets
%                         associated with multiple MNX IDs are repeated
%                         that number of times.
%              metNames   List of model met names.
%           metFormulas   List of model met formulas.
%                metIDs   List of MNX IDs associated with each met.
%     FormulaMatchExact   Logical vector indicating whether the met formula
%                         in the model structure ('metFormulas') is an 
%                         EXACT match to the formula in the MNX database
%                         corresponding to the associated MNX ID.
%    FormulaMatchNoProt   Logical vector indicating whether the met formula
%                         matches the MNX formula when protons are removed.
%             NameMatch   Logical vector indicating whether any of the
%                         names associated with met in the model (in
%                         'metNames' or 'metNamesAlt') match with any of
%                         the names corresponding to the associated MNX ID.
%          mismatchMets   List of model mets that were originally mapped to
%                         one or more MNX IDs, but failed to match the
%                         formula or name associated with any of those IDs.
%


% append model metNames field with metNamesAlt field
if isfield(model,'metNamesAlt')
    model.metNames = [model.metNames,model.metNamesAlt];
end

% convert model metMNXIDs field into column vector (flatten cell array)
metIDs = model.metMNXID';
empty_inds = cellfun(@isempty,model.metMNXID);
metIDs = metIDs(~empty_inds');

% now flatten other relevant model met fields to align with metIDs
n = length(model.mets);
mets = arrayfun(@(i) repmat(model.mets(i),sum(~empty_inds(i,:),2),1),[1:n]','UniformOutput',false);
mets = vertcat(mets{:});
metNames = arrayfun(@(i) repmat(model.metNames(i,:),sum(~empty_inds(i,:),2),1),[1:n]','UniformOutput',false);
metNames = vertcat(metNames{:});
metFormulas = arrayfun(@(i) repmat(model.metFormulas(i),sum(~empty_inds(i,:),2),1),[1:n]','UniformOutput',false);
metFormulas = vertcat(metFormulas{:});


% get formulas without protons
metFormulasNoProt = regexprep(metFormulas,'H\d*','');

% initialize results structure
results.mets = mets;
results.metNames = metNames;
results.metFormulas = metFormulas;
results.metIDs = metIDs;
results.FormulaMatchExact = false(size(mets));
results.FormulaMatchNoProt = false(size(mets));
results.NameMatch = false(size(mets));


% retrieve metabolite formula information from mnx structure
[~,mnx_form_inds] = ismember(metIDs,mnx.metMNXID);
mnx_metFormulas = mnx.metFormulas(mnx_form_inds);
mnx_metFormulasNoProt = regexprep(mnx_metFormulas,'H\d*','');  % get formulas without protons

% compare formulas
noFormula = cellfun(@isempty,metFormulas) | cellfun(@isempty,mnx_metFormulas);
results.FormulaMatchExact(~noFormula) = strcmp(metFormulas(~noFormula),mnx_metFormulas(~noFormula));
results.FormulaMatchNoProt(~noFormula) = strcmp(metFormulasNoProt(~noFormula),mnx_metFormulasNoProt(~noFormula));


% get MNXID-name pairs from MNX data structure
mnxID2name = mnx.mnxID2name;

% remove entries that don't match any IDs
ind = ismember(mnxID2name(:,1),unique(metIDs));
mnxID2name(~ind,:) = [];

% make lowercase and ignore special characters in metabolite names
mnxID2name(:,2) = lower(regexprep(mnxID2name(:,2),'[^a-zA-Z0-9]',''));
metNames = lower(regexprep(metNames,'[^a-zA-Z0-9]',''));

% retrieve names corresponding to each of the MNX IDs
mnx_metNames = cellfun(@(id) mnxID2name(ismember(mnxID2name(:,1),id),2),metIDs,'UniformOutput',false);

% compress the metNames cell array into a column vector
metNames = nestCell(metNames,true);

% compare metabolite names
results.NameMatch = arrayfun(@(i) any(ismember(metNames{i},mnx_metNames{i})),[1:length(mets)]');


% determine which metabolites do not match formulas or names among any of
% their associated MNXIDs
metList = model.mets(~empty_inds(:,1));  % ignore mets that had no MNXIDs to begin with
results.mismatchMets = {};
for i = 1:length(metList)
    ind = ismember(mets,metList(i));
    if ~any(results.FormulaMatchExact(ind) | results.FormulaMatchNoProt(ind) | results.NameMatch(ind))
        results.mismatchMets = [results.mismatchMets;metList(i)];
    end
end















