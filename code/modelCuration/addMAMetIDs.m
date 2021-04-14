

% load data
metIDs = importTsvFile('metabolites.tsv');

% get compartment id for non-HMR mets that will be used later
compID = cellfun(@(x) regexprep(x, '^.+_(\w)$', '$1'), metIDs.mets, 'UniformOutput', false);
matche2NonHMR = regexp(metIDs.mets, '_\w$', 'match');
emptymask = cellfun('isempty', matche2NonHMR);
compID(emptymask) = {''};    % empty HMR met ids

% generate MA met ids from HMR ids only with "MA-" prefix
metMAID = cellfun(@(x) regexprep(x, '^m(\d\d\d\d\d\w)$', 'MA-M$1'), metIDs.mets, 'UniformOutput', false);

% check consistency between HMR ids and empty elemnets of non-HMR ids
matche2HMR = regexp(metMAID, '^MA-M\d\d\d\d\d\w$', 'match');
ind2HMR = ~cellfun('isempty', matche2HMR);
check = find(ind2HMR ~= emptymask);
metIDs.mets(check)
%    {'temp001c'}
%    {'temp001s'}
% these two mets are from HMR but in non-standard format, will be treated as
% non-HMR ids
compID(check) = {'c';'s'};

% kepp HMR ids in metMAID and clean others
ind_MAID = startsWith(metMAID, 'MA-M');
metMAID(~ind_MAID) = {''};
metMAIDNoComp = cellfun(@(x) regexprep(x, '^(.+)\w$', '$1'), metMAID, 'UniformOutput', false);

% get index of nonHMR ids from compID and remove their comppart id
nonHMRidInd = ~cellfun('isempty', compID);
nonHMRidNoComp = cellfun(@(x) regexprep(x, '^(.+)_\w$', '$1'), metIDs.mets(nonHMRidInd), 'UniformOutput', false);

% get unique list of nonHMR ids
[uniqueID, ia, ic] = unique(nonHMRidNoComp);
%isequal(nonHMRidNoComp, uniqueID(ic))
%isequal(uniqueID, nonHMRidNoComp(ia))

% prepare standard format MA ids for nonHMR mets
idsToAdd = cell(10000,1);
for i=1:10000
    idsToAdd{i} = strcat('MA-M',sprintf('%05d', i));
end
idsToAdd = setdiff(idsToAdd, metMAIDNoComp);     % exclude existing ones
idsToAdd = idsToAdd(1:length(uniqueID));

% get ids to fill in
idsToFill = idsToAdd(ic);

% append compartment id to idsToFill
idsToFill = strcat(idsToFill,compID(nonHMRidInd));

% complete MA met ids
metMAID(nonHMRidInd) = idsToFill;

% format check
all(cell2mat(regexp(metMAID, '^MA-M\d\d\d\d\d\w$')))
% unique check
isequal(length(metMAID), length(sort(metMAID)))
% empty check
B = cellfun('isempty', metMAID);
all(B(:) == 0)

% add new field and output
metIDs.metMAID = metMAID;
exportTsvFile(metIDs,'metabolites.tsv');



