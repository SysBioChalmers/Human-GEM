% this script is to provide complete reaction identifiers (#174) that will be used for Met Atlas web-portal.

% load reaction ids
rxnIDs = importTsvFile('reactions.tsv');

% replace underscore with `R`
rxnIDs.rxnMAID = strrep(rxnIDs.rxnMAID,'_','R');
indToFill = getNonEmptyList(rxnIDs.rxnMAID,false);

% prepare a list to fill
idsToAdd = cell(19999,1);
for i=1:19999
    idsToAdd{i} = strcat('MAR',sprintf('%05d', i));
end
idsToAdd = setdiff(idsToAdd, rxnIDs.rxnMAID);

% sequentially fill in blank MA ids
rxnIDs.rxnMAID(indToFill) = idsToAdd(1:length(indToFill));
if isequal(length(rxnIDs.rxnMAID), length(unique(rxnIDs.rxnMAID)))
    fprintf('the filling looks okay.\n');
end

% insert dash between MA and RDDDDD
rxnIDs.rxnMAID = strrep(rxnIDs.rxnMAID,'MA','MA-');


exportTsvFile(rxnIDs,'reactions.tsv'); 

