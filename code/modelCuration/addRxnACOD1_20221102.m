% load model and new reaction info
ihuman = importYaml('../../model/Human-GEM.yml');
rxnsToAdd = importTsvFile('../../data/modelCuration/addRxnACOD1_20221102.tsv');

% add new genes to Human-GEM
newGEM = ihuman;
newGEM.genes = [newGEM.genes; rxnsToAdd.grRules];
newGEM.rxnGeneMat(:, end+1) = 0;

% reformat subsystem 
if ~iscell(rxnsToAdd.subSystems{1})
    rxnsToAdd.subSystems = cellfun(@(s) {{s}}, rxnsToAdd.subSystems);
end
% add new reaction
newGEM = addRxns(newGEM, rxnsToAdd, 3);

% add reaction annotation
rxnAssoc = importTsvFile('../../model/reactions.tsv');
annoNames = fieldnames(rxnAssoc);
for i=1:length(annoNames)
    if ismember(annoNames{i}, fieldnames(rxnsToAdd))
        rxnAssoc.(annoNames{i}) = [rxnAssoc.(annoNames{i}); rxnsToAdd.(annoNames{i})];
    elseif isequal(annoNames{i},'spontaneous')
        rxnAssoc.(annoNames{i}) = [rxnAssoc.(annoNames{i}); 0];
    else
        rxnAssoc.(annoNames{i}) = [rxnAssoc.(annoNames{i}); {''}];
    end
end

% update yaml model and reaction association file
exportYaml(newGEM, '../../model/Human-GEM.yml');
exportTsvFile(rxnAssoc,'../../model/reactions.tsv');

