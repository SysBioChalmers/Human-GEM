

% load Human-GEM
model = importYaml('../../model/Human-GEM.yml');


% specify pairs of duplicate reactions, where reactions in first column are
% those that will be kept and those in the second column removed
rxns = {
    'MAR08971','MAR04306'
    'MAR00121','MAR00607'
    'MAR04411','MAR04911'
    'MAR00079','MAR06324'
    'MAR00101','MAR02199'
    'MAR02339','MAR00445'
    'MAR01728','MAR07441'
    'MAR01827','MAR07742'
    'MAR11237','MAR04455'
    'MAR11251','MAR04463'
    'MAR11252','MAR04478'
    'MAR11253','MAR04502'
    'MAR11254','MAR04517'
    'MAR11247','MAR04576'
    'MAR11249','MAR04609'
    'MAR11248','MAR04622'
    'MAR11255','MAR04626'
    'MAR11256','MAR04653'
    'MAR11257','MAR04659'
    'MAR11258','MAR04669'
    'MAR11260','MAR04713'
    'MAR11261','MAR04724'
    'MAR11262','MAR04747'
    'MAR11263','MAR04751'
    'MAR11264','MAR04753'
    'MAR11265','MAR04761'
    'MAR11238','MAR04794'
    'MAR11266','MAR04795'
    'MAR11267','MAR04798'
    'MAR11239','MAR04800'
    'MAR11268','MAR04801'
    'MAR11240','MAR06494'
    'MAR11241','MAR06497'
    'MAR11242','MAR06504'
    'MAR11243','MAR06541'
    'MAR11244','MAR06582'
    'MAR11246','MAR06590'
    'MAR11271','MAR06594'
    'MAR11245','MAR06706'
    'MAR11280','MAR06853'
    'MAR11278','MAR06856'
    'MAR11279','MAR06857'
    'MAR11259','MAR06858'
    'MAR09057','MAR07191'};

% verify that the reaction pairs have identical stoichiometry (either forward or reverse)
for i = 1:size(rxns, 1)
    indx = getIndexes(model, rxns(i,:), 'rxns');
    coeffs = full(model.S(any(model.S(:, indx) ~= 0, 2), indx));
    coeffs_norm = coeffs ./ max(abs(coeffs), [], 1);
    if ~isequal(coeffs_norm(:,1), coeffs_norm(:,2)) && ...
            ~isequal(-coeffs_norm(:,1), coeffs_norm(:,2))
        error('Reactions in row %u are not duplicates!', i);
    end
end

% normalize the stoichiometric coefficients of MAR08971
%  currently: 3 glucose-6-phosphate[c] + 3 NADP+[c] => 3 glucono-1,5-lactone-6-phosphate[c] + 3 H+[c] + 3 NADPH[c]
% normalized:   glucose-6-phosphate[c] +   NADP+[c] =>   glucono-1,5-lactone-6-phosphate[c] +   H+[c] +   NADPH[c]
rxn_indx = getIndexes(model, 'MAR08971', 'rxns');
model.S(:, rxn_indx) = model.S(:, rxn_indx) / 3;

% merge rxn annotations from the deleted duplicate rxn into the kept rxn
rxnAssocFile = '../../model/reactions.tsv';
rxnAssoc = importTsvFile(rxnAssocFile);
[~, rxn_indx] = ismember(rxns, rxnAssoc.rxns);
mergeFields = setdiff(fieldnames(rxnAssoc), {'rxns', 'rxnRetired', 'spontaneous'});
for i = 1:size(rxns, 1)
    for f = 1:numel(mergeFields)
        entry1 = strsplit(rxnAssoc.(mergeFields{f}){rxn_indx(i,1)}, ';');
        entry2 = strsplit(rxnAssoc.(mergeFields{f}){rxn_indx(i,2)}, ';');
        merged_entry = setdiff(union(entry1, entry2), {''});
        if ~isempty(merged_entry)
            rxnAssoc.(mergeFields{f}){rxn_indx(i,1)} = strjoin(merged_entry, ';');
        end
    end
end

% save annotations for to-be-deleted reactions to the deprecated rxns file
depRxnFile = '../../data/deprecatedIdentifiers/deprecatedReactions.tsv';
deprecRxns = importTsvFile(depRxnFile);
deprecRxnsTable = struct2table(deprecRxns);
rxnAssocTable = struct2table(rxnAssoc);
exportTsvFile([deprecRxnsTable; rxnAssocTable(rxn_indx(:,2), :)], depRxnFile);

% merge additional model fields from duplicate rxns with the kept rxns
mergeFields = {'eccodes'; 'rxnReferences'};
for i = 1:size(rxns, 1)
    for f = 1:numel(mergeFields)
        entry1 = strsplit(model.(mergeFields{f}){rxn_indx(i,1)}, ';');
        entry2 = strsplit(model.(mergeFields{f}){rxn_indx(i,2)}, ';');
        merged_entry = setdiff(union(entry1, entry2), {''});
        if ~isempty(merged_entry)
            model.(mergeFields{f}){rxn_indx(i,1)} = strjoin(merged_entry, ';');
        end
    end
end

% delete reactions from model and annotation file
model = removeReactions(model, rxns(:,2));
exportYaml(model, '../../model/Human-GEM.yml');
rxnAssocTable(rxn_indx(:,2), :) = [];
exportTsvFile(rxnAssocTable, rxnAssocFile);


