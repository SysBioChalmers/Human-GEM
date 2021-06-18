%
% FILE NAME:    implementeMAID.m
% 
% PURPOSE: This script is to implement MA reaction and metabolite ids to
%          Human-GEM, as planned in 265#.
%


%% Load model and tsv annotaiton files

% load HumanGEM
ihuman = importYaml('Human-GEM.yml');
ihuman_orig = ihuman;  % to track changes

% load tsv files
rxnAssoc = importTsvFile('reactions.tsv');
metAssoc = importTsvFile('metabolites.tsv');

rxnAssoc_orig = rxnAssoc;
metAssoc_orig = metAssoc;


%% swap id columns

% sanity check and remove dash from MA ids
if isequal(rxnAssoc.rxns, ihuman.rxns) && isequal(metAssoc.mets, ihuman.mets)
    fprintf('sanity check passed and move on.\n');
    metMAID = replace(metAssoc.metMAID, 'MA-M', 'MAM');
    rxnMAID = replace(rxnAssoc.rxnMAID, 'MA-R', 'MAR');
end

% back up retired ids
rxnAssoc.rxnRetired = rxnAssoc.rxns;
metAssoc.metRetired = metAssoc.mets;

% implement new ids
ihuman.rxns   = rxnMAID;
rxnAssoc.rxns = rxnMAID;

ihuman.mets   = metMAID;
metAssoc.mets = metMAID;

% regenerate metsNoComp column by striping off the tailing compartment id
% from mets column
metAssoc.metsNoComp = cellfun(@(x) regexprep(x, '^(.+)\w$', '$1'),...
    metAssoc.mets, 'UniformOutput', false);


%% sanity checks

if isequal(rxnAssoc.rxns, ihuman.rxns) && isequal(metAssoc.mets, ihuman.mets)
    fprintf('sanity check passed.\n');
end

compareArrayStructure(ihuman, ihuman_orig);
compareArrayStructure(rxnAssoc, rxnAssoc_orig);
compareArrayStructure(metAssoc, metAssoc_orig);
% everything looks good


%% save model and annotaiton files
exportYaml(ihuman,'../../../model/Human-GEM.yml');

rxnAssoc = rmfield(rxnAssoc, {'rxnMAID'});
exportTsvFile(rxnAssoc, '../../../model/reactions.tsv');

metAssoc = rmfield(metAssoc, {'metMAID'});
exportTsvFile(metAssoc, '../../../model/metabolites.tsv');

