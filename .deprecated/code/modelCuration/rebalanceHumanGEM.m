%
% FILE NAME:    rebalanceHumanGEM.m
% 
% PURPOSE: This script curates the formulas and annotation of many
%          metabolites and the equations and annotation of many reactions
%          to achieve full stoichiomietric consistency and nearly complete
%          mass and charge balance.
%


%% Load Human-GEM and annotation structures

% load Human-GEM model
load('HumanGEM.mat');
ihuman_orig = ihuman;  % keep copy of original version

% load metabolite and reaction annotation data
metAssoc = jsondecode(fileread('../../ComplementaryData/annotation/humanGEMMetAssoc.JSON'));
rxnAssoc = jsondecode(fileread('../../ComplementaryData/annotation/humanGEMRxnAssoc.JSON'));
changeNotes = {};

% verify that HumanGEM and annotation structures are aligned
if ~isequal(ihuman.mets, metAssoc.mets) || ~isequal(ihuman.rxns, rxnAssoc.rxns)
    error('HumanGEM is not synced with the metAssoc and/or rxnAssoc structure!');
end


%% Update rev and lb fields to be consistent

% There are 11 reactions with rev=1 but lb=0, these reactions should be
% updated so that their rev=0. This change to reversibility was made to
% these 11 reactions in curateATPmetabolism (implemented in PR#113), but
% their "rev" field was not properly updated at that time.
ihuman.rev(ihuman.lb == 0 & ihuman.ub ~= 0) = 0;


%% Add new field to rxnAssoc structure

% create a new field for HMR reaction IDs in the rxnAssoc structure
if ~isfield(rxnAssoc,'rxnHMRID')
    % include all IDs of the form "HMR_####". Note: ignore "HMR_#####"
    % reactions (5 numbers), as these are new to HumanGEM.
    rxnAssoc.rxnHMRID = repmat({''},size(rxnAssoc.rxns));
    hmr_rxns = startsWith(rxnAssoc.rxns,'HMR_') & (cellfun(@numel,rxnAssoc.rxns) == 8);
    rxnAssoc.rxnHMRID(hmr_rxns) = rxnAssoc.rxns(hmr_rxns);
end


%% Update four metabolite names to avoid parsing errors
% Four metabolite names begin with a number or number with commas, followed
% by a space, which some functions can confuse with stoichiometric
% coefficients when parsing reaction equations. To fix this, the space is
% replaced with a dash (-).
nameArray = {'1 Acyl Phosphoglycerol', '1-Acyl Phosphoglycerol'
             '15, 31-O-Didesmethyl-tacrolimus', '15,31-O-Didesmethyl-tacrolimus'
             '2,6 Dimethylheptanoyl Coenzyme A', '2,6-Dimethylheptanoyl Coenzyme A'
             '4,8 Dimethylnonanoyl Coenzyme A', '4,8-Dimethylnonanoyl Coenzyme A'};
[hasMatch,nameInd] = ismember(ihuman.metNames, nameArray(:,1));
if any(hasMatch)
    ihuman.metNames(hasMatch) = nameArray(nameInd(hasMatch),2);
    changeNotes = [changeNotes; [ihuman.mets(hasMatch),...
    repmat({'Added dash to name to avoid confusion with stoich coeffs'},sum(hasMatch),1)]];
end


%% Add new metabolites to the model

% load new metabolite information from file
fid = fopen('../../ComplementaryData/modelCuration/fullRebalance/rebalance_mets_new.tsv');
metData = textscan(fid,'%s%s%s%s%d%s%s%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% verify that none of the metabolite IDs (ignoring compartment) exist in the current model
if any(startsWith(ihuman.mets, regexprep(metData{1},'.$','')))
    error('One or more metabolite IDs to be added already exist in the model.');
end

% add new metabolites
metsToAdd = {};
metsToAdd.mets = metData{1};
metsToAdd.metNames = metData{2};
metsToAdd.compartments =  metData{3};
metsToAdd.metFormulas =  metData{4};
metsToAdd.metCharges =  metData{5};
ihuman = addMets(ihuman, metsToAdd);

% add metabolites to the metAssoc structure
numOrigMets = numel(metAssoc.mets);
numNewMets = numel(metData{1});
newMetInd = (numOrigMets+1:numOrigMets+numNewMets)';
f = fieldnames(metAssoc);
for i = 1:numel(f)
    % initialize structure with empty entries
    metAssoc.(f{i})(newMetInd) = {''};
end
metAssoc.mets(newMetInd) = metData{1};
metAssoc.metsNoComp(newMetInd) = regexprep(metData{1},'.$','');
metAssoc.metKEGGID(newMetInd) = metData{6};
metAssoc.metBiGGID(newMetInd) = metData{7};
metAssoc.metChEBIID(newMetInd) = metData{8};
metAssoc.metMetaNetXID(newMetInd) = metData{9};


%% Add existing metabolites to new compartments

% load metabolite information from file
fid = fopen('../../ComplementaryData/modelCuration/fullRebalance/rebalance_mets_newCompVersions.tsv');
metData = textscan(fid,'%s%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% verify that none of the metabolites exist in the current model
if any(ismember(metData{1},ihuman.mets))
    error('One or more metabolites to be added already exist in the model.');
end

% add metabolites to the model
metsToAdd = {};
metsToAdd.mets = metData{1};
metsToAdd.metNames = metData{2};
metsToAdd.compartments = metData{3};
ihuman = addMets(ihuman, metsToAdd);

% add metabolites to the metAssoc structure
numOrigMets = numel(metAssoc.mets);
numNewMets = numel(metData{1});
newMetInd = (numOrigMets+1:numOrigMets+numNewMets)';
newMetsNoComp = regexprep(metData{1},'.$','');

% copy association data from other compartment versions of metabolites
f = fieldnames(metAssoc);
for i = 1:numel(f)
    for j = 1:numNewMets
        [~,ind] = ismember(newMetsNoComp(j), metAssoc.metsNoComp);
        if ind == 0
            error('Metabolite does not exist in any model compartments!');
        end
        metAssoc.(f{i})(newMetInd(j)) = metAssoc.(f{i})(ind);
    end
end
metAssoc.mets(newMetInd) = metData{1};
metAssoc.metsNoComp(newMetInd) = newMetsNoComp;


%% Add new reactions to the model

% load new reaction information from file
fid = fopen('../../ComplementaryData/modelCuration/fullRebalance/rebalance_rxns_new.tsv');
rxnData = textscan(fid,'%s%s%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% verify that none of the reactions exist in the current model
if any(ismember(rxnData{1},ihuman.rxns))
    error('One or more reactions to be added already exist in the model.');
end

% add reactions to the model
rxnsToAdd = {};
rxnsToAdd.rxns = rxnData{1};
rxnsToAdd.equations = rxnData{2};
rxnsToAdd.rxnNames = rxnData{3};
rxnsToAdd.subSystems = cellfun(@(s) {{s}},rxnData{4});
ihuman = addRxns(ihuman, rxnsToAdd, 3);
ihuman.priorCombiningGrRules(end+1:end+numel(rxnsToAdd.rxns)) = {''};

% add reactions to the annotation structure
numOrigRxns = numel(rxnAssoc.rxns);
numNewRxns = numel(rxnData{1});
newRxnInd = (numOrigRxns+1:numOrigRxns+numNewRxns)';
f = fieldnames(rxnAssoc);
for i = 1:numel(f)
    % initialize structure with empty entries
    rxnAssoc.(f{i})(newRxnInd) = {''};
end
rxnAssoc.rxns(newRxnInd) = rxnData{1};


%% Update metabolite formulas and charges

% load metabolite information from file
fid = fopen('../../ComplementaryData/modelCuration/fullRebalance/rebalance_mets_updatedFormula.tsv');
metData = textscan(fid,'%s%s%d','Delimiter','\t','Headerlines',1);
fclose(fid);

% extract data
metNames = metData{1};
metFormulas = metData{2};
metCharges = metData{3};

% update metabolite formulas and/or charges
for i = 1:numel(metNames)
    met_ind = ismember(ihuman.metNames, metNames(i));
    if ~any(met_ind)
        error('Metabolite "%s" not found in model.',metNames{i});
    end
    ihuman.metFormulas(met_ind) = metFormulas(i);
    ihuman.metCharges(met_ind) = metCharges(i);
end


%% Update reaction equations

% load reaction information from file
fid = fopen('../../ComplementaryData/modelCuration/fullRebalance/rebalance_rxns_updatedEqn.tsv');
rxnData = textscan(fid,'%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% update reaction equations
ihuman = changeRxns(ihuman, rxnData{1}, rxnData{2}, 3);


%% Delete duplicated reactions

% load information on reactions identified as duplicates
fid = fopen('../../ComplementaryData/modelCuration/fullRebalance/rebalance_rxns_duplicated.tsv');
rxnData = textscan(fid,'%s%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% extract data
rxnIDdup = rxnData{1};  % duplicated rxns to remove
rxnIDmatch = rxnData{2};  % matching rxns that will be kept
rxnDupNotes = rxnData{3};  % notes

% add data to rxnAssoc structure
for i = 1:numel(rxnIDdup)
    rmatch = strsplit(rxnIDmatch{i},'; ');
    for j = 1:numel(rmatch)
        % find rxn in rxnAssoc structure
        [~,rmatch_ind] = ismember(rmatch(j), rxnAssoc.rxns);
        if any(rmatch_ind == 0)
            error('Reaction "%s" not found in rxnAssoc structure!',rmatch{j});
        end
        
        % determine if duplicated rxn is from HMR or Recon3D
        if startsWith(rxnIDdup(i),'HMR_')
            idtype = 'rxnHMRID';
        else
            idtype = 'rxnRecon3DID';
        end
        
        % add rxn ID to rxnAssoc structure
        if isempty(rxnAssoc.(idtype){rmatch_ind})
            rxnAssoc.(idtype)(rmatch_ind) = rxnIDdup(i);
        else
            rxnAssoc.(idtype)(rmatch_ind) = join([rxnAssoc.(idtype)(rmatch_ind), rxnIDdup(i)],'; ');
        end
    end
end

% delete reactions from model
ihuman = removeReactionsFull(ihuman,rxnIDdup);

% delete reactions from rxnAssoc
[~,remInd] = ismember(rxnIDdup,rxnAssoc.rxns);
f = fieldnames(rxnAssoc);
for i = 1:numel(f)
    rxnAssoc.(f{i})(remInd) = [];
end
changeNotes = [changeNotes; [rxnIDdup, rxnDupNotes]];


%% Inactivate 196 invalid reactions
% 196 reactions that are imbalanced and/or unsupported by any literature or
% databases will be inactivated (lb = ub = 0). These reactions will
% eventually be fully deleted from the model if they cannot be properly
% revised or repaired.

% load information on reactions to be inactivated
fid = fopen('../../ComplementaryData/modelCuration/fullRebalance/rebalance_rxns_inactivate.tsv');
rxnData = textscan(fid,'%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% inactivate reactions
ihuman = setParam(ihuman,'eq',rxnData{1},0);
changeNotes = [changeNotes; [rxnData{1}, rxnData{2}]];


%% Reactivate 10 previously inactivated reactions
% Ten reactions were repaired or are no longer in violation of mass
% balances or other problems that led to their prior inactivation. These
% reactions will be reactivated here (UB set to 1000; it was confirmed that
% none of these reactions were previously reversible, so the LB will remain
% at zero).

% load information on reactions to be reactivated
fid = fopen('../../ComplementaryData/modelCuration/fullRebalance/rebalance_rxns_reactivate.tsv');
rxnData = textscan(fid,'%s','Headerlines',1);
fclose(fid);

% reactivate reactions
ihuman = setParam(ihuman,'ub',rxnData{1},1000);
changeNotes = [changeNotes; [rxnData{1},...
repmat({'reaction is no longer invalid/inconsistent and was reactivated'},numel(rxnData{1}),1)]];


%% Merge and delete duplicated metabolites

% load information on duplicated metabolites
fid = fopen('../../ComplementaryData/modelCuration/fullRebalance/rebalance_mets_duplicated.tsv');
metData = textscan(fid,'%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% extract data
dupMetName = metData{1};
keepMetName = metData{2};

% merge duplicated metabolites with their matches in stoich matrix
for i = 1:numel(dupMetName)
    metComps = unique(ihuman.metComps(ismember(ihuman.metNames, dupMetName(i))));
    for j = 1:numel(metComps)
        dupMetInd = find(ismember(ihuman.metNames,dupMetName(i)) & ihuman.metComps == metComps(j));
        keepMetInd = find(ismember(ihuman.metNames,keepMetName(i)) & ihuman.metComps == metComps(j));
        ihuman.S(keepMetInd,:) = ihuman.S(keepMetInd,:) + ihuman.S(dupMetInd,:);
    end
end

% remove duplicated metabolites from the model
remMetInd = find(ismember(ihuman.metNames, dupMetName)); % needed for editing metAssoc below
remMet = ihuman.mets(remMetInd); % needed for editing metAssoc below
ihuman = removeMets(ihuman,dupMetName,true);

% remove duplicated metabolites from the metAssoc structure
f = fieldnames(metAssoc);
for i = 1:numel(f)
    metAssoc.(f{i})(remMetInd) = [];
end
changeNotes = [changeNotes; [remMet,...
repmat({'metabolite is a duplicate and was therefore removed'},numel(remMet),1)]];


%% Update annotation information for 27 unique metabolites

% load updated metabolite annotation information
fid = fopen('../../ComplementaryData/modelCuration/fullRebalance/rebalance_mets_updatedAnnotation.tsv');
metData = textscan(fid,'%s%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% extract data
metName = metData{1};
annotField = metData{2};
annotValue = metData{3};

% update annotation information
for i = 1:numel(metName)
    metInd = find(ismember(ihuman.metNames,metName(i)));
    metAssoc.(annotField{i})(metInd) = annotValue(i);
end


%% Update MetaNetX, BiGG, and KEGG ids for 15 reactions

% load updated reaction annotation information
fid = fopen('../../ComplementaryData/modelCuration/fullRebalance/rebalance_rxns_updatedAnnotation.tsv');
rxnData = textscan(fid,'%s%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% extract data
rxn = rxnData{1};
annotField = rxnData{2};
annotValue = rxnData{3};

% update annotation information
for i = 1:numel(rxn)
    [~,rxnInd] = ismember(rxn(i),ihuman.rxns);
    rxnAssoc.(annotField{i})(rxnInd) = annotValue(i);
end


%% Update 38 reaction grRules
% These new gene associations were found through some of the reaction
% duplication cases

ihuman = updateGrRules('fullRebalance/rebalance_rxns_updatedGrRules.tsv',1,2,false,ihuman);


%% Finalize model changes

% remove unused metabolites and/or genes from the model
metsOrig = ihuman.mets;
ihuman = removeReactions(ihuman,[],true,true);
metsRemoved = setdiff(metsOrig,ihuman.mets);

% if any metabolites were removed, also remove them from metAssoc
if ~isempty(metsRemoved)
    [~,remInd] = ismember(metsRemoved,metAssoc.mets);
    f = fieldnames(metAssoc);
    for i = 1:numel(f)
        metAssoc.(f{i})(remInd) = [];
    end
    changeNotes = [changeNotes; [metsRemoved,...
    repmat({'metabolite no longer used after removing reactions'},numel(metsRemoved),1)]];
end

% update bounds and reversibility field to be consistent with each other
ihuman.rev(ihuman.lb == ihuman.ub) = 0;  % inactivated rxns are not reversible
ihuman.lb(ihuman.rev == 1) = -1000;
ihuman.lb(ihuman.rev == 0) = 0;
ihuman.rev = double(ihuman.lb < 0);

% update unconstrained field
ihuman.unconstrained = double(ihuman.metComps == 9);


%% Document changes and export files

% verify that HumanGEM and annotation structures are aligned
if ~isequal(ihuman.mets, metAssoc.mets) || ~isequal(ihuman.rxns, rxnAssoc.rxns)
    error('HumanGEM is not synced with the metAssoc and/or rxnAssoc structure!');
end

% export reaction and metabolite annotation structures to JSON
jsonStr = jsonencode(rxnAssoc);
fid = fopen('../../ComplementaryData/annotation/humanGEMRxnAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);
jsonStr = jsonencode(metAssoc);
fid = fopen('../../ComplementaryData/annotation/humanGEMMetAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);

% determine and document model changes
modelChanges = docModelChanges(ihuman_orig,ihuman,changeNotes);
writeModelChanges(modelChanges,'../../ComplementaryData/modelCuration/fullRebalance/rebalance_modelChanges.tsv');

% export HumanGEM
exportHumanGEM(ihuman,'HumanGEM','../../',{'mat','yml'},false,false);

% clear unneeded variables
clearvars -except ihuman_orig ihuman modelChanges


%% Test and report model balance stats

% remove inactivated reactions before performing tests
model_orig = simplifyModel(ihuman_orig,false,false,true);
model = simplifyModel(ihuman,false,false,true);
num_rxns_orig = numel(model_orig.rxns);
num_rxns = numel(model.rxns);

% check reaction mass balances before and after changes
bal_orig = getElementalBalance(model_orig);
num_unbal_orig = sum(bal_orig.balanceStatus ~= 1);
bal = getElementalBalance(model);
num_unbal = sum(bal.balanceStatus ~= 1);

fprintf('Number of mass-imbalanced reactions before changes: %u (%.2f%%)\n', num_unbal_orig, num_unbal_orig/num_rxns_orig*100);
fprintf('Number of mass-imbalanced reactions after changes: %u (%.2f%%)\n\n', num_unbal, num_unbal/num_rxns*100);

% check reaction charge balances
bal_orig = model_orig.S' * double(model_orig.metCharges);
num_unbal_orig = sum(bal_orig ~= 0);
bal = model.S' * double(model.metCharges);
num_unbal = sum(bal ~= 0);

fprintf('Number of charge-imbalanced reactions before changes: %u (%.2f%%)\n', num_unbal_orig, num_unbal_orig/num_rxns_orig*100);
fprintf('Number of charge-imbalanced reactions after changes: %u (%.2f%%)\n\n', num_unbal, num_unbal/num_rxns*100);

% check stoichiometric consistency (requires Cobra)
[~,m_orig] = checkStoichiometricConsistency(model_orig);
[~,m] = checkStoichiometricConsistency(model);

fprintf('Model consistency (conserved metabolites) before changes: %.2f%%\n', sum(m_orig ~= 0)/numel(model_orig.mets)*100);
fprintf('Model consistency (conserved metabolites) after changes: %.2f%%\n', sum(m ~= 0)/numel(model.mets)*100);




%% Old code for intermediate analyses during model rebalancing
% 
% 
% %% load MNX data
% 
% % load MetaNetX metabolite annotations (loads as "MNXMets" structure)
% load('../../ComplementaryData/MetaNetX/MNXMets.mat');
% 
% % remove unneccessary MNXMets entries to reduce structure size
% allModelMNXids = cellfun(@(x) strsplit(x,'; '), metAssoc.metMNXID, 'UniformOutput', false);
% allModelMNXids = unique([allModelMNXids{:}]');
% allModelMNXids(cellfun(@isempty, allModelMNXids)) = [];  % remove empty ID
% rem_ind = ~ismember(MNXMets.mets, allModelMNXids);
% f = fieldnames(MNXMets);
% for i = 1:numel(f)
%     MNXMets.(f{i})(rem_ind) = [];
% end
% 
% 
% %% Retrieve metabolite formulas and charges from corresponding MNX IDs
% 
% % specify metabolite indices for which information will be retreived
% % met_ind = find(contains(lower(ihuman.metNames),{''}));
% 
% % ignore mets repeated in different compartments
% [~,uniq_ind] = unique(ihuman.metNames);
% met_ind = intersect(met_ind, uniq_ind);
% 
% % for each metabolite, check if there is another formula available from MNX
% mnxIDs = repmat({''},size(met_ind));
% mnxFormulas = repmat({''},size(met_ind));
% mnxCharges = NaN(size(met_ind));
% for i = 1:numel(met_ind)
%     
%     % get MNX ID(s) for the metabolite
%     if isempty(metAssoc.metMNXID{met_ind(i)})
%         continue
%     else
%         mnxIDs{i} = metAssoc.metMNXID{met_ind(i)};
%         mnx_id = strsplit(mnxIDs{i}, '; ');
%     end
%     
%     % retrieve formula and charge from MNXMets structure
%     [~,mnxid_ind] = ismember(mnx_id, MNXMets.mets);
%     mnx_formula = unique(MNXMets.metFormulas(mnxid_ind));
%     mnx_formula(cellfun(@isempty, mnx_formula)) = [];
%     if ~isempty(mnx_formula)
%         mnxFormulas(i) = join(mnx_formula, '; ');
%     end
%     
%     mnx_charge = unique(MNXMets.metCharges(mnxid_ind));
%     mnx_charge(arrayfun(@isnan, mnx_charge)) = [];
%     if ~isempty(mnx_charge)
%         mnxCharges(i) = mnx_charge;
%     end
%     
% end
% 
% % organize information
% x = [{'mets','metNames','metFormulas','metCharges','MNXIDs','MNXformulas','MNXcharges'};
%     [ihuman.mets(met_ind), ihuman.metNames(met_ind), ihuman.metFormulas(met_ind), ...
%     num2cell(ihuman.metCharges(met_ind)), mnxIDs, mnxFormulas, num2cell(mnxCharges)]];
% 
% 
% 
% %% Analyze effect of model changes on reaction balance status
% 
% % initialize variables
% metNames = {};
% metFormulas = {};
% metCharges = int64([]);
% rxns = {};
% rxnEqns = {};
% 
% % change metabolite formulas and/or charges
% for i = 1:numel(metNames)
%     ind = ismember(ihuman.metNames, metNames(i));
%     if ~any(ind)
%         error('Metabolite "%s" not found in model.',metNames{i});
%     end
%     ihuman.metFormulas(ind) = metFormulas(i);
%     ihuman.metCharges(ind) = metCharges(i);
% end
% 
% % change reaction equations
% ihuman = changeRxns(ihuman, rxns, rxnEqns, 3);
% 
% % perform mass and charge balance analysis
% [massImbal1, imBalMass1, imBalCharge1, imBalRxnBool1, Elements1, missFormulaeBool1, balMetBool1] = ...
%     checkMassChargeBalance(ihuman_orig);
% [massImbal2, imBalMass2, imBalCharge2, imBalRxnBool2, Elements2, missFormulaeBool2, balMetBool2] = ...
%     checkMassChargeBalance(ihuman);
% 
% % check only charge balances
% imBalCharge = ihuman.S' * double(ihuman.metCharges);
% 
% % determine reactions that became balanced
% newBalMass = find(imBalRxnBool1 & ~imBalRxnBool2);
% newBalCharge = find(imBalCharge1 ~= 0 & imBalCharge2 == 0);
% 
% % determine reactions that became unbalanced
% unBalMass = find(~imBalRxnBool1 & imBalRxnBool2);
% unBalCharge = find(imBalCharge1 == 0 & imBalCharge2 ~= 0);
% 
% 
% % check balance status of a selected set of reactions
% bal = getElementalBalance(ihuman, r_ind);
% elDiff = elementalMatrixToFormulae(bal.rightComp - bal.leftComp, bal.elements.abbrevs);
% 


