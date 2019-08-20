%
% FILE NAME:    curateMetFormulas2.m
% 
% DATE CREATED: 2019-07-22
%     MODIFIED: 2019-07-22
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: This script curates the formulas of several metabolites to
%          improve model accuracy and reaction mass balances.
%

%% Load Human-GEM and metAssoc structure

% load Human-GEM model
load('humanGEM.mat');
ihuman_orig = ihuman;  % keep copy of original version

% load metabolite annotation data
metAssoc = jsondecode(fileread('../../ComplementaryData/annotation/humanGEMMetAssoc.JSON'));

% verify that HumanGEM and metAssoc structure are aligned
if ~isequal(ihuman.mets, metAssoc.mets)
    error('HumanGEM metabolites are not synced with the metAssoc structure!');
end

%% Load MetaNetX metabolite annotation data

% load MetaNetX metabolite annotations (loads as "MNXMets" structure)
load('../../ComplementaryData/MetaNetX/MNXMets.mat');

% remove unneccessary MNXMets entries to reduce structure size
allModelMNXids = cellfun(@(x) strsplit(x,'; '), metAssoc.metMNXID, 'UniformOutput', false);
allModelMNXids = unique([allModelMNXids{:}]');
allModelMNXids(cellfun(@isempty, allModelMNXids)) = [];  % remove empty ID
rem_ind = ~ismember(MNXMets.mets, allModelMNXids);
f = fieldnames(MNXMets);
for i = 1:numel(f)
    MNXMets.(f{i})(rem_ind) = [];
end


%% Retrieve metabolite formulas and charges from corresponding MNX IDs

% specify metabolite indices for which information will be retreived
% met_ind = find(contains(lower(ihuman.metNames),{'dolichol';'dolichyl'}));

% ignore mets repeated in different compartments
[~,uniq_ind] = unique(ihuman.metNames);
met_ind = intersect(met_ind, uniq_ind);

% for each metabolite, check if there is another formula available from MNX
mnxIDs = repmat({''},size(met_ind));
mnxFormulas = repmat({''},size(met_ind));
mnxCharges = NaN(size(met_ind));
for i = 1:numel(met_ind)
    
    % get MNX ID(s) for the metabolite
    if isempty(metAssoc.metMNXID{met_ind(i)})
        continue
    else
        mnxIDs{i} = metAssoc.metMNXID{met_ind(i)};
        mnx_id = strsplit(mnxIDs{i}, '; ');
    end
    
    % retrieve formula and charge from MNXMets structure
    [~,mnxid_ind] = ismember(mnx_id, MNXMets.mets);
    mnx_formula = unique(MNXMets.metFormulas(mnxid_ind));
    mnx_formula(cellfun(@isempty, mnx_formula)) = [];
    if ~isempty(mnx_formula)
        mnxFormulas(i) = join(mnx_formula, '; ');
    end
    
    mnx_charge = unique(MNXMets.metCharges(mnxid_ind));
    mnx_charge(arrayfun(@isnan, mnx_charge)) = [];
    if ~isempty(mnx_charge)
        mnxCharges(i) = mnx_charge;
    end
    
end

% % remove empty entries, or those with no change in ID
% sameID = strcmp(standardizeMetFormulas(ihuman.metFormulas(X_met_ind)), ...
%                 standardizeMetFormulas(altFormulas));
% emptyID = cellfun(@isempty, altFormulas);
% X_met_ind(sameID | emptyID) = [];
% altFormulas(sameID | emptyID) = [];

% organize information
x = [{'mets','metNames','metFormulas','metCharges','MNXIDs','MNXformulas','MNXcharges'};
    [ihuman.mets(met_ind), ihuman.metNames(met_ind), ihuman.metFormulas(met_ind), ...
    num2cell(ihuman.metCharges(met_ind)), mnxIDs, mnxFormulas, num2cell(mnxCharges)]];


%% Analyze effect of model changes on reaction balance status

% initialize variables
metNames = {};
metFormulas = {};
metCharges = int64([]);
rxns = {};
rxnEqns = {};

% change metabolite formulas and/or charges
for i = 1:numel(metNames)
    ind = ismember(ihuman.metNames, metNames(i));
    ihuman.metFormulas(ind) = metFormulas(i);
    ihuman.metCharges(ind) = metCharges(i);
end

% change reaction equations
ihuman = changeRxns(ihuman, rxns, rxnEqns, 3);

% perform mass and charge balance analysis
[massImbal1, imBalMass1, imBalCharge1, imBalRxnBool1, Elements1, missFormulaeBool1, balMetBool1] = ...
    checkMassChargeBalance(ihuman_orig);
[massImbal2, imBalMass2, imBalCharge2, imBalRxnBool2, Elements2, missFormulaeBool2, balMetBool2] = ...
    checkMassChargeBalance(ihuman);

% determine reactions that became balanced
newBalMass = find(imBalRxnBool1 & ~imBalRxnBool2);
newBalCharge = find(imBalCharge1 ~= 0 & imBalCharge2 == 0);

% determine reactions that became unbalanced
unBalMass = find(~imBalRxnBool1 & imBalRxnBool2);
unBalCharge = find(imBalCharge1 == 0 & imBalCharge2 ~= 0);


% check balance status of a selected set of reactions
bal = getElementalBalance(ihuman, r_ind);
elDiff = elementalMatrixToFormulae(bal.rightComp - bal.leftComp, bal.elements.abbrevs);



