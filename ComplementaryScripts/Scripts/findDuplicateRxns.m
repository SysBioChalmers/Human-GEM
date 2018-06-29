% Script for finding duplicate reactions in Recon4.
% Uses a modified form of the Cobra "checkDuplicateRxn" function.


%% load latest Recon4 model
x = load('ComplementaryScripts/ModelIntegration/HMR3_02.mat');
f = fields(x);
model = x.(f{1});  % ugly way to rename model without knowing variable name in advance



%% Add new field to model: rxnRecon3DID
% This is to keep track of which reactions from Recon3D were "merged" with
% each of the HMR reactions when generating the combined model (Recon4).

% load Recon3D-HMR rxn mapping information (loads "Recon3D" variable)
load('ComplementaryScripts/RxnAssociation/Recon3Rxns2HMR.mat');
model.rxnRecon3DID = repmat({''},size(model.rxns));

% get unique list of all HMR rxn IDs mapped to the Recon3D rxns
HMRrxns = unique(flattenCell(cellfun(@(r) strsplit(r,';'),Recon3D.rxnHMRID,'UniformOutput',false),true));
HMRrxns(ismember(HMRrxns,{''})) = [];  % remove empty cell if present

% associate Recon3D rxnIDs to each of the reactions in Recon4
for i = 1:length(HMRrxns)
    r3ind = contains(Recon3D.rxnHMRID,HMRrxns{i});
    r4ind = ismember(model.rxns,HMRrxns{i});
    model.rxnRecon3DID{r4ind} = strjoin(Recon3D.rxns(r3ind),'; ');
end


%% Identify duplicated reactions in Recon4

% specify names of metabolites to ignore (e.g., protons)
ignoreMets = {'H+'};
ignoreMetInds = ismember(model.metNames,ignoreMets);
if any(~ismember(ignoreMets,model.metNames))
    error('At least one of the specified ignoreMets was not found in the model.');
end

% Set the stoich matrix entries for ignored mets to zero, except in
% reactions where the ignored met appears in multiple compartments. This is
% to avoid flagging reactions involving, for example, proton transport, as
% the same as reactions without proton transport.
model_temp = model;
ignoreRxnInds = true(size(model.rxns));
for i = 1:length(ignoreMets)
    ind = ismember(model.metNames,ignoreMets(i));
    ignoreRxnInds(sum(model.S(ind,:) ~= 0,1) > 1) = false;
end
model_temp.S(ignoreMetInds,ignoreRxnInds) = 0;

% Check for duplicated reactions based on S-matrix.
% Note that reactions with stoichiometry that differs only by scalar
% multiplication will be treated as identical, and reactions that are
% identical but in reverse direction from one another will NOT be treated
% as identical.
[~,~,~,duplicateRxnSets] = checkDuplicateRxn_mod(model_temp,'S',false);

% generate reaction equations for manual inspection of reactions
rxnEqns = constructEquations(model);

% organize duplicateRxnSets into easily comparable rxn eqn format to
% facilitate manual curation
duplicateRxnInfo = {};
for i = 1:length(duplicateRxnSets)
    duplicateRxnInfo = [duplicateRxnInfo;[{''},{''}];[model.rxns(duplicateRxnSets{i}),rxnEqns(duplicateRxnSets{i})]];
end








