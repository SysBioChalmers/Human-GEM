%
% FILE NAME:    curateATPmetabolism.m
% 
% PURPOSE: This script makes a number of improvements to the modeling of 
%          ATP metabolism (electron transport chain) in HumanGEM.
%
%          In addition, this script removes the now-outdated "rxnRecon3DID"
%          model field. Such annotations are now stored in the separate
%          humanGEMRxnAssoc.JSON file.
%
%          Finally, the "version" model field is cleared, as it is now
%          assigned/documented in an alternative manner.


%% Load model and initialize variables

% load HumanGEM
load('humanGEM.mat');
ihuman_orig = ihuman;  % to track changes

% initialize reaction change notes
rxnNotes = {};


%% Remove rxnRecon3DID model field and clear version field
% this information is identical to the association information stored in
% the humanGEMRxnAssoc.JSON file, and is therefore unnecessary
ihuman = rmfield(ihuman,'rxnRecon3DID');
ihuman.version = '';


%% Remove duplicate reaction
% A proton transporter is currently duplicated in HumanGEM v1.0.3:
%
%   HMR_7638:  H+[i] => H+[m]
%       Htmi:  H+[i] => H+[m]
%
% The original form of HMR_7638 (from HMR2) involved transport of protons
% from the cytoplasm to the mitochondria: H+[c] => H+[m], but was modified
% to the current form in HumanGEM v0.5.0. Since this reaction is necessary
% to restore protons that are lost through the adenine nucleotide
% transporter (ATP [m] => ATP[c]), the HMR_7638 reaction equation will be
% reverted to its original form (H+[c] => H+[m]).

% get relevant reaction and metabolite
rxn_ind = getIndexes(ihuman,{'HMR_7638';'Htmi'},'rxns');
Hi_ind = getIndexes(ihuman,'H+[i]','metscomps');
Hc_ind = getIndexes(ihuman,'H+[c]','metscomps');

% update stoichiometry of HMR_7638
ihuman.S(Hi_ind,rxn_ind(1)) = 0;
ihuman.S(Hc_ind,rxn_ind(1)) = -1;

% need to update some associations related to these reactions
rxnAssoc = jsondecode(fileread('../../ComplementaryData/annotation/humanGEMRxnAssoc.JSON'));
if ~isequal(ihuman.rxns,rxnAssoc.rxns)
    % model reactions must be consistent with the annotation file
    error('Model reactions are inconsistent with humanGEMRxnAssoc.JSON file!');
end
rxnAssoc.rxnRecon3DID(rxn_ind(1)) = {'Htm'};  % update the Recon3D ID association for HMR_7638
rxnAssoc.rxnMNXID(rxn_ind) = {'MNXR100765'};  % MNXIDs for HMR_7638 and Htmi are the same, because compartments are ignored

rxnNotes = [rxnNotes; [{'HMR_7638'} {'Reaction was identical to Htmi, and proton transport is needed from cytosol to mitochondria, so the reaction was reverted to its original form.'}]];


%% Update the stoichiometry of ATP synthase
% [Issue #98]
% The stoichiometry for protons in ATP synthase (HMR_6916) should be 3
% instead of 4. Also, it needs to be balanced by removing one additional
% proton from the products.
%   Current: ADP[m] + Pi[m] + 4 H+[i] => ATP[m] + 4 H+[m] + H2O[m]
%   Revised: ADP[m] + Pi[m] + 3 H+[i] => ATP[m] + 2 H+[m] + H2O[m]
rxn_ind = getIndexes(ihuman,'HMR_6916','rxns');
Hi = getIndexes(ihuman,'H+[i]','metscomps');
Hm = getIndexes(ihuman,'H+[m]','metscomps');
ihuman.S([Hi;Hm],rxn_ind) = [-3;2];

rxnNotes = [rxnNotes; [{'HMR_6916'} {'Balanced reaction mass and charge, and changed to 3 protons pumped per ATP produced (PMID: 15620362)'}]];


%% Prevent free transport of Pi from [c] to [m]
% [Issue #99]
% The following reactions should only move Pi from [m] to [c]
%   HMR_3971    fumarate[m] + Pi[c] <=> fumarate[c] + Pi[m]
%   HMR_4862    Pi[c] + succinate[m] <=> Pi[m] + succinate[c]
%   HMR_4865    malate[m] + Pi[c] <=> malate[c] + Pi[m]
%   HMR_4870    malonate[m] + Pi[c] <=> malonate[c] + Pi[m]
%   HMR_4940    GSH[c] + Pi[m] <=> GSH[m] + Pi[c]
%   HMR_6330    AKG[c] + Pi[m] <=> AKG[m] + Pi[c]
%   HMR_6331    oxalate[m] + Pi[c] <=> oxalate[c] + Pi[m]
rxns = {'HMR_3971';'HMR_4862';'HMR_4865';'HMR_4870';'HMR_4940';'HMR_6330';'HMR_6331'};
rxn_ind = getIndexes(ihuman,rxns,'rxns');
Pi_m = getIndexes(ihuman,'Pi[m]','metscomps');

% reactions written with Pi[c] -> Pi[m] need to be turned around
flip_ind = full(ihuman.S(Pi_m,rxn_ind)') > 0;
ihuman.S(:,rxn_ind(flip_ind)) = -ihuman.S(:,rxn_ind(flip_ind));

% now make all reactions irreversible
ihuman.lb(rxn_ind) = 0;

rxnNotes = [rxnNotes; [rxns, repmat({'Reaction reversibility updated to prevent free transport of Pi from cytosol to mitochondria, which should require 1 proton co-transported (PMID:15620362;22733773;21706682)'},numel(rxns),1)]];


%% ATP synthesis from propanoate fermentation
% [Issue #100]
% In the presence of glucose and absence of oxygen the model synthesizes
% propanoate instead of lactate with an ATP yield of ~5 ATP/glucose instead
% of 2. This is caused by the reversibility of the following reactions:
%
%   HMR_0153    AMP[c] + PPi[c] + propanoyl-CoA[c] <=> ATP[c] + CoA[c] + propanoate[c]
%   HMR_4459    ATP[c] + H+[c] + propanoate[c] <=> PPi[c] + propinol adenylate[c]
%
% Both reactions should be irreversible, where HMR_0153 should only proceed
% in the backward direction (consuming propanoate).
rxns = {'HMR_0153';'HMR_4459'};
rxn_ind = getIndexes(ihuman,rxns,'rxns');
ihuman.S(:,rxn_ind(1)) = -ihuman.S(:,rxn_ind(1));

% make reactions irreversible
ihuman.lb(rxn_ind) = 0;

rxnNotes = [rxnNotes; [rxns, repmat({'Reaction reversibility updated to prevent synthesis of propanoate instead of lactate in the absence of oxygen, which leads to an ATP yield of ~5 ATP/glucose instead of 2'},numel(rxns),1)]];


%% Flux through complex I under anaerobic conditions
% [Issue #101]
% In the absence of oxygen and presence of glucose the model has flux
% through complex I and synthesizes 3-Methyl-Glutaconate and
% 2-Methyl-3-Hydroxy-Valerate. This results in ATP yield > 4 per glucose 
% instead of 2. This occurs because the model is able to use these as
% electron acceptors, via ubiquinol and FAD. This can be resolved by making
% RE1519X and HMR_3212 irreversible.
%   RE1519X     4-cis-decenoyl-CoA[p] + FAD[p] <=> 2-trans-4-cis-decadienoyl-CoA[p] + FADH2[p]
%   HMR_3212    FAD[m] + propanoyl-CoA[m] <=> acrylyl-CoA[m] + FADH2[m]
rxns = {'RE1519X';'HMR_3212'};
rxn_ind = getIndexes(ihuman,rxns,'rxns');
ihuman.lb(rxn_ind) = 0;

rxnNotes = [rxnNotes; [rxns, repmat({'Reaction reversibility updated to prevent use of 3-Methyl-Glutaconate and 2-Methyl-3-Hydroxy-Valerate as electron acceptors.'},numel(rxns),1)]];


%% ATP and carbon from Pi and O2
%[Issue #102]
% The model is able to produce infinite ATP and CO2 from Pi and O2. This
% can be prevented by blocking the flux through the pool reaction HMR_0686:
% HMR_0686  fatty acid-LD-TG2 pool[c] => 0.0001 (10Z)-heptadecenoic acid[c] + 0.0001 (11Z,14Z)-eicosadienoic acid[c] + 0.0001 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0001 (13Z)-eicosenoic acid[c] + 0.0001 (13Z)-octadecenoic acid[c] + 0.0001 (13Z,16Z)-docosadienoic acid[c] + 0.0016 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0001 (6Z,9Z)-octadecadienoic acid[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0001 (7Z)-octadecenoic acid[c] + 0.0001 (7Z)-tetradecenoic acid[c] + 0.0001 (9E)-tetradecenoic acid[c] + 0.0001 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0001 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0001 10,13,16,19-docosatetraenoic acid[c] + 0.0001 10,13,16-docosatriynoic acid[c] + 0.0001 12,15,18,21-tetracosatetraenoic acid[c] + 0.0001 13,16,19-docosatrienoic acid[c] + 0.0001 7-palmitoleic acid[c] + 0.0001 8,11-eicosadienoic acid[c] + 0.0001 9-eicosenoic acid[c] + 0.0001 9-heptadecylenic acid[c] + 0.0032 adrenic acid[c] + 0.0521 arachidonate[c] + 0.0001 behenic acid[c] + 0.0001 cerotic acid[c] + 0.0001 cis-cetoleic acid[c] + 0.0001 cis-erucic acid[c] + 0.0001 cis-gondoic acid[c] + 0.0217 cis-vaccenic acid[c] + 0.0081 DHA[c] + 0.0023 dihomo-gamma-linolenate[c] + 0.0018 DPA[c] + 0.0001 eicosanoate[c] + 0.0001 elaidate[c] + 0.0003 EPA[c] + 0.003 gamma-linolenate[c] + 0.0001 henicosanoic acid[c] + 0.0001 lauric acid[c] + 0.0001 lignocerate[c] + 0.281 linoleate[c] + 0.0119 linolenate[c] + 0.0001 margaric acid[c] + 0.0001 mead acid[c] + 0.0024 myristic acid[c] + 0.0001 nervonic acid[c] + 0.0001 nonadecylic acid[c] + 0.4226 oleate[c] + 0.0001 omega-3-arachidonic acid[c] + 0.12 palmitate[c] + 0.0229 palmitolate[c] + 0.0001 pentadecylic acid[c] + 0.0001 physeteric acid[c] + 0.0384 stearate[c] + 0.0025 stearidonic acid[c] + 0.0001 tricosanoic acid[c] + 0.0001 tridecylic acid[c] + 0.0001 ximenic acid[c]
ihuman = setParam(ihuman,'eq','HMR_0686',0);
rxnNotes = [rxnNotes; [{'HMR_0686'} {'This pool reaction enables infinite ATP and CO production from Pi and O2, and therefore should be inactivated until it can be properly reformulated or removed entirely.'}]];


%% Verify model changes (OPTIONAL)

% load metabolic tasks that assess issues addressed in this and previous
% model curation processes
taskStruct = parseTaskList('../../ComplementaryData/metabolicTasks/metabolicTasks_VerifyModel.xls');

% check original HumanGEM task performance
checkTasks(ihuman_orig,[],true,false,false,taskStruct);

% check curated HumanGEM task performance
checkTasks(ihuman,[],true,false,false,taskStruct);


%% Write updated reaction annotation file

jsonStr = jsonencode(rxnAssoc);
fid = fopen('../../ComplementaryData/annotation/humanGEMRxnAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);


%% Document reaction changes

rxnChanges = docRxnChanges(ihuman_orig,ihuman,rxnNotes);
writeRxnChanges(rxnChanges,'../../ComplementaryData/modelCuration/curateATPmetabolism_rxnChanges.tsv');


%% Export updated HumanGEM

exportHumanGEM(ihuman,'humanGEM','../../',{'mat','yml'},false,false);




