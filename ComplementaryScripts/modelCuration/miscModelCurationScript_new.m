%
% FILE NAME:    miscModelCurationScript_new.m
% 
% DATE CREATED: 2018-09-21
%     MODIFIED: 2018-09-21
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script for performing a number of different curations, updates,
%          and corrections to HumanGEM.
%


%% Update of GPRs

% The following reaction:
%   HMR_2116: H+[c] + vitamin D3[c] <=> H+[m] + vitamin D3[m]
% Is associated with (only) the following gene:
%   ENSG00000175592 (FOSL1)
% However, this FOSL1 protein is described as a transcription factor, is
% not associated with any catalytic activity, and UniProt shows it as
% localized to the nucleus. It therefore has no apparent relationship with
% the above reaction, and should be removed from the associated grRule.
[~,rxn_ind] = ismember('HMR_2116', ihuman.rxns);
if strcmp(ihuman.grRules(rxn_ind),'ENSG00000175592')
    ihuman.grRules(rxn_ind) = {''};
    ihuman.rxnGeneMat(rxn_ind,:) = 0;
end


%% Removal of invalid/unsupported reactions

% This Recon3D reaction does not exist in the literature or any databases:
%     r1453: proline[m] + ubiquinol[m] <=> 1-pyrroline-5-carboxylate[m] + 5 H+[m] + ubiquinone[m]
% It is also not present on the BiGG Database. It should thus be removed.
del_rxns = [del_rxns; {'r1453'}];


% These reactions are strange in that they involve an additional
% metabolite "Adenosine-5'-Triphosphate-Energy". This metabolite was not
% found anywhere else in the model, and therefore these reactions are
% dead-end. They should be removed from the model:
%   r1319: ATP[c] + H2O[c] => ADP[c] + Pi[c] + Adenosine-5'-Triphosphate-Energy[c]
%   r1320: ATP[m] + H2O[m] => ADP[m] + Pi[m] + Adenosine-5'-Triphosphate-Energy[m]
del_rxns = [del_rxns; {'r1319'; 'r1320'}];



%% Add a reaction allowing transport of protons from I --> C
% This reaction does not generate energy, and is primarily to allow mass
% balancing of protons in the model.
rxnsToAdd = {};
rxnsToAdd.rxns = {'Htransport_I_C'};  % this needs to be replaced with a systematic name
rxnsToAdd.equations = {'H+[i] => H+[c]'};
rxnsToAdd.rxnNames = rxnsToAdd.rxns;
rxnsToAdd.lb = 0;
rxnsToAdd.ub = 1000;

% this is an ugly way to get around the automatic "standardization" of
% grRules which is enforced by the "addRxns" function, which throws a
% warning which is not relevant for us at the moment.
grRules_orig = ihuman.grRules;
ihuman.grRules(:) = {''};

% add the new reactions
ihuman = addRxns(ihuman,rxnsToAdd,3,[],false);
ihuman.grRules(1:length(grRules_orig)) = grRules_orig;  % restore original rules

% update other non-standard fields
ind = length(ihuman.rxns);
ihuman.rxnKEGGID(ind) = {''};
ihuman.rxnEHMNID(ind) = {''};
ihuman.rxnBiGGID(ind) = {''};
ihuman.rxnHepatoNET1ID(ind) = {''};
ihuman.rxnREACTOMEID(ind) = {''};
ihuman.rxnRecon3DID(ind) = {''};
ihuman.prRules(ind) = {''};
ihuman.rxnProtMat(ind,:) = 0;
ihuman.priorCombiningGrRules(ind) = {''};




%% Treatment of ubiquinone and FAD+
% The model essentially treats FAD/FADH2 the same as ubiquinone/ubiquinol,
% as exhibited by the presence of the following reaction:
%  HMR_6911: FADH2[m] + ubiquinone[m] <=> FAD[m] + ubiquinol[m]
%
% Some reactions in the model are duplicated such that one version of the
% reaction uses FAD, whereas the other uses ubiquinone. These are
% effectively identical reactions, and one should be removed (unless
% evidence suggests otherwise). 

% (Not yet implemented)

