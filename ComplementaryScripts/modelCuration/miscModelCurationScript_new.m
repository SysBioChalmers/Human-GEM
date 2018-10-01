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

