%
% FILE NAME:    curateMitochMembraneComp.m
% 
% DATE CREATED: 2018-09-17
%     MODIFIED: 2018-09-17
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script for analyzing/curating the inner mitochondrial matrix
%          compartment.
%


%% Analysis of reactions transporting protons into/out of the mitochondria

% load HumanGEM model
load('humanGEM.mat');

% construct rxn equations
eqns = constructEquations(ihuman);

% get indices of protons in relevant compartments
[~,comp_inds] = ismember({'cytosol','mitochondria','inner mitochondria'},lower(ihuman.compNames));
if any(comp_inds == 0)
    error('one or more compartment names could not be found.');
end
Hc = find( ismember(ihuman.metNames,'H+') & (ihuman.metComps == comp_inds(1)) );
Hm = find( ismember(ihuman.metNames,'H+') & (ihuman.metComps == comp_inds(2)) );
Hi = find( ismember(ihuman.metNames,'H+') & (ihuman.metComps == comp_inds(3)) );


% find all rxns that can transport protons from [m] to [c]
fwd_rxn_inds = find((ihuman.S(Hm,:) < 0) & (ihuman.S(Hc,:) > 0) & (ihuman.ub > 0)')';  % written in forward direction
rev_rxn_inds = find((ihuman.S(Hm,:) > 0) & (ihuman.S(Hc,:) < 0) & (ihuman.lb < 0)')';  % written in reverse direction
rxn_inds = [fwd_rxn_inds; rev_rxn_inds];

% This resulted in 8 fwd rxns, and 21 rev rxns, for a total of 29 rxns


% determine which of the reactions can carry flux in the M --> C direction
% model = simplifyModel(ihuman);  % make model simulation-ready by deleting boundary mets
% model.lb(startsWith(model.rxns,'sink_')) = 0;  % don't allow sink reactions to operate in reverse
% max_flux = [];
% for i = 1:length(rxn_inds)
%     fprintf('%u\n',i);
%     m = model;
%     m.c(:) = 0;
%     if ismember(rxn_inds(i),fwd_rxn_inds)
%         m.c(rxn_inds(i)) = 1;
%     else
%         m.c(rxn_inds(i)) = -1;
%     end
%     sol = solveLP(m);
%     max_flux(i,1) = sol.x(rxn_inds(i));
% end


%% Update rxn bounds to prevent M --> C proton pumping

% Three reactions are allowed to pump protons from [m] to [c]:
% Complex IV (HMR_6914), Complex III (HMR_6918), and Complex I (HMR_6921).
[~,allowed_inds] = ismember({'HMR_6914';'HMR_6918';'HMR_6921'},ihuman.rxns);
rxn_inds(ismember(rxn_inds,allowed_inds)) = [];  % don't constrain these reactions

% constrain the remaining reactions such that protons can only be moved
% from [c] to [m]
ihuman.ub(rxn_inds(ismember(rxn_inds,fwd_rxn_inds))) = 0;
ihuman.lb(rxn_inds(ismember(rxn_inds,rev_rxn_inds))) = 0;
ihuman.rev = double(ihuman.lb < 0 & ihuman.ub > 0);  % update ".rev" field

% For reactions that are now irreversible, but only in the direction that
% is backward from the way it is written, turn those reactions around, so
% they are written in the direction they proceed.
flip_ind = ihuman.lb < 0 & ihuman.ub == 0;
ihuman.ub(flip_ind) = -ihuman.lb(flip_ind);
ihuman.lb(flip_ind) = 0;
ihuman.S(:,flip_ind) = -ihuman.S(:,flip_ind);
fprintf('The following %u rxns have been re-written in the reverse direction:\n',sum(flip_ind));
fprintf('\t%s\n',ihuman.rxns{flip_ind});


% This resulted in 2 rxns with lb and ub equal to zero: 
%   HMR_6455:  5-carboxy-gamma-chromanol[m] + H+[m] => 5-carboxy-gamma-chromanol[c] + H+[c]
%      r1330:  H+[m] => H+[c] + Proton-Gradient[m]
%
% Both of these rxns are dead-end (cannot carry flux), so they should be
% deleted.
del_rxns = {'HMR_6455';'r1330'};



%% Update rxns involving C <--> M proton transport to be I <--> M
% protons in the mitochondrial compartment (m) will now transport to/from 
% the inner mitochrondrial membrane compartment (i)

% identify all rxns involving C <--> M proton transport (either direction)
rxn_inds = find( (ihuman.S(Hm,:) ~= 0) & (sign(ihuman.S(Hm,:)) == -sign(ihuman.S(Hc,:))) );

% change the compartment of all cytoplasmic protons (C) to inner mitochondrial membrane (I)
ihuman.S(Hi,rxn_inds) = ihuman.S(Hc,rxn_inds);
ihuman.S(Hc,rxn_inds) = 0;



%% Ensure that all rxns involving ATP-driven transport cannot generate ATP
% For example, met[c] + ATP[c] + H2O[c] --> met[p] + ADP[c] + Pi[c] should
% not be able to operate in the reverse direction.

% identify all rxns involving transport of any met(s) between compartments
reactants = arrayfun(@(i) ihuman.metNames(ihuman.S(:,i) < 0),(1:length(ihuman.rxns))','UniformOutput',false);
products  = arrayfun(@(i) ihuman.metNames(ihuman.S(:,i) > 0),(1:length(ihuman.rxns))','UniformOutput',false);
trans_inds = arrayfun(@(i) any(ismember(reactants{i},products{i})),(1:length(ihuman.rxns))');

% identify all rxns involving ATP <--> ADP conversion
atp_inds = find(ismember(ihuman.metNames,'ATP'));
adp_inds = find(ismember(ihuman.metNames,'ADP'));
atp_consume = (any(ihuman.S(atp_inds,:) < 0) & any(ihuman.S(adp_inds,:) > 0))';
atp_produce = (any(ihuman.S(atp_inds,:) > 0) & any(ihuman.S(adp_inds,:) < 0))';

% identify transport rxns that can produce ATP
atp_trans_inds = find(trans_inds & (atp_consume & (ihuman.lb < 0)) | ...
                      trans_inds & (atp_produce & (ihuman.ub > 0)) );
% Results in 9 reactions.

% some reactions are allowed (ATP synthase, and ADP/ADP translocases)
allowed_rxns = {'HMR_6916';'HMR_6328';'HMR_4908';'HMR_4907';'HMR_4906';'ATPtg';'ATPS4mi'};
atp_trans_inds(ismember(atp_trans_inds,find(ismember(ihuman.rxns,allowed_rxns)))) = [];

% Results in 2 remaining reactions: 
%        HMR_1696: 3alpha,7alpha-dihydroxy-5beta-cholestanate[c] + ATP[c] + H2O[c] <=> 3alpha,7alpha-dihydroxy-5beta-cholestanate[p] + ADP[c] + Pi[c]
%  4GLU56DIHDINDt: ATP[c] + H2O[c] + 4-S-Glutathionyl-5,6-Dihydroxyindoline[c] <=> ADP[c] + H+[c] + Pi[c] + 4-S-Glutathionyl-5,6-Dihydroxyindoline[s]

% restrict the bounds on the remaining reactions to prevent ATP generation
ihuman.lb(atp_trans_inds(ismember(atp_trans_inds,find(atp_consume)))) = 0;
ihuman.ub(atp_trans_inds(ismember(atp_trans_inds,find(atp_produce)))) = 0;

% update ".rev" field in the model, if necessary
ihuman.rev = double(ihuman.lb < 0 & ihuman.ub > 0);



%% Fix other problematic reactions

% This reaction is similar to an existing HMR rxn, but is reversible and
% missing the FAD(H2) metabolite:
%   RE1573M: cis,cis-3,6-dodecadienoyl-CoA[m] <=> trans,cis-lauro-2,6-dienoyl-CoA[m]
%  HMR_3288: cis,cis-3,6-dodecadienoyl-CoA[m] + FAD[m] => FADH2[m] + trans,cis-lauro-2,6-dienoyl-CoA[m]
% Therefore, the Recon3D version of the reaction should be deleted.
del_rxns = {'RE1573M'};


% Reactions: r1330  'H+[m] => H+[c] + Proton-Gradient[m]'
%            r1331  'H+[c] => H+[s] + Proton-Gradient[c]'
% These reactions originate from Recon3D. They appear to be an attempt to
% fix ATP leakage, but will not be necessary/compatible with the updated
% handling of mitochondrial proton transport.
del_rxns = [del_rxns; {'r1330'; 'r1331'}];


% This Recon3D reaction does not exist in the literature or any databases:
%     r1453: proline[m] + ubiquinol[m] <=> 1-pyrroline-5-carboxylate[m] + 5 H+[m] + ubiquinone[m]
% It is also not present on the BiGG Database. It should thus be removed.
del_rxns = [del_rxns; {'r1453'}];




% These rxns should not be allowed to operate in the reverse direction:
%  RE3238C: (11Z)-eicosenoyl-CoA[c] + H2O[c] <=> cis-gondoic acid[c] + CoA[c] + H+[c]
%  RE3239C: (13Z)-docosenoyl-CoA[c] + H2O[c] <=> cis-erucic acid[c] + CoA[c] + H+[c]
%  RE2649C: H2O[c] + propanoyl-CoA[c] <=> CoA[c] + H+[c] + propanoate[c]
ihuman.lb(ismember(ihuman.rxns,{'RE3238C';'RE3239C';'RE2649C'})) = 0;



% These reactions concern the glycerol phosphate shuttle:
%  HMR_0483: DHAP[c] + ubiquinol[m] => sn-glycerol-3-phosphate[c] + ubiquinone[m]
%  HMR_0482: DHAP[c] + FADH2[c] => FAD[c] + sn-glycerol-3-phosphate[c]
%    r0202m: NAD+[m] + sn-glycerol-3-phosphate[m] => DHAP[m] + H+[m] + NADH[m]
% The HMR reactions are written in the wrong direction, and should be
% reversed. The Recon3D reaction (r0202m) is also in the wrong direction,
% but should also not take place in the mitochondria, and should therefore
% be deleted.
ihuman.S(:,ismember(ihuman.rxns,{'HMR_0483';'HMR_0482'})) = -ihuman.S(:,ismember(ihuman.rxns,{'HMR_0483';'HMR_0482'}));
del_rxns = [del_rxns; {'r0202m'}];





% get unique list of rxns to delete
del_rxns = unique(del_rxns);


% constrain both bounds of problematic reactions to zero
del_rxn_ind = ismember(ihuman.rxns,del_rxns);
ihuman.lb(del_rxn_ind) = 0;
ihuman.ub(del_rxn_ind) = 0;



%% Perform flux analysis to assess infinite ATP generation

model = simplifyModel(ihuman);
exch_inds = sum(model.S ~= 0) == 1;
model.ub(exch_inds) = 0;
model.lb(exch_inds) = 0;
model.c(:) = 0;

[~,ATPhydr_ind] = ismember('HMR_3964',ihuman.rxns);
model.c(ATPhydr_ind) = 1;







