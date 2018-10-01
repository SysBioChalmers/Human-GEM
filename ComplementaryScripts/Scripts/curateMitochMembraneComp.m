%
% FILE NAME:    curateMitochMembraneComp.m
% 
% DATE CREATED: 2018-09-17
%     MODIFIED: 2018-09-19
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script for analyzing/curating the inner mitochondrial matrix
%          compartment.
%


%% Load model

% load HumanGEM model if not already loaded
if ~exist('ihuman','var')
    load('humanGEM.mat');
end


%% Remove Recon3D ATP synthase, complex I,  reaction

% the HMR ATP synthase reaction will be kept, but we need to document the
% reaction association in the model and the rxnAssoc.mat file
%
%  *ATP synthase
%     HMR_6916: ADP[m] + 4 H+[c] + Pi[m] => ATP[m] + 4 H+[m] + H2O[m]
%      ATPS4mi: ADP[m] + 4 H+[i] + Pi[m] => ATP[m] + 3 H+[m] + H2O[m]
%
%  *Complex I
%     HMR_6921: 5 H+[m] + NADH[m] + ubiquinone[m] => NAD+[m] + ubiquinol[m] + 4 H+[c]
%  NADH2_u10mi: 5 H+[m] + NADH[m] + ubiquinone[m] => NAD+[m] + ubiquinol[m] + 4 H+[i]
%
%  *Complex III
%     HMR_6918: 2 ferricytochrome C[m] + 2 H+[m] + ubiquinol[m] => 2 ferrocytochrome C[m] + ubiquinone[m] + 4 H+[c]
%   CYOR_u10mi: 2 ferricytochrome C[m] + 2 H+[m] + ubiquinol[m] => 2 ferrocytochrome C[m] + ubiquinone[m] + 4 H+[i]
%
%  *Complex IV
%     HMR_6914: 4 ferrocytochrome C[m] + 8 H+[m] + O2[m] => 4 ferricytochrome C[m] + 2 H2O[m] + 4 H+[c]
%      CYOOm2i: 4 ferrocytochrome C[m] + 8 H+[m] + O2[m] => 4 ferricytochrome C[m] + 2 H2O[m] + 4 H+[i]
%

hmr_rxns = {'HMR_6916';'HMR_6921';'HMR_6918';'HMR_6914'};
[~,hmr_ind] = ismember(hmr_rxns,ihuman.rxns);
r3_rxns = {'ATPS4mi';'NADH2_u10mi';'CYOR_u10mi';'CYOOm2i'};
[~,r3_ind] = ismember(r3_rxns,ihuman.rxns);

% if some of the Recon3D reactions have already been removed from the
% model, then exclude those from the process
remove_ind = (r3_ind == 0);
hmr_rxns(remove_ind) = [];
hmr_ind(remove_ind) = [];
r3_rxns(remove_ind) = [];
r3_ind(remove_ind) = [];

if ~isempty(r3_ind)
    
    fprintf('%u reactions (%s) will be removed from the model\n',length(r3_ind),strjoin(r3_rxns,', '));
    fprintf('because they are duplicates of existing HMR reactions.\n\n');
    
    % update "rxnRecon3DID" field
    for i = 1:length(hmr_rxns)
        if isempty(ihuman.rxnRecon3DID{hmr_ind(i)})
            ihuman.rxnRecon3DID(hmr_ind(i)) = r3_rxns(i);
        else
            ihuman.rxnRecon3DID(hmr_ind(i)) = join(unique([strsplit(ihuman.rxnRecon3DID{hmr_ind(i)},';'), r3_rxns(i)]), ';');
        end
    end
    
    % load rxnAssoc.mat file
    load('ComplementaryScripts/modelIntegration/rxnAssoc.mat');
    
    % ignore any associations that already exist in rxnAssoc.mat
    remove_ind = ismember(join([hmr_rxns,r3_rxns],'***'),join([rxnAssoc.rxnHMRID,rxnAssoc.rxnRecon3DID],'***'));
    add_hmr_rxns = hmr_rxns(~remove_ind);
    add_hmr_ind = hmr_ind(~remove_ind);
    add_r3_rxns = r3_rxns(~remove_ind);
    add_r3_ind = r3_ind(~remove_ind);
    
    % add associations to rxnAssoc structure
    rxnAssoc.rxnHMRID = [rxnAssoc.rxnHMRID; add_hmr_rxns];
    rxnAssoc.rxnRecon3DID = [rxnAssoc.rxnRecon3DID; add_r3_rxns];
    rxnAssoc.lbHMR = [rxnAssoc.lbHMR; ihuman.lb(add_hmr_ind)];
    rxnAssoc.ubHMR = [rxnAssoc.ubHMR; ihuman.ub(add_hmr_ind)];
    rxnAssoc.lbRecon3D = [rxnAssoc.lbRecon3D; ihuman.lb(add_r3_ind)];
    rxnAssoc.ubRecon3D = [rxnAssoc.ubRecon3D; ihuman.ub(add_r3_ind)];

    % save new rxnAssoc.mat file
    fprintf('The rxnAssoc.mat file has been updated with rxns related to the electron transport chain.\n\n');
    save('rxnAssoc_new.mat','rxnAssoc');
    
    % delete Recon3D reactions from model
    ihuman = removeReactionsFull(ihuman,r3_rxns);
    
end



%% Analysis of reactions transporting protons into/out of the mitochondria
% Identify all reactions that involve proton transport between the
% cytoplasm [c] and mitochondria [m].

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

fprintf('A total of %u reactions are able to transport protons from [m] to [c].\n\n',length(rxn_inds));

% This resulted in 9 fwd rxns, and 21 rev rxns, for a total of 30 rxns


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


%% Update rxn bounds to prevent M --> C proton pumping for electrogenic rxns
% Reactions that involve proton transport between the mitochondria [m] and
% the cytoplasm [c] will not be modified as long as they are
% non-electrogenic (i.e., do not involve a net change in either compartment
% charge). If they ARE electrogenic, then proton pumping from [m] to [c]
% should only be permitted for the electron transport chain, but restricted
% otherwise.

% determine which of the M --> C proton pump reactions are electrogenic
charge_diff = crossCompChargeDiff(ihuman,rxn_inds,false);
non_elec = all(charge_diff.mat == 0,2);
fprintf('%u of the %u rxns able to transport protons from [m] to [c] are non-electrogenic,\n',sum(non_elec),length(rxn_inds));
fprintf('and will therefore be permitted.\n\n');
rxn_inds(non_elec) = [];  % don't constrain non-electrogenic transport rxns

% Three reactions are allowed to pump protons from [m] to [c]:
% Complex IV (HMR_6914), Complex III (HMR_6918), Complex I (HMR_6921).
[~,allowed_inds] = ismember({'HMR_6914';'HMR_6918';'HMR_6921'},ihuman.rxns);
non_elec(ismember(rxn_inds,allowed_inds)) = [];
rxn_inds(ismember(rxn_inds,allowed_inds)) = [];  % don't constrain these reactions
fprintf('ETC rxns involving Complex IV (HMR_6914), Complex III (HMR_6918), and Complex I (HMR_6921)\n');
fprintf('are also permitted to transport protons from [m] to [c].\n\n');


% constrain the remaining electrogenic reactions such that protons can only
% be moved from [c] to [m]
fprintf('The remaining %u electrogenic, proton-transporting rxns will be constrained\n',length(rxn_inds));
fprintf('to prevent their proton transport from [m] to [c].\n\n');
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
fprintf('\n');


% This resulted in one rxn with lb and ub equal to zero: 
%   
%      r1330:  H+[m] => H+[c] + Proton-Gradient[m]
%
% This rxn is anyway dead-end (cannot carry flux), so it should be deleted
% from the model.
del_rxns = {'r1330'};



%% Update electrogenic rxns involving C <--> M proton transport to be I <--> M
% Protons in the mitochondrial compartment [m] will now transport to/from 
% the inner mitochrondrial membrane compartment [i] for reactions that are
% electrogenic.
% NOTE: Reactions that transport protons from [c] to [m] (and vice versa)
% but are NOT electrogenic will remain in the original form (i.e., the
% cytoplasmic proton will NOT be changed to the [i] compartment).

% identify all electrogenic rxns involving C <--> M proton transport (either direction)
charge_diff = crossCompChargeDiff(ihuman,[],false);
electrogenic = any(charge_diff.mat ~= 0, 2);
rxn_inds = find( (ihuman.S(Hm,:) ~= 0) & (sign(ihuman.S(Hm,:)) == -sign(ihuman.S(Hc,:))) & electrogenic' );
num_non_elec = full(sum( (ihuman.S(Hm,:) ~= 0) & (sign(ihuman.S(Hm,:)) == -sign(ihuman.S(Hc,:))) & ~electrogenic' ));

% change the compartment of all cytoplasmic protons (C) to inner mitochondrial membrane (I)
ihuman.S(Hi,rxn_inds) = ihuman.S(Hc,rxn_inds);
ihuman.S(Hc,rxn_inds) = 0;
fprintf('%u electrogenic rxns involving proton transport between [m] and [c] (either direction) were modified\n',length(rxn_inds));
fprintf('to replace cytoplasmic protons (H+[c]) with inner mitochondrial membrane protons (H+[i]).\n\n');
fprintf('The remaining %u non-electrogenic rxns involving proton transport between [m] and [c] will not be\n',num_non_elec);
fprintf('modified, and their cytoplasmic protons (H+[c]) will be left as is.\n\n');



%% Add a reaction allowing transport of protons from I --> C
% This reaction does not generate energy, and is primarily to allow mass
% balancing of protons in the model.
rxnsToAdd = {};
rxnsToAdd.rxns = {'Htransport_I_C'};  % this needs to be replaced with a systematic name
rxnsToAdd.equations = {'H+[i] => H+[c]'};
rxnsToAdd.rxnNames = rxnsToAdd.rxns;
rxnsToAdd.lb = 0;
rxnsToAdd.ub = 1000;
ihuman = addRxns(ihuman,rxnsToAdd,3,[],false);




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
allowed_rxns = {'HMR_6916';'HMR_6328';'HMR_4908';'HMR_4907';'HMR_4906';'ATPtg'};
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


% The following HMR reaction involves a proton from a compartment other
% than the one containing the other metabolites in the reaction:
%   HMR_8616: ATP[m] + cob(I)alamin[m] + H+[c] <=> cobamide-coenzyme[m] + triphosphate[m]
% However, there is no evidence in the literature for such transport, and
% it is therefore suspected to be a compartment labeling error. The proton
% compartment should therefore be corrected to mitochondria [m] to be
% consistent with the other compounds in the reaction.
[~,comp_inds] = ismember({'cytosol','mitochondria'},lower(ihuman.compNames));
Hc = find( ismember(ihuman.metNames,'H+') & (ihuman.metComps == comp_inds(1)) );
Hm = find( ismember(ihuman.metNames,'H+') & (ihuman.metComps == comp_inds(2)) );
[~,rxn_ind] = ismember('HMR_8616',ihuman.rxns);
ihuman.S(Hc,rxn_ind) = 0;
ihuman.S(Hm,rxn_ind) = -1;


% This reaction is similar to an existing HMR rxn, but is reversible and
% missing the FAD(H2) metabolite:
%    RE1573M: cis,cis-3,6-dodecadienoyl-CoA[m] <=> trans,cis-lauro-2,6-dienoyl-CoA[m]
%   HMR_3288: cis,cis-3,6-dodecadienoyl-CoA[m] + FAD[m] => FADH2[m] + trans,cis-lauro-2,6-dienoyl-CoA[m]
% Therefore, the Recon3D version of the reaction should be deleted.
del_rxns = [del_rxns; {'RE1573M'}];


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


% The following reactions are part of the melatonin degradation pathway:
%   HMR_4551: formyl-N-acetyl-5-methoxykynurenamine[c] + 2 H+[c] + H2O2[c] <=> formate[c] + H2O[c] + N-acetyl-5-methoxykynuramine[c]'
%    RE2440C: 2 formyl-N-acetyl-5-methoxykynurenamine[c] + H2O2[c] <=> CO2[c] + formate[c] + H+[c] + 2 N-acetyl-5-methoxykynuramine[c]
% Although they are written as reversible, the reverse reaction should not
% be possible (see, e.g., PMID: 19573038). Therefore, these reactions
% should be constrained to proceed only in the forward direction.
ihuman.lb(ismember(ihuman.rxns,{'HMR_4551';'RE2440C'})) = 0;



% The following reaction from Recon3D:
%   r1466: gamma-linolenoyl-CoA[r] + 4 H+[r] + 2 malonyl-CoA[r] + 2 NADPH[r] + 2 O2[r] <=> 4 CO2[r] + 2 CoA[r] + dihomo-gamma-linolenoyl-CoA[r] + 2 H2O[r] + 2 NADP+[r]
% should NOT be reversible. Only the forward direction is possible.
% Therefore, the lower bound of this reaction will be set to zero.
ihuman.lb(ismember(ihuman.rxns,{'r1466'})) = 0;



% The following reaction from Recon3D:
%   r0698: chenodeoxycholoyl-CoA[p] + 4 H+[p] + propanoyl-CoA[p] => 25(R)DHCA-CoA[p] + CoA[p] + H2O[p]
% has no evidence supporting its existence, and is charge-imbalanced (there
% is no source providing the electrons to reduce the four protons). It will
% therefore be removed from the model.
del_rxns = [del_rxns; {'r0698'}];





% The following reactions from Recon3D:
%
%   ARTFR61: FADH2[m] + (2E)-hexadecenoyl-CoA[c] => FAD[m] + R Group 6 Coenzyme A[c]
%     RTOT6: R Group 6 Coenzyme A[c] => R Total Coenzyme A[c]
%   ARTPLM1: R Total Coenzyme A[c] => palmitoyl-CoA[c]
%
% effectively sum to the combined reaction:
%
%            FADH2[m] + (2E)-hexadecenoyl-CoA[c] => FAD[m] + palmitoyl-CoA[c]
%
% However, this reaction requires NADPH to proceed, and there even exists a
% reaction (from Recon3D) in the model that does just this:
%
%     r0309: NADP+[m] + palmitoyl-CoA[m] <=> (2E)-hexadecenoyl-CoA[m] + H+[m] + NADPH[m]
%
% Therefore, the cofactor in the first rxn above (ARTFR61) should be
% changed to NADPH instead of FADH2. The reaction should probably be
% deleted entirely due to it being redundant, but this fix is sufficient
% for now. The same situation was found for the reactions:
%
%   ARTFR46: (2E)-octadecenoyl-CoA[c] + FADH2[m] => FAD[m] + 1.125 R Group 4 Coenzyme A[c]
%   ARTFR32: FADH2[m] + oleoyl-CoA[c] => FAD[m] + 1.125 R Group 3 Coenzyme A[c]
%   ARTFR41: FADH2[m] + palmitoleoyl-CoA[c] => FAD[m] + R Group 4 Coenzyme A[c]
%   ARTFR42: FADH2[m] + oleoyl-CoA[c] => FAD[m] + 1.125 R Group 4 Coenzyme A[c]
%   ARTFR12: FADH2[m] + palmitoleoyl-CoA[c] => FAD[m] + R Group 1 Coenzyme A[c]
%   ARTFR33: FADH2[m] + 11-Octadecenoyl Coenzyme A[c] => FAD[m] + 1.125 R Group 3 Coenzyme A[c]
%   ARTFR34: (6Z,9Z)-octadecadienoyl-CoA[c] + 2 FADH2[m] => 2 FAD[m] + 1.125 R Group 3 Coenzyme A[c]
%   ARTFR43: FADH2[m] + 11-Octadecenoyl Coenzyme A[c] => FAD[m] + 1.125 R Group 4 Coenzyme A[c]
%   ARTFR44: (6Z,9Z)-octadecadienoyl-CoA[c] + 2 FADH2[m] => 2 FAD[m] + 1.125 R Group 4 Coenzyme A[c]
%   ARTFR45: (15Z)-tetracosenoyl-CoA[c] + FADH2[m] => FAD[m] + 1.5 R Group 4 Coenzyme A[c]

[~,comp_num] = ismember('mitochondria',lower(ihuman.compNames));
fadh2_ind = find(ismember(ihuman.metNames,'FADH2') & (ihuman.metComps == comp_num));
fad_ind = find(ismember(ihuman.metNames,'FAD') & (ihuman.metComps == comp_num));
nadp_ind = find(ismember(ihuman.metNames,'NADP+') & (ihuman.metComps == comp_num));
nadph_ind = find(ismember(ihuman.metNames,'NADPH') & (ihuman.metComps == comp_num));
h_ind = find(ismember(ihuman.metNames,'H+') & (ihuman.metComps == comp_num));

[~,rxn_ind] = ismember({'ARTFR61';'ARTFR46';'ARTFR32';'ARTFR41';'ARTFR42';'ARTFR12';'ARTFR33';'ARTFR34';'ARTFR43';'ARTFR44';'ARTFR45'},ihuman.rxns);
for i = 1:length(rxn_ind)
    if ihuman.S(fadh2_ind,rxn_ind(i)) ~= 0  % check that the rxn hasn't been changed already
        ihuman.S([fadh2_ind; fad_ind], rxn_ind(i)) = 0;
        ihuman.S([nadp_ind; nadph_ind; h_ind], rxn_ind(i)) = [1; -1; -1];
    end
end



% The model contains the following two reactions from Recon3D:
%   DHCR241r: FADH2[r] + zymosterol[r] => FAD[r] + cholestenol[r]
%      r1380: H+[r] + NADPH[r] + zymosterol[r] <=> NADP+[r] + cholestenol[r]
% The reactions are identical, except one uses FADH2, whereas the other
% uses NADPH. Based on the available databases, the reaction should use
% NADPH. Therefore, the first reaction should be deleted, since it will be
% identical to the second after replacing its cofactor with NADPH.
% The same is true for the following reactions:
%   DHCR242r: FADH2[r] + 5alpha-cholesta-7,24-dien-3beta-ol[r] => FAD[r] + lathosterol[r]
del_rxns = [del_rxns; {'DHCR241r'; 'DHCR242r'}];


% The following reaction from Recon3D:
%   r1479: CoA[m] + 3-Oxolaur-Cis-5-Enoyl Coenzyme A[m] => (3Z)-dodecenoyl-CoA[m] + acetyl-CoA[m]
% was found as part of a reaction loop creating energy and mass. This
% reaction is imbalanced (difference of C2H2), and was not found on the
% BiGG database, suggesting that it has been removed or updated recently.
% Therefore, this reaction should be removed from the model.
del_rxns = [del_rxns; {'r1479'}];


% The following reaction from Recon3D:
%   DOPACCL: dopamine-O-quinone[c] => H+[c] + leukoaminochrome[c]
% has no gene associated with the reaction, and no sources, and it could
% not be found in the literature. A paper (PMID: 20600874) shows part of
% the pathway of dopamine-derived quinone metabolism, where it is clear
% that the above reaction would not take place. Therefore it should be
% removed from the model.
del_rxns = [del_rxns; {'DOPACCL'}];


% The following reaction from Recon3D:
%   ALPA_HSx: acylglycerone-phosphate[p] + H+[p] + NADP+[p] => NADPH[p] + Lysophosphatidic Acid[p]
% should consume NADPH, not produce it (see KEGG rxn R02756, where
% "Lysophosphatidic Acid" is also known as 1-Acyl-sn-glycerol 3-phosphate.
[~,comp_num] = ismember('peroxisome',lower(ihuman.compNames));
nadp_ind = find(ismember(ihuman.metNames,'NADP+') & (ihuman.metComps == comp_num));
nadph_ind = find(ismember(ihuman.metNames,'NADPH') & (ihuman.metComps == comp_num));
[~,rxn_ind] = ismember('ALPA_HSx',ihuman.rxns);
if ihuman.S(nadp_ind,rxn_ind) == -1  % check that the rxn hasn't been changed already
    ihuman.S([nadp_ind, nadph_ind], rxn_ind) = [1,-1];
end



% The following reaction from Recon3D:
%
%   DOLGPP_Ler: 0.1 dolichyl-D-glucosyl-phosphate[r] + H2O[r] => 0.1 dolichyl-phosphate[r] + glucose[r] + H+[r]
%
% Is not properly formulated, as it is creating 1 equivalent of glucose
% from 0.1 equivalents. It is nearly the same as the following reaction:
%
%   HMR_8692: dolichyl-D-glucosyl-phosphate[r] + H2O[r] => dolichyl-phosphate[r] + glucose[r]
%
% except for the coefficients and the additional proton.
%
% *NOTE: This appears to be a problem with the difference in formula for
% dolichyl-phosphate between Recon3D and HMR:
%
%   Recon3D formula: C1080H1758O40P10
%       HMR formula: C20H37O4P(C5H8)n
%
% Other pairs of reactions that share this same problem are:
%
%   DOLASNT_Ler: 0.1 (Glc)3 (GlcNAc)2 (Man)9 (PP-Dol)1[r] + [protein]-L-asparagine[r] => (alpha-D-Glucosyl)3-(alpha-D-mannosyl)8-beta-D-mannosyl-diacetylchitobiosyl-L-asparagine (protein)[r] + 0.1 dolichyl-diphosphate[r] + H+[r]
%      HMR_7285:     (Glc)3 (GlcNAc)2 (Man)9 (PP-Dol)1[r] + [protein]-L-asparagine[r] => (alpha-D-Glucosyl)3-(alpha-D-mannosyl)8-beta-D-mannosyl-diacetylchitobiosyl-L-asparagine (protein)[r] +     dolichyl-diphosphate[r]
%
%    DOLDPP_Ler: 0.1 dolichyl-diphosphate[r] + H2O[r] => 0.1 dolichyl-phosphate[r] + H+[r] + Pi[r]
%      HMR_8691:     dolichyl-diphosphate[r] + H2O[r] =>     dolichyl-phosphate[r] + Pi[r]
%
%        DOLK_L: CTP[c] + 0.1 dolichol[c] => CDP[c] + 0.1 dolichyl-phosphate[c] + H+[c]
%      HMR_7263: CTP[c] +     dolichol[c] => CDP[c] +     dolichyl-phosphate[c]
% 	
%  DOLMANP_Lter: 0.1 dolichyl-phosphate-D-mannose[c] => 0.1 dolichyl-phosphate-D-mannose[r]
%      HMR_7272:     dolichyl-phosphate-D-mannose[c] =>     dolichyl-phosphate-D-mannose[r]
% 	
%   DOLPMT3_Ler: 0.1 dolichyl-phosphate[c] + GDP-mannose[c] => 0.1 dolichyl-phosphate-D-mannose[c] + GDP[c]
%      HMR_7271:     dolichyl-phosphate[c] + GDP-mannose[c] =>    dolichyl-phosphate-D-mannose[c] + GDP[c]
% 	
%     GPIMTer_L: 0.1 dolichyl-phosphate-D-mannose[r] + glucosaminyl-acylphosphatidylinositol[r] => 0.1 dolichyl-phosphate[r] + H+[r] + mgacpail heparan sulfate[r]
%      HMR_8383:     dolichyl-phosphate-D-mannose[r] + glucosaminyl-acylphosphatidylinositol[r] =>     dolichyl-phosphate[r] + mgacpail heparan sulfate[r]
% 	
%    GLCNACPT_L: 0.1 dolichyl-phosphate[c] + UDP-N-acetylglucosamine[c] => 0.1 N-acetyl-D-glucosaminyldiphosphodolichol[c] + UMP[c]
%      HMR_7264:     dolichyl-phosphate[c] + UDP-N-acetylglucosamine[c] =>     N-acetyl-D-glucosaminyldiphosphodolichol[c] + UMP[c]
%
%   DOLPGT3_Ler: 0.1 dolichyl-phosphate[r] + H2O[r] => 0.1 dolichol[r] + Pi[r]
%      HMR_7261:     dolichyl-phosphate[r] + H2O[r] =>     dolichol[r] + Pi[r]
%
% DOLICHOL_Lter: 0.1 dolichol[r] <=> 0.1 dolichol[c]
%      HMR_7262:     dolichol[c] <=>     dolichol[r]
%
%      DEDOLR_L: 0.1 dehydrodolichol[c] + H+[c] + NADPH[c] => 0.1 dolichol[c] + NADP+[c]
%      HMR_7260:     dehydrodolichol[c] + H+[c] + NADPH[c] =>     dolichol[c] + NADP+[c]
%
% where again the Recon3D version is using 0.1 equivalents of the dolichyl 
% component, which creates mass when used together with the HMR reactions.
% Therefore, these Recon3D reactions should be removed from the model.
%
% In addition, there were a few other reactions that did not have an HMR
% equivalent, but should have their stoich coeffs adjusted from 0.1 to 1
%
%      H8MTer_L: 0.1 dolichyl-phosphate-D-mannose[r] + gpi heparan sulfate[r] => 0.1 dolichyl-phosphate[r] + H+[r] + Mannosyl-3-(Phosphoethanolaminyl-Mannosyl)-Glucosaminyl-Acylphosphatidylinositol (M4A)[r]
%      H8MTer_U: 0.1 dolichyl-phosphate-D-mannose[r] + gpi heparan sulfate[r] => 0.1 dolichyl-phosphate[r] + H+[r] + HMA[r]
%    UDPDOLPT_L: 0.1 dolichyl-phosphate[c] + UDP-glucose[c] => UDP[c] + 0.1 dolichyl-D-glucosyl-phosphate[c]
%
del_rxns = [del_rxns; {'DOLGPP_Ler';'DOLASNT_Ler';'DOLDPP_Ler';'DOLK_L';'DOLMANP_Lter';...
                       'DOLPMT3_Ler';'GPIMTer_L';'GLCNACPT_L';'DOLPGT3_Ler';'DOLICHOL_Lter';'DEDOLR_L'}];
rxn_ind = getIndexes(ihuman,{'H8MTer_L';'H8MTer_U';'UDPDOLPT_L'},'rxns');
ihuman.S(:,rxn_ind) = sign(ihuman.S(:,rxn_ind));  % convert all nonzero values to +/- 1


% The following reaction from Recon3D:
%   TAG_HSad_E: (4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoyl-CoA[c] + (5Z,8Z,11Z,14Z,17Z)-eicosapentaenoyl-CoA[c] + arachidonyl-CoA[c] + 2 H2O[c] + linolenoyl-CoA[c] + linoleoyl-CoA[c] + 2 sn-glycerol-3-phosphate[c] => 6 CoA[c] + 2 Pi[c] + 2 TAG-VLDL pool[c]
% is has one extra CoA in the products compared to reactants. The
% stoichiometric coefficient of the product CoA should therefore be changed
% to 5 to correct this imbalance.
met_ind = getIndexes(ihuman,'CoA[c]','metscomps');
rxn_ind = getIndexes(ihuman,'TAG_HSad_E','rxns');
if ihuman.S(met_ind,rxn_ind) == 6
    ihuman.S(met_ind,rxn_ind) = 5;
end



% The following reaction from HMR:
%
%   HMR_5238: LDL remnant[l] => 25 2-lysolecithin pool[l] + 110 CDP-diacylglycerol-LD-PI pool[l] + 425 PC-LD pool[l] + 30 PE-LD pool[l] + 160 SM pool[l] + apoB100[l] + 1515 cholesterol-ester pool[l] + 680 cholesterol[l]
%
% is nearly identical to the reverse of another reaction:
%
%   HMR_5239: 25 2-lysolecithin pool[r] + apoB100[r] + 110 CDP-diacylglycerol-LD-PI pool[r] + 680 cholesterol[r] + 1515 cholesterol-ester pool[r] + 425 PC-LD pool[r] + 30 PE-LD pool[r] + 160 SM pool[r] => LDL[r]
%
% However, the first reaction consumes LDL remnant, whereas the second
% reaction produces LDL, which is a problem because the following reaction
% also exists:
%
%   HMR_0014: LDL[s] => 680 cholesterol[s] + 1515 cholesterol-ester pool[s] + LDL remnant[s]
%
% Therefore, these three reactions can be used to create cholesterol and
% cholesterol-ester pool metabolites from nothing. To fix this, the LDL
% remnant should not produce cholesterol and cholesterol-ester pool
% (i.e., HMR_5238 should have these mets removed from its stoichiometry).
% This exact same situtaion was also found for HDL and HDL remnant,
% requiring the modification of the following rxn as well:
%
%   HMR_5233: HDL remnant[l] => 90 PC-LD pool[l] + 25 PE-LD pool[l] + 30 PS-LD pool[l] + 75 SM pool[l] + 2 apoA1[l] + 160 cholesterol-ester pool[l] + 20 cholesterol[l]
met_ind = getIndexes(ihuman,{'cholesterol-ester pool[l]';'cholesterol[l]'},'metscomps');
rxn_ind = getIndexes(ihuman,{'HMR_5238';'HMR_5233'},'rxns');
ihuman.S(met_ind,rxn_ind) = 0;


% The following reactions from Recon3D:
%   AGPAT1: 2 H+[c] + R Total 2 Coenzyme A[c] + Lysophosphatidic Acid[c] => CoA[c] + phosphatidate-LD-TAG pool[c]
%   AGPAT2: palmitoyl-CoA[c] + Lysophosphatidic Acid[c] => CoA[c] + phosphatidate-LD-TAG pool[c]
%   AGPAT3: oleoyl-CoA[c] + Lysophosphatidic Acid[c] => CoA[c] + phosphatidate-LD-TAG pool[c]
%   AGPAT4: linoleoyl-CoA[c] + Lysophosphatidic Acid[c] => CoA[c] + phosphatidate-LD-TAG pool[c]
% are inconsistent in their formulas, leading to a product which now is
% ambiguous in mass. Therefore, these reactions should be removed from the
% model (or properly balanced in the future).
del_rxns = [del_rxns; {'AGPAT1';'AGPAT2';'AGPAT3';'AGPAT4'}];





% get unique list of rxns to delete
del_rxns = unique(del_rxns);


% constrain both bounds of problematic reactions to zero
del_rxn_ind = ismember(ihuman.rxns,del_rxns);
ihuman.lb(del_rxn_ind) = 0;
ihuman.ub(del_rxn_ind) = 0;


%% Perform flux analysis to assess infinite ATP generation

model = simplifyModel(ihuman);
exch_inds = sum(model.S ~= 0) == 1;
model.S(:,exch_inds) = -abs(model.S(:,exch_inds));  % make all exch rxns the same direction (met --> 0)
model.ub(exch_inds) = 1000;  % allow only export
model.lb(exch_inds) = 0;


% Set the new objective
% [~,obj_ind] = ismember('HMR_3964',model.rxns);  % ATP hydrolysis
% [~,obj_ind] = ismember('HMR_6916',model.rxns);  % ATP synthase
% [~,obj_ind] = ismember('HMR_9034',model.rxns);  % glucose export
[~,obj_ind] = ismember('HMR_9058',model.rxns);  % CO2 export
model.c(:) = 0;
model.c(obj_ind) = 1;



% temporarily constrain some reactions for investigation
constrain_rxns = {'CYOOm3i';    % 4 ferrocytochrome C[m] + 7.92 H+[m] + O2[m] => 4 ferricytochrome C[m] + 1.96 H2O[m] + 0.02 O2-[m] + 4 H+[i]
                  'HMR_0685';   % fatty acid-LD-TG1 pool[c] => 0.0001 (10Z)-heptadecenoic acid[c] + 0.0001 (11Z,14Z)-eicosadienoic acid[c] + 0.0001 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0001 (13Z)-eicosenoic acid[c] + 0.0001 (13Z)-octadecenoic acid[c] + 0.0001 (13Z,16Z)-docosadienoic acid[c] + 0.0001 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0001 (6Z,9Z)-octadecadienoic acid[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0001 (7Z)-octadecenoic acid[c] + 0.0001 (7Z)-tetradecenoic acid[c] + 0.0001 (9E)-tetradecenoic acid[c] + 0.0001 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0001 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0001 10,13,16,19-docosatetraenoic acid[c] + 0.0001 10,13,16-docosatriynoic acid[c] + 0.0001 12,15,18,21-tetracosatetraenoic acid[c] + 0.0001 13,16,19-docosatrienoic acid[c] + 0.0001 7-palmitoleic acid[c] + 0.0001 8,11-eicosadienoic acid[c] + 0.0001 9-eicosenoic acid[c] + 0.0001 9-heptadecylenic acid[c] + 0.0005 DHA[c] + 0.0017 DPA[c] + 0.0013 EPA[c] + 0.0001 adrenic acid[c] + 0.0024 arachidonate[c] + 0.0001 behenic acid[c] + 0.0001 cerotic acid[c] + 0.0001 cis-cetoleic acid[c] + 0.0001 cis-erucic acid[c] + 0.0001 cis-gondoic acid[c] + 0.0208 cis-vaccenic acid[c] + 0.0013 dihomo-gamma-linolenate[c] + 0.0001 eicosanoate[c] + 0.0001 elaidate[c] + 0.0048 gamma-linolenate[c] + 0.0001 henicosanoic acid[c] + 0.0001 lauric acid[c] + 0.0001 lignocerate[c] + 0.0741 linoleate[c] + 0.0037 linolenate[c] + 0.0001 margaric acid[c] + 0.0001 mead acid[c] + 0.0096 myristic acid[c] + 0.0001 nervonic acid[c] + 0.0001 nonadecylic acid[c] + 0.122 oleate[c] + 0.0001 omega-3-arachidonic acid[c] + 0.64 palmitate[c] + 0.0286 palmitolate[c] + 0.0001 pentadecylic acid[c] + 0.0001 physeteric acid[c] + 0.082 stearate[c] + 0.0028 stearidonic acid[c] + 0.0001 tricosanoic acid[c] + 0.0001 tridecylic acid[c] + 0.0001 ximenic acid[c]
                  'HMR_0686';   % fatty acid-LD-TG2 pool[c] => 0.0001 (10Z)-heptadecenoic acid[c] + 0.0001 (11Z,14Z)-eicosadienoic acid[c] + 0.0001 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0001 (13Z)-eicosenoic acid[c] + 0.0001 (13Z)-octadecenoic acid[c] + 0.0001 (13Z,16Z)-docosadienoic acid[c] + 0.0016 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0001 (6Z,9Z)-octadecadienoic acid[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0001 (7Z)-octadecenoic acid[c] + 0.0001 (7Z)-tetradecenoic acid[c] + 0.0001 (9E)-tetradecenoic acid[c] + 0.0001 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0001 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0001 10,13,16,19-docosatetraenoic acid[c] + 0.0001 10,13,16-docosatriynoic acid[c] + 0.0001 12,15,18,21-tetracosatetraenoic acid[c] + 0.0001 13,16,19-docosatrienoic acid[c] + 0.0001 7-palmitoleic acid[c] + 0.0001 8,11-eicosadienoic acid[c] + 0.0001 9-eicosenoic acid[c] + 0.0001 9-heptadecylenic acid[c] + 0.0081 DHA[c] + 0.0018 DPA[c] + 0.0003 EPA[c] + 0.0032 adrenic acid[c] + 0.0521 arachidonate[c] + 0.0001 behenic acid[c] + 0.0001 cerotic acid[c] + 0.0001 cis-cetoleic acid[c] + 0.0001 cis-erucic acid[c] + 0.0001 cis-gondoic acid[c] + 0.0217 cis-vaccenic acid[c] + 0.0023 dihomo-gamma-linolenate[c] + 0.0001 eicosanoate[c] + 0.0001 elaidate[c] + 0.003 gamma-linolenate[c] + 0.0001 henicosanoic acid[c] + 0.0001 lauric acid[c] + 0.0001 lignocerate[c] + 0.281 linoleate[c] + 0.0119 linolenate[c] + 0.0001 margaric acid[c] + 0.0001 mead acid[c] + 0.0024 myristic acid[c] + 0.0001 nervonic acid[c] + 0.0001 nonadecylic acid[c] + 0.4226 oleate[c] + 0.0001 omega-3-arachidonic acid[c] + 0.12 palmitate[c] + 0.0229 palmitolate[c] + 0.0001 pentadecylic acid[c] + 0.0001 physeteric acid[c] + 0.0384 stearate[c] + 0.0025 stearidonic acid[c] + 0.0001 tricosanoic acid[c] + 0.0001 tridecylic acid[c] + 0.0001 ximenic acid[c]
                  'HMR_0687';   % fatty acid-LD-TG3 pool[c] => 0.0001 (10Z)-heptadecenoic acid[c] + 0.0001 (11Z,14Z)-eicosadienoic acid[c] + 0.0001 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0001 (13Z)-eicosenoic acid[c] + 0.0001 (13Z)-octadecenoic acid[c] + 0.0001 (13Z,16Z)-docosadienoic acid[c] + 0.0011 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0001 (6Z,9Z)-octadecadienoic acid[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0001 (7Z)-octadecenoic acid[c] + 0.0001 (7Z)-tetradecenoic acid[c] + 0.0001 (9E)-tetradecenoic acid[c] + 0.0001 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0001 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0001 10,13,16,19-docosatetraenoic acid[c] + 0.0001 10,13,16-docosatriynoic acid[c] + 0.0001 12,15,18,21-tetracosatetraenoic acid[c] + 0.0001 13,16,19-docosatrienoic acid[c] + 0.0001 7-palmitoleic acid[c] + 0.0001 8,11-eicosadienoic acid[c] + 0.0001 9-eicosenoic acid[c] + 0.0001 9-heptadecylenic acid[c] + 0.0221 DHA[c] + 0.0038 DPA[c] + 0.0035 EPA[c] + 0.0022 adrenic acid[c] + 0.0264 arachidonate[c] + 0.0001 behenic acid[c] + 0.0001 cerotic acid[c] + 0.0001 cis-cetoleic acid[c] + 0.0001 cis-erucic acid[c] + 0.0001 cis-gondoic acid[c] + 0.0194 cis-vaccenic acid[c] + 0.0141 dihomo-gamma-linolenate[c] + 0.0001 eicosanoate[c] + 0.0001 elaidate[c] + 0.0081 gamma-linolenate[c] + 0.0001 henicosanoic acid[c] + 0.0001 lauric acid[c] + 0.0001 lignocerate[c] + 0.278 linoleate[c] + 0.0023 linolenate[c] + 0.0001 margaric acid[c] + 0.0001 mead acid[c] + 0.003 myristic acid[c] + 0.0001 nervonic acid[c] + 0.0001 nonadecylic acid[c] + 0.4519 oleate[c] + 0.0001 omega-3-arachidonic acid[c] + 0.0849 palmitate[c] + 0.014 palmitolate[c] + 0.0001 pentadecylic acid[c] + 0.0001 physeteric acid[c] + 0.0584 stearate[c] + 0.0026 stearidonic acid[c] + 0.0001 tricosanoic acid[c] + 0.0001 tridecylic acid[c] + 0.0001 ximenic acid[c]
                  'HMR_0688';   % fatty acid-LD-PC pool[c] <=> 0.0005 (10Z)-heptadecenoic acid[c] + 0.0005 (11Z,14Z)-eicosadienoic acid[c] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0005 (13Z)-eicosenoic acid[c] + 0.0005 (13Z)-octadecenoic acid[c] + 0.0005 (13Z,16Z)-docosadienoic acid[c] + 0.0039 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0005 (6Z,9Z)-octadecadienoic acid[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0005 (7Z)-octadecenoic acid[c] + 0.0005 (7Z)-tetradecenoic acid[c] + 0.0005 (9E)-tetradecenoic acid[c] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0005 10,13,16,19-docosatetraenoic acid[c] + 0.0005 10,13,16-docosatriynoic acid[c] + 0.0005 12,15,18,21-tetracosatetraenoic acid[c] + 0.0005 13,16,19-docosatrienoic acid[c] + 0.0005 7-palmitoleic acid[c] + 0.0005 8,11-eicosadienoic acid[c] + 0.0005 9-eicosenoic acid[c] + 0.0005 9-heptadecylenic acid[c] + 0.0252 DHA[c] + 0.0056 DPA[c] + 0.0109 EPA[c] + 0.0011 adrenic acid[c] + 0.0709 arachidonate[c] + 0.0005 behenic acid[c] + 0.0005 cerotic acid[c] + 0.0005 cis-cetoleic acid[c] + 0.0005 cis-erucic acid[c] + 0.0005 cis-gondoic acid[c] + 0.0271 cis-vaccenic acid[c] + 0.0228 dihomo-gamma-linolenate[c] + 0.0005 eicosanoate[c] + 0.0005 elaidate[c] + 0.0036 gamma-linolenate[c] + 0.0005 henicosanoic acid[c] + 0.0005 lauric acid[c] + 0.0005 lignocerate[c] + 0.2479 linoleate[c] + 0.0047 linolenate[c] + 0.0005 margaric acid[c] + 0.0005 mead acid[c] + 0.0054 myristic acid[c] + 0.0005 nervonic acid[c] + 0.0005 nonadecylic acid[c] + 0.1143 oleate[c] + 0.0178 omega-3-arachidonic acid[c] + 0.2781 palmitate[c] + 0.0105 palmitolate[c] + 0.0005 pentadecylic acid[c] + 0.0005 physeteric acid[c] + 0.1271 stearate[c] + 0.0026 stearidonic acid[c] + 0.0005 tricosanoic acid[c] + 0.0005 tridecylic acid[c] + 0.0005 ximenic acid[c]
                  'HMR_0689';   % fatty acid-LD-PE pool[c] => 0.0005 (10Z)-heptadecenoic acid[c] + 0.0005 (11Z,14Z)-eicosadienoic acid[c] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0005 (13Z)-eicosenoic acid[c] + 0.0005 (13Z)-octadecenoic acid[c] + 0.0005 (13Z,16Z)-docosadienoic acid[c] + 0.0052 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0005 (6Z,9Z)-octadecadienoic acid[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0005 (7Z)-octadecenoic acid[c] + 0.0005 (7Z)-tetradecenoic acid[c] + 0.0005 (9E)-tetradecenoic acid[c] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0005 10,13,16,19-docosatetraenoic acid[c] + 0.0005 10,13,16-docosatriynoic acid[c] + 0.0005 12,15,18,21-tetracosatetraenoic acid[c] + 0.0005 13,16,19-docosatrienoic acid[c] + 0.0005 7-palmitoleic acid[c] + 0.0005 8,11-eicosadienoic acid[c] + 0.0005 9-eicosenoic acid[c] + 0.0005 9-heptadecylenic acid[c] + 0.0459 DHA[c] + 0.0067 DPA[c] + 0.0221 EPA[c] + 0.0007 adrenic acid[c] + 0.2125 arachidonate[c] + 0.0005 behenic acid[c] + 0.0005 cerotic acid[c] + 0.0005 cis-cetoleic acid[c] + 0.0005 cis-erucic acid[c] + 0.0005 cis-gondoic acid[c] + 0.0285 cis-vaccenic acid[c] + 0.0411 dihomo-gamma-linolenate[c] + 0.0005 eicosanoate[c] + 0.0005 elaidate[c] + 0.0011 gamma-linolenate[c] + 0.0005 henicosanoic acid[c] + 0.0005 lauric acid[c] + 0.0005 lignocerate[c] + 0.1312 linoleate[c] + 0.0166 linolenate[c] + 0.0005 margaric acid[c] + 0.0005 mead acid[c] + 0.0319 myristic acid[c] + 0.0005 nervonic acid[c] + 0.0005 nonadecylic acid[c] + 0.0619 oleate[c] + 0.0163 omega-3-arachidonic acid[c] + 0.1243 palmitate[c] + 0.0327 palmitolate[c] + 0.0005 pentadecylic acid[c] + 0.0005 physeteric acid[c] + 0.1983 stearate[c] + 0.0025 stearidonic acid[c] + 0.0005 tricosanoic acid[c] + 0.0005 tridecylic acid[c] + 0.0005 ximenic acid[c]
                  'HMR_0690';   % fatty acid-LD-PS pool[c] => 0.0005 (10Z)-heptadecenoic acid[c] + 0.0005 (11Z,14Z)-eicosadienoic acid[c] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0005 (13Z)-eicosenoic acid[c] + 0.0005 (13Z)-octadecenoic acid[c] + 0.0005 (13Z,16Z)-docosadienoic acid[c] + 0.0055 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0005 (6Z,9Z)-octadecadienoic acid[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0005 (7Z)-octadecenoic acid[c] + 0.0005 (7Z)-tetradecenoic acid[c] + 0.0005 (9E)-tetradecenoic acid[c] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0005 10,13,16,19-docosatetraenoic acid[c] + 0.0005 10,13,16-docosatriynoic acid[c] + 0.0005 12,15,18,21-tetracosatetraenoic acid[c] + 0.0005 13,16,19-docosatrienoic acid[c] + 0.0005 7-palmitoleic acid[c] + 0.0005 8,11-eicosadienoic acid[c] + 0.0005 9-eicosenoic acid[c] + 0.0005 9-heptadecylenic acid[c] + 0.0931 DHA[c] + 0.0329 DPA[c] + 0.0286 EPA[c] + 0.0043 adrenic acid[c] + 0.2005 arachidonate[c] + 0.0005 behenic acid[c] + 0.0005 cerotic acid[c] + 0.0005 cis-cetoleic acid[c] + 0.0005 cis-erucic acid[c] + 0.0005 cis-gondoic acid[c] + 0.0131 cis-vaccenic acid[c] + 0.0089 dihomo-gamma-linolenate[c] + 0.0005 eicosanoate[c] + 0.0005 elaidate[c] + 0.0011 gamma-linolenate[c] + 0.0005 henicosanoic acid[c] + 0.0005 lauric acid[c] + 0.0005 lignocerate[c] + 0.0177 linoleate[c] + 0.0035 linolenate[c] + 0.0005 margaric acid[c] + 0.0005 mead acid[c] + 0.0055 myristic acid[c] + 0.0005 nervonic acid[c] + 0.0005 nonadecylic acid[c] + 0.0323 oleate[c] + 0.0174 omega-3-arachidonic acid[c] + 0.0337 palmitate[c] + 0.0058 palmitolate[c] + 0.0005 pentadecylic acid[c] + 0.0005 physeteric acid[c] + 0.4731 stearate[c] + 0.0025 stearidonic acid[c] + 0.0005 tricosanoic acid[c] + 0.0005 tridecylic acid[c] + 0.0005 ximenic acid[c]
                  'HMR_0691';   % fatty acid-LD-PI pool[c] => 0.0005 (10Z)-heptadecenoic acid[c] + 0.0005 (11Z,14Z)-eicosadienoic acid[c] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0005 (13Z)-eicosenoic acid[c] + 0.0005 (13Z)-octadecenoic acid[c] + 0.0005 (13Z,16Z)-docosadienoic acid[c] + 0.0124 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0005 (6Z,9Z)-octadecadienoic acid[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0005 (7Z)-octadecenoic acid[c] + 0.0005 (7Z)-tetradecenoic acid[c] + 0.0005 (9E)-tetradecenoic acid[c] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0005 10,13,16,19-docosatetraenoic acid[c] + 0.0005 10,13,16-docosatriynoic acid[c] + 0.0005 12,15,18,21-tetracosatetraenoic acid[c] + 0.0005 13,16,19-docosatrienoic acid[c] + 0.0005 7-palmitoleic acid[c] + 0.0005 8,11-eicosadienoic acid[c] + 0.0005 9-eicosenoic acid[c] + 0.0005 9-heptadecylenic acid[c] + 0.0269 DHA[c] + 0.0101 DPA[c] + 0.0043 EPA[c] + 0.0031 adrenic acid[c] + 0.2118 arachidonate[c] + 0.0005 behenic acid[c] + 0.0005 cerotic acid[c] + 0.0005 cis-cetoleic acid[c] + 0.0005 cis-erucic acid[c] + 0.0005 cis-gondoic acid[c] + 0.0266 cis-vaccenic acid[c] + 0.0228 dihomo-gamma-linolenate[c] + 0.0005 eicosanoate[c] + 0.0005 elaidate[c] + 0.0012 gamma-linolenate[c] + 0.0005 henicosanoic acid[c] + 0.0005 lauric acid[c] + 0.0005 lignocerate[c] + 0.0678 linoleate[c] + 0.0019 linolenate[c] + 0.0005 margaric acid[c] + 0.0005 mead acid[c] + 0.0056 myristic acid[c] + 0.0005 nervonic acid[c] + 0.0005 nonadecylic acid[c] + 0.136 oleate[c] + 0.0179 omega-3-arachidonic acid[c] + 0.0678 palmitate[c] + 0.0053 palmitolate[c] + 0.0005 pentadecylic acid[c] + 0.0005 physeteric acid[c] + 0.3553 stearate[c] + 0.0027 stearidonic acid[c] + 0.0005 tricosanoic acid[c] + 0.0005 tridecylic acid[c] + 0.0005 ximenic acid[c]
                  'HMR_0692';   % fatty acid-LD-SM pool[c] <=> 0.002 (10Z)-heptadecenoic acid[c] + 0.002 (11Z,14Z)-eicosadienoic acid[c] + 0.002 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.002 (13Z)-eicosenoic acid[c] + 0.002 (13Z)-octadecenoic acid[c] + 0.002 (13Z,16Z)-docosadienoic acid[c] + 0.002 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.002 (6Z,9Z)-octadecadienoic acid[c] + 0.002 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.002 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.002 (7Z)-octadecenoic acid[c] + 0.002 (7Z)-tetradecenoic acid[c] + 0.002 (9E)-tetradecenoic acid[c] + 0.002 (9Z,12Z,15Z,18Z)-TTA[c] + 0.002 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.002 10,13,16,19-docosatetraenoic acid[c] + 0.002 10,13,16-docosatriynoic acid[c] + 0.002 12,15,18,21-tetracosatetraenoic acid[c] + 0.002 13,16,19-docosatrienoic acid[c] + 0.002 7-palmitoleic acid[c] + 0.002 8,11-eicosadienoic acid[c] + 0.002 9-eicosenoic acid[c] + 0.002 9-heptadecylenic acid[c] + 0.003 DHA[c] + 0.002 DPA[c] + 0.0068 EPA[c] + 0.002 adrenic acid[c] + 0.0136 arachidonate[c] + 0.002 behenic acid[c] + 0.002 cerotic acid[c] + 0.002 cis-cetoleic acid[c] + 0.002 cis-erucic acid[c] + 0.002 cis-gondoic acid[c] + 0.0453 cis-vaccenic acid[c] + 0.002 dihomo-gamma-linolenate[c] + 0.002 eicosanoate[c] + 0.002 elaidate[c] + 0.002 gamma-linolenate[c] + 0.002 henicosanoic acid[c] + 0.002 lauric acid[c] + 0.002 lignocerate[c] + 0.025 linoleate[c] + 0.0075 linolenate[c] + 0.002 margaric acid[c] + 0.002 mead acid[c] + 0.011 myristic acid[c] + 0.002 nervonic acid[c] + 0.002 nonadecylic acid[c] + 0.061 oleate[c] + 0.002 omega-3-arachidonic acid[c] + 0.557 palmitate[c] + 0.0358 palmitolate[c] + 0.002 pentadecylic acid[c] + 0.002 physeteric acid[c] + 0.138 stearate[c] + 0.002 stearidonic acid[c] + 0.002 tricosanoic acid[c] + 0.002 tridecylic acid[c] + 0.002 ximenic acid[c]
                  'HMR_3537'};  % cholesterol-ester pool[l] => 0.0001 cholesterol-ester-10,13,16,19-docosa[l] + 0.0001 cholesterol-ester-10,13,16-docosa[l] + 0.0001 cholesterol-ester-10-hepta[l] + 0.0001 cholesterol-ester-11,14,17-eico[l] + 0.0001 cholesterol-ester-11,14-eicosa[l] + 0.0001 cholesterol-ester-11-docose[l] + 0.0001 cholesterol-ester-11-eico[l] + 0.0001 cholesterol-ester-12,15,18,21-tetracosa[l] + 0.0001 cholesterol-ester-13,16,19-doco[l] + 0.0001 cholesterol-ester-13,16-docosa[l] + 0.0001 cholesterol-ester-13-docose[l] + 0.0001 cholesterol-ester-13-eicose[l] + 0.0001 cholesterol-ester-13-octade[l] + 0.0001 cholesterol-ester-15-tetra[l] + 0.0041 cholesterol-ester-4,7,10,13,16,19-doco[l] + 0.0001 cholesterol-ester-4,7,10,13,16-docosa[l] + 0.0071 cholesterol-ester-5,8,11,14,17-eico[l] + 0.0001 cholesterol-ester-5,8,11-eico[l] + 0.0001 cholesterol-ester-5-tetra[l] + 0.0001 cholesterol-ester-6,9,12,15,18,21-tetra[l] + 0.0001 cholesterol-ester-6,9,12,15,18-tetraco[l] + 0.0001 cholesterol-ester-6,9,12,15-octa[l] + 0.0001 cholesterol-ester-6,9-octa[l] + 0.0001 cholesterol-ester-7,10,13,16,19-docosa[l] + 0.0001 cholesterol-ester-7,10,13,16-docosa[l] + 0.0406 cholesterol-ester-7-hexa[l] + 0.0001 cholesterol-ester-7-octade[l] + 0.0001 cholesterol-ester-7-tetrade[l] + 0.0001 cholesterol-ester-8,11,14,17-eico[l] + 0.0001 cholesterol-ester-8,11-eico[l] + 0.0001 cholesterol-ester-9,12,15,18,21-tetra[l] + 0.0001 cholesterol-ester-9,12,15,18-tetraco[l] + 0.0001 cholesterol-ester-9-eicose[l] + 0.0001 cholesterol-ester-9-heptade[l] + 0.0001 cholesterol-ester-9-octa[l] + 0.0001 cholesterol-ester-9-tetrade[l] + 0.0518 cholesterol-ester-arach[l] + 0.0001 cholesterol-ester-cis-vac[l] + 0.005 cholesterol-ester-dihomo-gamma[l] + 0.0001 cholesterol-ester-docosa[l] + 0.0001 cholesterol-ester-eico[l] + 0.0001 cholesterol-ester-gamma-lin[l] + 0.0001 cholesterol-ester-heneico[l] + 0.0001 cholesterol-ester-hepta[l] + 0.0001 cholesterol-ester-hexacosa[l] + 0.0001 cholesterol-ester-hexecose[l] + 0.0001 cholesterol-ester-laur[l] + 0.5254 cholesterol-ester-lin[l] + 0.0061 cholesterol-ester-linolen[l] + 0.0081 cholesterol-ester-myrist[l] + 0.0001 cholesterol-ester-nanode[l] + 0.1958 cholesterol-ester-ol[l] + 0.138 cholesterol-ester-palm[l] + 0.0001 cholesterol-ester-palmn[l] + 0.0001 cholesterol-ester-penta[l] + 0.0132 cholesterol-ester-stea[l] + 0.0001 cholesterol-ester-tetraco[l] + 0.0001 cholesterol-ester-trico[l] + 0.0001 cholesterol-ester-tridec[l]
model.lb(ismember(model.rxns,constrain_rxns)) = 0;
model.ub(ismember(model.rxns,constrain_rxns)) = 0;



%% clear intermediate variables

clearvars -except ihuman model del_rxns
