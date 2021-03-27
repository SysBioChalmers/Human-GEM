%
% FILE NAME:    repairModelLeaks.m
% 
% PURPOSE: Script to identify and update/constrain/remove reactions that 
%          allow the creation of mass and/or energy (which results in a
%          "leaky" model). The script is divided into two main sections:
%
%          1. Changes to reaction bounds/direction/reversibility
%             - Only the bounds (lb or ub) or reaction direction is
%               modified in the reaction.
%
%          2. Changes to reaction stoichiometry/metabolites
%             - The metabolite(s) involved in a reaction are
%               changed/added/removed, and/or the stoichiometric
%               coefficients in a reaction are modified.
%
%          3. Reaction inactivation
%             - A total of 236 reactions (12 were also modified in above two
%               steps) need to be constrained (i.e. ub=lb=0), in order to
%               achieve all the metabolic tasks listed in
%               `metabolicTasks_LeakCheck.xls`. These inactivate reactions
%               are archieved in `inactivationRxns.tsv` and should be further
%               assessed to determine which, if any, should be completely
%               removed from the model.
%
% Note: This script is partitioned into sections, each includes changes aim
% to achieve certain task(s), whose ID numbers are marked and correspond to
% the task ID numbers listed in the "metabolicTasks_LeakCheck.xls" task list.
% Since these changes indirectly affect each other and cannot be fully separated.
% The ID number(s) listed in each section indicate the primary task(s) that
% the changes are designed to address but may not fully fix them. However, 
% all targeted tasks can be achived by collective implemention of the changes
% in this script.


%% Load model and initialize some variables

% load HumanGEM model (if not already loaded)
if ~exist('ihuman','var')
    load('humanGEM.mat');  % version 0.6.0
end
ihuman_orig = ihuman;  % to keep track of changes made

% initialize vars
rxnNotes = {};


%% Changes to reaction bounds/direction/reversibility

% The following reaction from Recon3D:
%
%   NMNATr: ATP[c] + H+[c] + nicotinamide D-ribonucleotide[c] <=> NAD+[c] + PPi[c]
%
% should not be reversible (HumanCyc, rxn 2.7.7.1).
% Relevant metabolic task IDs: 2-9, 17
rxn_ind = ismember(ihuman.rxns,'NMNATr');
ihuman.lb(rxn_ind) = 0;
rxnNotes = [rxnNotes; {'NMNATr','rxn should not be reversible (HumanCyc rxn 2.7.7.1)'}];


% The following reaction pairs are present in the model:
%
%   HMR_1358: 9-peroxy-(5Z,7E,11Z,14Z)-eicosatetraenoate[c] + H+[c] <=> arachidonate[c] + O2-[c]
%    RE3449C: 3 arachidonate[c] + 5 H+[c] + 5 O2-[c] <=> 3 9-peroxy-(5Z,7E,11Z,14Z)-eicosatetraenoate[c] + 4 H2O[c]
%
%   HMR_1364: 12-peroxy-(5Z,8Z,10E,14Z)-eicosatetraenoate[c] + H+[c] <=> arachidonate[c] + O2-[c]
%    RE3452C: 3 arachidonate[c] + 5 H+[c] + 5 O2-[c] <=> 3 12-peroxy-(5Z,8Z,10E,14Z)-eicosatetraenoate[c] + 4 H2O[c]
%
%   HMR_1361: 11-peroxy-(5Z,8Z,12E,14Z)-eicosatetraenoate[c] + H+[c] <=> arachidonate[c] + O2-[c]
%    RE3458C: 3 arachidonate[c] + 5 H+[c] + 5 O2-[c] <=> 3 11-peroxy-(5Z,8Z,12E,14Z)-eicosatetraenoate[c] + 4 H2O[c]
%
%   HMR_1405: 10-peroxy-docosahexaenoate[c] + H+[c] <=> DHA[c] + O2-[c]
%    RE2852C: 3 DHA[c] + 5 H+[c] + 5 O2-[c] <=> 3 10-peroxy-docosahexaenoate[c] + 4 H2O[c]
%
% The problem with these reaction pairs is that, together, they can
% generate superoxide (and eventually oxygen) and protons from water, which
% should not be possible. The Recon3D versions appear to be an attempt to
% charge-balance the reactions. To prevent the model from splitting water
% into protons and oxygen, these reactions should be constrained such that
% they are irreversible, and in the direction of peroxidation (consumption
% of O2-).
% Relevant metabolic task IDs: 1
rxn_ind = ismember(ihuman.rxns,{'HMR_1358';'RE3449C';'HMR_1361';'RE3458C';'HMR_1405';'RE2852C';'HMR_1364';'RE3452C'});
o2s_ind = getIndexes(ihuman,'O2-[c]','metscomps');

% turn reactions around so the forward direction consumes O2-
ihuman.S(:,rxn_ind) = -(ihuman.S(:,rxn_ind) .* sign(ihuman.S(o2s_ind,rxn_ind)));
ihuman.lb(rxn_ind) = 0;  % constrain lower bound to zero
rxnNotes = [rxnNotes; [ihuman.rxns(rxn_ind), repmat({'this rxn should not be able to produce superoxide'},sum(rxn_ind),1)]];


% The following pairs of reactions:
%
%   HMR_1357: 5,9-cyclo-6,8,12-trihydroxy-(10E,14Z)-eicosadienoic acid[c] + dehydroascorbic acid[c] + H2O[c] <=> 5,9-cyclo-6,8-cycloperoxy-12-hydroperoxy-(10E,14Z)-eicosadienoate[c] + ascorbate[c] + 2 H+[c]
%    RE3456C: 5,9-cyclo-6,8-cycloperoxy-12-hydroperoxy-(10E,14Z)-eicosadienoate[c] + 2 ascorbate[c] + 2 H+[c] => 5,9-cyclo-6,8,12-trihydroxy-(10E,14Z)-eicosadienoic acid[c] + 2 dehydroascorbic acid[c] + H2O[c]
%
%   HMR_1360: 5,9,11-trihydroxyprosta-(6E,14Z)-dien-1-oate[c] + dehydroascorbic acid[c] + H2O[c] <=> 9,11-cycloperoxy-5-hydroperoxy-(6E,14Z)-eicosadienoate[c] + ascorbate[c] + 2 H+[c]
%    RE3450C: 9,11-cycloperoxy-5-hydroperoxy-(6E,14Z)-eicosadienoate[c] + 2 ascorbate[c] + 2 H+[c] => 5,9,11-trihydroxyprosta-(6E,14Z)-dien-1-oate[c] + 2 dehydroascorbic acid[c] + H2O[c]
%
% have different stoichiometries for ascorbate, and therefore when run
% together can oxidize ascorbate to dehydroascorbic acid without any other
% metabolites. These reactions should not be reversible (i.e., water cannot
% be used to peroxidate those compounds). The HMR versions will therefore
% be made irreversible such that they produce water, as is the case with
% the Recon3D versions.
% Relevant metabolic task IDs: 1
rxn_ind = ismember(ihuman.rxns,{'HMR_1357';'HMR_1360'});
h2o_ind = getIndexes(ihuman,'H2O[c]','metscomps');

% turn reactions around so the forward direction produces H2O
ihuman.S(:,rxn_ind) = ihuman.S(:,rxn_ind) .* sign(ihuman.S(h2o_ind,rxn_ind));
ihuman.lb(rxn_ind) = 0;  % constrain lower bound to zero
rxnNotes = [rxnNotes; [ihuman.rxns(rxn_ind), repmat({'H2O should not be able to peroxidate compound'},sum(rxn_ind),1)]];


% These rxns should not be allowed to operate in the reverse direction 
% (i.e., they should require ATP consumption if run in reverse):
%
%  RE3238C: (11Z)-eicosenoyl-CoA[c] + H2O[c] <=> cis-gondoic acid[c] + CoA[c] + H+[c]
%  RE3239C: (13Z)-docosenoyl-CoA[c] + H2O[c] <=> cis-erucic acid[c] + CoA[c] + H+[c]
%  RE2649C: H2O[c] + propanoyl-CoA[c] <=> CoA[c] + H+[c] + propanoate[c]
%
% The lower bound of these reactions will therefore be constrained to zero.
% Relevant metabolic task IDs: 2-9, 17
rxn_ind = ismember(ihuman.rxns,{'RE3238C';'RE3239C';'RE2649C'});
ihuman.lb(rxn_ind) = 0;
rxnNotes = [rxnNotes; [ihuman.rxns(rxn_ind), repmat({'reverse rxn requires ATP, therefore rxn was made irreversible'},sum(rxn_ind),1)]];


% These reactions concern the glycerol phosphate shuttle:
%
%  HMR_0483: DHAP[c] + ubiquinol[m] => sn-glycerol-3-phosphate[c] + ubiquinone[m]
%  HMR_0482: DHAP[c] + FADH2[c] => FAD[c] + sn-glycerol-3-phosphate[c]
%    r0202m: NAD+[m] + sn-glycerol-3-phosphate[m] => DHAP[m] + H+[m] + NADH[m]
%
% The HMR reactions are written in the wrong direction, and should be
% reversed. The Recon3D reaction (r0202m) is also in the wrong direction,
% but should also not take place in the mitochondria, and should therefore
% be deleted.
% Relevant metabolic task IDs: 10-12
[~,rxn_ind] = ismember({'HMR_0483';'HMR_0482'},ihuman.rxns);
DHAP_ind = getIndexes(ihuman,'DHAP[c]','metscomps');
rxn_ind(ihuman.S(DHAP_ind,rxn_ind) > 0) = [];  % exclude rxns that have already been fixed
if ~isempty(rxn_ind)
    ihuman.S(:,rxn_ind) = -ihuman.S(:,rxn_ind);  % turn rxns around
    rxnNotes = [rxnNotes; [ihuman.rxns(rxn_ind), repmat({'directionality changed to reflect proper activity of glycerol phosphate shuttle'},length(rxn_ind),1)]];
end


% The following reactions are part of the melatonin degradation pathway:
%
%   HMR_4551: formyl-N-acetyl-5-methoxykynurenamine[c] + 2 H+[c] + H2O2[c] <=> formate[c] + H2O[c] + N-acetyl-5-methoxykynuramine[c]'
%    RE2440C: 2 formyl-N-acetyl-5-methoxykynurenamine[c] + H2O2[c] <=> CO2[c] + formate[c] + H+[c] + 2 N-acetyl-5-methoxykynuramine[c]
%
% Although they are written as reversible, the reverse reaction should not
% be possible (see, e.g., PMID: 19573038). Therefore, these reactions
% should be constrained to proceed only in the forward direction.
% Relevant metabolic task IDs: 1
[~,rxn_ind] = ismember({'HMR_4551';'RE2440C'},ihuman.rxns);
ihuman.lb(rxn_ind) = 0;
rxnNotes = [rxnNotes; [ihuman.rxns(rxn_ind), repmat({'the reverse reaction (formate consumption) should not be possible (PMID: 19573038)'},length(rxn_ind),1)]];


% The following reaction from Recon3D:
%
%   r1466: gamma-linolenoyl-CoA[r] + 4 H+[r] + 2 malonyl-CoA[r] + 2 NADPH[r] + 2 O2[r] <=> 4 CO2[r] + 2 CoA[r] + dihomo-gamma-linolenoyl-CoA[r] + 2 H2O[r] + 2 NADP+[r]
%
% should NOT be reversible. Only the forward direction is possible.
% Therefore, the lower bound of this reaction will be set to zero.
% Relevant metabolic task IDs: 10-12
ihuman.lb(ismember(ihuman.rxns,{'r1466'})) = 0;
rxnNotes = [rxnNotes; {'r1466', 'reaction should not be reversible'}];


% The following Recon3D reactions:
%
%   RE1448N: 1-phosphatidyl-1D-myo-inositol-5-phosphate[n] + H2O[n] <=> PI pool[n] + Pi[n]
%   RE3273C: H2O[c] + PI pool[c] <=> H+[c] + inositol[c] + phosphatidate-LD-TAG pool[c]
%   RE3273G: H2O[g] + PI pool[g] <=> H+[g] + inositol[g] + phosphatidate-LD-TAG pool[g]
%   RE3273R: H2O[r] + PI pool[r] <=> H+[r] + inositol[r] + phosphatidate-LD-TAG pool[r]
%
% should not be reversible. They should only proceed in the forward direction.
% Relevant metabolic task IDs: 2-9, 17
rxn_ind = ismember(ihuman.rxns,{'RE1448N';'RE3273C';'RE3273G';'RE3273R'});
ihuman.lb(rxn_ind) = 0;
ihuman.rev(rxn_ind) = 0;
rxnNotes = [rxnNotes; [ihuman.rxns(rxn_ind), repmat({'reaction should not be reversible'},sum(rxn_ind),1)]];


%% Changes to reaction stoichiometry/metabolites

% The following reaction from HMR:
%
%   HMR_7625: [protein]-L-arginine[c] + NAD+[c] => N(omega)-(ADP-D-ribosyl)-L-arginine[c] + nicotinamide[c]
%
% involves a protein-bound amino acid being converted to an unbound
% metabolite, which should require an additional water as a reactant.
% Relevant metabolic task IDs: 1
rxn_ind = ismember(ihuman.rxns,'HMR_7625');
h2o_ind = getIndexes(ihuman,'H2O[c]','metscomps');
ihuman.S(h2o_ind,rxn_ind) = -1;
rxnNotes = [rxnNotes; {'HMR_7625', 'water is required as a reactant'}];


% The following HMR reaction involves a proton from a compartment other
% than the one containing the other metabolites in the reaction:
%
%   HMR_8616: ATP[m] + cob(I)alamin[m] + H+[c] <=> cobamide-coenzyme[m] + triphosphate[m]
%
% However, there is no evidence in the literature for such transport, and
% it is therefore suspected to be a compartment labeling error. The proton
% compartment should therefore be corrected to mitochondria [m] to be
% consistent with the other compounds in the reaction.
% Relevant metabolic task IDs: 11
Hc = getIndexes(ihuman,'H+[c]','metscomps');
Hm = getIndexes(ihuman,'H+[m]','metscomps');
[~,rxn_ind] = ismember('HMR_8616',ihuman.rxns);
ihuman.S(Hc,rxn_ind) = 0;
ihuman.S(Hm,rxn_ind) = -1;
rxnNotes = [rxnNotes; {'HMR_8616', 'corrected suspected mistake in proton compartment'}];


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
% However, this reaction should use NADPH to proceed, and there even exists
% a reaction (from Recon3D) in the model that does just this:
%
%     r0309: NADP+[m] + palmitoyl-CoA[m] <=> (2E)-hexadecenoyl-CoA[m] + H+[m] + NADPH[m]
%
% Therefore, the cofactor in the first rxn above (ARTFR61) should be
% changed to NADPH instead of FADH2. The reaction should probably be
% deleted entirely due to it being redundant, but this fix is sufficient
% for now. The same situation was found for the following reactions:
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
%
% Relevant metabolic task IDs: 2-15, 17
fadh2_ind = getIndexes(ihuman,'FADH2[m]','metscomps');
fad_ind = getIndexes(ihuman,'FAD[m]','metscomps');
nadp_ind = getIndexes(ihuman,'NADP+[m]','metscomps');
nadph_ind = getIndexes(ihuman,'NADPH[m]','metscomps');
h_ind = getIndexes(ihuman,'H+[m]','metscomps');

[~,rxn_ind] = ismember({'ARTFR61';'ARTFR46';'ARTFR32';'ARTFR41';'ARTFR42';'ARTFR12';'ARTFR33';'ARTFR34';'ARTFR43';'ARTFR44';'ARTFR45'},ihuman.rxns);
rxn_ind(ihuman.S(fadh2_ind,rxn_ind) == 0) = [];  % don't change rxns that have been fixed already
for i = 1:length(rxn_ind)
    ihuman.S([fadh2_ind; fad_ind], rxn_ind(i)) = 0;
    ihuman.S([nadp_ind; nadph_ind; h_ind], rxn_ind(i)) = [1; -1; -1];
end
rxnNotes = [rxnNotes; [ihuman.rxns(rxn_ind), repmat({'cofactor changed from FADH2 to NADPH'},length(rxn_ind),1)]];


% The following reaction from Recon3D:
%
%   ALPA_HSx: acylglycerone-phosphate[p] + H+[p] + NADP+[p] => NADPH[p] + Lysophosphatidic Acid[p]
%
% should consume NADPH, not produce it (see KEGG rxn R02756, where
% "Lysophosphatidic Acid" is also known as 1-Acyl-sn-glycerol 3-phosphate.
% Relevant metabolic task IDs: 10-12
nadp_ind = getIndexes(ihuman,'NADP+[p]','metscomps');
nadph_ind = getIndexes(ihuman,'NADPH[p]','metscomps');
[~,rxn_ind] = ismember('ALPA_HSx',ihuman.rxns);
if ihuman.S(nadp_ind,rxn_ind) == -1  % check that the rxn hasn't been changed already
    ihuman.S([nadp_ind, nadph_ind], rxn_ind) = [1,-1];
    rxnNotes = [rxnNotes; {'ALPA_HSx', 'reaction should consume NADPH, not produce it (KEGG R02756)'}];
end


% The following reaction from Recon3D:
%
%   TAG_HSad_E: (4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoyl-CoA[c] + (5Z,8Z,11Z,14Z,17Z)-eicosapentaenoyl-CoA[c] + arachidonyl-CoA[c] + 2 H2O[c] + linolenoyl-CoA[c] + linoleoyl-CoA[c] + 2 sn-glycerol-3-phosphate[c] => 6 CoA[c] + 2 Pi[c] + 2 TAG-VLDL pool[c]
%
% is has one extra CoA in the products compared to reactants. The
% stoichiometric coefficient of the product CoA should therefore be changed
% to 5 to correct this imbalance.
% Relevant metabolic task IDs: 16
met_ind = getIndexes(ihuman,'CoA[c]','metscomps');
rxn_ind = getIndexes(ihuman,'TAG_HSad_E','rxns');
if ihuman.S(met_ind,rxn_ind) == 6
    ihuman.S(met_ind,rxn_ind) = 5;
    rxnNotes = [rxnNotes; {'TAG_HSad_E', 'adjusted the stoichiometric coefficient of CoA in the products to balance the reaction'}];
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
%
% Relevant metabolic task IDs: 16
met_ind = getIndexes(ihuman,{'cholesterol-ester pool[l]';'cholesterol[l]'},'metscomps');
rxn_ind = getIndexes(ihuman,{'HMR_5238';'HMR_5233'},'rxns');
ihuman.S(met_ind,rxn_ind) = 0;
rxnNotes = [rxnNotes; [ihuman.rxns(rxn_ind), repmat({'balanced mass by removing cholesterol and cholesterol-ester pool from products'},length(rxn_ind),1)]];


% The following reaction from Recon3D:
%
%   DOLGPP_Ler: 0.1 dolichyl-D-glucosyl-phosphate[r] + H2O[r] => 0.1 dolichyl-phosphate[r] + glucose[r] + H+[r]
%
% Is not properly formulated, as it is creating 1 equivalent of glucose
% from 0.1 equivalents. It is nearly the same as the following reaction:
%
%   HMR_8692: dolichyl-D-glucosyl-phosphate[r] + H2O[r] => dolichyl-phosphate[r] + glucose[r]
%
% except for the coefficients and the additional proton. This appears to be
% a problem with the difference in formula for dolichyl-phosphate between
% Recon3D and HMR:
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
% Therefore, these Recon3D reactions should be removed from the model and
% dealt with in script constrainReactions.m.
%
% In addition, there were a few other reactions that did not have an HMR
% equivalent, but need to have their stoich coeffs adjusted from 0.1 to 1, to
% be consistent with how other reactions in the model are treating the mass
% of these dolichol compounds.
%
%      H8MTer_L: 0.1 dolichyl-phosphate-D-mannose[r] + gpi heparan sulfate[r] => 0.1 dolichyl-phosphate[r] + H+[r] + Mannosyl-3-(Phosphoethanolaminyl-Mannosyl)-Glucosaminyl-Acylphosphatidylinositol (M4A)[r]
%      H8MTer_U: 0.1 dolichyl-phosphate-D-mannose[r] + gpi heparan sulfate[r] => 0.1 dolichyl-phosphate[r] + H+[r] + HMA[r]
%    UDPDOLPT_L: 0.1 dolichyl-phosphate[c] + UDP-glucose[c] => UDP[c] + 0.1 dolichyl-D-glucosyl-phosphate[c]
%
% Relevant metabolic task IDs: 16
rxn_ind = getIndexes(ihuman,{'H8MTer_L';'H8MTer_U';'UDPDOLPT_L'},'rxns');
ihuman.S(:,rxn_ind) = sign(ihuman.S(:,rxn_ind));  % convert all nonzero values to +/- 1
rxnNotes = [rxnNotes; [ihuman.rxns(rxn_ind), repmat({'corrected dolichol-related mets stoich coeffs to be consistent with its effective mass elsewhere in the model'},length(rxn_ind),1)]];


%% Apply reaction constraining and update rxnNotes

% Load the inactivation reaction list and constrain them
fid = fopen('inactivationRxns.tsv','r');
input = textscan(fid,'%s %s','Delimiter','\t','Headerlines',1);
fclose(fid);
constrainRxnNotes = [input{1}(1:end-2), input{2}(1:end-2)];
ihuman = setParam(ihuman, 'eq', constrainRxnNotes(:,1), 0);

% Update rxnNotes because there are overlap between reaction sets with
% bound/coefficient adjustment and constraining
[a, b] = ismember(constrainRxnNotes(:,1), rxnNotes(:,1));
overlapInd = b(a);
rxnNotes(overlapInd,2) = strcat(rxnNotes(overlapInd,2), '; ', constrainRxnNotes(a,2));
% Add non-overlap reaction sets to notes
[~, nonOverlapInd]= setdiff(constrainRxnNotes(:,1), rxnNotes(:,1));
rxnNotes = [rxnNotes; constrainRxnNotes(nonOverlapInd,:)];


%% Generate model change report
rxnChanges = docRxnChanges(ihuman_orig,ihuman,rxnNotes);
writeRxnChanges(rxnChanges,'repairModelLeaks_rxnChanges',true);


%% Clear intermediate vars and save model file
clearvars -except ihuman

save('../../model/Human-GEM.mat','ihuman');
movefile('repairModelLeaks_rxnChanges.tsv','../../ComplementaryData/modelCuration/');

