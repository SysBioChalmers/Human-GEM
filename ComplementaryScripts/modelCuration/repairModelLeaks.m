%
% FILE NAME:    repairModelLeaks.m
% 
% DATE CREATED: 2018-10-01
%     MODIFIED: 2018-10-03
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script to identify and update/constrain/remove reactions that 
%          allow the creation of mass or energy (which results in a "leaky"
%          model).
%

% initialize del_rxns list
del_rxns = {};

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



% The following reaction from Recon3D:
%
%   r1453: proline[m] + ubiquinol[m] <=> 1-pyrroline-5-carboxylate[m] + 5 H+[m] + ubiquinone[m]
%
% is incorrect. The ubiquinol and ubiquinone should be on opposite sides of
% the reaction, as it is in the HMR version:
%
%   HMR_3838: 1-pyrroline-5-carboxylate[m] + H+[m] + ubiquinol[m] <=> proline[m] + ubiquinone[m]
%
% Therefore the Recon3D reaction should be removed from the model.
del_rxns = [del_rxns; {'r1453'}];



% These rxns should not be allowed to operate in the reverse direction 
% (i.e., they should require ATP consumption if run in reverse):
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


% The following reactions from Recon3D deal with the generation and
% breakdown of various types of lipoproteins:
% 
%     VLDL_HSDEG: 5 H2O[s] + Very Low Density Lipoprotein[s] => 2 cholesterol[s] + 5 glycerol[s] + 2 PC-LD pool[s] + 5 R Total[s] + 5 R Total 2 Position[s] + 5 R Total 3 Position[s] + 0.2 apoB100[s] + 0.2 apoC1[s] + 0.2 apoC2[s] + 0.2 apoC3[s]
%      IDL_HSDEG: 4 H2O[s] + Intermediate Density Lipoprotein[s] => 4 cholesterol[s] + 4 glycerol[s] + 4 R Total[s] + 4 R Total 2 Position[s] + 4 R Total 3 Position[s] + 0.5 apoB100[s] + 0.5 apoE[s]
%      LDL_HSDEG: H2O[s] + Low Density Lipoprotein[s] => 5 cholesterol[s] + glycerol[s] + 2 PC-LD pool[s] + R Total[s] + R Total 2 Position[s] + R Total 3 Position[s] + 2 apoB100[s]
%      HDL_HSDEG: H2O[s] + High Density Lipoprotein[s] => apoA1[s] + 2 cholesterol[s] + glycerol[s] + 2 PC-LD pool[s] + R Total[s] + R Total 2 Position[s] + R Total 3 Position[s] + apoC1[s] + apoC2[s] + apoC3[s] + apoE[s]
%    CHYLO_HSDEG: H2O[s] + Chylomicron Lipoprotein[s] => apoA1[s] + glycerol[s] + R Total[s] + R Total 2 Position[s] + R Total 3 Position[s] + apoB100[s] + apoC1[s] + apoC2[s] + apoC3[s] + apoE[s]
% 
%     VLDL_HSSYN: 0.2 apoB100[c] + 0.2 apoC1[c] + 0.2 apoC2[c] + 0.2 apoC3[c] + 2 cholesterol[c] + 2 PC-LD pool[c] + 5 TAG-VLDL pool[c] => Very Low Density Lipoprotein[c]
%      IDL_HSSYN: 4 cholesterol[s] + 4 TAG-VLDL pool[s] + 0.5 apoB100[s] + 0.5 apoE[s] => Intermediate Density Lipoprotein[s]
%      LDL_HSSYN: 5 cholesterol[s] + 2 PC-LD pool[s] + TAG-VLDL pool[s] + 2 apoB100[s] => Low Density Lipoprotein[s]
%      HDL_HSSYN: apoA1[s] + 2 cholesterol[s] + 2 PC-LD pool[s] + TAG-VLDL pool[s] + apoC1[s] + apoC2[s] + apoC3[s] + apoE[s] => High Density Lipoprotein[s]
%   MYELIN_HSSYN: cholesterol[c] + PC-LD pool[c] + PE-LD pool[c] + PI pool[c] + PS-LD pool[c] + SM pool[c] + sulfatide galactocerebroside[c] => Myelin Sheath
%    CHYLO_HSSYN: apoA1[c] + apoB100[c] + apoC1[c] + apoC2[c] + apoC3[c] + apoE[c] + TAG-VLDL pool[c] => Chylomicron Lipoprotein[c]
%
% The problem with these reactions is that the synthesis reactions involve
% the "TAG-VLDL pool" metabolite, whereas their degradation does not
% include this metabolite, but instead forms other components. This causes
% a problem because the makeup of this metabolite seems to be different
% between HMR and Recon3D, so using these reactions leads to
% inconsistencies in its mass. Therefore, these reactions should be
% constrained until they can be properly re-balanced/verified.
del_rxns = [del_rxns; {'VLDL_HSDEG';'IDL_HSDEG';'LDL_HSDEG';'HDL_HSDEG';'CHYLO_HSDEG';
                       'VLDL_HSSYN';'IDL_HSSYN';'LDL_HSSYN';'HDL_HSSYN';'MYELIN_HSSYN';'CHYLO_HSSYN'}];

%


% The following two reactions from HMR:
%   HMR_2128: calcitriol[m] + H+[m] + NADPH[m] + O2[m] => calcitroic acid[m] + H2O[m] + NADP+[m]
%   HMR_2141: calcitriol[m] + H+[m] + NADPH[m] + O2[m] => calcitetrol[m] + H2O[m] + NADP+[m]
% are identical except for one of the products. An investigation of these
% metabolites revealed that the first reaction is imbalanced, as calcitroic
% acid has 4 fewer carbons than calcitriol (or calcitetrol). Therefore, the
% first reaction (HMR_2128) should be removed from the model.
del_rxns = [del_rxns; {'HMR_2128'}];


% The following reactions from Recon3D:
%    DHAPA: DHAP[c] + R Total Coenzyme A[c] => acylglycerone-phosphate[c] + CoA[c]
%   DHAPAx: DHAP[p] + R Total Coenzyme A[p] => acylglycerone-phosphate[p] + CoA[p]
% are imbalanced, because Recon3D treats "R Total Coenzyme A" as equivalent
% to palmitoyl-CoA. This leads to reactions that involve 3 + 37 = 40
% carbons consumed to produce 4 + 21 = 25 carbons. They should therefore be
% removed from the model.
del_rxns = [del_rxns; {'DHAPA';'DHAPAx'}];



% The following Recon3D reaction:
%
%   DSAT: sphinganine[c] + R Total Coenzyme A[c] => CoA[c] + dihydroceramide pool[c] + H+[c]
%
% produces the "dihydroceramide pool" metabolite, but its mass is treated
% differently here than in HMR reactions, e.g.:
%
%  HMR_0753:  dihydroceramide pool[c] + H+[c] + H2O[c] <=> fatty acid-LD-SM pool[c] + sphinganine[c]
%  HMR_0692:  fatty acid-LD-SM pool[c] <=> 0.002 (10Z)-heptadecenoic acid[c] + 0.002 (11Z,14Z)-eicosadienoic acid[c] + 0.002 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.002 (13Z)-eicosenoic acid[c] + 0.002 (13Z)-octadecenoic acid[c] + 0.002 (13Z,16Z)-docosadienoic acid[c] + 0.002 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.002 (6Z,9Z)-octadecadienoic acid[c] + 0.002 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.002 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.002 (7Z)-octadecenoic acid[c] + 0.002 (7Z)-tetradecenoic acid[c] + 0.002 (9E)-tetradecenoic acid[c] + 0.002 (9Z,12Z,15Z,18Z)-TTA[c] + 0.002 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.002 10,13,16,19-docosatetraenoic acid[c] + 0.002 10,13,16-docosatriynoic acid[c] + 0.002 12,15,18,21-tetracosatetraenoic acid[c] + 0.002 13,16,19-docosatrienoic acid[c] + 0.002 7-palmitoleic acid[c] + 0.002 8,11-eicosadienoic acid[c] + 0.002 9-eicosenoic acid[c] + 0.002 9-heptadecylenic acid[c] + 0.003 DHA[c] + 0.002 DPA[c] + 0.0068 EPA[c] + 0.002 adrenic acid[c] + 0.0136 arachidonate[c] + 0.002 behenic acid[c] + 0.002 cerotic acid[c] + 0.002 cis-cetoleic acid[c] + 0.002 cis-erucic acid[c] + 0.002 cis-gondoic acid[c] + 0.0453 cis-vaccenic acid[c] + 0.002 dihomo-gamma-linolenate[c] + 0.002 eicosanoate[c] + 0.002 elaidate[c] + 0.002 gamma-linolenate[c] + 0.002 henicosanoic acid[c] + 0.002 lauric acid[c] + 0.002 lignocerate[c] + 0.025 linoleate[c] + 0.0075 linolenate[c] + 0.002 margaric acid[c] + 0.002 mead acid[c] + 0.011 myristic acid[c] + 0.002 nervonic acid[c] + 0.002 nonadecylic acid[c] + 0.061 oleate[c] + 0.002 omega-3-arachidonic acid[c] + 0.557 palmitate[c] + 0.0358 palmitolate[c] + 0.002 pentadecylic acid[c] + 0.002 physeteric acid[c] + 0.138 stearate[c] + 0.002 stearidonic acid[c] + 0.002 tricosanoic acid[c] + 0.002 tridecylic acid[c] + 0.002 ximenic acid[c]
%
% Therefore, the Recon3D reaction should be deleted/constrained until this
% mass imbalance is resolved.
del_rxns = [del_rxns; {'DSAT'}];


% The following Recon3D reactions:
%
% RE3301C: H2O[c] + PS-LD pool[c] <=> H+[c] + phosphatidate-LD-TAG pool[c] + serine[c]
% RE3301G: H2O[g] + PS-LD pool[g] <=> H+[g] + phosphatidate-LD-TAG pool[g] + serine[g]
% RE3301R: H2O[r] + PS-LD pool[r] <=> H+[r] + phosphatidate-LD-TAG pool[r] + serine[r]
%
% differ from the HMR reaction:
%
% HMR_0660: H2O[c] + PS-LD pool[c] => phosphatidate-LD-PS pool[c] + serine[c]
%
% suggesting a different treatment of these pools from the different
% models. Therefore, the Recon3D reactions should be removed/constrained
% until they can be properly re-balanced/integrated.
del_rxns = [del_rxns; {'RE3301C';'RE3301G';'RE3301R'}];


% The following Recon3D reactions:
%   RE1448N: 1-phosphatidyl-1D-myo-inositol-5-phosphate[n] + H2O[n] <=> PI pool[n] + Pi[n]
%   RE3273C: H2O[c] + PI pool[c] <=> H+[c] + inositol[c] + phosphatidate-LD-TAG pool[c]
%   RE3273G: H2O[g] + PI pool[g] <=> H+[g] + inositol[g] + phosphatidate-LD-TAG pool[g]
%   RE3273R: H2O[r] + PI pool[r] <=> H+[r] + inositol[r] + phosphatidate-LD-TAG pool[r]
% should not be reversible. They should only proceed in the forward direction.
rxn_ind = ismember(ihuman.rxns,{'RE1448N';'RE3273C';'RE3273G';'RE3273R'});
ihuman.lb(rxn_ind) = 0;
ihuman.rev(rxn_ind) = 0;



% The following Recon3D reactions:
%
%        CDS: CTP[c] + H+[c] + phosphatidate-LD-TAG pool[c] => CDP-diacylglycerol-LD-PI pool[c] + PPi[c]
%       CDSm: CTP[m] + H+[m] + phosphatidate-LD-TAG pool[m] => CDP-diacylglycerol-LD-PI pool[m] + PPi[m]
%
% is inconsistent with a similar HMR reaction:
%
%   HMR_0607: CDP-diacylglycerol-LD-PI pool[c] + PPi[c] <=> CTP[c] + phosphatidate-LD-PI pool[c]
%
% Therefore, the Recon3D reaction should be removed.
del_rxns = [del_rxns; {'CDS';'CDSm'}];



% The following Recon3D reaction:
%
%   CHOLESACATc: cholesterol[c] + R Total Coenzyme A[c] => cholesterol-ester pool[c] + CoA[c]
%
% generates the metabolite "cholesterol-ester pool". However, this
% metabolite is treated differently by the HMR reactions, so this reaction
% creates a mass imbalance with existing reactions. Therefore, this
% reaction should be removed.
del_rxns = [del_rxns; {'CHOLESACATc'}];


% The following Recon3D reactions:
%
%       LPS: H2O[c] + TAG-VLDL pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + H+[c] + R Total 3 Position[c]
%      LPSe: H2O[s] + TAG-VLDL pool[s] => 1,2-diacylglycerol-LD-TAG pool[s] + H+[s] + R Total 3 Position[s]
%      DGAT: 1,2-diacylglycerol-LD-TAG pool[c] + R Total 3 Coenzyme A[c] => CoA[c] + TAG-VLDL pool[c]
% 
% are inconsistent with existing HMR reactions, e.g.:
%
%  HMR_0007: H2O[s] + TAG-VLDL pool[s] => 1,2-diacylglycerol-VLDL pool[s] + fatty acid-VLDL pool[s]
%
% This leads to a mass imbalance, and therefore the Recon3D reactions
% should be removed.
del_rxns = [del_rxns; {'LPS';'LPSe';'DGAT'}];



% The following reactions from HMR:
%
% HMR_3537: cholesterol-ester pool[l] => 0.0001 cholesterol-ester-10,13,16,19-docosa[l] + 0.0001 cholesterol-ester-10,13,16-docosa[l] + 0.0001 cholesterol-ester-10-hepta[l] + 0.0001 cholesterol-ester-11,14,17-eico[l] + 0.0001 cholesterol-ester-11,14-eicosa[l] + 0.0001 cholesterol-ester-11-docose[l] + 0.0001 cholesterol-ester-11-eico[l] + 0.0001 cholesterol-ester-12,15,18,21-tetracosa[l] + 0.0001 cholesterol-ester-13,16,19-doco[l] + 0.0001 cholesterol-ester-13,16-docosa[l] + 0.0001 cholesterol-ester-13-docose[l] + 0.0001 cholesterol-ester-13-eicose[l] + 0.0001 cholesterol-ester-13-octade[l] + 0.0001 cholesterol-ester-15-tetra[l] + 0.0041 cholesterol-ester-4,7,10,13,16,19-doco[l] + 0.0001 cholesterol-ester-4,7,10,13,16-docosa[l] + 0.0071 cholesterol-ester-5,8,11,14,17-eico[l] + 0.0001 cholesterol-ester-5,8,11-eico[l] + 0.0001 cholesterol-ester-5-tetra[l] + 0.0001 cholesterol-ester-6,9,12,15,18,21-tetra[l] + 0.0001 cholesterol-ester-6,9,12,15,18-tetraco[l] + 0.0001 cholesterol-ester-6,9,12,15-octa[l] + 0.0001 cholesterol-ester-6,9-octa[l] + 0.0001 cholesterol-ester-7,10,13,16,19-docosa[l] + 0.0001 cholesterol-ester-7,10,13,16-docosa[l] + 0.0406 cholesterol-ester-7-hexa[l] + 0.0001 cholesterol-ester-7-octade[l] + 0.0001 cholesterol-ester-7-tetrade[l] + 0.0001 cholesterol-ester-8,11,14,17-eico[l] + 0.0001 cholesterol-ester-8,11-eico[l] + 0.0001 cholesterol-ester-9,12,15,18,21-tetra[l] + 0.0001 cholesterol-ester-9,12,15,18-tetraco[l] + 0.0001 cholesterol-ester-9-eicose[l] + 0.0001 cholesterol-ester-9-heptade[l] + 0.0001 cholesterol-ester-9-octa[l] + 0.0001 cholesterol-ester-9-tetrade[l] + 0.0518 cholesterol-ester-arach[l] + 0.0001 cholesterol-ester-cis-vac[l] + 0.005 cholesterol-ester-dihomo-gamma[l] + 0.0001 cholesterol-ester-docosa[l] + 0.0001 cholesterol-ester-eico[l] + 0.0001 cholesterol-ester-gamma-lin[l] + 0.0001 cholesterol-ester-heneico[l] + 0.0001 cholesterol-ester-hepta[l] + 0.0001 cholesterol-ester-hexacosa[l] + 0.0001 cholesterol-ester-hexecose[l] + 0.0001 cholesterol-ester-laur[l] + 0.5254 cholesterol-ester-lin[l] + 0.0061 cholesterol-ester-linolen[l] + 0.0081 cholesterol-ester-myrist[l] + 0.0001 cholesterol-ester-nanode[l] + 0.1958 cholesterol-ester-ol[l] + 0.138 cholesterol-ester-palm[l] + 0.0001 cholesterol-ester-palmn[l] + 0.0001 cholesterol-ester-penta[l] + 0.0132 cholesterol-ester-stea[l] + 0.0001 cholesterol-ester-tetraco[l] + 0.0001 cholesterol-ester-trico[l] + 0.0001 cholesterol-ester-tridec[l]
% HMR_3622: cholesterol-ester pool[r] <=> 0.0001 cholesterol-ester-10,13,16,19-docosa[r] + 0.0001 cholesterol-ester-10,13,16-docosa[r] + 0.0001 cholesterol-ester-10-hepta[r] + 0.0001 cholesterol-ester-11,14,17-eico[r] + 0.0001 cholesterol-ester-11,14-eicosa[r] + 0.0001 cholesterol-ester-11-docose[r] + 0.0001 cholesterol-ester-11-eico[r] + 0.0001 cholesterol-ester-12,15,18,21-tetracosa[r] + 0.0001 cholesterol-ester-13,16,19-doco[r] + 0.0001 cholesterol-ester-13,16-docosa[r] + 0.0001 cholesterol-ester-13-docose[r] + 0.0001 cholesterol-ester-13-eicose[r] + 0.0001 cholesterol-ester-13-octade[r] + 0.0001 cholesterol-ester-15-tetra[r] + 0.0041 cholesterol-ester-4,7,10,13,16,19-doco[r] + 0.0001 cholesterol-ester-4,7,10,13,16-docosa[r] + 0.0071 cholesterol-ester-5,8,11,14,17-eico[r] + 0.0001 cholesterol-ester-5,8,11-eico[r] + 0.0001 cholesterol-ester-5-tetra[r] + 0.0001 cholesterol-ester-6,9,12,15,18,21-tetra[r] + 0.0001 cholesterol-ester-6,9,12,15,18-tetraco[r] + 0.0001 cholesterol-ester-6,9,12,15-octa[r] + 0.0001 cholesterol-ester-6,9-octa[r] + 0.0001 cholesterol-ester-7,10,13,16,19-docosa[r] + 0.0001 cholesterol-ester-7,10,13,16-docosa[r] + 0.0406 cholesterol-ester-7-hexa[r] + 0.0001 cholesterol-ester-7-octade[r] + 0.0001 cholesterol-ester-7-tetrade[r] + 0.0001 cholesterol-ester-8,11,14,17-eico[r] + 0.0001 cholesterol-ester-8,11-eico[r] + 0.0001 cholesterol-ester-9,12,15,18,21-tetra[r] + 0.0001 cholesterol-ester-9,12,15,18-tetraco[r] + 0.0001 cholesterol-ester-9-eicose[r] + 0.0001 cholesterol-ester-9-heptade[r] + 0.0001 cholesterol-ester-9-octa[r] + 0.0001 cholesterol-ester-9-tetrade[r] + 0.0518 cholesterol-ester-arach[r] + 0.0001 cholesterol-ester-cis-vac[r] + 0.005 cholesterol-ester-dihomo-gamma[r] + 0.0001 cholesterol-ester-docosa[r] + 0.0001 cholesterol-ester-eico[r] + 0.0001 cholesterol-ester-gamma-lin[r] + 0.0001 cholesterol-ester-heneico[r] + 0.0001 cholesterol-ester-hepta[r] + 0.0001 cholesterol-ester-hexacosa[r] + 0.0001 cholesterol-ester-hexecose[r] + 0.0001 cholesterol-ester-laur[r] + 0.5254 cholesterol-ester-lin[r] + 0.0061 cholesterol-ester-linolen[r] + 0.0081 cholesterol-ester-myrist[r] + 0.0001 cholesterol-ester-nanode[r] + 0.1958 cholesterol-ester-ol[r] + 0.138 cholesterol-ester-palm[r] + 0.0001 cholesterol-ester-palmn[r] + 0.0001 cholesterol-ester-penta[r] + 0.0132 cholesterol-ester-stea[r] + 0.0001 cholesterol-ester-tetraco[r] + 0.0001 cholesterol-ester-trico[r] + 0.0001 cholesterol-ester-tridec[r]
% HMR_0012: fatty acid-uptake pool[s] => 0.0001 (10Z)-heptadecenoic acid[s] + 0.0001 (11Z,14Z)-eicosadienoic acid[s] + 0.0001 (11Z,14Z,17Z)-eicosatrienoic acid[s] + 0.0001 (13Z)-eicosenoic acid[s] + 0.0001 (13Z)-octadecenoic acid[s] + 0.0001 (13Z,16Z)-docosadienoic acid[s] + 0.0001 (4Z,7Z,10Z,13Z,16Z)-DPA[s] + 0.0001 (6Z,9Z)-octadecadienoic acid[s] + 0.0001 (6Z,9Z,12Z,15Z,18Z)-TPA[s] + 0.0001 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[s] + 0.0001 (7Z)-octadecenoic acid[s] + 0.0001 (7Z)-tetradecenoic acid[s] + 0.0001 (9E)-tetradecenoic acid[s] + 0.0001 (9Z,12Z,15Z,18Z)-TTA[s] + 0.0001 (9Z,12Z,15Z,18Z,21Z)-TPA[s] + 0.0001 10,13,16,19-docosatetraenoic acid[s] + 0.0001 10,13,16-docosatriynoic acid[s] + 0.0001 12,15,18,21-tetracosatetraenoic acid[s] + 0.0001 13,16,19-docosatrienoic acid[s] + 0.0001 7-palmitoleic acid[s] + 0.0001 8,11-eicosadienoic acid[s] + 0.0001 9-eicosenoic acid[s] + 0.0001 9-heptadecylenic acid[s] + 0.0001 adrenic acid[s] + 0.0082 arachidonate[s] + 0.0001 behenic acid[s] + 0.0001 cerotic acid[s] + 0.0001 cis-cetoleic acid[s] + 0.0001 cis-erucic acid[s] + 0.0001 cis-gondoic acid[s] + 0.0001 cis-vaccenic acid[s] + 0.0041 DHA[s] + 0.002 dihomo-gamma-linolenate[s] + 0.0001 DPA[s] + 0.0001 eicosanoate[s] + 0.0001 elaidate[s] + 0.001 EPA[s] + 0.0001 gamma-linolenate[s] + 0.0001 henicosanoic acid[s] + 0.0001 lauric acid[s] + 0.0001 lignocerate[s] + 0.1535 linoleate[s] + 0.0092 linolenate[s] + 0.0001 margaric acid[s] + 0.0001 mead acid[s] + 0.0338 myristic acid[s] + 0.0001 nervonic acid[s] + 0.0001 nonadecylic acid[s] + 0.3837 oleate[s] + 0.0001 omega-3-arachidonic acid[s] + 0.3015 palmitate[s] + 0.0522 palmitolate[s] + 0.0001 pentadecylic acid[s] + 0.0001 physeteric acid[s] + 0.046 stearate[s] + 0.0001 stearidonic acid[s] + 0.0001 tricosanoic acid[s] + 0.0001 tridecylic acid[s] + 0.0001 ximenic acid[s]
% HMR_0685: fatty acid-LD-TG1 pool[c] => 0.0001 (10Z)-heptadecenoic acid[c] + 0.0001 (11Z,14Z)-eicosadienoic acid[c] + 0.0001 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0001 (13Z)-eicosenoic acid[c] + 0.0001 (13Z)-octadecenoic acid[c] + 0.0001 (13Z,16Z)-docosadienoic acid[c] + 0.0001 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0001 (6Z,9Z)-octadecadienoic acid[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0001 (7Z)-octadecenoic acid[c] + 0.0001 (7Z)-tetradecenoic acid[c] + 0.0001 (9E)-tetradecenoic acid[c] + 0.0001 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0001 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0001 10,13,16,19-docosatetraenoic acid[c] + 0.0001 10,13,16-docosatriynoic acid[c] + 0.0001 12,15,18,21-tetracosatetraenoic acid[c] + 0.0001 13,16,19-docosatrienoic acid[c] + 0.0001 7-palmitoleic acid[c] + 0.0001 8,11-eicosadienoic acid[c] + 0.0001 9-eicosenoic acid[c] + 0.0001 9-heptadecylenic acid[c] + 0.0001 adrenic acid[c] + 0.0024 arachidonate[c] + 0.0001 behenic acid[c] + 0.0001 cerotic acid[c] + 0.0001 cis-cetoleic acid[c] + 0.0001 cis-erucic acid[c] + 0.0001 cis-gondoic acid[c] + 0.0208 cis-vaccenic acid[c] + 0.0005 DHA[c] + 0.0013 dihomo-gamma-linolenate[c] + 0.0017 DPA[c] + 0.0001 eicosanoate[c] + 0.0001 elaidate[c] + 0.0013 EPA[c] + 0.0048 gamma-linolenate[c] + 0.0001 henicosanoic acid[c] + 0.0001 lauric acid[c] + 0.0001 lignocerate[c] + 0.0741 linoleate[c] + 0.0037 linolenate[c] + 0.0001 margaric acid[c] + 0.0001 mead acid[c] + 0.0096 myristic acid[c] + 0.0001 nervonic acid[c] + 0.0001 nonadecylic acid[c] + 0.122 oleate[c] + 0.0001 omega-3-arachidonic acid[c] + 0.64 palmitate[c] + 0.0286 palmitolate[c] + 0.0001 pentadecylic acid[c] + 0.0001 physeteric acid[c] + 0.082 stearate[c] + 0.0028 stearidonic acid[c] + 0.0001 tricosanoic acid[c] + 0.0001 tridecylic acid[c] + 0.0001 ximenic acid[c]
% HMR_0686: fatty acid-LD-TG2 pool[c] => 0.0001 (10Z)-heptadecenoic acid[c] + 0.0001 (11Z,14Z)-eicosadienoic acid[c] + 0.0001 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0001 (13Z)-eicosenoic acid[c] + 0.0001 (13Z)-octadecenoic acid[c] + 0.0001 (13Z,16Z)-docosadienoic acid[c] + 0.0016 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0001 (6Z,9Z)-octadecadienoic acid[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0001 (7Z)-octadecenoic acid[c] + 0.0001 (7Z)-tetradecenoic acid[c] + 0.0001 (9E)-tetradecenoic acid[c] + 0.0001 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0001 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0001 10,13,16,19-docosatetraenoic acid[c] + 0.0001 10,13,16-docosatriynoic acid[c] + 0.0001 12,15,18,21-tetracosatetraenoic acid[c] + 0.0001 13,16,19-docosatrienoic acid[c] + 0.0001 7-palmitoleic acid[c] + 0.0001 8,11-eicosadienoic acid[c] + 0.0001 9-eicosenoic acid[c] + 0.0001 9-heptadecylenic acid[c] + 0.0032 adrenic acid[c] + 0.0521 arachidonate[c] + 0.0001 behenic acid[c] + 0.0001 cerotic acid[c] + 0.0001 cis-cetoleic acid[c] + 0.0001 cis-erucic acid[c] + 0.0001 cis-gondoic acid[c] + 0.0217 cis-vaccenic acid[c] + 0.0081 DHA[c] + 0.0023 dihomo-gamma-linolenate[c] + 0.0018 DPA[c] + 0.0001 eicosanoate[c] + 0.0001 elaidate[c] + 0.0003 EPA[c] + 0.003 gamma-linolenate[c] + 0.0001 henicosanoic acid[c] + 0.0001 lauric acid[c] + 0.0001 lignocerate[c] + 0.281 linoleate[c] + 0.0119 linolenate[c] + 0.0001 margaric acid[c] + 0.0001 mead acid[c] + 0.0024 myristic acid[c] + 0.0001 nervonic acid[c] + 0.0001 nonadecylic acid[c] + 0.4226 oleate[c] + 0.0001 omega-3-arachidonic acid[c] + 0.12 palmitate[c] + 0.0229 palmitolate[c] + 0.0001 pentadecylic acid[c] + 0.0001 physeteric acid[c] + 0.0384 stearate[c] + 0.0025 stearidonic acid[c] + 0.0001 tricosanoic acid[c] + 0.0001 tridecylic acid[c] + 0.0001 ximenic acid[c]
% HMR_0687: fatty acid-LD-TG3 pool[c] => 0.0001 (10Z)-heptadecenoic acid[c] + 0.0001 (11Z,14Z)-eicosadienoic acid[c] + 0.0001 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0001 (13Z)-eicosenoic acid[c] + 0.0001 (13Z)-octadecenoic acid[c] + 0.0001 (13Z,16Z)-docosadienoic acid[c] + 0.0011 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0001 (6Z,9Z)-octadecadienoic acid[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0001 (7Z)-octadecenoic acid[c] + 0.0001 (7Z)-tetradecenoic acid[c] + 0.0001 (9E)-tetradecenoic acid[c] + 0.0001 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0001 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0001 10,13,16,19-docosatetraenoic acid[c] + 0.0001 10,13,16-docosatriynoic acid[c] + 0.0001 12,15,18,21-tetracosatetraenoic acid[c] + 0.0001 13,16,19-docosatrienoic acid[c] + 0.0001 7-palmitoleic acid[c] + 0.0001 8,11-eicosadienoic acid[c] + 0.0001 9-eicosenoic acid[c] + 0.0001 9-heptadecylenic acid[c] + 0.0022 adrenic acid[c] + 0.0264 arachidonate[c] + 0.0001 behenic acid[c] + 0.0001 cerotic acid[c] + 0.0001 cis-cetoleic acid[c] + 0.0001 cis-erucic acid[c] + 0.0001 cis-gondoic acid[c] + 0.0194 cis-vaccenic acid[c] + 0.0221 DHA[c] + 0.0141 dihomo-gamma-linolenate[c] + 0.0038 DPA[c] + 0.0001 eicosanoate[c] + 0.0001 elaidate[c] + 0.0035 EPA[c] + 0.0081 gamma-linolenate[c] + 0.0001 henicosanoic acid[c] + 0.0001 lauric acid[c] + 0.0001 lignocerate[c] + 0.278 linoleate[c] + 0.0023 linolenate[c] + 0.0001 margaric acid[c] + 0.0001 mead acid[c] + 0.003 myristic acid[c] + 0.0001 nervonic acid[c] + 0.0001 nonadecylic acid[c] + 0.4519 oleate[c] + 0.0001 omega-3-arachidonic acid[c] + 0.0849 palmitate[c] + 0.014 palmitolate[c] + 0.0001 pentadecylic acid[c] + 0.0001 physeteric acid[c] + 0.0584 stearate[c] + 0.0026 stearidonic acid[c] + 0.0001 tricosanoic acid[c] + 0.0001 tridecylic acid[c] + 0.0001 ximenic acid[c]
% HMR_0688: fatty acid-LD-PC pool[c] <=> 0.0005 (10Z)-heptadecenoic acid[c] + 0.0005 (11Z,14Z)-eicosadienoic acid[c] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0005 (13Z)-eicosenoic acid[c] + 0.0005 (13Z)-octadecenoic acid[c] + 0.0005 (13Z,16Z)-docosadienoic acid[c] + 0.0039 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0005 (6Z,9Z)-octadecadienoic acid[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0005 (7Z)-octadecenoic acid[c] + 0.0005 (7Z)-tetradecenoic acid[c] + 0.0005 (9E)-tetradecenoic acid[c] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0005 10,13,16,19-docosatetraenoic acid[c] + 0.0005 10,13,16-docosatriynoic acid[c] + 0.0005 12,15,18,21-tetracosatetraenoic acid[c] + 0.0005 13,16,19-docosatrienoic acid[c] + 0.0005 7-palmitoleic acid[c] + 0.0005 8,11-eicosadienoic acid[c] + 0.0005 9-eicosenoic acid[c] + 0.0005 9-heptadecylenic acid[c] + 0.0011 adrenic acid[c] + 0.0709 arachidonate[c] + 0.0005 behenic acid[c] + 0.0005 cerotic acid[c] + 0.0005 cis-cetoleic acid[c] + 0.0005 cis-erucic acid[c] + 0.0005 cis-gondoic acid[c] + 0.0271 cis-vaccenic acid[c] + 0.0252 DHA[c] + 0.0228 dihomo-gamma-linolenate[c] + 0.0056 DPA[c] + 0.0005 eicosanoate[c] + 0.0005 elaidate[c] + 0.0109 EPA[c] + 0.0036 gamma-linolenate[c] + 0.0005 henicosanoic acid[c] + 0.0005 lauric acid[c] + 0.0005 lignocerate[c] + 0.2479 linoleate[c] + 0.0047 linolenate[c] + 0.0005 margaric acid[c] + 0.0005 mead acid[c] + 0.0054 myristic acid[c] + 0.0005 nervonic acid[c] + 0.0005 nonadecylic acid[c] + 0.1143 oleate[c] + 0.0178 omega-3-arachidonic acid[c] + 0.2781 palmitate[c] + 0.0105 palmitolate[c] + 0.0005 pentadecylic acid[c] + 0.0005 physeteric acid[c] + 0.1271 stearate[c] + 0.0026 stearidonic acid[c] + 0.0005 tricosanoic acid[c] + 0.0005 tridecylic acid[c] + 0.0005 ximenic acid[c]
% HMR_0689: fatty acid-LD-PE pool[c] => 0.0005 (10Z)-heptadecenoic acid[c] + 0.0005 (11Z,14Z)-eicosadienoic acid[c] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0005 (13Z)-eicosenoic acid[c] + 0.0005 (13Z)-octadecenoic acid[c] + 0.0005 (13Z,16Z)-docosadienoic acid[c] + 0.0052 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0005 (6Z,9Z)-octadecadienoic acid[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0005 (7Z)-octadecenoic acid[c] + 0.0005 (7Z)-tetradecenoic acid[c] + 0.0005 (9E)-tetradecenoic acid[c] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0005 10,13,16,19-docosatetraenoic acid[c] + 0.0005 10,13,16-docosatriynoic acid[c] + 0.0005 12,15,18,21-tetracosatetraenoic acid[c] + 0.0005 13,16,19-docosatrienoic acid[c] + 0.0005 7-palmitoleic acid[c] + 0.0005 8,11-eicosadienoic acid[c] + 0.0005 9-eicosenoic acid[c] + 0.0005 9-heptadecylenic acid[c] + 0.0007 adrenic acid[c] + 0.2125 arachidonate[c] + 0.0005 behenic acid[c] + 0.0005 cerotic acid[c] + 0.0005 cis-cetoleic acid[c] + 0.0005 cis-erucic acid[c] + 0.0005 cis-gondoic acid[c] + 0.0285 cis-vaccenic acid[c] + 0.0459 DHA[c] + 0.0411 dihomo-gamma-linolenate[c] + 0.0067 DPA[c] + 0.0005 eicosanoate[c] + 0.0005 elaidate[c] + 0.0221 EPA[c] + 0.0011 gamma-linolenate[c] + 0.0005 henicosanoic acid[c] + 0.0005 lauric acid[c] + 0.0005 lignocerate[c] + 0.1312 linoleate[c] + 0.0166 linolenate[c] + 0.0005 margaric acid[c] + 0.0005 mead acid[c] + 0.0319 myristic acid[c] + 0.0005 nervonic acid[c] + 0.0005 nonadecylic acid[c] + 0.0619 oleate[c] + 0.0163 omega-3-arachidonic acid[c] + 0.1243 palmitate[c] + 0.0327 palmitolate[c] + 0.0005 pentadecylic acid[c] + 0.0005 physeteric acid[c] + 0.1983 stearate[c] + 0.0025 stearidonic acid[c] + 0.0005 tricosanoic acid[c] + 0.0005 tridecylic acid[c] + 0.0005 ximenic acid[c]
% HMR_0690: fatty acid-LD-PS pool[c] => 0.0005 (10Z)-heptadecenoic acid[c] + 0.0005 (11Z,14Z)-eicosadienoic acid[c] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0005 (13Z)-eicosenoic acid[c] + 0.0005 (13Z)-octadecenoic acid[c] + 0.0005 (13Z,16Z)-docosadienoic acid[c] + 0.0055 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0005 (6Z,9Z)-octadecadienoic acid[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0005 (7Z)-octadecenoic acid[c] + 0.0005 (7Z)-tetradecenoic acid[c] + 0.0005 (9E)-tetradecenoic acid[c] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0005 10,13,16,19-docosatetraenoic acid[c] + 0.0005 10,13,16-docosatriynoic acid[c] + 0.0005 12,15,18,21-tetracosatetraenoic acid[c] + 0.0005 13,16,19-docosatrienoic acid[c] + 0.0005 7-palmitoleic acid[c] + 0.0005 8,11-eicosadienoic acid[c] + 0.0005 9-eicosenoic acid[c] + 0.0005 9-heptadecylenic acid[c] + 0.0043 adrenic acid[c] + 0.2005 arachidonate[c] + 0.0005 behenic acid[c] + 0.0005 cerotic acid[c] + 0.0005 cis-cetoleic acid[c] + 0.0005 cis-erucic acid[c] + 0.0005 cis-gondoic acid[c] + 0.0131 cis-vaccenic acid[c] + 0.0931 DHA[c] + 0.0089 dihomo-gamma-linolenate[c] + 0.0329 DPA[c] + 0.0005 eicosanoate[c] + 0.0005 elaidate[c] + 0.0286 EPA[c] + 0.0011 gamma-linolenate[c] + 0.0005 henicosanoic acid[c] + 0.0005 lauric acid[c] + 0.0005 lignocerate[c] + 0.0177 linoleate[c] + 0.0035 linolenate[c] + 0.0005 margaric acid[c] + 0.0005 mead acid[c] + 0.0055 myristic acid[c] + 0.0005 nervonic acid[c] + 0.0005 nonadecylic acid[c] + 0.0323 oleate[c] + 0.0174 omega-3-arachidonic acid[c] + 0.0337 palmitate[c] + 0.0058 palmitolate[c] + 0.0005 pentadecylic acid[c] + 0.0005 physeteric acid[c] + 0.4731 stearate[c] + 0.0025 stearidonic acid[c] + 0.0005 tricosanoic acid[c] + 0.0005 tridecylic acid[c] + 0.0005 ximenic acid[c]
% HMR_0691: fatty acid-LD-PI pool[c] => 0.0005 (10Z)-heptadecenoic acid[c] + 0.0005 (11Z,14Z)-eicosadienoic acid[c] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0005 (13Z)-eicosenoic acid[c] + 0.0005 (13Z)-octadecenoic acid[c] + 0.0005 (13Z,16Z)-docosadienoic acid[c] + 0.0124 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0005 (6Z,9Z)-octadecadienoic acid[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0005 (7Z)-octadecenoic acid[c] + 0.0005 (7Z)-tetradecenoic acid[c] + 0.0005 (9E)-tetradecenoic acid[c] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0005 10,13,16,19-docosatetraenoic acid[c] + 0.0005 10,13,16-docosatriynoic acid[c] + 0.0005 12,15,18,21-tetracosatetraenoic acid[c] + 0.0005 13,16,19-docosatrienoic acid[c] + 0.0005 7-palmitoleic acid[c] + 0.0005 8,11-eicosadienoic acid[c] + 0.0005 9-eicosenoic acid[c] + 0.0005 9-heptadecylenic acid[c] + 0.0031 adrenic acid[c] + 0.2118 arachidonate[c] + 0.0005 behenic acid[c] + 0.0005 cerotic acid[c] + 0.0005 cis-cetoleic acid[c] + 0.0005 cis-erucic acid[c] + 0.0005 cis-gondoic acid[c] + 0.0266 cis-vaccenic acid[c] + 0.0269 DHA[c] + 0.0228 dihomo-gamma-linolenate[c] + 0.0101 DPA[c] + 0.0005 eicosanoate[c] + 0.0005 elaidate[c] + 0.0043 EPA[c] + 0.0012 gamma-linolenate[c] + 0.0005 henicosanoic acid[c] + 0.0005 lauric acid[c] + 0.0005 lignocerate[c] + 0.0678 linoleate[c] + 0.0019 linolenate[c] + 0.0005 margaric acid[c] + 0.0005 mead acid[c] + 0.0056 myristic acid[c] + 0.0005 nervonic acid[c] + 0.0005 nonadecylic acid[c] + 0.136 oleate[c] + 0.0179 omega-3-arachidonic acid[c] + 0.0678 palmitate[c] + 0.0053 palmitolate[c] + 0.0005 pentadecylic acid[c] + 0.0005 physeteric acid[c] + 0.3553 stearate[c] + 0.0027 stearidonic acid[c] + 0.0005 tricosanoic acid[c] + 0.0005 tridecylic acid[c] + 0.0005 ximenic acid[c]
% HMR_0692: fatty acid-LD-SM pool[c] <=> 0.002 (10Z)-heptadecenoic acid[c] + 0.002 (11Z,14Z)-eicosadienoic acid[c] + 0.002 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.002 (13Z)-eicosenoic acid[c] + 0.002 (13Z)-octadecenoic acid[c] + 0.002 (13Z,16Z)-docosadienoic acid[c] + 0.002 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.002 (6Z,9Z)-octadecadienoic acid[c] + 0.002 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.002 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.002 (7Z)-octadecenoic acid[c] + 0.002 (7Z)-tetradecenoic acid[c] + 0.002 (9E)-tetradecenoic acid[c] + 0.002 (9Z,12Z,15Z,18Z)-TTA[c] + 0.002 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.002 10,13,16,19-docosatetraenoic acid[c] + 0.002 10,13,16-docosatriynoic acid[c] + 0.002 12,15,18,21-tetracosatetraenoic acid[c] + 0.002 13,16,19-docosatrienoic acid[c] + 0.002 7-palmitoleic acid[c] + 0.002 8,11-eicosadienoic acid[c] + 0.002 9-eicosenoic acid[c] + 0.002 9-heptadecylenic acid[c] + 0.002 adrenic acid[c] + 0.0136 arachidonate[c] + 0.002 behenic acid[c] + 0.002 cerotic acid[c] + 0.002 cis-cetoleic acid[c] + 0.002 cis-erucic acid[c] + 0.002 cis-gondoic acid[c] + 0.0453 cis-vaccenic acid[c] + 0.003 DHA[c] + 0.002 dihomo-gamma-linolenate[c] + 0.002 DPA[c] + 0.002 eicosanoate[c] + 0.002 elaidate[c] + 0.0068 EPA[c] + 0.002 gamma-linolenate[c] + 0.002 henicosanoic acid[c] + 0.002 lauric acid[c] + 0.002 lignocerate[c] + 0.025 linoleate[c] + 0.0075 linolenate[c] + 0.002 margaric acid[c] + 0.002 mead acid[c] + 0.011 myristic acid[c] + 0.002 nervonic acid[c] + 0.002 nonadecylic acid[c] + 0.061 oleate[c] + 0.002 omega-3-arachidonic acid[c] + 0.557 palmitate[c] + 0.0358 palmitolate[c] + 0.002 pentadecylic acid[c] + 0.002 physeteric acid[c] + 0.138 stearate[c] + 0.002 stearidonic acid[c] + 0.002 tricosanoic acid[c] + 0.002 tridecylic acid[c] + 0.002 ximenic acid[c]
% HMR_5257: fatty acid-LD-PC pool[r] => 0.0005 (10Z)-heptadecenoic acid[r] + 0.0005 (11Z,14Z)-eicosadienoic acid[r] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[r] + 0.0005 (13Z)-eicosenoic acid[r] + 0.0005 (13Z)-octadecenoic acid[r] + 0.0005 (13Z,16Z)-docosadienoic acid[r] + 0.0039 (4Z,7Z,10Z,13Z,16Z)-DPA[r] + 0.0005 (6Z,9Z)-octadecadienoic acid[r] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[r] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[r] + 0.0005 (7Z)-octadecenoic acid[r] + 0.0005 (7Z)-tetradecenoic acid[r] + 0.0005 (9E)-tetradecenoic acid[r] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[r] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[r] + 0.0005 10,13,16,19-docosatetraenoic acid[r] + 0.0005 10,13,16-docosatriynoic acid[r] + 0.0005 12,15,18,21-tetracosatetraenoic acid[r] + 0.0005 13,16,19-docosatrienoic acid[r] + 0.0005 7-palmitoleic acid[r] + 0.0005 8,11-eicosadienoic acid[r] + 0.0005 9-eicosenoic acid[r] + 0.0005 9-heptadecylenic acid[r] + 0.0011 adrenic acid[r] + 0.0709 arachidonate[r] + 0.0005 behenic acid[r] + 0.0005 cerotic acid[r] + 0.0005 cis-cetoleic acid[r] + 0.0005 cis-erucic acid[r] + 0.0005 cis-gondoic acid[r] + 0.0271 cis-vaccenic acid[r] + 0.0252 DHA[r] + 0.0228 dihomo-gamma-linolenate[r] + 0.0056 DPA[r] + 0.0005 eicosanoate[r] + 0.0005 elaidate[r] + 0.0109 EPA[r] + 0.0036 gamma-linolenate[r] + 0.0005 henicosanoic acid[r] + 0.0005 lauric acid[r] + 0.0005 lignocerate[r] + 0.2479 linoleate[r] + 0.0047 linolenate[r] + 0.0005 margaric acid[r] + 0.0005 mead acid[r] + 0.0054 myristic acid[r] + 0.0005 nervonic acid[r] + 0.0005 nonadecylic acid[r] + 0.1143 oleate[r] + 0.0178 omega-3-arachidonic acid[r] + 0.2781 palmitate[r] + 0.0105 palmitolate[r] + 0.0005 pentadecylic acid[r] + 0.0005 physeteric acid[r] + 0.1271 stearate[r] + 0.0026 stearidonic acid[r] + 0.0005 tricosanoic acid[r] + 0.0005 tridecylic acid[r] + 0.0005 ximenic acid[r]
%
% involve the breakdown of a pool into its components. These particular
% reactions do not appear to be fully or clearly balanced, or at least
% create problems when merged with Recon3D reactions. Therefore, these
% reactions will be restricted to only allow the GENERATION of the pool
% (i.e., the reverse direction in which it is written), but NOT the 
% breakdown into components.
rxn_ind = ismember(ihuman.rxns,{'HMR_3537';'HMR_3622';'HMR_0012';'HMR_0685';
                                'HMR_0686';'HMR_0687';'HMR_0688';'HMR_0689';
                                'HMR_0690';'HMR_0691';'HMR_0692';'HMR_5257'});
ihuman.lb(rxn_ind) = -1000;
ihuman.ub(rxn_ind) = 0;






%% Delete reactions

del_ind = ismember(ihuman.rxns,del_rxns);
ihuman.ub(del_ind) = 0;
ihuman.lb(del_ind) = 0;
ihuman.rev(del_ind) = 0;




%% Perform flux analysis to assess infinite ATP generation

model = simplifyModel(ihuman);
exch_inds = sum(model.S ~= 0) == 1;
model.S(:,exch_inds) = -abs(model.S(:,exch_inds));  % make all exch rxns the same direction (met --> 0)
model.ub(exch_inds) = 1000;   % specify export bound
model.lb(exch_inds) = 0;  % specify import bound

% allow import of specified mets
model.lb(ismember(model.rxns,{'HMR_9034'})) = -1000;  % glucose
model.lb(ismember(model.rxns,{'HMR_9048'})) = -1000;  % oxygen


% Set the new objective
% [~,obj_ind] = ismember('biomass_reaction',model.rxns);         % 0.50563 alanine[c] + 0.35926 arginine[c] + 0.27942 asparagine[c] + 0.35261 aspartate[c] + 20.7045 ATP[c] + 0.020401 cholesterol[c] + 0.011658 CL pool[c] + 0.039036 CTP[c] + 0.046571 cysteine[c] + 0.013183 dATP[n] + 0.009442 dCTP[n] + 0.009898 dGTP[n] + 0.013091 dTTP[n] + 0.27519 glucose-6-phosphate[c] + 0.38587 glutamate[c] + 0.326 glutamine[c] + 0.53889 glycine[c] + 0.036117 GTP[c] + 20.6508 H2O[c] + 0.12641 histidine[c] + 0.28608 isoleucine[c] + 0.54554 leucine[c] + 0.59211 lysine[c] + 0.15302 methionine[c] + 0.15446 PC-LD pool[c] + 0.055374 PE-LD pool[c] + 0.25947 phenylalanine[c] + 0.023315 PI pool[c] + 0.41248 proline[c] + 0.005829 PS-LD pool[c] + 0.39253 serine[c] + 0.017486 SM pool[c] + 0.31269 threonine[c] + 0.013306 tryptophan[c] + 0.15967 tyrosine[c] + 0.053446 UTP[c] + 0.35261 valine[c] + 0.002914 PG-CL pool[c] => 20.6508 ADP[c] + 20.6508 H+[c] + 20.6508 Pi[c]
[~,obj_ind] = ismember('biomass_components',model.rxns);       % alanine[c] + arginine[c] + asparagine[c] + aspartate[c] + cholesterol[c] + cholesterol-ester pool[r] + CL pool[c] + cysteine[c] + DNA[n] + DNA-5-methylcytosine[n] + glutamate[c] + glutamine[c] + glycine[c] + histidine[c] + isoleucine[c] + leucine[c] + lipid droplet[c] + lysine[c] + methionine[c] + phenylalanine[c] + phosphatidate-LD-TAG pool[c] + PI pool[c] + proline[c] + RNA[c] + serine[c] + SM pool[c] + threonine[c] + tryptophan[c] + tyrosine[c] + valine[c] + glycogen[c] + cofactors and vitamins[c] => biomass[x]
% [~,obj_ind] = ismember('HMR_biomass_Renalcancer',model.rxns);  % 0.0031466 1,2-diacylglycerol-LD-TAG pool[c] + 0.0038601 1-acylglycerol-3P-LD-TG1 pool[c] + 0.0012642 2-lysolecithin pool[c] + 0.11871 alanine[c] + 0.0009274 arginine[c] + 0.012056 asparagine[c] + 0.69277 aspartate[c] + 0.0082708 cholesterol[c] + 0.0065713 cholesterol-ester pool[r] + 0.00057356 CL pool[c] + 0.0012056 cysteine[c] + 0.007581 DNA[n] + 0.16786 glutamate[c] + 0.19104 glutamine[c] + 0.13911 glycine[c] + 0.028749 histidine[c] + 0.0037096 isoleucine[c] + 0.011129 leucine[c] + 0.009274 lysine[c] + 0.0018548 methionine[c] + 0.017171 PC-LD pool[c] + 0.020051 PE-LD pool[c] + 0.0037096 phenylalanine[c] + 0.016414 phosphatidate-LD-TAG pool[c] + 0.0030972 PI pool[c] + 0.017992 proline[c] + 0.005891 PS-LD pool[c] + 0.030815 RNA[c] + 0.037096 serine[c] + 0.0024143 SM pool[c] + 0.027482 TAG-LD pool[c] + 0.020403 threonine[c] + 0.0009274 tryptophan[c] + 0.0055644 tyrosine[c] + 0.012056 valine[c] + 0.26551 glycogen[c] + 0.001 cofactors and vitamins[c] => biomass[c]
% [~,obj_ind] = ismember('cofactors_vitamins',model.rxns);         % [protein]-N6-(lipoyl)lysine[c] + biotin[c] + CoA[c] + cobamide-coenzyme[m] + cytochrome-C[m] + FADH2[c] + L-carnitine[c] + NADH[c] + NADPH[c] + riboflavin[c] + tetrahydrobiopterin[c] + THF[c] + ubiquinol[m] + 0.1 vitamin A derivatives[c] + 0.1 vitamin D derivatives[c] + 0.1 vitamin E derivatives[c] => cofactors and vitamins[c]
% [~,obj_ind] = ismember('HMR_3964',model.rxns);  % ATP hydrolysis
% [~,obj_ind] = ismember('HMR_6916',model.rxns);  % ATP synthase
% [~,obj_ind] = ismember('HMR_9034',model.rxns);  % glucose export
% [~,obj_ind] = ismember('HMR_9058',model.rxns);  % CO2 export
model.c(:) = 0;
model.c(obj_ind) = 1;




% delete all sink and demand reactions
sink_dm = find(startsWith(ihuman.rxns,{'sink_','DM_'}));
model = removeReactionsFull(ihuman,sink_dm);

% simplify model
model = simplifyModel(model);

% specify metUptake
model = specifyMetUptake(model,mets);

