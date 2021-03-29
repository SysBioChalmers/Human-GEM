%
% FILE NAME:    constrainReactions.m
% 
% PURPOSE: Script to prepare a defined list of reactions for constraining:
%
%          1. reactions allow the creation of mass and/or energy (which
%          results in a "leaky" model).
%
%          2. reactions involve an identical set of metabolites except for 
%          one, and that one different metabolite does not have the same 
%          mass in each reaction. For example, the following reactions were
%          identified as a set of "mass variable reactions":
%
%           LCAT39e: cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosahexenoylglycerophosphocholine (Delta 4, 7, 10, 13, 16, 19), Sn1-Lpc (22:6)[s]
%            LCAT5e: cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosadienoylglycerophosphocholine (Delta 11,14)[s]
%           LCAT31e: cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Octadeca-Trienoylglycerophosphocholine, Sn1-Lpc (18:3, Delta 6, 9, 12)[s]
%
%
% NOTE: These reactions will be constrained to zero for now ("inactivated")
% and scheduled for potential future "hard deletion" from the model.
%


%% Initialize some variables

rxnNotes = {};
del_rxns = {};


%% Group1: reactions allow the creation of mass and/or energy and results in a "leaky" model

% These reactions concern the glycerol phosphate shuttle:
%
%  HMR_0483: DHAP[c] + ubiquinol[m] => sn-glycerol-3-phosphate[c] + ubiquinone[m]
%  HMR_0482: DHAP[c] + FADH2[c] => FAD[c] + sn-glycerol-3-phosphate[c]
%    r0202m: NAD+[m] + sn-glycerol-3-phosphate[m] => DHAP[m] + H+[m] + NADH[m]
%
% The HMR reactions are written in the wrong direction, and would be
% reversed in script repairModelLeaks.m. The Recon3D reaction (r0202m) is 
% also in the wrong direction, but should not take place in the
% mitochondria, and therefore be deleted here.
del_rxns = [del_rxns; {'r0202m'}];
rxnNotes = [rxnNotes; {'r0202m', 'reaction should not take place in mitochondria, should be DELETED'}];


% This reaction is similar to an existing HMR rxn, but is reversible and
% missing the FAD(H2) metabolite:
%
%    RE1573M: cis,cis-3,6-dodecadienoyl-CoA[m] <=> trans,cis-lauro-2,6-dienoyl-CoA[m]
%   HMR_3288: cis,cis-3,6-dodecadienoyl-CoA[m] + FAD[m] => FADH2[m] + trans,cis-lauro-2,6-dienoyl-CoA[m]
%
% Therefore, the Recon3D version of the reaction should be deleted.
del_rxns = [del_rxns; {'RE1573M'}];
rxnNotes = [rxnNotes; {'RE1573M', 'reaction is missing FADH2 and would be identical to HMR_3288, so it should be DELETED'}];


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
rxnNotes = [rxnNotes; {'r1453', 'reaction is same as HMR_3838 but ubiquinol is on wrong side of equation; should be DELETED'}];


% The following reaction from Recon3D:
%
%   r0698: chenodeoxycholoyl-CoA[p] + 4 H+[p] + propanoyl-CoA[p] => 25(R)DHCA-CoA[p] + CoA[p] + H2O[p]
%
% has no evidence supporting its existence, and is charge-imbalanced (there
% is no source providing the electrons to reduce the four protons). It will
% therefore be removed from the model.
del_rxns = [del_rxns; {'r0698'}];
rxnNotes = [rxnNotes; {'r0698', 'no evidence supporting such a reaction, and missing electron source; should be DELETED'}];


% The model contains the following two reactions from Recon3D:
%
%   DHCR241r: FADH2[r] + zymosterol[r] => FAD[r] + cholestenol[r]
%      r1380: H+[r] + NADPH[r] + zymosterol[r] <=> NADP+[r] + cholestenol[r]
%
% The reactions are identical, except one uses FADH2, whereas the other
% uses NADPH. Based on the available databases, the reaction should use
% NADPH. Therefore, the first reaction should be removed, since it will be
% identical to the second after replacing its cofactor with NADPH.
% A similar situation was found for the following reaction pair:
%
%   DHCR242r: 5alpha-cholesta-7,24-dien-3beta-ol[r] + FADH2[r] => lathosterol[r] + FAD[r]
%   HMR_1533: 5alpha-cholesta-7,24-dien-3beta-ol[c] + H+[c] + NADPH[c] => lathosterol[c] + NADP+[c]
%
del_rxns = [del_rxns; {'DHCR241r'; 'DHCR242r'}];
rxnNotes = [rxnNotes; {'DHCR241r', 'reaction is identical to r1380, but uses incorrect cofactor (FADH2); should be DELETED'}];
rxnNotes = [rxnNotes; {'DHCR242r', 'reaction is identical to HMR_1533, but uses incorrect cofactor (FADH2); should be DELETED'}];


% The following reaction from Recon3D:
%
%   r1479: CoA[m] + 3-Oxolaur-Cis-5-Enoyl Coenzyme A[m] => (3Z)-dodecenoyl-CoA[m] + acetyl-CoA[m]
%
% was found as part of a reaction loop creating energy and mass. This
% reaction is imbalanced (generates C2H2), and was not found on the
% BiGG database, suggesting that it has been removed or updated recently.
% Therefore, this reaction should be removed from the model.
del_rxns = [del_rxns; {'r1479'}];
rxnNotes = [rxnNotes; {'r1479', 'reaction is mass-imbalanced and generates carbon, and should therefore be DELETED'}];


% The following reactions from Recon3D:
%
%   FAOXC2251836m: (4Z,7Z,10Z,13Z,16Z)-docosapentaenoyl-CoA[m] + 2 CoA[m] + 2 H2O[m] + 2 NAD+[m] => 2 acetyl-CoA[m] + gamma-linolenoyl-CoA[m] + 2 H+[m] + 2 NADH[m]
%   FAOXC2251836x: (4Z,7Z,10Z,13Z,16Z)-docosapentaenoyl-CoA[p] + 2 CoA[p] + 2 H2O[p] + 2 NAD+[p] => 2 acetyl-CoA[p] + gamma-linolenoyl-CoA[p] + 2 H+[p] + 2 NADH[p]
%
% enable the biosynthesis of linolenate, which is an essential fatty acid
% for humans (i.e., humans are unable to synthesize this compound). 
% Furthermore, there are no references associated with these reactions, and
% databases (KEGG, METACYC, etc.) do not support the existence of such a 
% reaction. Therefore, these reactions should be removed.
del_rxns = [del_rxns; {'FAOXC2251836m';'FAOXC2251836x'}];
rxnNotes = [rxnNotes; [{'FAOXC2251836m';'FAOXC2251836x'}, repmat({'reaction enables biosynthesis of an essential fatty acid (linoleic acid), which is physiologically incorrect. Reaction should be DELETED, unless sufficient evidence supporting its inclusion can be provided.'},2,1)]];


% The following reactions from HMR:
%
%   HMR_3451: 3-oxo-dihomo-gamma-linolenoyl-CoA[m] + CoA[m] => acetyl-CoA[m] + gamma-linolenoyl-CoA[m]
%   HMR_3465: 3-oxo-dihomo-gamma-linolenoyl-CoA[p] + CoA[p] => acetyl-CoA[p] + gamma-linolenoyl-CoA[p]
%
% enable the biosynthesis of linolenate, which is an essential fatty acid
% for humans that cannot be produced; only the reverse of these reactions
% are supported by literature (though with a slight difference):
%
%   HMR_2371: gamma-linolenoyl-CoA[c] + H+[c] + malonyl-CoA[c] => 3-oxo-dihomo-gamma-linolenoyl-CoA[c] + CO2[c] + CoA[c]
%    RE3103R: gamma-linolenoyl-CoA[r] + H+[r] + malonyl-CoA[r] => 3-oxo-dihomo-gamma-linolenoyl-CoA[r] + CO2[r] + CoA[r]
%
% The references associated with the problematic reactions (HMR_3451 and 
% HMR_3465) do not contain any evidence supporting such a reaction.
% Therefore, these reactions should be removed from the model.
del_rxns = [del_rxns; {'HMR_3451';'HMR_3465'}];
rxnNotes = [rxnNotes; [{'HMR_3451';'HMR_3465'}, repmat({'reaction enables biosynthesis of an essential fatty acid (linoleic acid), which is physiologically incorrect. Reaction should be DELETED, unless sufficient evidence supporting its inclusion can be provided.'},2,1)]];


% The following reaction from Recon3D:
%
%   r1169: cholesterol[r] + gamma-linolenoyl-CoA[r] => cholesterol-ester-linolen[r] + CoA[r]
%
% treats gamma-linolenoyl-CoA as equivalent to linolenoyl-CoA, which is not
% the case; see for example the following reactions:
%
%   HMR_3720: cholesterol-ester-linolen[r] + H2O[r] => cholesterol[r] + H+[r] + linolenate[r]
%   HMR_3674: cholesterol[r] + gamma-linolenoyl-CoA[r] => cholesterol-ester-gamma-lin[r] + CoA[r]
%
% Therefore, the Recon-derived reaction should be removed from the model.
del_rxns = [del_rxns; {'r1169'}];
rxnNotes = [rxnNotes; {'r1169', 'reaction incorrectly treats gamma-linolenoyl-CoA and linolenoyl-CoA as equivalent (see e.g. HMR_3720 and HMR_3674), and should therefore be DELETED'}];


% The following reaction from Recon3D:
%
%   DOPACCL: dopamine-O-quinone[c] => H+[c] + leukoaminochrome[c]
%
% has no gene associated with the reaction, and no sources, and it could
% not be found in the literature. A paper (PMID: 20600874) shows part of
% the pathway of dopamine-derived quinone metabolism, where it is clear
% that the above reaction would not take place. Therefore it should be
% removed from the model.
del_rxns = [del_rxns; {'DOPACCL'}];
rxnNotes = [rxnNotes; {'DOPACCL', 'no sources supporting this reaction, and literature (PMID: 20600874) suggests that it could not occur; reaction should be DELETED'}];


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
% Therefore, these Recon3D reactions should be removed from the model.
%
% In addition, there were a few other reactions that did not have an HMR
% equivalent, but would have their stoich coeffs adjusted from 0.1 to 1, to
% be consistent with how other reactions in the model are treating the mass
% of these dolichol compounds, in script repairModelLeaks.m.
%
%      H8MTer_L: 0.1 dolichyl-phosphate-D-mannose[r] + gpi heparan sulfate[r] => 0.1 dolichyl-phosphate[r] + H+[r] + Mannosyl-3-(Phosphoethanolaminyl-Mannosyl)-Glucosaminyl-Acylphosphatidylinositol (M4A)[r]
%      H8MTer_U: 0.1 dolichyl-phosphate-D-mannose[r] + gpi heparan sulfate[r] => 0.1 dolichyl-phosphate[r] + H+[r] + HMA[r]
%    UDPDOLPT_L: 0.1 dolichyl-phosphate[c] + UDP-glucose[c] => UDP[c] + 0.1 dolichyl-D-glucosyl-phosphate[c]
%
rxns = {'DOLGPP_Ler';'DOLASNT_Ler';'DOLDPP_Ler';'DOLK_L';'DOLMANP_Lter';...
                       'DOLPMT3_Ler';'GPIMTer_L';'GLCNACPT_L';'DOLPGT3_Ler';'DOLICHOL_Lter';'DEDOLR_L'};
del_rxns = [del_rxns; rxns];
rxnNotes = [rxnNotes; [rxns, repmat({'reaction assumes dolichol-related mets have greater mass than other rxns in model, thus creating mass imbalances; rxn should be DELETED'},length(rxns),1)]];



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
rxns = {'VLDL_HSDEG';'IDL_HSDEG';'LDL_HSDEG';'HDL_HSDEG';'CHYLO_HSDEG';
                       'VLDL_HSSYN';'IDL_HSSYN';'LDL_HSSYN';'HDL_HSSYN';'MYELIN_HSSYN';'CHYLO_HSSYN'};
del_rxns = [del_rxns; rxns];
rxnNotes = [rxnNotes; [rxns, repmat({'treatment of TAG-VLDL pool here is different than other (HMR) reactions in the model, resulting in mass imbalances; rxn should be DELETED'},length(rxns),1)]];


% The following two reactions from HMR:
%
%   HMR_2128: calcitriol[m] + H+[m] + NADPH[m] + O2[m] => calcitroic acid[m] + H2O[m] + NADP+[m]
%   HMR_2141: calcitriol[m] + H+[m] + NADPH[m] + O2[m] => calcitetrol[m] + H2O[m] + NADP+[m]
%
% are identical except for one of the products. An investigation of these
% metabolites revealed that the first reaction is imbalanced, as calcitroic
% acid has 4 fewer carbons than calcitriol (or calcitetrol). Therefore, the
% first reaction (HMR_2128) should be removed from the model.
del_rxns = [del_rxns; {'HMR_2128'}];
rxnNotes = [rxnNotes; {'HMR_2128', 'rxn is mass imbalanced, but correction would result in identical reaction as HMR_2141; should therefore be DELETED'}];


% The following reactions from Recon3D:
%
%   RE3273C:  H2O[c] + PI pool[c] <=> H+[c] + inositol[c] + phosphatidate-LD-TAG pool[c]
%   RE3273G:  H2O[g] + PI pool[g] <=> H+[g] + inositol[g] + phosphatidate-LD-TAG pool[g]
%   RE3273R:  H2O[r] + PI pool[r] <=> H+[r] + inositol[r] + phosphatidate-LD-TAG pool[r]
%
% Treat the mass of the PI pool and/or phosphatidate-LD-TAG pool different
% from other reactions in the model (e.g., HMR_0610), and are also
% inconsistent with the treatment among many of the reactions from Recon3D.
% Therefore, these reactions should be constrained to zero until they can
% be properly re-balanced, or removed entirely.
rxns = {'RE3273C';'RE3273G';'RE3273R'};
del_rxns = [del_rxns; rxns];
rxnNotes = [rxnNotes; [rxns, repmat({'treatment of PI pool and/or phosphatidate-LD-TAG pool here is different than other (HMR) reactions in the model, resulting in mass imbalances; rxn should therefore be constrained until imbalances can be addressed, otherwise DELETED'},length(rxns),1)]];


% The following reactions from Recon3D:
%
%    DHAPA: DHAP[c] + R Total Coenzyme A[c] => acylglycerone-phosphate[c] + CoA[c]
%   DHAPAx: DHAP[p] + R Total Coenzyme A[p] => acylglycerone-phosphate[p] + CoA[p]
%
% are imbalanced, because Recon3D treats "R Total Coenzyme A" as equivalent
% to palmitoyl-CoA. This leads to reactions that involve 3 + 37 = 40
% carbons consumed to produce 4 + 21 = 25 carbons. They should therefore be
% removed from the model.
del_rxns = [del_rxns; {'DHAPA';'DHAPAx'}];
rxnNotes = [rxnNotes; [{'DHAPA';'DHAPAx'}, repmat({'R Total Coenzyme A is effectively equal to palmitoyl-CoA, so the reaction is therefore mass-imbalanced; should be DELETED'},2,1)]];


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
rxnNotes = [rxnNotes; {'DSAT', 'rxn treats mass of dihydroceramide pool metabolite differently than others in model (e.g. HMR_0753, HMR_0692), resulting in mass imbalances; rxn should therefore be constrained until imbalances can be addressed, otherwise DELETED'}];


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
rxns = {'RE3301C';'RE3301G';'RE3301R'};
del_rxns = [del_rxns; rxns];
rxnNotes = [rxnNotes; [rxns, repmat({'rxn treats PS-LD/phosphatidate-LD-TAG pool differently than others in model (e.g., HMR_0660), resulting in mass imbalances; rxn should therefore be constrained until imbalances can be addressed, otherwise DELETED'},length(rxns),1)]];


% The following Recon3D reactions:
%
%        CDS: CTP[c] + H+[c] + phosphatidate-LD-TAG pool[c] => CDP-diacylglycerol-LD-PI pool[c] + PPi[c]
%       CDSm: CTP[m] + H+[m] + phosphatidate-LD-TAG pool[m] => CDP-diacylglycerol-LD-PI pool[m] + PPi[m]
%
% are inconsistent with a similar HMR reaction:
%
%   HMR_0607: CDP-diacylglycerol-LD-PI pool[c] + PPi[c] <=> CTP[c] + phosphatidate-LD-PI pool[c]
%
% Therefore, the Recon3D reaction should be constrained until properly 
% re-balanced/integrated.
del_rxns = [del_rxns; {'CDS';'CDSm'}];
rxnNotes = [rxnNotes; [{'CDS';'CDSm'}, repmat({'rxn treats phosphatidate-LD-TAG pool and/or CDP-diacylglycerol-LD-PI pool differently than others in the model (e.g., HMR_0607), resulting in mass imbalances; rxn should therefore be constrained until imbalances can be addressed, otherwise DELETED'},2,1)]];


% The following Recon3D reaction:
%
%   CHOLESACATc: cholesterol[c] + R Total Coenzyme A[c] => cholesterol-ester pool[c] + CoA[c]
%
% generates the metabolite "cholesterol-ester pool". However, this
% metabolite is treated differently by the HMR reactions, so this reaction
% creates a mass imbalance with existing reactions. Therefore, this
% reaction should be constrained until it can be properly integrated.
del_rxns = [del_rxns; {'CHOLESACATc'}];
rxnNotes = [rxnNotes; {'CHOLESACATc', 'rxn treats mass of cholesterol-ester pool metabolite differently than others in model, resulting in mass imbalances; rxn should therefore be constrained until imbalances can be addressed, otherwise DELETED'}];


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
rxns = {'LPS';'LPSe';'DGAT'};
del_rxns = [del_rxns; rxns];
rxnNotes = [rxnNotes; [rxns, repmat({'rxn treats TAG-VLDL pool and/or 1,2-diacylglycerol-LD-TAG pool differently than others in model (e.g., HMR_0007), resulting in mass imbalances; rxn should therefore be constrained until imbalances can be addressed, otherwise DELETED'},length(rxns),1)]];


% The following Recon3D reaction:
%
%   r0626: NAD+[c] + 3alpha,7alpha-dihydroxy-5beta-cholest-24-enoyl-CoA[c] => 3alpha,7alpha,12alpha-trihydroxy-5beta-cholestan-26-al[c] + NADH[c]
%
% is mass-imbalanced (for CoA and additional elements), and should
% therefore be constrained or deleted unless it can be properly
% re-balanced.
del_rxns = [del_rxns; {'r0626'}];
rxnNotes = [rxnNotes; {'r0626', 'reaction is mass-imbalanced and should therefore be constrained until imbalances can be addressed, otherwise DELETED'}];


% The following Recon3D reaction:
%
%   r1386: lysine[c] => procollagen-L-lysine[c]
%
% is mass-imbalanced, as procollagen-L-lysine has more carbons that lysine.
% The reaction should therefore be deleted, unless it can be corrected.
del_rxns = [del_rxns; {'r1386'}];
rxnNotes = [rxnNotes; {'r1386', 'reaction is mass-imbalanced and should therefore be constrained until imbalances can be addressed, otherwise DELETED'}];


% The following Recon3D reaction:
%
%   HC02191c: H2O[c] + NADP+[c] + 3beta-hydroxy-5-cholestenal[c] => 2 H+[c] + lithocholate[c] + NADPH[c]
%
% is carbon-imbalanced, and the associated references (PMIDs) do not
% provide any evidence for such a reaction. It should therefore be deleted.
del_rxns = [del_rxns; {'HC02191c'}];
rxnNotes = [rxnNotes; {'HC02191c', 'reaction is carbon-imbalanced and should therefore be constrained until imbalances can be addressed, otherwise DELETED'}];


% The following Recon3D reaction:
%
%   r1254: ATP[c] + CoA[c] + 8 H+[c] + stearidonic acid[c] => AMP[c] + PPi[c] + stearoyl-CoA[c]
%
% should be producing stearidonoyl-CoA, NOT stearoyl-CoA. There exists a
% similar HMR reaction:
%
%   HMR_0353: (6Z,9Z,12Z,15Z)-octadecatetraenoyl-CoA[c] + AMP[c] + PPi[c] <=> ATP[c] + CoA[c] + stearidonic acid[c]
%
% where (6Z,9Z,12Z,15Z)-octadecatetraenoyl-CoA is stearidonoyl-CoA.
% Therefore, the HMR reaction should be kept, while the Recon3D reaction
% shoud be removed.
del_rxns = [del_rxns; {'r1254'}];
rxnNotes = [rxnNotes; {'r1254', 'rxn should produce stearidonoyl-CoA, NOT stearoyl-CoA (see HMR_0353); rxn should therefore be DELETED'}];


% The following Recon3D reactions:
%
%   RE3267E: CDP-ethanolamine[s] + 1,2-Diacyl-Sn-Glycerol (Didodecanoyl, N-C12:0)[s] => CMP[s] + H+[s] + PE-LD pool[s]
%   RE3267G: CDP-ethanolamine[g] + 1,2-Diacyl-Sn-Glycerol (Didodecanoyl, N-C12:0)[g] => CMP[g] + H+[g] + PE-LD pool[g]
%   RE3267M: CDP-ethanolamine[m] + 1,2-Diacyl-Sn-Glycerol (Didodecanoyl, N-C12:0)[m] => CMP[m] + H+[m] + PE-LD pool[m]
%   RE3267N: CDP-ethanolamine[n] + 1,2-Diacyl-Sn-Glycerol (Didodecanoyl, N-C12:0)[n] => CMP[n] + H+[n] + PE-LD pool[n]
%   RE3267R: CDP-ethanolamine[r] + 1,2-Diacyl-Sn-Glycerol (Didodecanoyl, N-C12:0)[r] => CMP[r] + H+[r] + PE-LD pool[r]
%
% are similar to the HMR reaction:
%
%   HMR_0614: 1,2-diacylglycerol-LD-PE pool[c] + CDP-ethanolamine[c] => CMP[c] + PE-LD pool[c]
%
% However, the Recon3D reactions involve the reaction of two non-pool
% metabolites with set masses (CDP-ethanolamine and 1,2-Diacyl-Sn-Glycerol)
% to form a pool metabolite (PE-LD pool), which results in a different
% assumed elemental composition of PE-LD pool than the HMR reaction.
% Therefore, these Recon3D reactions should be constrained/removed.
rxns = {'RE3267E';'RE3267G';'RE3267M';'RE3267N';'RE3267R'};
del_rxns = [del_rxns; rxns];
rxnNotes = [rxnNotes; [rxns, repmat({'rxn treats PE-LD pool differently than others in model (e.g., HMR_0614), resulting in mass imbalances; rxn should therefore be constrained until imbalances can be addressed, otherwise DELETED'},length(rxns),1)]];


% The following reactions from Recon3D:
%
%   r0001: S-adenosylmethioninamine[c] => 5-methylthioadenosine[c] + Adenosylmethioninamine-Potential[c]
%   r1319: ATP[c] + H2O[c] => ADP[c] + Pi[c] + Adenosine-5'-Triphosphate-Energy[c]
%   r1320: ATP[m] + H2O[m] => ADP[m] + Pi[m] + Adenosine-5'-Triphosphate-Energy[m]
%   r1321: NADH[r] => NAD+[r] + Nadh-Redox-Potential[r]
%   r1322: NADH[c] => NAD+[c] + Nadh-Redox-Potential[c]
%   r1323: NADH[m] => NAD+[m] + Nadh-Redox-Potential[m]
%   r1324: NADH[p] => NAD+[p] + Nadh-Redox-Potential[p]
%   r1325: NADPH[r] => NADP+[r] + Nadph-Redox-Potential[r]
%   r1326: NADPH[c] => NADP+[c] + Nadph-Redox-Potential[c]
%   r1327: NADPH[m] => NADP+[m] + Nadph-Redox-Potential[m]
%   r1328: NADPH[p] => NADP+[p] + Nadph-Redox-Potential[p]
%   r1329: FADH2[c] => FAD[c] + Fadh-Redox-Potential[c]
%
% all involve the generation of an artificial "Potential" or "Energy"
% metabolite, which cannot be balanced (i.e., they are all dead-end
% reactions). They should therefore be removed from the model.
rxns = {'r0001';'r1319';'r1320';'r1321';'r1322';'r1323';'r1324';'r1325';'r1326';'r1327';'r1328';'r1329'};
del_rxns = [del_rxns; rxns];
rxnNotes = [rxnNotes; [rxns, repmat({'reaction is artificial and involves dead-end artificial metabolite with no apparent purpose, and should therefore be DELETED'},length(rxns),1)]];


% The following reactions from HMR:
%
% HMR_0689: fatty acid-LD-PE pool[c] => 0.0005 (10Z)-heptadecenoic acid[c] + 0.0005 (11Z,14Z)-eicosadienoic acid[c] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0005 (13Z)-eicosenoic acid[c] + 0.0005 (13Z)-octadecenoic acid[c] + 0.0005 (13Z,16Z)-docosadienoic acid[c] + 0.0052 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0005 (6Z,9Z)-octadecadienoic acid[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0005 (7Z)-octadecenoic acid[c] + 0.0005 (7Z)-tetradecenoic acid[c] + 0.0005 (9E)-tetradecenoic acid[c] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0005 10,13,16,19-docosatetraenoic acid[c] + 0.0005 10,13,16-docosatriynoic acid[c] + 0.0005 12,15,18,21-tetracosatetraenoic acid[c] + 0.0005 13,16,19-docosatrienoic acid[c] + 0.0005 7-palmitoleic acid[c] + 0.0005 8,11-eicosadienoic acid[c] + 0.0005 9-eicosenoic acid[c] + 0.0005 9-heptadecylenic acid[c] + 0.0007 adrenic acid[c] + 0.2125 arachidonate[c] + 0.0005 behenic acid[c] + 0.0005 cerotic acid[c] + 0.0005 cis-cetoleic acid[c] + 0.0005 cis-erucic acid[c] + 0.0005 cis-gondoic acid[c] + 0.0285 cis-vaccenic acid[c] + 0.0459 DHA[c] + 0.0411 dihomo-gamma-linolenate[c] + 0.0067 DPA[c] + 0.0005 eicosanoate[c] + 0.0005 elaidate[c] + 0.0221 EPA[c] + 0.0011 gamma-linolenate[c] + 0.0005 henicosanoic acid[c] + 0.0005 lauric acid[c] + 0.0005 lignocerate[c] + 0.1312 linoleate[c] + 0.0166 linolenate[c] + 0.0005 margaric acid[c] + 0.0005 mead acid[c] + 0.0319 myristic acid[c] + 0.0005 nervonic acid[c] + 0.0005 nonadecylic acid[c] + 0.0619 oleate[c] + 0.0163 omega-3-arachidonic acid[c] + 0.1243 palmitate[c] + 0.0327 palmitolate[c] + 0.0005 pentadecylic acid[c] + 0.0005 physeteric acid[c] + 0.1983 stearate[c] + 0.0025 stearidonic acid[c] + 0.0005 tricosanoic acid[c] + 0.0005 tridecylic acid[c] + 0.0005 ximenic acid[c]
% HMR_0690: fatty acid-LD-PS pool[c] => 0.0005 (10Z)-heptadecenoic acid[c] + 0.0005 (11Z,14Z)-eicosadienoic acid[c] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0005 (13Z)-eicosenoic acid[c] + 0.0005 (13Z)-octadecenoic acid[c] + 0.0005 (13Z,16Z)-docosadienoic acid[c] + 0.0055 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0005 (6Z,9Z)-octadecadienoic acid[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0005 (7Z)-octadecenoic acid[c] + 0.0005 (7Z)-tetradecenoic acid[c] + 0.0005 (9E)-tetradecenoic acid[c] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0005 10,13,16,19-docosatetraenoic acid[c] + 0.0005 10,13,16-docosatriynoic acid[c] + 0.0005 12,15,18,21-tetracosatetraenoic acid[c] + 0.0005 13,16,19-docosatrienoic acid[c] + 0.0005 7-palmitoleic acid[c] + 0.0005 8,11-eicosadienoic acid[c] + 0.0005 9-eicosenoic acid[c] + 0.0005 9-heptadecylenic acid[c] + 0.0043 adrenic acid[c] + 0.2005 arachidonate[c] + 0.0005 behenic acid[c] + 0.0005 cerotic acid[c] + 0.0005 cis-cetoleic acid[c] + 0.0005 cis-erucic acid[c] + 0.0005 cis-gondoic acid[c] + 0.0131 cis-vaccenic acid[c] + 0.0931 DHA[c] + 0.0089 dihomo-gamma-linolenate[c] + 0.0329 DPA[c] + 0.0005 eicosanoate[c] + 0.0005 elaidate[c] + 0.0286 EPA[c] + 0.0011 gamma-linolenate[c] + 0.0005 henicosanoic acid[c] + 0.0005 lauric acid[c] + 0.0005 lignocerate[c] + 0.0177 linoleate[c] + 0.0035 linolenate[c] + 0.0005 margaric acid[c] + 0.0005 mead acid[c] + 0.0055 myristic acid[c] + 0.0005 nervonic acid[c] + 0.0005 nonadecylic acid[c] + 0.0323 oleate[c] + 0.0174 omega-3-arachidonic acid[c] + 0.0337 palmitate[c] + 0.0058 palmitolate[c] + 0.0005 pentadecylic acid[c] + 0.0005 physeteric acid[c] + 0.4731 stearate[c] + 0.0025 stearidonic acid[c] + 0.0005 tricosanoic acid[c] + 0.0005 tridecylic acid[c] + 0.0005 ximenic acid[c]
%
% involve the breakdown of a pool into its components. These particular
% reactions do not appear to be fully or clearly balanced, and are able to
% create mass when coupled with other reactions in the model. This problem
% persists even when all Recon3D-derived reactions are constrained,
% suggesting that it did not result from merging the models. These
% reactions should therefore be constrained until they can be properly
% re-formulated in such a way as to prevent mass imbalances.
%
% NOTE: Other, similar, pool reactions are also present in the model:
%
%   HMR_0012: fatty acid-uptake pool[s] => 0.0001 (10Z)-heptadecenoic acid[s] + 0.0001 (11Z,14Z)-eicosadienoic acid[s] + 0.0001 (11Z,14Z,17Z)-eicosatrienoic acid[s] + 0.0001 (13Z)-eicosenoic acid[s] + 0.0001 (13Z)-octadecenoic acid[s] + 0.0001 (13Z,16Z)-docosadienoic acid[s] + 0.0001 (4Z,7Z,10Z,13Z,16Z)-DPA[s] + 0.0001 (6Z,9Z)-octadecadienoic acid[s] + 0.0001 (6Z,9Z,12Z,15Z,18Z)-TPA[s] + 0.0001 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[s] + 0.0001 (7Z)-octadecenoic acid[s] + 0.0001 (7Z)-tetradecenoic acid[s] + 0.0001 (9E)-tetradecenoic acid[s] + 0.0001 (9Z,12Z,15Z,18Z)-TTA[s] + 0.0001 (9Z,12Z,15Z,18Z,21Z)-TPA[s] + 0.0001 10,13,16,19-docosatetraenoic acid[s] + 0.0001 10,13,16-docosatriynoic acid[s] + 0.0001 12,15,18,21-tetracosatetraenoic acid[s] + 0.0001 13,16,19-docosatrienoic acid[s] + 0.0001 7-palmitoleic acid[s] + 0.0001 8,11-eicosadienoic acid[s] + 0.0001 9-eicosenoic acid[s] + 0.0001 9-heptadecylenic acid[s] + 0.0001 adrenic acid[s] + 0.0082 arachidonate[s] + 0.0001 behenic acid[s] + 0.0001 cerotic acid[s] + 0.0001 cis-cetoleic acid[s] + 0.0001 cis-erucic acid[s] + 0.0001 cis-gondoic acid[s] + 0.0001 cis-vaccenic acid[s] + 0.0041 DHA[s] + 0.002 dihomo-gamma-linolenate[s] + 0.0001 DPA[s] + 0.0001 eicosanoate[s] + 0.0001 elaidate[s] + 0.001 EPA[s] + 0.0001 gamma-linolenate[s] + 0.0001 henicosanoic acid[s] + 0.0001 lauric acid[s] + 0.0001 lignocerate[s] + 0.1535 linoleate[s] + 0.0092 linolenate[s] + 0.0001 margaric acid[s] + 0.0001 mead acid[s] + 0.0338 myristic acid[s] + 0.0001 nervonic acid[s] + 0.0001 nonadecylic acid[s] + 0.3837 oleate[s] + 0.0001 omega-3-arachidonic acid[s] + 0.3015 palmitate[s] + 0.0522 palmitolate[s] + 0.0001 pentadecylic acid[s] + 0.0001 physeteric acid[s] + 0.046 stearate[s] + 0.0001 stearidonic acid[s] + 0.0001 tricosanoic acid[s] + 0.0001 tridecylic acid[s] + 0.0001 ximenic acid[s]
%   HMR_3537: cholesterol-ester pool[l] => 0.0001 cholesterol-ester-10,13,16,19-docosa[l] + 0.0001 cholesterol-ester-10,13,16-docosa[l] + 0.0001 cholesterol-ester-10-hepta[l] + 0.0001 cholesterol-ester-11,14,17-eico[l] + 0.0001 cholesterol-ester-11,14-eicosa[l] + 0.0001 cholesterol-ester-11-docose[l] + 0.0001 cholesterol-ester-11-eico[l] + 0.0001 cholesterol-ester-12,15,18,21-tetracosa[l] + 0.0001 cholesterol-ester-13,16,19-doco[l] + 0.0001 cholesterol-ester-13,16-docosa[l] + 0.0001 cholesterol-ester-13-docose[l] + 0.0001 cholesterol-ester-13-eicose[l] + 0.0001 cholesterol-ester-13-octade[l] + 0.0001 cholesterol-ester-15-tetra[l] + 0.0041 cholesterol-ester-4,7,10,13,16,19-doco[l] + 0.0001 cholesterol-ester-4,7,10,13,16-docosa[l] + 0.0071 cholesterol-ester-5,8,11,14,17-eico[l] + 0.0001 cholesterol-ester-5,8,11-eico[l] + 0.0001 cholesterol-ester-5-tetra[l] + 0.0001 cholesterol-ester-6,9,12,15,18,21-tetra[l] + 0.0001 cholesterol-ester-6,9,12,15,18-tetraco[l] + 0.0001 cholesterol-ester-6,9,12,15-octa[l] + 0.0001 cholesterol-ester-6,9-octa[l] + 0.0001 cholesterol-ester-7,10,13,16,19-docosa[l] + 0.0001 cholesterol-ester-7,10,13,16-docosa[l] + 0.0406 cholesterol-ester-7-hexa[l] + 0.0001 cholesterol-ester-7-octade[l] + 0.0001 cholesterol-ester-7-tetrade[l] + 0.0001 cholesterol-ester-8,11,14,17-eico[l] + 0.0001 cholesterol-ester-8,11-eico[l] + 0.0001 cholesterol-ester-9,12,15,18,21-tetra[l] + 0.0001 cholesterol-ester-9,12,15,18-tetraco[l] + 0.0001 cholesterol-ester-9-eicose[l] + 0.0001 cholesterol-ester-9-heptade[l] + 0.0001 cholesterol-ester-9-octa[l] + 0.0001 cholesterol-ester-9-tetrade[l] + 0.0518 cholesterol-ester-arach[l] + 0.0001 cholesterol-ester-cis-vac[l] + 0.005 cholesterol-ester-dihomo-gamma[l] + 0.0001 cholesterol-ester-docosa[l] + 0.0001 cholesterol-ester-eico[l] + 0.0001 cholesterol-ester-gamma-lin[l] + 0.0001 cholesterol-ester-heneico[l] + 0.0001 cholesterol-ester-hepta[l] + 0.0001 cholesterol-ester-hexacosa[l] + 0.0001 cholesterol-ester-hexecose[l] + 0.0001 cholesterol-ester-laur[l] + 0.5254 cholesterol-ester-lin[l] + 0.0061 cholesterol-ester-linolen[l] + 0.0081 cholesterol-ester-myrist[l] + 0.0001 cholesterol-ester-nanode[l] + 0.1958 cholesterol-ester-ol[l] + 0.138 cholesterol-ester-palm[l] + 0.0001 cholesterol-ester-palmn[l] + 0.0001 cholesterol-ester-penta[l] + 0.0132 cholesterol-ester-stea[l] + 0.0001 cholesterol-ester-tetraco[l] + 0.0001 cholesterol-ester-trico[l] + 0.0001 cholesterol-ester-tridec[l]
%   HMR_3622: cholesterol-ester pool[r] <=> 0.0001 cholesterol-ester-10,13,16,19-docosa[r] + 0.0001 cholesterol-ester-10,13,16-docosa[r] + 0.0001 cholesterol-ester-10-hepta[r] + 0.0001 cholesterol-ester-11,14,17-eico[r] + 0.0001 cholesterol-ester-11,14-eicosa[r] + 0.0001 cholesterol-ester-11-docose[r] + 0.0001 cholesterol-ester-11-eico[r] + 0.0001 cholesterol-ester-12,15,18,21-tetracosa[r] + 0.0001 cholesterol-ester-13,16,19-doco[r] + 0.0001 cholesterol-ester-13,16-docosa[r] + 0.0001 cholesterol-ester-13-docose[r] + 0.0001 cholesterol-ester-13-eicose[r] + 0.0001 cholesterol-ester-13-octade[r] + 0.0001 cholesterol-ester-15-tetra[r] + 0.0041 cholesterol-ester-4,7,10,13,16,19-doco[r] + 0.0001 cholesterol-ester-4,7,10,13,16-docosa[r] + 0.0071 cholesterol-ester-5,8,11,14,17-eico[r] + 0.0001 cholesterol-ester-5,8,11-eico[r] + 0.0001 cholesterol-ester-5-tetra[r] + 0.0001 cholesterol-ester-6,9,12,15,18,21-tetra[r] + 0.0001 cholesterol-ester-6,9,12,15,18-tetraco[r] + 0.0001 cholesterol-ester-6,9,12,15-octa[r] + 0.0001 cholesterol-ester-6,9-octa[r] + 0.0001 cholesterol-ester-7,10,13,16,19-docosa[r] + 0.0001 cholesterol-ester-7,10,13,16-docosa[r] + 0.0406 cholesterol-ester-7-hexa[r] + 0.0001 cholesterol-ester-7-octade[r] + 0.0001 cholesterol-ester-7-tetrade[r] + 0.0001 cholesterol-ester-8,11,14,17-eico[r] + 0.0001 cholesterol-ester-8,11-eico[r] + 0.0001 cholesterol-ester-9,12,15,18,21-tetra[r] + 0.0001 cholesterol-ester-9,12,15,18-tetraco[r] + 0.0001 cholesterol-ester-9-eicose[r] + 0.0001 cholesterol-ester-9-heptade[r] + 0.0001 cholesterol-ester-9-octa[r] + 0.0001 cholesterol-ester-9-tetrade[r] + 0.0518 cholesterol-ester-arach[r] + 0.0001 cholesterol-ester-cis-vac[r] + 0.005 cholesterol-ester-dihomo-gamma[r] + 0.0001 cholesterol-ester-docosa[r] + 0.0001 cholesterol-ester-eico[r] + 0.0001 cholesterol-ester-gamma-lin[r] + 0.0001 cholesterol-ester-heneico[r] + 0.0001 cholesterol-ester-hepta[r] + 0.0001 cholesterol-ester-hexacosa[r] + 0.0001 cholesterol-ester-hexecose[r] + 0.0001 cholesterol-ester-laur[r] + 0.5254 cholesterol-ester-lin[r] + 0.0061 cholesterol-ester-linolen[r] + 0.0081 cholesterol-ester-myrist[r] + 0.0001 cholesterol-ester-nanode[r] + 0.1958 cholesterol-ester-ol[r] + 0.138 cholesterol-ester-palm[r] + 0.0001 cholesterol-ester-palmn[r] + 0.0001 cholesterol-ester-penta[r] + 0.0132 cholesterol-ester-stea[r] + 0.0001 cholesterol-ester-tetraco[r] + 0.0001 cholesterol-ester-trico[r] + 0.0001 cholesterol-ester-tridec[r]
%   HMR_0685: fatty acid-LD-TG1 pool[c] => 0.0001 (10Z)-heptadecenoic acid[c] + 0.0001 (11Z,14Z)-eicosadienoic acid[c] + 0.0001 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0001 (13Z)-eicosenoic acid[c] + 0.0001 (13Z)-octadecenoic acid[c] + 0.0001 (13Z,16Z)-docosadienoic acid[c] + 0.0001 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0001 (6Z,9Z)-octadecadienoic acid[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0001 (7Z)-octadecenoic acid[c] + 0.0001 (7Z)-tetradecenoic acid[c] + 0.0001 (9E)-tetradecenoic acid[c] + 0.0001 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0001 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0001 10,13,16,19-docosatetraenoic acid[c] + 0.0001 10,13,16-docosatriynoic acid[c] + 0.0001 12,15,18,21-tetracosatetraenoic acid[c] + 0.0001 13,16,19-docosatrienoic acid[c] + 0.0001 7-palmitoleic acid[c] + 0.0001 8,11-eicosadienoic acid[c] + 0.0001 9-eicosenoic acid[c] + 0.0001 9-heptadecylenic acid[c] + 0.0001 adrenic acid[c] + 0.0024 arachidonate[c] + 0.0001 behenic acid[c] + 0.0001 cerotic acid[c] + 0.0001 cis-cetoleic acid[c] + 0.0001 cis-erucic acid[c] + 0.0001 cis-gondoic acid[c] + 0.0208 cis-vaccenic acid[c] + 0.0005 DHA[c] + 0.0013 dihomo-gamma-linolenate[c] + 0.0017 DPA[c] + 0.0001 eicosanoate[c] + 0.0001 elaidate[c] + 0.0013 EPA[c] + 0.0048 gamma-linolenate[c] + 0.0001 henicosanoic acid[c] + 0.0001 lauric acid[c] + 0.0001 lignocerate[c] + 0.0741 linoleate[c] + 0.0037 linolenate[c] + 0.0001 margaric acid[c] + 0.0001 mead acid[c] + 0.0096 myristic acid[c] + 0.0001 nervonic acid[c] + 0.0001 nonadecylic acid[c] + 0.122 oleate[c] + 0.0001 omega-3-arachidonic acid[c] + 0.64 palmitate[c] + 0.0286 palmitolate[c] + 0.0001 pentadecylic acid[c] + 0.0001 physeteric acid[c] + 0.082 stearate[c] + 0.0028 stearidonic acid[c] + 0.0001 tricosanoic acid[c] + 0.0001 tridecylic acid[c] + 0.0001 ximenic acid[c]
%   HMR_0686: fatty acid-LD-TG2 pool[c] => 0.0001 (10Z)-heptadecenoic acid[c] + 0.0001 (11Z,14Z)-eicosadienoic acid[c] + 0.0001 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0001 (13Z)-eicosenoic acid[c] + 0.0001 (13Z)-octadecenoic acid[c] + 0.0001 (13Z,16Z)-docosadienoic acid[c] + 0.0016 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0001 (6Z,9Z)-octadecadienoic acid[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0001 (7Z)-octadecenoic acid[c] + 0.0001 (7Z)-tetradecenoic acid[c] + 0.0001 (9E)-tetradecenoic acid[c] + 0.0001 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0001 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0001 10,13,16,19-docosatetraenoic acid[c] + 0.0001 10,13,16-docosatriynoic acid[c] + 0.0001 12,15,18,21-tetracosatetraenoic acid[c] + 0.0001 13,16,19-docosatrienoic acid[c] + 0.0001 7-palmitoleic acid[c] + 0.0001 8,11-eicosadienoic acid[c] + 0.0001 9-eicosenoic acid[c] + 0.0001 9-heptadecylenic acid[c] + 0.0032 adrenic acid[c] + 0.0521 arachidonate[c] + 0.0001 behenic acid[c] + 0.0001 cerotic acid[c] + 0.0001 cis-cetoleic acid[c] + 0.0001 cis-erucic acid[c] + 0.0001 cis-gondoic acid[c] + 0.0217 cis-vaccenic acid[c] + 0.0081 DHA[c] + 0.0023 dihomo-gamma-linolenate[c] + 0.0018 DPA[c] + 0.0001 eicosanoate[c] + 0.0001 elaidate[c] + 0.0003 EPA[c] + 0.003 gamma-linolenate[c] + 0.0001 henicosanoic acid[c] + 0.0001 lauric acid[c] + 0.0001 lignocerate[c] + 0.281 linoleate[c] + 0.0119 linolenate[c] + 0.0001 margaric acid[c] + 0.0001 mead acid[c] + 0.0024 myristic acid[c] + 0.0001 nervonic acid[c] + 0.0001 nonadecylic acid[c] + 0.4226 oleate[c] + 0.0001 omega-3-arachidonic acid[c] + 0.12 palmitate[c] + 0.0229 palmitolate[c] + 0.0001 pentadecylic acid[c] + 0.0001 physeteric acid[c] + 0.0384 stearate[c] + 0.0025 stearidonic acid[c] + 0.0001 tricosanoic acid[c] + 0.0001 tridecylic acid[c] + 0.0001 ximenic acid[c]
%   HMR_0687: fatty acid-LD-TG3 pool[c] => 0.0001 (10Z)-heptadecenoic acid[c] + 0.0001 (11Z,14Z)-eicosadienoic acid[c] + 0.0001 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0001 (13Z)-eicosenoic acid[c] + 0.0001 (13Z)-octadecenoic acid[c] + 0.0001 (13Z,16Z)-docosadienoic acid[c] + 0.0011 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0001 (6Z,9Z)-octadecadienoic acid[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0001 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0001 (7Z)-octadecenoic acid[c] + 0.0001 (7Z)-tetradecenoic acid[c] + 0.0001 (9E)-tetradecenoic acid[c] + 0.0001 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0001 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0001 10,13,16,19-docosatetraenoic acid[c] + 0.0001 10,13,16-docosatriynoic acid[c] + 0.0001 12,15,18,21-tetracosatetraenoic acid[c] + 0.0001 13,16,19-docosatrienoic acid[c] + 0.0001 7-palmitoleic acid[c] + 0.0001 8,11-eicosadienoic acid[c] + 0.0001 9-eicosenoic acid[c] + 0.0001 9-heptadecylenic acid[c] + 0.0022 adrenic acid[c] + 0.0264 arachidonate[c] + 0.0001 behenic acid[c] + 0.0001 cerotic acid[c] + 0.0001 cis-cetoleic acid[c] + 0.0001 cis-erucic acid[c] + 0.0001 cis-gondoic acid[c] + 0.0194 cis-vaccenic acid[c] + 0.0221 DHA[c] + 0.0141 dihomo-gamma-linolenate[c] + 0.0038 DPA[c] + 0.0001 eicosanoate[c] + 0.0001 elaidate[c] + 0.0035 EPA[c] + 0.0081 gamma-linolenate[c] + 0.0001 henicosanoic acid[c] + 0.0001 lauric acid[c] + 0.0001 lignocerate[c] + 0.278 linoleate[c] + 0.0023 linolenate[c] + 0.0001 margaric acid[c] + 0.0001 mead acid[c] + 0.003 myristic acid[c] + 0.0001 nervonic acid[c] + 0.0001 nonadecylic acid[c] + 0.4519 oleate[c] + 0.0001 omega-3-arachidonic acid[c] + 0.0849 palmitate[c] + 0.014 palmitolate[c] + 0.0001 pentadecylic acid[c] + 0.0001 physeteric acid[c] + 0.0584 stearate[c] + 0.0026 stearidonic acid[c] + 0.0001 tricosanoic acid[c] + 0.0001 tridecylic acid[c] + 0.0001 ximenic acid[c]
%   HMR_0688: fatty acid-LD-PC pool[c] <=> 0.0005 (10Z)-heptadecenoic acid[c] + 0.0005 (11Z,14Z)-eicosadienoic acid[c] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0005 (13Z)-eicosenoic acid[c] + 0.0005 (13Z)-octadecenoic acid[c] + 0.0005 (13Z,16Z)-docosadienoic acid[c] + 0.0039 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0005 (6Z,9Z)-octadecadienoic acid[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0005 (7Z)-octadecenoic acid[c] + 0.0005 (7Z)-tetradecenoic acid[c] + 0.0005 (9E)-tetradecenoic acid[c] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0005 10,13,16,19-docosatetraenoic acid[c] + 0.0005 10,13,16-docosatriynoic acid[c] + 0.0005 12,15,18,21-tetracosatetraenoic acid[c] + 0.0005 13,16,19-docosatrienoic acid[c] + 0.0005 7-palmitoleic acid[c] + 0.0005 8,11-eicosadienoic acid[c] + 0.0005 9-eicosenoic acid[c] + 0.0005 9-heptadecylenic acid[c] + 0.0011 adrenic acid[c] + 0.0709 arachidonate[c] + 0.0005 behenic acid[c] + 0.0005 cerotic acid[c] + 0.0005 cis-cetoleic acid[c] + 0.0005 cis-erucic acid[c] + 0.0005 cis-gondoic acid[c] + 0.0271 cis-vaccenic acid[c] + 0.0252 DHA[c] + 0.0228 dihomo-gamma-linolenate[c] + 0.0056 DPA[c] + 0.0005 eicosanoate[c] + 0.0005 elaidate[c] + 0.0109 EPA[c] + 0.0036 gamma-linolenate[c] + 0.0005 henicosanoic acid[c] + 0.0005 lauric acid[c] + 0.0005 lignocerate[c] + 0.2479 linoleate[c] + 0.0047 linolenate[c] + 0.0005 margaric acid[c] + 0.0005 mead acid[c] + 0.0054 myristic acid[c] + 0.0005 nervonic acid[c] + 0.0005 nonadecylic acid[c] + 0.1143 oleate[c] + 0.0178 omega-3-arachidonic acid[c] + 0.2781 palmitate[c] + 0.0105 palmitolate[c] + 0.0005 pentadecylic acid[c] + 0.0005 physeteric acid[c] + 0.1271 stearate[c] + 0.0026 stearidonic acid[c] + 0.0005 tricosanoic acid[c] + 0.0005 tridecylic acid[c] + 0.0005 ximenic acid[c]
%   HMR_0691: fatty acid-LD-PI pool[c] => 0.0005 (10Z)-heptadecenoic acid[c] + 0.0005 (11Z,14Z)-eicosadienoic acid[c] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.0005 (13Z)-eicosenoic acid[c] + 0.0005 (13Z)-octadecenoic acid[c] + 0.0005 (13Z,16Z)-docosadienoic acid[c] + 0.0124 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.0005 (6Z,9Z)-octadecadienoic acid[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.0005 (7Z)-octadecenoic acid[c] + 0.0005 (7Z)-tetradecenoic acid[c] + 0.0005 (9E)-tetradecenoic acid[c] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[c] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.0005 10,13,16,19-docosatetraenoic acid[c] + 0.0005 10,13,16-docosatriynoic acid[c] + 0.0005 12,15,18,21-tetracosatetraenoic acid[c] + 0.0005 13,16,19-docosatrienoic acid[c] + 0.0005 7-palmitoleic acid[c] + 0.0005 8,11-eicosadienoic acid[c] + 0.0005 9-eicosenoic acid[c] + 0.0005 9-heptadecylenic acid[c] + 0.0031 adrenic acid[c] + 0.2118 arachidonate[c] + 0.0005 behenic acid[c] + 0.0005 cerotic acid[c] + 0.0005 cis-cetoleic acid[c] + 0.0005 cis-erucic acid[c] + 0.0005 cis-gondoic acid[c] + 0.0266 cis-vaccenic acid[c] + 0.0269 DHA[c] + 0.0228 dihomo-gamma-linolenate[c] + 0.0101 DPA[c] + 0.0005 eicosanoate[c] + 0.0005 elaidate[c] + 0.0043 EPA[c] + 0.0012 gamma-linolenate[c] + 0.0005 henicosanoic acid[c] + 0.0005 lauric acid[c] + 0.0005 lignocerate[c] + 0.0678 linoleate[c] + 0.0019 linolenate[c] + 0.0005 margaric acid[c] + 0.0005 mead acid[c] + 0.0056 myristic acid[c] + 0.0005 nervonic acid[c] + 0.0005 nonadecylic acid[c] + 0.136 oleate[c] + 0.0179 omega-3-arachidonic acid[c] + 0.0678 palmitate[c] + 0.0053 palmitolate[c] + 0.0005 pentadecylic acid[c] + 0.0005 physeteric acid[c] + 0.3553 stearate[c] + 0.0027 stearidonic acid[c] + 0.0005 tricosanoic acid[c] + 0.0005 tridecylic acid[c] + 0.0005 ximenic acid[c]
%   HMR_0692: fatty acid-LD-SM pool[c] <=> 0.002 (10Z)-heptadecenoic acid[c] + 0.002 (11Z,14Z)-eicosadienoic acid[c] + 0.002 (11Z,14Z,17Z)-eicosatrienoic acid[c] + 0.002 (13Z)-eicosenoic acid[c] + 0.002 (13Z)-octadecenoic acid[c] + 0.002 (13Z,16Z)-docosadienoic acid[c] + 0.002 (4Z,7Z,10Z,13Z,16Z)-DPA[c] + 0.002 (6Z,9Z)-octadecadienoic acid[c] + 0.002 (6Z,9Z,12Z,15Z,18Z)-TPA[c] + 0.002 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[c] + 0.002 (7Z)-octadecenoic acid[c] + 0.002 (7Z)-tetradecenoic acid[c] + 0.002 (9E)-tetradecenoic acid[c] + 0.002 (9Z,12Z,15Z,18Z)-TTA[c] + 0.002 (9Z,12Z,15Z,18Z,21Z)-TPA[c] + 0.002 10,13,16,19-docosatetraenoic acid[c] + 0.002 10,13,16-docosatriynoic acid[c] + 0.002 12,15,18,21-tetracosatetraenoic acid[c] + 0.002 13,16,19-docosatrienoic acid[c] + 0.002 7-palmitoleic acid[c] + 0.002 8,11-eicosadienoic acid[c] + 0.002 9-eicosenoic acid[c] + 0.002 9-heptadecylenic acid[c] + 0.002 adrenic acid[c] + 0.0136 arachidonate[c] + 0.002 behenic acid[c] + 0.002 cerotic acid[c] + 0.002 cis-cetoleic acid[c] + 0.002 cis-erucic acid[c] + 0.002 cis-gondoic acid[c] + 0.0453 cis-vaccenic acid[c] + 0.003 DHA[c] + 0.002 dihomo-gamma-linolenate[c] + 0.002 DPA[c] + 0.002 eicosanoate[c] + 0.002 elaidate[c] + 0.0068 EPA[c] + 0.002 gamma-linolenate[c] + 0.002 henicosanoic acid[c] + 0.002 lauric acid[c] + 0.002 lignocerate[c] + 0.025 linoleate[c] + 0.0075 linolenate[c] + 0.002 margaric acid[c] + 0.002 mead acid[c] + 0.011 myristic acid[c] + 0.002 nervonic acid[c] + 0.002 nonadecylic acid[c] + 0.061 oleate[c] + 0.002 omega-3-arachidonic acid[c] + 0.557 palmitate[c] + 0.0358 palmitolate[c] + 0.002 pentadecylic acid[c] + 0.002 physeteric acid[c] + 0.138 stearate[c] + 0.002 stearidonic acid[c] + 0.002 tricosanoic acid[c] + 0.002 tridecylic acid[c] + 0.002 ximenic acid[c]
%   HMR_5257: fatty acid-LD-PC pool[r] => 0.0005 (10Z)-heptadecenoic acid[r] + 0.0005 (11Z,14Z)-eicosadienoic acid[r] + 0.0005 (11Z,14Z,17Z)-eicosatrienoic acid[r] + 0.0005 (13Z)-eicosenoic acid[r] + 0.0005 (13Z)-octadecenoic acid[r] + 0.0005 (13Z,16Z)-docosadienoic acid[r] + 0.0039 (4Z,7Z,10Z,13Z,16Z)-DPA[r] + 0.0005 (6Z,9Z)-octadecadienoic acid[r] + 0.0005 (6Z,9Z,12Z,15Z,18Z)-TPA[r] + 0.0005 (6Z,9Z,12Z,15Z,18Z,21Z)-THA[r] + 0.0005 (7Z)-octadecenoic acid[r] + 0.0005 (7Z)-tetradecenoic acid[r] + 0.0005 (9E)-tetradecenoic acid[r] + 0.0005 (9Z,12Z,15Z,18Z)-TTA[r] + 0.0005 (9Z,12Z,15Z,18Z,21Z)-TPA[r] + 0.0005 10,13,16,19-docosatetraenoic acid[r] + 0.0005 10,13,16-docosatriynoic acid[r] + 0.0005 12,15,18,21-tetracosatetraenoic acid[r] + 0.0005 13,16,19-docosatrienoic acid[r] + 0.0005 7-palmitoleic acid[r] + 0.0005 8,11-eicosadienoic acid[r] + 0.0005 9-eicosenoic acid[r] + 0.0005 9-heptadecylenic acid[r] + 0.0011 adrenic acid[r] + 0.0709 arachidonate[r] + 0.0005 behenic acid[r] + 0.0005 cerotic acid[r] + 0.0005 cis-cetoleic acid[r] + 0.0005 cis-erucic acid[r] + 0.0005 cis-gondoic acid[r] + 0.0271 cis-vaccenic acid[r] + 0.0252 DHA[r] + 0.0228 dihomo-gamma-linolenate[r] + 0.0056 DPA[r] + 0.0005 eicosanoate[r] + 0.0005 elaidate[r] + 0.0109 EPA[r] + 0.0036 gamma-linolenate[r] + 0.0005 henicosanoic acid[r] + 0.0005 lauric acid[r] + 0.0005 lignocerate[r] + 0.2479 linoleate[r] + 0.0047 linolenate[r] + 0.0005 margaric acid[r] + 0.0005 mead acid[r] + 0.0054 myristic acid[r] + 0.0005 nervonic acid[r] + 0.0005 nonadecylic acid[r] + 0.1143 oleate[r] + 0.0178 omega-3-arachidonic acid[r] + 0.2781 palmitate[r] + 0.0105 palmitolate[r] + 0.0005 pentadecylic acid[r] + 0.0005 physeteric acid[r] + 0.1271 stearate[r] + 0.0026 stearidonic acid[r] + 0.0005 tricosanoic acid[r] + 0.0005 tridecylic acid[r] + 0.0005 ximenic acid[r]
%
% However, analyses did not indicate that the presence of these reactions
% led to mass/energy imbalances, so they will not be constrained or
% modified at this time.
rxns = {'HMR_0689';'HMR_0690'};
del_rxns = [del_rxns; rxns];
rxnNotes = [rxnNotes; [rxns, repmat({'mass balance analyses show that this reaction contributes to a flux solution that can generate mass; the rxn should therefore be constrained until imbalances can be addressed'},length(rxns),1)]];


%% Group2 reactions:
% These reactions have the same reactants (or products) but differ in products
% (or reactants), which results in imbalanced mass

% RXN ID    RXN EQUATION
%
% AGPAT1    2 H+[c] + R Total 2 Coenzyme A[c] + Lysophosphatidic Acid[c] => CoA[c] + phosphatidate-LD-TAG pool[c]
% AGPAT2    palmitoyl-CoA[c] + Lysophosphatidic Acid[c] => CoA[c] + phosphatidate-LD-TAG pool[c]
% AGPAT3    oleoyl-CoA[c] + Lysophosphatidic Acid[c] => CoA[c] + phosphatidate-LD-TAG pool[c]
% AGPAT4    linoleoyl-CoA[c] + Lysophosphatidic Acid[c] => CoA[c] + phosphatidate-LD-TAG pool[c]
%  
% LCAT55e   cholesterol[s] + PI pool[s] => cholesterol-ester pool[s] + 1-Palmitoylglycerophosphoinositol[s]
% LCAT17e   cholesterol[s] + PI pool[s] => cholesterol-ester pool[s] + 1-Stearoylglycerophosphoinositol[s]
% LCAT4e    cholesterol[s] + PI pool[s] => cholesterol-ester pool[s] + 1-Arachidonoylglycerophosphoinositol[s]
%  
% LCAT12e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Oleoylglycerophosphoethanolamine (Delta 9)[s]
% LCAT54e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Palmitoylglycerophosphoethanolamine[s]
% LCAT16e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Stearoylglycerophosphoethanolamine[s]
% LCAT19e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 2-Linoleoylglycerophosphoethanolamine[s]
% LCAT3e    cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Arachidonoyl-Sn-Glycero-3-Phosphoethanolamine[s]
% LCAT40e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosatrienoylglycerophosphoethanolamine (Delta 11, 14, 17), Lpe (20:3)[s]
% LCAT41e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Docosahexenoylglyceroethanolamine (Delta 4, 7, 10, 13, 16, 19), Lpe (22:6)[s]
% LCAT42e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Docosatetraenoyglycerophosphoethanolamine (22:4, Delta 7, 10, 13, 16)[s]
% LCAT43e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Dihomo-Linolenoylglycerophosphoethanolamine (20:3, Delta 8, 11, 14)[s]
% LCAT44e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Didecanoylglycerophosphoethanolamine (C12:0 Pe)[s]
% LCAT45e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Myristoylglycerophosphoethanolamine (C14:0 Pe)[s]
% LCAT56e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Hexadecenoylglycerophosphoethanolamine (C16:1 Pe, Delta 9)[s]
% LCAT46e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Tridecanoylglycerophosphoethanolamine (C13:0 Pe)[s]
% LCAT47e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Pentadecanoylglycerophosphoethanolamine (C15:0 Pe)[s]
% LCAT48e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Heptadecanoylglycerophosphoethanolamine (C17:0 Pe)[s]
% LCAT9e    cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Linoleoylglycerophosphoethanolamine (Delta 9,12)[s]
% 
% LCAT10e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Myristoylglycerophosphocholine[s]
% LCAT11e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Oleoylglycerophosphocholine (Delta 9)[s]
% LCAT13e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Palmitoleoylglycerophosphocholine (Delta 9)[s]
% LCAT14e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Palmitoylglycerophosphocholine[s]
% LCAT15e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Stearoylglycerophosphocholine[s]
% LCAT18e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 2-Linoleoylglycerophosphocholine[s]
% LCAT20e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 2-Oleoylglycerophosphocholine[s]
% LCAT21e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 2-Palmitoylglycerophosphocholine[s]
% LCAT22e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 2-Stearoylglycerophosphocholine[s]
% LCAT23e   cholesterol[s] + PC-LD pool[s] => 2-lysolecithin pool[s] + 1-Gamma-Linolenoyl-Cholesterol, Cholesterol-Ester (18:3, Delta 6, 9, 12)[s]
% LCAT25e   cholesterol[s] + PC-LD pool[s] => 2-lysolecithin pool[s] + 1-Vaccenoyl-Cholesterol, Cholesterol-Ester (18:1, Delta 11)[s]
% LCAT26e   cholesterol[s] + PC-LD pool[s] => 2-lysolecithin pool[s] + 1-Timnodnoyl-Cholesterol, Cholesterol-Ester (20:5, Delta 5,8,11,14,17)[s]
% LCAT27e   cholesterol[s] + PC-LD pool[s] => 2-lysolecithin pool[s] + Cholesteryl Arachidonate, Cholesterol-Ester (20:4, Delta 5,8,11,14)[s]
% LCAT28e   cholesterol[s] + PC-LD pool[s] => 2-lysolecithin pool[s] + Cholesteryl Docosahexanoate, Cholesterol-Ester (22:6, Delta 4,7,10,13,16,19)[s]
% LCAT29e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Pentadecanoylglycerophosphocholine, Sn1-Lpc (15:0)[s]
% LCAT2e    cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Arachidonoyl-Glycero-3-Phosphocholine[s]
% LCAT30e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Octadeca-Trienoylglycerophosphocholine, Sn1-Lpc (18:3, Delta 9, 12, 15)[s]
% LCAT31e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Octadeca-Trienoylglycerophosphocholine, Sn1-Lpc (18:3, Delta 6, 9, 12)[s]
% LCAT32e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Nonadecanoylglycerophosphocholine, Sn1-Lpc (19:0)[s]
% LCAT33e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosenoylglycerophosphocholine (Delta 11) ,Sn1-Lpc (20:1)[s]
% LCAT34e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosatetraenoylglycerophosphocholine (Delta 8, 11, 14, 17), Sn1-Lpc (20:4)[s]
% LCAT35e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosapentenoylglycerophosphocholine (Delta 5, 8, 11, 14, 17), Sn1-Lpc (20:5)[s]
% LCAT36e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosatetraenoylglycerophosphocholine (Delta 7, 10, 13, 16), Sn1-Lpc (22:4)[s]
% LCAT37e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosapentenoylglycerophosphocholine (Delta 7, 10, 13, 16, 19), Sn1-Lpc (22:5)-W3[s]
% LCAT38e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosapentenoylglycerophosphocholine (Delta 4, 7, 10, 13, 16), Sn1-Lpc (22:5)-W6[s]
% LCAT39e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosahexenoylglycerophosphocholine (Delta 4, 7, 10, 13, 16, 19), Sn1-Lpc (22:6)[s]
% LCAT49e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Dihomo-Linolenoylglycerophosphocholine (20:3, Delta 8, 11, 14), Lysopc A C20:3[s]
% LCAT50e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Lignocericylglycerophosphocholine (24:0), Lysopc A C24[s]
% LCAT51e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + Lysopc A C26:1 (Delta 5)[s]
% LCAT52e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + Lysopc A C28:1 (Delta 5)[s]
% LCAT53e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + Lysopc A C28:0[s]
% LCAT57e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosahexaenoylglycerophosphocholine[s]
% LCAT5e    cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosadienoylglycerophosphocholine (Delta 11,14)[s]
% LCAT6e    cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosatrienoylglycerophosphocholine (Delta 11, 14, 17)[s]
% LCAT7e    cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Heptadecanoylglycerophosphocholine[s]
% LCAT8e    cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Linoleoylglycerophosphocholine (Delta 9,12)[s]
% 
% SMS1      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/14:0), Sphingomyelin[c]
% SMS10     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/21:0), Sphingomyelin[c]
% SMS11     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/22:1), Sphingomyelin[c]
% SMS12     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/22:0), Sphingomyelin[c]
% SMS16     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/23:0), Sphingomyelin[c]
% SMS13     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:0/24:1), Sphingomyelin[c]
% SMS14     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:0/24:0), Sphingomyelin[c]
% SMS15     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:0/25:0), Sphingomyelin[c]
% SMS2      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/15:0), Sphingomyelin[c]
% SMS3      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/16:1), Sphingomyelin[c]
% SMS4      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/16:0), Sphingomyelin[c]
% SMS5      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/17:0), Sphingomyelin[c]
% SMS6      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/18:0), Sphingomyelin[c]
% SMS7      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/18:1), Sphingomyelin[c]
% SMS8      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/20:1), Sphingomyelin[c]
% SMS9      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/20:0), Sphingomyelin[c]
% 
% PEOLE_HSPLA2      H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Oleoylglycerophosphoethanolamine (Delta 9)[c]
% PEPALM_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Palmitoylglycerophosphoethanolamine[c]
% PE2LINL_HSPLA2    H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 2-Linoleoylglycerophosphoethanolamine[c]
% PEAR_HSPLA2       H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Arachidonoyl-Sn-Glycero-3-Phosphoethanolamine[c]
% PE203_HSPLA2      H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Eicosatrienoylglycerophosphoethanolamine (Delta 11, 14, 17), Lpe (20:3)[c]
% PE226_HSPLA2      H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosahexenoylglyceroethanolamine (Delta 4, 7, 10, 13, 16, 19), Lpe (22:6)[c]
% PE224_HSPLA2      H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosatetraenoyglycerophosphoethanolamine (22:4, Delta 7, 10, 13, 16)[c]
% PEDH203_HSPLA2    H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Dihomo-Linolenoylglycerophosphoethanolamine (20:3, Delta 8, 11, 14)[c]
% PEDH12_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Didecanoylglycerophosphoethanolamine (C12:0 Pe)[c]
% PEDH14_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Myristoylglycerophosphoethanolamine (C14:0 Pe)[c]
% PEDH161_HSPLA2    H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Hexadecenoylglycerophosphoethanolamine (C16:1 Pe, Delta 9)[c]
% PEDH13_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Tridecanoylglycerophosphoethanolamine (C13:0 Pe)[c]
% PEDH15_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Pentadecanoylglycerophosphoethanolamine (C15:0 Pe)[c]
% PEDH17_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Heptadecanoylglycerophosphoethanolamine (C17:0 Pe)[c]
% PELINL_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Linoleoylglycerophosphoethanolamine (Delta 9,12)[c]
%  
% PLA2_2            H2O[c] + PC-LD pool[c] => 2-lysolecithin pool[c] + H+[c] + R Total 2 Position[c]
% PLA2_2e           H2O[s] + PC-LD pool[s] => H+[s] + 2-lysolecithin pool[s] + R Total 2 Position[s]
% PCHOLMYR_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Myristoylglycerophosphocholine[c]
% PCHOLOLE_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Oleoylglycerophosphocholine (Delta 9)[c]
% PCHOLPALME_HSPLA2 H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Palmitoleoylglycerophosphocholine (Delta 9)[c]
% PCHOLPALM_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Palmitoylglycerophosphocholine[c]
% PCHOLSTE_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Stearoylglycerophosphocholine[c]
% PCHOL2LINL_HSPLA2 H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 2-Linoleoylglycerophosphocholine[c]
% PCHOL2OLE_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 2-Oleoylglycerophosphocholine[c]
% PCHOL2PALM_HSPLA2 H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 2-Palmitoylglycerophosphocholine[c]
% PCHOL2STE_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 2-Stearoylglycerophosphocholine[c]
% PCHOLN15_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Pentadecanoylglycerophosphocholine, Sn1-Lpc (15:0)[c]
% PCHOLAR_HSPLA2    H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Arachidonoyl-Glycero-3-Phosphocholine[c]
% PCHOLN183_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Octadeca-Trienoylglycerophosphocholine, Sn1-Lpc (18:3, Delta 9, 12, 15)[c]
% PCHOLN1836_HSPLA2 H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Octadeca-Trienoylglycerophosphocholine, Sn1-Lpc (18:3, Delta 6, 9, 12)[c]
% PCHOLN19_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Nonadecanoylglycerophosphocholine, Sn1-Lpc (19:0)[c]
% PCHOLN201_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Eicosenoylglycerophosphocholine (Delta 11) ,Sn1-Lpc (20:1)[c]
% PCHOLN204_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Eicosatetraenoylglycerophosphocholine (Delta 8, 11, 14, 17), Sn1-Lpc (20:4)[c]
% PCHOLN205_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Eicosapentenoylglycerophosphocholine (Delta 5, 8, 11, 14, 17), Sn1-Lpc (20:5)[c]
% PCHOLN224_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosatetraenoylglycerophosphocholine (Delta 7, 10, 13, 16), Sn1-Lpc (22:4)[c]
% PCHOLN225_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosapentenoylglycerophosphocholine (Delta 7, 10, 13, 16, 19), Sn1-Lpc (22:5)-W3[c]
% PCHOLN2254_HSPLA2 H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosapentenoylglycerophosphocholine (Delta 4, 7, 10, 13, 16), Sn1-Lpc (22:5)-W6[c]
% PCHOLN226_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosahexenoylglycerophosphocholine (Delta 4, 7, 10, 13, 16, 19), Sn1-Lpc (22:6)[c]
% PCHOLN203_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Dihomo-Linolenoylglycerophosphocholine (20:3, Delta 8, 11, 14), Lysopc A C20:3[c]
% PCHOLN24_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Lignocericylglycerophosphocholine (24:0), Lysopc A C24[c]
% PCHOLN261_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + Lysopc A C26:1 (Delta 5)[c]
% PCHOLN281_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + Lysopc A C28:1 (Delta 5)[c]
% PCHOLN28_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + Lysopc A C28:0[c]
% PCHOLDOC_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosahexaenoylglycerophosphocholine[c]
% PCHOLDEIC_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Eicosadienoylglycerophosphocholine (Delta 11,14)[c]
% PCHOLDET_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Eicosatrienoylglycerophosphocholine (Delta 11, 14, 17)[c]
% PCHOLHEP_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Heptadecanoylglycerophosphocholine[c]
% PCHOLLINL_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Linoleoylglycerophosphocholine (Delta 9,12)[c]
%  
% LPS2e         H2O[s] + 1,2-diacylglycerol-LD-TAG pool[s] => H+[s] + 1-acylglycerol-3P-LD-TG1 pool[s] + R Total[s]
% MAGLINL_HSe   H2O[s] + 1,2-diacylglycerol-LD-TAG pool[s] => H+[s] + R Total[s] + 1-Linoleoylglycerol[s]
% MAGOLE_HSe    H2O[s] + 1,2-diacylglycerol-LD-TAG pool[s] => H+[s] + R Total[s] + 1-Oleoylglycerol[s]
% LPS5e         H2O[s] + 1,2-diacylglycerol-LD-TAG pool[s] => H+[s] + R Total[s] + 1-Palmitoylglycerol[s]
% LPS6e         H2O[s] + 1,2-diacylglycerol-LD-TAG pool[s] => H+[s] + R Total[s] + 1-Stearoylglycerol[s]
% LPS7e         H2O[s] + 1,2-diacylglycerol-LD-TAG pool[s] => H+[s] + R Total[s] + 1-Arachidonoyl Glycerol[s]

% ARTFR13   myristoyl-CoA[c] => 0.875 R Group 1 Coenzyme A[c]
% ARTFR202  2 FADH2[m] + H+[m] + linolenoyl-CoA[c] + NADPH[m] => 2 FAD[m] + NADP+[m] + 1.125 R Group 2 Coenzyme A[c]
% ARTFR203  2 FADH2[m] + gamma-linolenoyl-CoA[c] + H+[m] + NADPH[m] => 2 FAD[m] + NADP+[m] + 1.125 R Group 2 Coenzyme A[c]
% ARTFR204  (6Z,9Z,12Z,15Z)-octadecatetraenoyl-CoA[c] + 2 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 2 FAD[m] + 2 NADP+[m] + 1.125 R Group 2 Coenzyme A[c]
% ARTFR205  dihomo-gamma-linolenoyl-CoA[c] + 2 FADH2[m] + H+[m] + NADPH[m] => 2 FAD[m] + NADP+[m] + 1.25 R Group 2 Coenzyme A[c]
% ARTFR206  arachidonyl-CoA[c] + 2 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 2 FAD[m] + 2 NADP+[m] + 1.25 R Group 2 Coenzyme A[c]
% ARTFR207  eicosanoyl-CoA[c] => 1.25 R Group 2 Coenzyme A[c]
% ARTFR208  (5Z,8Z,11Z,14Z,17Z)-eicosapentaenoyl-CoA[c] + 3 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 3 FAD[m] + 2 NADP+[m] + 1.25 R Group 2 Coenzyme A[c]
% ARTFR209  (7Z,10Z,13Z,16Z)-docosatetraenoyl-CoA[c] + 2 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 2 FAD[m] + 2 NADP+[m] + 1.375 R Group 2 Coenzyme A[c]
% ARTFR210  FADH2[m] + H+[m] + linoleoyl-CoA[c] + NADPH[m] => FAD[m] + NADP+[m] + 1.125 R Group 2 Coenzyme A[c]
% ARTFR211  (7Z,10Z,13Z,16Z,19Z)-docosapentaenoyl-CoA[c] + 3 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 3 FAD[m] + 2 NADP+[m] + 1.375 R Group 2 Coenzyme A[c]
% ARTFR212  (4Z,7Z,10Z,13Z,16Z)-docosapentaenoyl-CoA[c] + 3 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 3 FAD[m] + 2 NADP+[m] + 1.375 R Group 2 Coenzyme A[c]
% ARTFR213  (4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoyl-CoA[c] + 3 FADH2[m] + 3 H+[m] + 3 NADPH[m] => 3 FAD[m] + 3 NADP+[m] + 1.375 R Group 2 Coenzyme A[c]
% ARTFR31   stearoyl-CoA[c] => 1.125 R Group 3 Coenzyme A[c]
% ARTFR32   FADH2[m] + oleoyl-CoA[c] => FAD[m] + 1.125 R Group 3 Coenzyme A[c]
% ARTFR33   FADH2[m] + 11-Octadecenoyl Coenzyme A[c] => FAD[m] + 1.125 R Group 3 Coenzyme A[c]
% ARTFR34   (6Z,9Z)-octadecadienoyl-CoA[c] + 2 FADH2[m] => 2 FAD[m] + 1.125 R Group 3 Coenzyme A[c]
% ARTFR42   FADH2[m] + oleoyl-CoA[c] => FAD[m] + 1.125 R Group 4 Coenzyme A[c]
% ARTFR43   FADH2[m] + 11-Octadecenoyl Coenzyme A[c] => FAD[m] + 1.125 R Group 4 Coenzyme A[c]
% ARTFR44   (6Z,9Z)-octadecadienoyl-CoA[c] + 2 FADH2[m] => 2 FAD[m] + 1.125 R Group 4 Coenzyme A[c]
% ARTFR45   (15Z)-tetracosenoyl-CoA[c] + FADH2[m] => FAD[m] + 1.5 R Group 4 Coenzyme A[c]
% ARTFR46   (2E)-octadecenoyl-CoA[c] + FADH2[m] => FAD[m] + 1.125 R Group 4 Coenzyme A[c]
% ARTFR51   tetracosanoyl-CoA[c] => 1.5 R Group 5 Coenzyme A[c]
% ARTFR52   hexacosanoyl-CoA[c] => 1.625 R Group 5 Coenzyme A[c]
% ARTFR53   (8Z,11Z,14Z,17Z)-eicosatetraenoyl-CoA[c] + 2 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 2 FAD[m] + 2 NADP+[m] + 1.25 R Group 5 Coenzyme A[c]
% ARTFR54   (6Z,9Z,12Z,15Z,18Z)-tetracosapentaenoyl-CoA[c] + 3 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 3 FAD[m] + 2 NADP+[m] + 1.5 R Group 5 Coenzyme A[c]
% ARTFR55   (9Z,12Z,15Z,18Z,21Z)-tetracosapentaenoyl-CoA[c] + 3 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 3 FAD[m] + 2 NADP+[m] + 1.5 R Group 5 Coenzyme A[c]
% ARTFR56   (9Z,12Z,15Z,18Z)-tetracosatetraenoyl-CoA[c] + 2 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 2 FAD[m] + 2 NADP+[m] + 1.5 R Group 5 Coenzyme A[c]
% ARTFR57   (6Z,9Z,12Z,15Z,18Z,21Z)-tetracosahexaenoyl-CoA[c] + 3 FADH2[m] + 3 H+[m] + 3 NADPH[m] => 3 FAD[m] + 3 NADP+[m] + 1.5 R Group 5 Coenzyme A[c]
%
% TAG_HSad      2 H2O[c] + 2 linoleoyl-CoA[c] + 2 oleoyl-CoA[c] + palmitoyl-CoA[c] + 2 sn-glycerol-3-phosphate[c] + stearoyl-CoA[c] => 6 CoA[c] + 2 Pi[c] + 2 TAG-VLDL pool[c]
% TAG_HSad_NE   2 H2O[c] + myristoyl-CoA[c] + 2 oleoyl-CoA[c] + palmitoleoyl-CoA[c] + palmitoyl-CoA[c] + 2 sn-glycerol-3-phosphate[c] + stearoyl-CoA[c] => 6 CoA[c] + 2 Pi[c] + 2 TAG-VLDL pool[c]
% TAG_HSad_E    (4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoyl-CoA[c] + (5Z,8Z,11Z,14Z,17Z)-eicosapentaenoyl-CoA[c] + arachidonyl-CoA[c] + 2 H2O[c] + linolenoyl-CoA[c] + linoleoyl-CoA[c] + 2 sn-glycerol-3-phosphate[c] => 5 CoA[c] + 2 Pi[c] + 2 TAG-VLDL pool[c]
%
% This reaction was found to exhibit variable mass with other reactions
% involving cholesterol and cholesterol-ester pools in the model:
% CHOLESTle     cholesterol-ester pool[s] + H2O[s] => cholesterol[s] + H+[s] + R Total[s]
%

% load list of the above reactions to constrain (stored in a tsv file)
constrain_rxns = importdata('../../ComplementaryData/modelCuration/variable_mass_rxns_to_constrain.tsv');

% update the rxnNotes array
rxnNotes = [rxnNotes; [constrain_rxns, repmat({'reaction is artificial and involves dead-end artificial metabolite with no apparent purpose, and should therefore be DELETED'},length(constrain_rxns),1)]];


%% log inactivation reactions into inactivationRxns.tsv
writecell2file(rxnNotes,'../../ComplementaryData/modelCuration/inactivationRxns.tsv',true,'\t','',true);

