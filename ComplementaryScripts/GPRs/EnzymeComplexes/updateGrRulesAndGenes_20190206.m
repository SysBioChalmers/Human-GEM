%
% FILE NAME:    updateGrRulesAndGenes_20190206.m
% 
% DATE CREATED: 2019-02-06
%     MODIFIED: 2019-02-07
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: This script updates humanGEM gene-reaction rules (grRules), such
%          that only primary assembly gene IDs are used. All gene IDs that
%          are associated with a non-primary assembly (i.e., allele
%          variants), will be removed from the model. In addition, the
%          grRules for a few reactions will be corrected based on
%          information from literature/databases.



%% Load humanGEM model

load('humanGEM.mat');  % version ???


%% Remove non-primary assembly Ensembl gene IDs from humanGEM

% get list of primary assmbly genes
tmp = readtable('../../../ComplementaryData/Ensembl/ensembl_ID_mapping_20190207.txt');
ensg_ids = unique(tmp.Gene_stable_ID);  % get list of primary assembly IDs
ensg_ids(cellfun(@isempty,ensg_ids)) = [];  % remove empty ID if it exists

% obtain list of humanGEM genes, with non-primary genes removed
new_ids = ihuman.genes;
new_ids(~ismember(new_ids,ensg_ids)) = {''};

% remove non-primary genes from model
[grRules,genes,rxnGeneMat] = translateGeneRules(ihuman.grRules,[ihuman.genes,new_ids]);
ihuman.grRules = grRules;
ihuman.genes = genes;
ihuman.rxnGeneMat = rxnGeneMat;


%% Update of miscellaneous grRules

% Incorporate the (DLD and DLST and OGDH) complex into two reactions
% HMR_4239: 2-oxoadipate[m] + CoA[m] + NAD+[m] => CO2[m] + glutaryl-CoA[m] + NADH[m]
% HMR_5297: AKG[m] + CoA[m] + NAD+[m] => CO2[m] + NADH[m] + succinyl-CoA[m]
% Ref: KEGG (R01933 and R08549)
[~,rxn_ind] = ismember({'HMR_4239';'HMR_5297'},ihuman.rxns);
ihuman.grRules(rxn_ind) = {'ENSG00000091140 and ENSG00000119689 and ENSG00000105953'};


% HMGCS1 is soluble form, whereas HMGCS2 is mitochondrial
% Ref: NCBI (ID: 3158)
% HMR_4604: acetoacetyl-CoA[p] + acetyl-CoA[p] + H2O[p] => CoA[p] + H+[p] + HMG-CoA[p]
% HMR_1437: acetoacetyl-CoA[c] + acetyl-CoA[c] + H2O[c] => CoA[c] + H+[c] + HMG-CoA[c]
[~,rxn_ind] = ismember({'HMR_4604';'HMR_1437'},ihuman.rxns);
ihuman.grRules(rxn_ind) = {'ENSG00000112972'};  % non-mitochondrial form (HMGCS1)

% HMR_1573: acetoacetyl-CoA[m] + acetyl-CoA[m] + H2O[m] => CoA[m] + H+[m] + HMG-CoA[m]
[~,rxn_ind] = ismember({'HMR_1573'},ihuman.rxns);
ihuman.grRules(rxn_ind) = {'ENSG00000134240'};  % mitochondrial form (HMGCS2)


% Several tRNA synthetase reactions in the model are associated with both
% the cytoplasmic and mitochondrial version of the gene, despite the
% reaction only taking place in the cytoplasm. Therefore, the mitochondrial
% version of the gene will be removed.
[~,rxn_ind] = ismember({'HMR_5130';'HMR_5131';'HMR_5132';'HMR_5133';'HMR_5134';...
                    'HMR_5135';'HMR_5137';'HMR_5139';'HMR_5140';'HMR_5141';...
                    'HMR_5143';'HMR_5145';'HMR_5146';'HMR_5147';'HMR_5148';...
                    'HMR_5149';'HMR_5150'},ihuman.rxns);
ihuman.grRules(rxn_ind) = {'ENSG00000134684';'ENSG00000090861';'ENSG00000113643';'ENSG00000134440';'ENSG00000115866';
                           'ENSG00000110619 or ENSG00000278191';'ENSG00000136628';'ENSG00000170445';'ENSG00000196305';'ENSG00000133706';
                           'ENSG00000166986';'ENSG00000116120 and ENSG00000179115';'ENSG00000136628';'ENSG00000031698';'ENSG00000113407';
                           'ENSG00000140105';'ENSG00000096171 or ENSG00000204394 or ENSG00000224264 or ENSG00000226589 or ENSG00000227686 or ENSG00000231116 or ENSG00000231945'};

                   
% The RNA polymerase reactions:
% HMR_7161: 0.18 ATP[c] + 0.3 CTP[c] + 0.34 GTP[c] + 0.18 UTP[c] => PPi[c] + RNA[c]
% HMR_7162: 0.18 ADP[c] + 0.3 CDP[c] + 0.34 GDP[c] + 0.18 UDP[c] => Pi[c] + RNA[c]
% do not properly represent the enzyme complexes in their grRules, and are
% associated with some genes (NQO1, REG3A, RIPOR1) that are unrelated.
% Therefore, the rule will be updated to incorporate complex information,
% and to remove erroneously associated genes.
[~,rxn_ind] = ismember({'HMR_7161';'HMR_7162'},ihuman.rxns);
ihuman.grRules(rxn_ind) = {'(ENSG00000113356 or ENSG00000121851) and (ENSG00000090060 or ENSG00000115421 or ENSG00000164329 or ENSG00000218823) and (ENSG00000083223 or ENSG00000134744 or ENSG00000149016) and (ENSG00000181222 or ENSG00000284832) and (ENSG00000058600 or ENSG00000284282) and (ENSG00000066379 or ENSG00000206502 or ENSG00000224859 or ENSG00000233795 or ENSG00000235176 or ENSG00000235443 or ENSG00000236808 or ENSG00000236949) and ENSG00000005075 and ENSG00000013503 and ENSG00000047315 and ENSG00000068654 and ENSG00000099817 and ENSG00000099821 and ENSG00000100142 and ENSG00000100413 and ENSG00000102978 and ENSG00000105258 and ENSG00000107951 and ENSG00000125630 and ENSG00000132664 and ENSG00000137054 and ENSG00000144231 and ENSG00000147669 and ENSG00000148606 and ENSG00000161980 and ENSG00000163882 and ENSG00000168002 and ENSG00000168495 and ENSG00000171453 and ENSG00000177700 and ENSG00000186141 and ENSG00000186184'};


% Need to account for enzyme complex and thioredoxin dependency in grRules
% of reactions involving conversion of ribonucleotides into deoxyribonucleotides
[~,rxn_ind] = ismember({'HMR_4611';'HMR_4612';'HMR_4614';'HMR_4615';'HMR_4617';...
                        'HMR_4618';'HMR_4619';'HMR_4621';'HMR_5415';'HMR_5416';...
                        'HMR_6621';'HMR_6622';'r1431';'r1432'},ihuman.rxns);
ihuman.grRules(rxn_ind) = {'(ENSG00000100348 or ENSG00000136810) and ENSG00000048392 and ENSG00000167325 and ENSG00000171848'};
                               

%% Update other gene- and protein-related fields

% update genes and rxnGeneMat fields
[genes,rxnGeneMat] = getGenesFromGrRules(ihuman.grRules);
ihuman.genes = genes;
ihuman.rxnGeneMat = rxnGeneMat;

% update protein fields
[prRules,proteins,rxnProtMat] = translateGeneRules(ihuman.grRules,'UniProt','ENSG');
ihuman.prRules = prRules;
ihuman.proteins = proteins;
ihuman.rxnProtMat = rxnProtMat;


%% Clear intermediate variables and export model

clear ensg_ids genes grRules new_ids proteins prRules rxn_ind rxnGeneMat
clear rxnProtMat tmp

save('../../../ModelFiles/mat/humanGEM.mat','ihuman');



