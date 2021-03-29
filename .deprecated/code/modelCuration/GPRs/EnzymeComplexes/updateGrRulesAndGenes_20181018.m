%
% FILE NAME:    updateGrRulesAndGenes_20181018.m
% 
% PURPOSE: This script corrects small errors found in two of the humanGEM
%          grRules, introduced when incorporating CORUM enzyme complex
%          information into the grRules (miscModelCurationScript_20181005):
%   
%          1. For reaction HMR_9579, the associated grRule contains a
%             quotes character ("). This will be removed from the rule.
%
%          2. For reaction HMR_6921, two genes were erroneously introduced
%             into the grRule: ENSG00000198764 and ENSG00000198765, neither
%             of which are associated with that reaction. These genes also
%             do not occur anywhere else in the model. Both of these genes
%             will be removed from the associated grRule.
%
%          In addition, this script will update the .genes and .rxnGeneMat
%          fields, which will remove any genes that no longer appear in any
%          of the grRules.
%
%          Finally, the .proteins, .rxnProtMat, and .prRules fields will
%          also be updated.
%


%% Load humanGEM model

load('humanGEM.mat');  % version 0.5.0


%% Remove erroneous genes from grRules

% specify problematic genes
rem_genes = {'"';'ENSG00000198764';'ENSG00000198765'};

% identify grRule containing ("), and remove it, along with preceeding "or"
rxn_ind = contains(ihuman.grRules,rem_genes(1));
ihuman.grRules(rxn_ind) = regexprep(ihuman.grRules(rxn_ind),' or "','');

% identify grRule containing two erroneous genes, and remove them, along
% with the preceeding "and".
rxn_ind = contains(ihuman.grRules,rem_genes(2:3));
ihuman.grRules(rxn_ind) = regexprep(ihuman.grRules(rxn_ind),' and ENSG00000198764','');
ihuman.grRules(rxn_ind) = regexprep(ihuman.grRules(rxn_ind),' and ENSG00000198765','');


%% Clean grRules, and regenerate the "genes" and "rxnGeneMat" fields

% update the gene-related fields
[ihuman.grRules,ihuman.genes,ihuman.rxnGeneMat] = translateGrRules(ihuman.grRules,'ENSG','ENSG');

% also update the protein-related fields
[ihuman.prRules,ihuman.proteins,ihuman.rxnProtMat] = translateGrRules(ihuman.grRules,'UniProt','ENSG');


%% Export model

save('../../../model/Human-GEM.mat','ihuman');



