%
% FILE NAME:    removeNonPrimaryGeneIDs.m
% 
% PURPOSE: This script updates humanGEM gene-reaction rules (grRules), such
%          that only primary assembly gene IDs are used. All gene IDs that
%          are associated with a non-primary assembly (i.e., allele
%          variants), will be removed from the model.


%% Load humanGEM model

load('humanGEM.mat');  % version 0.8.3


%% Remove non-primary assembly Ensembl gene IDs from humanGEM

% get list of primary assmbly genes
fid = fopen('../../ComplementaryData/Ensembl/ensembl_ID_mapping.tsv');
tmp = textscan(fid,'%s%s%s%s%s%s','Delimiter','\t','Headerlines',2);
fclose(fid);

% get list of primary assembly ENSG IDs
ensg_ids = unique(tmp{1});
ensg_ids(cellfun(@isempty,ensg_ids)) = [];  % remove empty ID if it exists

% obtain list of humanGEM genes, with non-primary genes removed
new_ids = ihuman.genes;
new_ids(~ismember(new_ids,ensg_ids)) = {''};

% remove non-primary genes from model
[grRules,genes,rxnGeneMat] = translateGrRules(ihuman.grRules,[ihuman.genes,new_ids]);
ihuman.grRules = grRules;
ihuman.genes = genes;
ihuman.rxnGeneMat = rxnGeneMat;


%% Update protein-related fields

% update protein fields
[prRules,proteins,rxnProtMat] = translateGrRules(ihuman.grRules,'UniProt','ENSG');
ihuman.prRules = prRules;
ihuman.proteins = proteins;
ihuman.rxnProtMat = rxnProtMat;


%% Clear intermediate variables and export model

clear ensg_ids genes grRules new_ids proteins prRules rxn_ind rxnGeneMat
clear rxnProtMat tmp fid

save('../../model/Human-GEM.mat','ihuman');



