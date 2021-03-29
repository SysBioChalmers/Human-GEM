%
% FILE NAME:    updateHMR2Genes.m
%
% PURPOSE: Update gene curation results into the model. The updated
%          fields include grRules, genes, rxnGeneMat and geneComps.
%


% Load the curated Ensembl gene ids with removal of confirmed pseudogenes
% and some updated ids
fid = fopen('../../ComplementaryData/modelCuration/pseudogeneCheck.tsv','r');
geneCuration=textscan(fid,'%s %s %s %s %s %s','Delimiter','\t','HeaderLines',4);
fclose(fid);

% Use the original and finalized ids for converstion
originalEnsemblID = geneCuration{1};
finalEnsemblID = geneCuration{5};

% Prepare a NX2 cell array, in which the first column contains original ids
% and the second column contains curated ones
geneConversion=cell(numel(originalEnsemblID),2);
geneConversion(:,:)={''};
geneConversion(:,1)=originalEnsemblID;
geneConversion(:,2)=finalEnsemblID;

% Re-generate the updated fields
load('HMRdatabase2_02.mat');
[newGrRules,newGenes,rxnGeneMat] = translateGrRules(ihuman.grRules,geneConversion);

% Update and save the model
ihuman.grRules=newGrRules;
ihuman.genes=newGenes;
ihuman.rxnGeneMat=rxnGeneMat;
ihuman.geneComps=ihuman.geneComps(1:numel(ihuman.genes));
save('../../model/Human-GEM.mat','ihuman');
