%
% FILE NAME:    geneAssoc.m
% 
% DATE CREATED: 2018-07-24
% 
% PROGRAMMER:   Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% 
% PURPOSE: This script aims to generate comprehensive gene associations
%          between HMR2 and Recon3D (i.e. relate genes from Entrenz ids
%          to Ensembl ids)
%

% Load Recon3D genes
load('Recon3DRaven.mat');    
geneAssoc.Recon3D=unique(regexprep(Recon3DRaven.genes, '\.\d+', ''));
geneAssoc.Recon3D(1)=[];    % Remove the mistake geneid '0'
geneAssoc.EnsemblID=cell(numel(geneAssoc.Recon3D),1);
geneAssoc.EnsemblID(:)={''};


% Load gene associations from Ensembl
load('Ensembl2NCBI.mat');
% Load the manual curation results of gene association from Ensembl to Entrez
T=readtable('/Users/haowa/Box Sync/HMR3/Curation Files/HMRdatabase2-20180724.xlsx',....
'Sheet','GENES','ReadVariableNames',1);
geneSheet=table2struct(T,'ToScalar',true);
geneSheet.EntrezID=strtrim(cellstr(num2str(geneSheet.GENEID2_INENTREZ)));
geneSheet.EnsemblID=geneSheet.GENEIDsToWorkWith_WhichWereMostlyBasedOnENTREZ_AndSecondlyENSEM;


% Try to get the one-to-one association first
for i=1:numel(geneAssoc.Recon3D)
		index1=find(strcmp(Ensembl2NCBI.NCBIGeneID,geneAssoc.Recon3D{i}));  % Based on Ensembl
		index2=find(strcmp(geneSheet.EntrezID,geneAssoc.Recon3D{i}));  % Based on manual curation
		if length(index1)==1
				geneAssoc.EnsemblID{i}=Ensembl2NCBI.genes{index1};
		elseif length(index2)==1
				geneAssoc.EnsemblID{i}=geneSheet.EnsemblID{index2};
		else
				% Get the Entrez ids with none or multiple associations to Ensembl
				disp([geneAssoc.Recon3D{i} ': ' num2str(length(index1)) ' - ' num2str(length(index2))]);
		end
end

%=====
%Entrez id: #assoc2Ensembl - #assoc2curation
%100507855: x - 0
%102724197: 0 - 0
%102724560: 2 - 0
%10554: 7 - 7
%10919: 7 - 4
%1589: 6 - 3
%160728: 2 - 2
%28: 2 - 2
%4758: 8 - 3
%5265: 2 - 0
%534: 7 - 2
%645740: x - 0
%7922: 6 - 3
%7923: 6 - 6
%79581: 2 - 0
%8041: x - 0
%8418: 0 - 0
%8705: 5 - 2
%8972: 2 - 2
%9374: 8 - 6
%=====

% The following assocaitons were manually conducted through searching the NCBI
% website, which gives unique association based on the primary assembly. The other
% Ensembl associations from alternate assemblies are ignored.
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'102724560'))}='ENSG00000274276';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'10554'))}='ENSG00000204310';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'10919'))}='ENSG00000204371';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'1589'))}='ENSG00000231852';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'160728'))}='ENSG00000256870';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'28'))}='ENSG00000175164';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'4758'))}='ENSG00000204386';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'5265'))}='ENSG00000197249';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'534'))}='ENSG00000213760';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'7922'))}='ENSG00000112473';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'7923'))}='ENSG00000204228';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'79581'))}='ENSG00000185803';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'8705'))}='ENSG00000235863';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'8972'))}='ENSG00000257335';
geneAssoc.EnsemblID{find(strcmp(geneAssoc.Recon3D,'9374'))}='ENSG00000221988';

save('geneAssoc.mat','geneAssoc');  % 2018-07-24
% This gene association results could be contineously modified/improved when
% mistakes/problems are spotted.
