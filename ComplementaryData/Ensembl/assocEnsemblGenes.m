%
%   FILE NAME: assocEnsemblGenes.m
% 
%   DATE CREATED: 2018-05-03
%        UPDATED: 2018-05-31
%	
%   PROGRAMMER:   Hao Wang
%                 Department of Biology and Biological Engineering
%                 Chalmers University of Technology
% 
% 
%   PURPOSE: Get exteranl associations of Ensembl genes
% 

% Move to Ensembl folder
cd('/Users/haowa/Box Sync/HMR3/Ensembl');

T=readtable('stableGeneID2NCBIgeneID.txt','ReadVariableNames',1);
T=table2struct(T,'ToScalar',true);

% Format NCBI Gene IDs
NCBIGeneID=num2cell(T.NCBIGeneID);   % convert vector to cell array
NCBIGeneID=cellfun(@num2str,NCBIGeneID,'UniformOutput',false);  % convert each element from double to string
% Clear elements without NCBI Gene ID associations
NaNidx=find(strcmp('NaN',NCBIGeneID));
NCBIGeneID(NaNidx)={''};

% Generate data structure for NCBIGeneID 
Ensembl2NCBI.genes=T.GeneStableID;
Ensembl2NCBI.NCBIGeneID=NCBIGeneID;
save('Ensembl2NCBI.mat','Ensembl2NCBI');

T=readtable('stableGeneID2UniprotID.txt','ReadVariableNames',1);
T=table2struct(T,'ToScalar',true);

% Generate data structure for UniprotID 
Ensembl2Uniprot.genes=T.GeneStableID;
Ensembl2Uniprot.ID=;
Ensembl2Uniprot.ID=;
save('Ensembl2Uniprot.mat','Ensembl2Uniprot');

