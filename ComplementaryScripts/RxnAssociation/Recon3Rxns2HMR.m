%
% FILE NAME:    Recon3Rxns2HMR.m
% 
% DATE CREATED: 2018-05-28
% 
% PROGRAMMER:   Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% 
% PURPOSE: This script attempts to assign proper HMR id(s) to each
%          Recon3D rxn, and output the association into a defined
%          array structure.
%

% Load rxn association from HMR2 and Recon3D to MNX
load('Recon3Rxns2MNX.mat');    %Recon3D rxn association to MNX
load('ihumanRxns2MNX.mat');    %HMR rxn association to MNX

% Add assocation based on BiGG database
Recon3D.BiGG2HMR=cell(numel(Recon3D.rxns),1);
Recon3D.BiGG2HMR(:)={''};
index=find(~cellfun(@isempty,Recon3D.rxnBiGGID));
[a, b]=ismember(Recon3D.rxnBiGGID(index),ihuman.HMR2BiGG);
I=find(a);
Recon3D.BiGG2HMR(index(I))=ihuman.rxns(b(I));
numel(find(~cellfun(@isempty,Recon3D.BiGG2HMR)))  % ans = 4522

% Check consistency
indexBiGG=find(~cellfun(@isempty,Recon3D.BiGG2HMR));
indexHMR=find(~cellfun(@isempty,Recon3D.rxnHMRID));
overlap=intersect(indexBiGG,indexHMR);
ind=find(~cellfun(@isequal,Recon3D.BiGG2HMR(overlap),Recon3D.rxnHMRID(overlap)));
fprintf('%s: %s-%s\n',num2str(overlap(ind)),Recon3D.BiGG2HMR{overlap(ind)},Recon3D.rxnHMRID{overlap(ind)});
% 11741: HMR_8671-HMR_6533
% Keep the direct association based on Recon3 rxn id and remove the other one
Recon3D.BiGG2HMR{11741}='';  

