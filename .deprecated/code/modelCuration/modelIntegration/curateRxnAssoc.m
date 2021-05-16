%
% FILE NAME:    curateRxnAssoc.m
% 
% PURPOSE: Curate rxnAssoc array structure by adding the reaction bounds
%          information extracted from HMR2 and Recon3D 
% 

% Load input
load('rxnAssoc.mat');       % rxnAssoc
load('humanGEM.mat');       % humanGEM v0.2.0
load('Recon3DRaven.mat');   % Recon3D

% Add fields for reaction bounds information
[~, indHMR]=ismember(rxnAssoc.rxnHMRID,ihuman.rxns);
rxnAssoc.lbHMR=ihuman.lb(indHMR);
rxnAssoc.ubHMR=ihuman.ub(indHMR);

[~, indRecon3D]=ismember(rxnAssoc.rxnRecon3DID,Recon3DRaven.rxns);
rxnAssoc.lbRecon3D=Recon3DRaven.lb(indRecon3D);
rxnAssoc.ubRecon3D=Recon3DRaven.ub(indRecon3D);

% Save to modelIntegration subfolder
save('rxnAssoc.mat','rxnAssoc');
