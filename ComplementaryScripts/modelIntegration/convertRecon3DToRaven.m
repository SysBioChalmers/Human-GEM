%
% FILE NAME:    convertRecon3DToRaven.m
% 
% DATE CREATED: 2018-04-25
% 
% PROGRAMMER:   Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% 
% PURPOSE: Convert Recon3D model to RAVEN format
%

% Load Recon3D
load('/Users/haowa/Box Sync/HMR3/Recon3D/Published/ModelFiles/Recon3D_301/Recon3D_301.mat');

% Load HMR2
load('HMRdatabase2_00.mat');

% Compartment information is missing, and is added here
Recon3D.comps=[ihuman.comps;'i'];
Recon3D.compNames=[ihuman.compNames;'Inner mitochondria'];

% Modify compartment system according to HMR2: e,x -> s,p
Recon3D.originalMets=Recon3D.mets;   % Make a backup
Recon3D.mets=regexprep(Recon3D.mets,'\[e\]','[s]');
Recon3D.mets=regexprep(Recon3D.mets,'\[x\]','[p]');

% The metPubChemID field is probablmatic and remove it
Recon3D=rmfield(Recon3D,'metPubChemID');

% Convert to RAVEN format
Recon3DRaven=ravenCobraWrapper(Recon3D);

% Add equations and orginal mets
Recon3DRaven.rxnEquations=constructEquations(Recon3DRaven);
Recon3DRaven.originalMets=Recon3D.originalMets;

save('Recon3DRaven.mat','Recon3DRaven');

