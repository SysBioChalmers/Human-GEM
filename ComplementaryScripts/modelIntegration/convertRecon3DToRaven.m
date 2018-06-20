%
% FILE NAME:    convertRecon3DToRaven.m
% 
% DATE CREATED: 2018-04-25
%     MODIFIED: 2018-06-20
% 
% PROGRAMMER:   Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% 
% PURPOSE: Convert Recon3D model to RAVEN format
%          replace met with HMR info
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

% Create data structure of one-to-one met assocation between HMR2 and Recon3D 
load('metAssocHMR2Recon3.mat');   %To get the manually curated met association info
% Directly save the unique associations
single_ind=find(cellfun(@numel,metAssocHMR2Recon3.metRecon3DID)==1);
metAssoc.metHMRID=metAssocHMR2Recon3.metHMRID(single_ind);
metAssoc.metRecon3DID=reformatElements(metAssocHMR2Recon3.metRecon3DID(single_ind),'cell2str');
metAssoc.metNames=metAssocHMR2Recon3.metNames(single_ind);
% Keep the associations between one HMR id to multiple Recon3D id
multi_ind=find(cellfun(@numel,metAssocHMR2Recon3.metRecon3DID)>1);
for i=1:length(multi_ind)
		m=multi_ind(i);
		num=numel(metAssocHMR2Recon3.metRecon3DID{m});
		temp=cell(num,1);
		temp(:)={metAssocHMR2Recon3.metHMRID{m}};
		names=cell(num,1);
		names(:)={metAssocHMR2Recon3.metNames{m}};
		metAssoc.metHMRID=[metAssoc.metHMRID;temp];
		metAssoc.metRecon3DID=[metAssoc.metRecon3DID;transpose(metAssocHMR2Recon3.metRecon3DID{m})];
		metAssoc.metNames=[metAssoc.metNames;names];
end
save('metAssoc.mat','metAssoc');         % 2018-06-20
% This file is designed for convinient model integration purpose,
% and should be updated once metAssocHMR2Recon3.mat is changed!

