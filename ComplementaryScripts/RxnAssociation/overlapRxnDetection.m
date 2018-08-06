%
% FILE NAME:    overlapRxnDetection.m
% 
% DATE CREATED: 2018-07-27
%
% 
% PROGRAMMER:   Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% 
% PURPOSE: detect duplication rxns between HMR2 and Recon3D after
%          comprehensive metaboite association
%


% Load HMR2 and Recon3D models
load('Recon3DRaven.mat');     %Recon3D
load('HMRdatabase2_00.mat');  %HMR2

% Reduce Recon3D model by removing the reactions that
% were directly taken from HMR2
[a, b]=ismember(Recon3DRaven.rxns,ihuman.rxns);
rxnToRemove=find(a);
reducedRecon3D=removeReactions(Recon3DRaven,rxnToRemove,1,1,1);
reducedRecon3D.id='reducedRecon3D';

% Merge two models and detect overlap reactions
HMR2Recon3D=mergeModels({ihuman reducedRecon3D});
overlap=detectDuplicateRxns(HMR2Recon3D,0);  %ignore the reaction direction

% Output overlap reaction pairs
index=find(cellfun(@numel, overlap.group)==2);
overlapRxns.rxnHMRID={};
overlapRxns.rxnRecon3DID={};
for i=1:length(index)
		m=index(i);
		if ismember(overlap.group{m}{1},ihuman.rxns) && ismember(overlap.group{m}{2},Recon3DRaven.rxns)
				overlapRxns.rxnHMRID=[overlapRxns.rxnHMRID;overlap.group{m}{1}];
				overlapRxns.rxnRecon3DID=[overlapRxns.rxnRecon3DID;overlap.group{m}{2}];
		end
end
save('overlapRxns.mat','overlapRxns');
