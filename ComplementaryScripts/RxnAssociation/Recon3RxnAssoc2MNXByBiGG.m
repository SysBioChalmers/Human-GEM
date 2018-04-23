%
%   FILE NAME: Recon3RxnAssoc2MNXByBiGG.m
% 
%   DATE CREATED: 2018-04-23
%
%	
%   PROGRAMMER:   Hao Wang
%                 Department of Biology and Biological Engineering
%                 Chalmers University of Technology
% 
% 
%   PURPOSE: Assocate Recon3 reactions through BiGG to MNX
%
%


% Load Recon3D
load('/Users/haowa/Box Sync/HMR3/Recon3D/Published/ModelFiles/Recon3D_301/Recon3D_301.mat');

% Load BiGGRxns database
load('BiGGRxns.mat');

% Associate Recon3D through BiGG to MNX
%===Comprehensive association based on bigg_id and oldids
% From BiGG to BiGG, start with bigg_id
Recon3D.BiGG2BiGG=cell(numel(Recon3D.rxns),1);
Recon3D.BiGG2BiGG(:)={''};
Recon3D.rxnMNXID=cell(numel(Recon3D.rxns),1);
Recon3D.rxnMNXID(:)={''};

% Direct association
[a, b]=ismember(Recon3D.rxns,BiGGRxns.bigg_id);
I=find(a);
Recon3D.BiGG2BiGG(I)=BiGGRxns.bigg_id(b(I));
Recon3D.rxnMNXID(I)=BiGGRxns.MNXrefid(b(I));
numel(find(~cellfun(@isempty,Recon3D.BiGG2BiGG)))  % ans = 5626
numel(find(~cellfun(@isempty,Recon3D.rxnMNXID)))  % ans = 5625

% Retrieve missing ids from old_bigg_ids
for i=1:numel(Recon3D.rxns)
		%Loop through for non-associated ids
		if isempty(Recon3D.BiGG2BiGG{i})
				for j=1:numel(BiGGRxns.oldids)
						if ismember(Recon3D.rxns{i},BiGGRxns.oldids{j})
								Recon3D.BiGG2BiGG{i}=BiGGRxns.bigg_id{j};
								Recon3D.rxnMNXID(i)=BiGGRxns.MNXrefid(j);
						end
				end
		end
end
numel(find(~cellfun(@isempty,Recon3D.BiGG2BiGG)))  % with BiGG association = 5840
numel(find(~cellfun(@isempty,Recon3D.rxnMNXID)))  % with MNX association = 5837


% Locate HMR rxns in Recon3D
Recon3D.rxnHMRID=cell(numel(Recon3D.rxns),1);
Recon3D.rxnHMRID(:)={''};

[a, b]=ismember(Recon3D.rxns,ihuman.rxns);
I=find(a);
Recon3D.rxnHMRID(I)=Recon3D.rxns(I);
numel(find(~cellfun(@isempty,Recon3D.rxnHMRID)))  % ans = 2486

save('Recon3Rxns2MNX.mat','Recon3D');   %2018-4-23

