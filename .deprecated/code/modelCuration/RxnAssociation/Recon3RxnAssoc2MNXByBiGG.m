%
%   FILE NAME: Recon3RxnAssoc2MNXByBiGG.m
% 
%   PURPOSE: Assocate Recon3 reactions through BiGG to MNX
%


% Load Recon3D
load('/Users/haowa/Box Sync/HMR3/Recon3D/Published/ModelFiles/Recon3D_301/Recon3D_301.mat');
% Load BiGGRxns database
load('BiGGRxns.mat');
% Load HMR2
load('HMRdatabase2_00.mat');


% Associate Recon3D through BiGG to MNX
%===Comprehensive association based on bigg_id and oldids
% From BiGG to BiGG, start with bigg_id
Recon3D.rxnBiGGID=cell(numel(Recon3D.rxns),1);
Recon3D.rxnBiGGID(:)={''};
Recon3D.rxnMNXID=cell(numel(Recon3D.rxns),1);
Recon3D.rxnMNXID(:)={''};

% Direct association
[a, b]=ismember(Recon3D.rxns,BiGGRxns.rxns);
I=find(a);
Recon3D.rxnBiGGID(I)=BiGGRxns.rxns(b(I));
Recon3D.rxnMNXID(I)=BiGGRxns.rxnMNXID(b(I));
numel(find(~cellfun(@isempty,Recon3D.rxnBiGGID)))  % ans = 9494
numel(find(~cellfun(@isempty,Recon3D.rxnMNXID)))   % ans = 5625

% Retrieve missing ids from old_bigg_ids
for i=1:numel(Recon3D.rxns)
		%Loop through for non-associated ids
		if isempty(Recon3D.rxnBiGGID{i})
				for j=1:numel(BiGGRxns.oldids)
						if ismember(Recon3D.rxns{i},BiGGRxns.oldids{j})
								Recon3D.rxnBiGGID{i}=BiGGRxns.rxns{j};
								Recon3D.rxnMNXID{i}=BiGGRxns.rxnMNXID{j};
						end
				end
		end
end
numel(find(~cellfun(@isempty,Recon3D.rxnBiGGID)))  % with BiGG association = 11755
numel(find(~cellfun(@isempty,Recon3D.rxnMNXID)))   % with MNX association = 6740


% Locate HMR rxns in Recon3D
Recon3D.rxnHMRID=cell(numel(Recon3D.rxns),1);
Recon3D.rxnHMRID(:)={''};

[a, b]=ismember(Recon3D.rxns,ihuman.rxns);
I=find(a);
Recon3D.rxnHMRID(I)=Recon3D.rxns(I);
numel(find(~cellfun(@isempty,Recon3D.rxnHMRID)))  % with HMR association = 2486

save('Recon3Rxns2MNX.mat','Recon3D');   % Save to rxnAssoc subfolder 2018-05-18

