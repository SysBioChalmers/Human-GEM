%
%   FILE NAME: Recon3MetAssoc2MNXByBiGG.m
% 
%   PURPOSE: Assocate Recon3 metabolites through BiGG DB to MNX
%


% Load Recon3Mets2MNX
load('Recon3Mets2MNX.mat');

% Load BiGGRxns database
load('BiGGMets.mat');

% Associate Recon3D mets through BiGG database to MNX
% Two new fields are added: metBiGGDB2BiGG and metBiGGDB2MNX
%===Comprehensive association based on bigg_id and oldids
% From BiGG to BiGG, start with bigg_id
Recon3D.metBiGGDB2BiGG=cell(numel(Recon3D.mets),1);
Recon3D.metBiGGDB2BiGG(:)={''};
Recon3D.metBiGGDB2MNX=cell(numel(Recon3D.mets),1);
Recon3D.metBiGGDB2MNX(:)={''};

% Direct association
Recon3D.mets=regexprep(Recon3D.mets,'[','_');
Recon3D.mets=regexprep(Recon3D.mets,']','');

[a, b]=ismember(Recon3D.mets,BiGGMets.mets);
I=find(a);
Recon3D.metBiGGDB2BiGG(I)=BiGGMets.mets(b(I));
Recon3D.metBiGGDB2MNX(I)=BiGGMets.metMNXID(b(I));
numel(find(~cellfun(@isempty,Recon3D.metBiGGDB2BiGG)))  % ans = 6812
numel(find(~cellfun(@isempty,Recon3D.metBiGGDB2MNX)))   % ans = 4731

% Retrieve missing ids from oldids
for i=1:numel(Recon3D.mets)
		%Loop through for non-associated ids
		if isempty(Recon3D.metBiGGDB2BiGG{i})
				for j=1:numel(BiGGMets.oldids)
						if ismember(Recon3D.mets{i},BiGGMets.oldids{j})
								Recon3D.metBiGGDB2BiGG{i}=BiGGMets.mets{j};
								Recon3D.metBiGGDB2MNX{i}=BiGGMets.metMNXID{j};
						end
				end
		end
end
numel(find(~cellfun(@isempty,Recon3D.metBiGGDB2BiGG)))  % with BiGG association = 6978
numel(find(~cellfun(@isempty,Recon3D.metBiGGDB2MNX)))   % with MNX association = 4895

save('Recon3Mets2MNX.mat','Recon3D');   % 2018-04-24

