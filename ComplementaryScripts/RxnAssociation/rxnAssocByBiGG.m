%
%   FILE NAME: rxnAssocByBiGG.m
% 
%   DATE CREATED: 2017-12-01
%        UPDATED: 2018-02-08
%	
%   PROGRAMMER:   Hao Wang
%                 Department of Biology and Biological Engineering
%                 Chalmers University of Technology
% 
% 
%   PURPOSE: Assocate HMR2 reactions to BiGG based on provided
%            external database identifiers
%


% Load HMR database ver 2.0.2
load('HMRdatabase2_02.mat');

num=numel(ihuman.rxns);

% Associate HMR2 to Recon2 reactions
% Import HMR2-Recon2 links provided by Adil
cd('/Users/haowa/Box Sync/HMR3/Recon2/');
T=readtable('Recon2_RAVEN.xlsx','Sheet','Reaction list','ReadVariableNames',1);
Recon2HMR2=table2struct(T,'ToScalar',true);
save('Recon2HMR2.mat','Recon2HMR2');

% Add Recon2 association
ihuman.rxnRecon2=cell(num,1);
ihuman.rxnRecon2(:,1)={''};
[a, b]=ismember(ihuman.rxns,Recon2HMR2.HMR_RXNS);
I=find(a);
ihuman.rxnRecon2(I)=Recon2HMR2.ReactionAbbreviation(b(I));
numel(find(~cellfun(@isempty,ihuman.rxnRecon2)))  % ans = 4845

% Compare BiGG/Recon1 and Recon2 association
x=0;  %compared pairs
y=0;  %matched pairs
for i=1:num
	if ~isempty(ihuman.rxnRecon2{i}) && ~isempty(ihuman.rxnBiGGID{i})
		x=x+1;
		if isequal(ihuman.rxnRecon2{i},ihuman.rxnBiGGID{i})
			y=y+1;
		end
	end
end
% Discard this list due to some mistakes, suggested by Adil

% Load BiGGRxns database (2018-01-17)
load('BiGGRxns.mat');

% A quick screening of matched rxn ids from HMR2 to BiGG
% BiGG id association
numel(find(ismember(ihuman.rxnBiGGID,BiGGRxns.bigg_id)))
numel(find(~cellfun(@isempty,ihuman.rxnBiGGID)))
% ans = 2375 of 3064 were mapped to BiGG database

% KEGG id association
numel(find(ismember(ihuman.rxnKEGGID,BiGGRxns.bigg_id)))
numel(find(~cellfun(@isempty,ihuman.rxnKEGGID)))
% ans = 0 of 1767 were mapped to BiGG database

% EHMN id association
numel(find(ismember(ihuman.rxnEHMNID,BiGGRxns.bigg_id)))
numel(find(~cellfun(@isempty,ihuman.rxnEHMNID)))
% ans = 719 of 1955 were mapped to BiGG database

% HepaoNET1 id association
numel(find(ismember(ihuman.rxnHepatoNET1ID,BiGGRxns.bigg_id)))
numel(find(~cellfun(@isempty,ihuman.rxnHepatoNET1ID)))
% ans = 1173 of 2383 were mapped to BiGG database

% Reactome id association
numel(find(ismember(ihuman.rxnREACTOMEID,BiGGRxns.bigg_id)))
numel(find(~cellfun(@isempty,ihuman.rxnREACTOMEID)))
% ans = 0 of 217 were mapped to BiGG database


%===Comprehensive association based on bigg_id and oldids
% From BiGG to BiGG, start with bigg_id
ihuman.BiGG2BiGG=cell(num,1);
ihuman.BiGG2BiGG(:,1)={''};
[a, b]=ismember(ihuman.rxnBiGGID,BiGGRxns.bigg_id);
I=find(a);
ihuman.BiGG2BiGG(I)=BiGGRxns.bigg_id(b(I));
numel(find(~cellfun(@isempty,ihuman.BiGG2BiGG)))  % ans = 2375
% Retrieve missing ids from old_bigg_ids
count=0;
for i=1:num
	%Loop through for non-associated ids
	if isempty(ihuman.BiGG2BiGG{i}) && ~isempty(ihuman.rxnBiGGID{i})
		for j=1:numel(BiGGRxns.oldids)
			if ismember(ihuman.rxnBiGGID{i},BiGGRxns.oldids{j})
				ihuman.BiGG2BiGG{i}=BiGGRxns.bigg_id{j};
				count=count+1;
			end
		end
	end
end
numel(find(~cellfun(@isempty,ihuman.BiGG2BiGG)))  % ans = 2484
%---

% From HepatoNET1 to BiGG
ihuman.HepatoNet12BiGG=cell(num,1);
ihuman.HepatoNet12BiGG(:,1)={''};
[a, b]=ismember(ihuman.rxnHepatoNET1ID,BiGGRxns.bigg_id);
I=find(a);
ihuman.HepatoNet12BiGG(I)=BiGGRxns.bigg_id(b(I));
numel(find(~cellfun(@isempty,ihuman.HepatoNet12BiGG)))  % ans = 1173
% Retrieve missing ids from old_bigg_ids
count=0;
for i=1:num
	%Loop through for non-associated ids
	if ~isempty(ihuman.rxnHepatoNET1ID{i}) && isempty(ihuman.HepatoNet12BiGG{i})
		for j=1:numel(BiGGRxns.oldids)
			if ismember(ihuman.rxnHepatoNET1ID{i},BiGGRxns.oldids{j})
				ihuman.HepatoNet12BiGG{i}=BiGGRxns.bigg_id{j};
				count=count+1;
			end
		end
	end
end
numel(find(~cellfun(@isempty,ihuman.HepatoNet12BiGG)))  % ans = 1242
%count=69

% From EHMN to BiGG
ihuman.EHMN2BiGG=cell(num,1);
ihuman.EHMN2BiGG(:,1)={''};
[a, b]=ismember(ihuman.rxnEHMNID,BiGGRxns.bigg_id);
I=find(a);
ihuman.EHMN2BiGG(I)=BiGGRxns.bigg_id(b(I));
numel(find(~cellfun(@isempty,ihuman.EHMN2BiGG)))  % ans = 719
% Retrieve missing ids from old_bigg_ids
count=0;
for i=1:num
	%Loop through for non-associated ids
	if ~isempty(ihuman.rxnEHMNID{i}) && isempty(ihuman.EHMN2BiGG{i})
		for j=1:numel(BiGGRxns.oldids)
			if ismember(ihuman.rxnEHMNID{i},BiGGRxns.oldids{j})
				ihuman.EHMN2BiGG{i}=BiGGRxns.bigg_id{j};
				count=count+1;
			end
		end		
	end
end
numel(find(~cellfun(@isempty,ihuman.EHMN2BiGG)))  % ans = 732
%count=13

% From KEGG to BiGG ids
% No KEGG id associated to BiGG ids

% From Reactome to BiGG ids
% No Reactome id associated to BiGG ids

%===Quick screening of matched rxn ids to BiGG name
index=find(~cellfun(@isempty,ihuman.rxnKEGGID));
numel(intersect(ihuman.rxnKEGGID(index),BiGGRxns.name))  % ans = 0  
index=find(~cellfun(@isempty,ihuman.rxnHepatoNET1ID));
numel(intersect(ihuman.rxnHepatoNET1ID(index),BiGGRxns.name))  % ans = 0
index=find(~cellfun(@isempty,ihuman.rxnREACTOMEID));
numel(intersect(ihuman.rxnREACTOMEID(index),BiGGRxns.name))  % ans = 0

index=find(~cellfun(@isempty,ihuman.rxnBiGGID));
numel(intersect(ihuman.rxnBiGGID(index),BiGGRxns.name))  % ans = 19
% Check the detail
[a, b]=ismember(ihuman.rxnBiGGID,BiGGRxns.name);
I=find(a);
isequal(ihuman.rxnBiGGID(I),BiGGRxns.name(b(I))) % ans = logical 1
count=0;
for i=1:numel(I)
	%Loop through BiGG ids matched to name
	if isempty(ihuman.BiGG2BiGG{I(i)}) && ~isempty(ihuman.rxnBiGGID{I(i)}) 
		count=count+1;
	end
end
%count=0, these ids had alrady been associated, so ignore them

index=find(~cellfun(@isempty,ihuman.rxnEHMNID));
numel(intersect(ihuman.rxnEHMNID(index),BiGGRxns.name))  % ans = 5
%    'RE0915'
%    'RE0926'
%    'RE0935'
%    'RE0944'
%    'RE0958'
% Check the detail
[a, b]=ismember(ihuman.rxnEHMNID,BiGGRxns.name);
I=find(a);
isequal(ihuman.rxnEHMNID(I),BiGGRxns.name(b(I))) % ans = logical 1
count=0;
for i=1:numel(I)
	%Loop through EHMN ids matched to name
	if ~isempty(ihuman.rxnEHMNID{I(i)}) && isempty(ihuman.EHMN2BiGG{I(i)})
		%disp([BiGGRxns.name{b(I(i))} ':' BiGGRxns.bigg_id{b(I(i))}]);
		ihuman.EHMN2BiGG{I(i)}=BiGGRxns.bigg_id{b(I(i))};
		count=count+1;
	end
end
%count=10, these rxns were associated accroding to BiGG reaction name
numel(find(~cellfun(@isempty,ihuman.EHMN2BiGG)))  % ans = 742

save('ihuman2BiGG.mat','ihuman');       %2018-2-9

%===Unify results
index=find(~cellfun(@isempty,ihuman.HepatoNet12BiGG));
numel(find(~cellfun(@isempty,ihuman.EHMN2BiGG(index))))  % ans = 0
% There is no overlap between EHMN and HepatoNet1, so direct merge
ihuman.HMR2BiGG=ihuman.EHMN2BiGG;
ihuman.HMR2BiGG(index)=ihuman.HepatoNet12BiGG(index);
numel(find(~cellfun(@isempty,ihuman.HMR2BiGG)))  % ans = 1984
% Combine with BiGG associations
index=find(~cellfun(@isempty,ihuman.BiGG2BiGG));
numel(find(~cellfun(@isempty,ihuman.HMR2BiGG(index))))  % ans = 7
% There are 7 overlaps, check detail
for i=1:numel(index)
	%Loop through BiGG id associations
	if isempty(ihuman.HMR2BiGG{index(i)})  % Direct merge non-overlap elements
		ihuman.HMR2BiGG{index(i)}=ihuman.BiGG2BiGG{index(i)};
	else
		% Leave conflict association for manual curation
		if ~isequal(ihuman.BiGG2BiGG{index(i)},ihuman.HMR2BiGG{index(i)})
			A=find(strcmp(ihuman.HMR2BiGG{index(i)},BiGGRxns.bigg_id));
			B=find(strcmp(ihuman.BiGG2BiGG{index(i)},BiGGRxns.bigg_id));
			fprintf('%d:\t%s(%s)-%s(%s)\n',index(i),ihuman.HMR2BiGG{index(i)},....
			BiGGRxns.MNXrefid{A},ihuman.BiGG2BiGG{index(i)},BiGGRxns.MNXrefid{B});
		else
			% This is consistent association, leave it and print out
			fprintf('%d:\t%s-%s\n',index(i),ihuman.HMR2BiGG{index(i)},ihuman.BiGG2BiGG{index(i)});
		end
	end
end
numel(find(~cellfun(@isempty,ihuman.HMR2BiGG)))  % ans = 4461
%---Manual curation following 7 associations
%669:	  SERD_L-SERD_L                            %leave it
%992:	  r0595(MNXR105354)-MCPST(MNXR101422)      %use MNX assoc to KEGG though
%1169:	RE3372C(MNXR103887)-FTHFCL(MNXR99668)
%3032:	RE2410C(MNXR97383)-DHCR71r(MNXR97383)    %use MNX assoc to KEGG though
%4089:	RE2626C(MNXR103704)-P45027A13m(MNXR102266)
%4417:	PPNCL2(MNXR103119)-PPNCL(MNXR103118)
%7594:	r0845(MNXR105086)-UGLCNACtg(MNXR105086)

% Curation results: an excel sheet was save in 'BiGG' subfolder
%992:		r0595(MNXR105354)			remove:MCPST(MNXR101422)
%1169:	RE3372C(MNXR103887)		remove:FTHFCL(MNXR99668)
%4089:	P45027A13m(MNXR102266)remove:RE2626C(MNXR103704)
%4417:	r0671/PPNCL2(MNXR103119)replace:PPNCL(MNXR103118) with PPNCL2

% Modify relevant external identifiers
ihuman.rxnBiGGID{992}='';
ihuman.rxnKEGGID{992}='R03105';
ihuman.rxnBiGGID{1169}='';
ihuman.rxnEHMNID{4089}='';
ihuman.rxnBiGGID{4417}='PPNCL2';

% Assign manually curated 5 ones
%ihuman.HMR2BiGG([992 1169 3032 4089 4417 7594])  % Have a check
ihuman.HMR2BiGG{4089}='P45027A13m';
ihuman.HMR2BiGG{4417}='PPNCL2';
numel(find(~cellfun(@isempty,ihuman.HMR2BiGG)))  % ans = 4461

save('ihuman2BiGG.mat','ihuman');       %2018-2-9

% Link HMR reactions through BiGG association to MNX
ihuman.HMR2BiGG2MNX=cell(num,1);
ihuman.HMR2BiGG2MNX(:,1)={''};
[a, b]=ismember(ihuman.HMR2BiGG,BiGGRxns.bigg_id);
I=find(a);
ihuman.HMR2BiGG2MNX(I)=BiGGRxns.MNXrefid(b(I));
numel(find(~cellfun(@isempty,ihuman.HMR2BiGG2MNX)))  % ans = 4460

% Find out the BiGG ids that miss MNX links
index=find(~cellfun(@isempty,ihuman.HMR2BiGG));
missingMNX=find(cellfun(@isempty,ihuman.HMR2BiGG2MNX(index)));
ihuman.HMR2BiGG(index(missingMNX))
% ans = 'PYK6'   (EHMN)

save('ihuman2BiGG.mat','ihuman');       %2018-2-9

% Rename mat file
load('ihuman2BiGG.mat');
save('ihumanRxns2BiGG.mat','ihuman');   %2018-4-11