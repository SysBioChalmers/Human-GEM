%
%   FILE NAME: rxnAssocByBiGG.m
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

% Load BiGGRxns database (2018-05-18)
load('BiGGRxns.mat');

% A quick screening of matched rxn ids from HMR2 to BiGG
% BiGG id association
numel(find(ismember(ihuman.rxnBiGGID,BiGGRxns.rxns)))
numel(find(~cellfun(@isempty,ihuman.rxnBiGGID)))
% ans = 2377 of 3064 were mapped to BiGG database

% KEGG id association
numel(find(ismember(ihuman.rxnKEGGID,BiGGRxns.rxns)))
numel(find(~cellfun(@isempty,ihuman.rxnKEGGID)))
% ans = 0 of 1767 were mapped to BiGG database

% EHMN id association
numel(find(ismember(ihuman.rxnEHMNID,BiGGRxns.rxns)))
numel(find(~cellfun(@isempty,ihuman.rxnEHMNID)))
% ans = 821 of 1955 were mapped to BiGG database

% HepaoNET1 id association
numel(find(ismember(ihuman.rxnHepatoNET1ID,BiGGRxns.rxns)))
numel(find(~cellfun(@isempty,ihuman.rxnHepatoNET1ID)))
% ans = 1218 of 2383 were mapped to BiGG database

% Reactome id association
numel(find(ismember(ihuman.rxnREACTOMEID,BiGGRxns.rxns)))
numel(find(~cellfun(@isempty,ihuman.rxnREACTOMEID)))
% ans = 0 of 217 were mapped to BiGG database


%===Comprehensive association based on bigg_id and oldids
% From BiGG to BiGG, start with bigg_id
ihuman.BiGG2BiGG=cell(num,1);
ihuman.BiGG2BiGG(:,1)={''};
[a, b]=ismember(ihuman.rxnBiGGID,BiGGRxns.rxns);
I=find(a);
ihuman.BiGG2BiGG(I)=BiGGRxns.rxns(b(I));
numel(find(~cellfun(@isempty,ihuman.BiGG2BiGG)))  % ans = 2377
% Retrieve missing ids from old_bigg_ids
count=0;
for i=1:num
	%Loop through for non-associated ids
	if isempty(ihuman.BiGG2BiGG{i}) && ~isempty(ihuman.rxnBiGGID{i})
		for j=1:numel(BiGGRxns.oldids)
			if ismember(ihuman.rxnBiGGID{i},BiGGRxns.oldids{j})
				ihuman.BiGG2BiGG{i}=BiGGRxns.rxns{j};
				count=count+1;
			end
		end
	end
end
numel(find(~cellfun(@isempty,ihuman.BiGG2BiGG)))  % ans = 2489
%count=112

% From HepatoNET1 to BiGG
ihuman.HepatoNet12BiGG=cell(num,1);
ihuman.HepatoNet12BiGG(:,1)={''};
[a, b]=ismember(ihuman.rxnHepatoNET1ID,BiGGRxns.rxns);
I=find(a);
ihuman.HepatoNet12BiGG(I)=BiGGRxns.rxns(b(I));
numel(find(~cellfun(@isempty,ihuman.HepatoNet12BiGG)))  % ans = 1218
% Retrieve missing ids from old_bigg_ids
count=0;
for i=1:num
	%Loop through for non-associated ids
	if ~isempty(ihuman.rxnHepatoNET1ID{i}) && isempty(ihuman.HepatoNet12BiGG{i})
		for j=1:numel(BiGGRxns.oldids)
			if ismember(ihuman.rxnHepatoNET1ID{i},BiGGRxns.oldids{j})
				ihuman.HepatoNet12BiGG{i}=BiGGRxns.rxns{j};
				count=count+1;
			end
		end
	end
end
numel(find(~cellfun(@isempty,ihuman.HepatoNet12BiGG)))  % ans = 1348
%count=130

% From EHMN to BiGG
ihuman.EHMN2BiGG=cell(num,1);
ihuman.EHMN2BiGG(:,1)={''};
[a, b]=ismember(ihuman.rxnEHMNID,BiGGRxns.rxns);
I=find(a);
ihuman.EHMN2BiGG(I)=BiGGRxns.rxns(b(I));
numel(find(~cellfun(@isempty,ihuman.EHMN2BiGG)))  % ans = 821
% Retrieve missing ids from old_bigg_ids
count=0;
for i=1:num
	%Loop through for non-associated ids
	if ~isempty(ihuman.rxnEHMNID{i}) && isempty(ihuman.EHMN2BiGG{i})
		for j=1:numel(BiGGRxns.oldids)
			if ismember(ihuman.rxnEHMNID{i},BiGGRxns.oldids{j})
				ihuman.EHMN2BiGG{i}=BiGGRxns.rxns{j};
				count=count+1;
			end
		end		
	end
end
numel(find(~cellfun(@isempty,ihuman.EHMN2BiGG)))  % ans = 847
%count=26

% From KEGG to BiGG ids
% No KEGG id associated to BiGG ids

% From Reactome to BiGG ids
% No Reactome id associated to BiGG ids

%===Quick screening of matched rxn ids to BiGG name
index=find(~cellfun(@isempty,ihuman.rxnKEGGID));
numel(intersect(ihuman.rxnKEGGID(index),BiGGRxns.rxnNames))  % ans = 0  
index=find(~cellfun(@isempty,ihuman.rxnHepatoNET1ID));
numel(intersect(ihuman.rxnHepatoNET1ID(index),BiGGRxns.rxnNames))  % ans = 0
index=find(~cellfun(@isempty,ihuman.rxnREACTOMEID));
numel(intersect(ihuman.rxnREACTOMEID(index),BiGGRxns.rxnNames))  % ans = 0

index=find(~cellfun(@isempty,ihuman.rxnBiGGID));
numel(intersect(ihuman.rxnBiGGID(index),BiGGRxns.rxnNames))  % ans = 14
% Check the detail
[a, b]=ismember(ihuman.rxnBiGGID,BiGGRxns.rxnNames);
I=find(a);
isequal(ihuman.rxnBiGGID(I),BiGGRxns.rxnNames(b(I))) % ans = logical 1
count=0;
for i=1:numel(I)
	%Loop through BiGG ids matched to name
	if isempty(ihuman.BiGG2BiGG{I(i)}) && ~isempty(ihuman.rxnBiGGID{I(i)}) 
		count=count+1;
	end
end
%count=0, these ids had alrady been associated, so ignore them

index=find(~cellfun(@isempty,ihuman.rxnEHMNID));
numel(intersect(ihuman.rxnEHMNID(index),BiGGRxns.rxnNames))  % ans = 5
%    'RE0915'
%    'RE0926'
%    'RE0935'
%    'RE0944'
%    'RE0958'
% Check the detail
[a, b]=ismember(ihuman.rxnEHMNID,BiGGRxns.rxnNames);
I=find(a);
isequal(ihuman.rxnEHMNID(I),BiGGRxns.rxnNames(b(I))) % ans = logical 1
count=0;
for i=1:numel(I)
	%Loop through EHMN ids matched to name
	if ~isempty(ihuman.rxnEHMNID{I(i)}) && isempty(ihuman.EHMN2BiGG{I(i)})
		ihuman.EHMN2BiGG{I(i)}=BiGGRxns.rxns{b(I(i))};
		count=count+1;
	end
end
%count=6, these rxns were associated accroding to BiGG reaction name
numel(find(~cellfun(@isempty,ihuman.EHMN2BiGG)))  % ans = 853


%===Unify results
% At first between HepatoNet1 and EHMN
index=find(~cellfun(@isempty,ihuman.HepatoNet12BiGG));
numel(find(~cellfun(@isempty,ihuman.EHMN2BiGG(index))))  % ans = 3
% There is 3 overlap between EHMN and HepatoNet1, so check them out:
overlapIdx=find(~cellfun(@isempty,ihuman.EHMN2BiGG(index)));
%ihuman.HepatoNet12BiGG(index(overlapIdx))
%ans = {'DHRT_ibcoa','r0706','r0706'}
%ihuman.EHMN2BiGG(index(overlapIdx))
%ans = {'DHRT_ibcoa','RE3247M','RE3247X'}

%directly combine them
ihuman.HMR2BiGG=ihuman.EHMN2BiGG;
ihuman.HMR2BiGG(index)=ihuman.HepatoNet12BiGG(index);
numel(find(~cellfun(@isempty,ihuman.HMR2BiGG)))  % ans = 2198 (1348+853-3)
for i=1:numel(overlapIdx)
		if ~isequal(ihuman.HMR2BiGG{index(overlapIdx(i))},ihuman.EHMN2BiGG{index(overlapIdx(i))})
				ihuman.HMR2BiGG{index(overlapIdx(i))}=strcat(ihuman.HMR2BiGG{index(overlapIdx(i))},';',ihuman.EHMN2BiGG{index(overlapIdx(i))});
		end
end

% Secondly combine with the associations from BiGG ids
indexBiGG=find(~cellfun(@isempty,ihuman.BiGG2BiGG));
indexOthers=find(~cellfun(@isempty,ihuman.HMR2BiGG));
overlapIdx=intersect(indexBiGG,indexOthers);  % ans = 15 overlaps, check later
nonOverlapIdx=setdiff(indexBiGG,indexOthers);
ihuman.HMR2BiGG(nonOverlapIdx)=ihuman.BiGG2BiGG(nonOverlapIdx);
numel(find(~cellfun(@isempty,ihuman.HMR2BiGG)))  % ans = 4672 (2198+2489-15)
% Resolve the conflicts
for i=1:numel(overlapIdx)
		if ~isequal(ihuman.HMR2BiGG{overlapIdx(i)},ihuman.BiGG2BiGG{overlapIdx(i)})
				%ihuman.HMR2BiGG{index(overlapIdx(i))}=strcat(ihuman.HMR2BiGG{index(overlapIdx(i))},';',ihuman.EHMN2BiGG{index(overlapIdx(i))});
				A=find(strcmp(ihuman.HMR2BiGG{overlapIdx(i)},BiGGRxns.rxns));
				B=find(strcmp(ihuman.BiGG2BiGG{overlapIdx(i)},BiGGRxns.rxns));
				fprintf('%s: %s(%s)-%s(%s)\n',num2str(overlapIdx(i)),ihuman.HMR2BiGG{overlapIdx(i)},....
				BiGGRxns.rxnMNXID{A},ihuman.BiGG2BiGG{overlapIdx(i)},BiGGRxns.rxnMNXID{B});
		end
end
%---Manual curation following 8 associations
%992: r0595(MNXR105354)-MCPST(MNXR101422)
%1076: r0669()-ECOAH2m(MNXR97886)
%1169: RE3372C(MNXR103887)-FTHFCL(MNXR99668)
%2954: STS4(MNXR104608)-STS4r(MNXR104608)
%3032: RE2410C(MNXR97383)-DHCR71r(MNXR97383)
%4089: RE2626C(MNXR103704)-P45027A13m(MNXR102266)  
%4417: PPNCL2(MNXR103119)-PPNCL(MNXR103118)
%7594: r0845(MNXR105086)-UGLCNACtg(MNXR105086)

% Curation results: an excel sheet was save in 'BiGG' subfolder
%992:		r0595(MNXR105354)			  remove:MCPST(MNXR101422)
%1076:  r0669()                 remove:ECOAH2m(MNXR97886)
%1169:	RE3372C(MNXR103887)		  remove:FTHFCL(MNXR99668)
%3032:	RE2410C(MNXR97383)		  remove:DHCR71r(MNXR97383)
%4089:	P45027A13m(MNXR102266)  remove:RE2626C(MNXR103704)
%4417:	PPNCL2(MNXR103119)      remove:PPNCL(MNXR103118)
%7594:  r0845(MNXR105086)       remove:UGLCNACtg(MNXR105086)

% Manually assign values
%ihuman.HMR2BiGG([992 1076 1169 2954 3032 4089 4090 4417 7594])  % Have a check
ihuman.HMR2BiGG{992}='r0595';
ihuman.HMR2BiGG{1076}='r0669';
ihuman.HMR2BiGG{1169}='RE3372C';
ihuman.HMR2BiGG{2954}='STS4r';
ihuman.HMR2BiGG{3032}='RE2410C';
ihuman.HMR2BiGG{4089}='';
ihuman.HMR2BiGG{4090}='P45027A13m';
ihuman.HMR2BiGG{4417}='PPNCL2';
ihuman.HMR2BiGG{7594}='r0845';
numel(find(~cellfun(@isempty,ihuman.HMR2BiGG)))  % ans = 4671

save('ihumanRxns2BiGG.mat','ihuman');   %2018-05-18

% Addtional manual curation
ihuman.BiGG2BiGG{992}='';
ihuman.BiGG2BiGG{1076}='';
ihuman.BiGG2BiGG{1169}='';
ihuman.EHMN2BiGG{2954}='';
ihuman.BiGG2BiGG{3032}='';
ihuman.EHMN2BiGG{4089}='';           %wrong cofactors NADP(H) and compartment
ihuman.BiGG2BiGG{4089}='';           %wrong cofactors NADP(H)
ihuman.HepatoNet12BiGG{4090}='';
ihuman.BiGG2BiGG{4417}='';
ihuman.BiGG2BiGG{7594}='';           %wrong compartment

% Now rxn associations to BiGG, EHMN and HepatoNet1 are unified to field HMR2BiGG
save('ihumanRxns2BiGG.mat','ihuman');   %2018-05-21
