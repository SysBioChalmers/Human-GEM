%
%   FILE NAME:  rxnAssocByMNX.m  
% 
%   PURPOSE: HMR2 reaction association to MNXref based on provided
%            external databasse identifiers
%


% Load HMR model with BiGG association
load('ihumanRxns2BiGG.mat');

% Some statistics
num=numel(ihuman.rxns);
rxnAssocNum=zeros(num,1);
for i=1:num
	% Go through each exteranl database and count
	count=0;
	if ~isempty(ihuman.rxnKEGGID{i}) 
		count=count+1;
	end
	if ~isempty(ihuman.rxnEHMNID{i}) 
		count=count+1;
	end
	if ~isempty(ihuman.rxnBiGGID{i}) 
		count=count+1;
	end
	if ~isempty(ihuman.rxnHepatoNET1ID{i}) 
		count=count+1;
	end
	if ~isempty(ihuman.rxnREACTOMEID{i}) 
		count=count+1;
	end
	rxnAssocNum(i,1)=count;
end
ihuman.rxnAssocNum=rxnAssocNum;
numel(find([rxnAssocNum(:,1)]==5))  %ans = 17 fully-associated

% The curation targets are ihuman.rxns(1:5127);
numel(find(~([rxnAssocNum(1:5127,1)]==0)))  %ans = 4340/5127 (85%)
numel(find(~cellfun(@isempty,ihuman.HMR2BiGG(1:5127))))  % ans = 2771 with BiGG association


%===Reaction mapping through MNXref database

% Load NNX reaction references version 3.0
load('MNXrefRxns.mat');

% From BiGG id (combined also from EHMN,HepatoNet1) to MNX id
ihuman.rxnBiGGDB2MNX=cell(num,1);
ihuman.rxnBiGGDB2MNX(:,1)={''};
[a, b]=ismember(ihuman.HMR2BiGG,MNXrefRxns.BiGGxref);
I=find(a);
ihuman.rxnBiGGDB2MNX(I)=MNXrefRxns.BiGGMNXid(b(I));
% Manually fix two elements:
ind=find(contains(ihuman.HMR2BiGG,'r0706'));
ihuman.rxnBiGGDB2MNX(ind)={'MNXR105369'};
numel(find(~cellfun(@isempty,ihuman.rxnBiGGDB2MNX)))  % ans = 4519/4671

% From KEGG to MNXref
ihuman.rxnKEGG2MNX=cell(num,1);
ihuman.rxnKEGG2MNX(:,1)={''};
[a, b]=ismember(ihuman.rxnKEGGID,MNXrefRxns.KEGGxref);
I=find(a);
ihuman.rxnKEGG2MNX(I)=MNXrefRxns.KEGGMNXid(b(I));
numel(find(~cellfun(@isempty,ihuman.rxnKEGG2MNX)))  % ans = 1752/1767

% From Reactome to MNXref
%=========================
% Add Reactome stable ids based on below cross-reference file that was
% downlaoded from Reactome website (https://reactome.org/download-data)
T=readtable('reactome_stable_ids.txt','HeaderLines',1,'Delimiter','tab');
ReactomID=table2struct(T,'ToScalar',true);
% Expand this structure (Stable_ID and old_identifier_s_) to 1-to-1 format
ReactomID.stableID={};
ReactomID.oldID={};
% This loop takes a lot of time
for i=1:numel(ReactomID.Stable_ID)
	if ~isempty(ReactomID.old_identifier_s_{i})
		oldID=transpose(strsplit(ReactomID.old_identifier_s_{i},','));
		stableID=cell(numel(oldID),1);
		stableID(:,1)={ReactomID.Stable_ID{i}};
		ReactomID.stableID=[ReactomID.stableID;stableID];
		ReactomID.oldID=[ReactomID.oldID;oldID];
	end
end
save('ReactomID.mat','ReactomID');  %2018-01-19
%=========================
load('ReactomID.mat','ReactomID');  %2018-2-9
% Map to new Stable IDs and add them to model
ihuman.rxnREACTOMEStableID=cell(num,1);
ihuman.rxnREACTOMEStableID(:,1)={''};
[a, b]=ismember(ihuman.rxnREACTOMEID,ReactomID.oldID);
I=find(a);
ihuman.rxnREACTOMEStableID(I)=ReactomID.stableID(b(I));
numel(find(~cellfun(@isempty,ihuman.rxnREACTOMEStableID)))  % ans = 216/217
% Because cannot associate REACT_22293 to stable Reactom id
% Then associate new Reactome ids to MNXref
ihuman.rxnReactome2MNX=cell(num,1);
ihuman.rxnReactome2MNX(:,1)={''};
[a, b]=ismember(ihuman.rxnREACTOMEStableID,MNXrefRxns.Reactomexref);
I=find(a);
ihuman.rxnReactome2MNX(I)=MNXrefRxns.ReactomeMNXid(b(I));
numel(find(~cellfun(@isempty,ihuman.rxnReactome2MNX)))  % ans = 210/216

%save('ihuman2MNX.mat','ihuman');  %2018-2-9

% Some statistics for unmatched MNX associations
% First conduct a comparison between rxnKEGG2MNX and rxnReactome2MNX
sharedIndex=intersect(find(~cellfun(@isempty,ihuman.rxnReactome2MNX)),find(~cellfun(@isempty,ihuman.rxnKEGG2MNX)));
isequal(ihuman.rxnReactome2MNX(sharedIndex),ihuman.rxnKEGG2MNX(sharedIndex))  %ans = 0
numel(find(~cellfun(@isequal,ihuman.rxnKEGG2MNX(sharedIndex),ihuman.rxnReactome2MNX(sharedIndex))))
% ans = 24, there are 24 conflicting pairs
% Then conduct a comparison between rxnKEGG2MNX and rxnBiGGDB2MNX
sharedIndex=intersect(find(~cellfun(@isempty,ihuman.rxnKEGG2MNX)),find(~cellfun(@isempty,ihuman.rxnBiGGDB2MNX)));
isequal(ihuman.rxnKEGG2MNX(sharedIndex),ihuman.rxnBiGGDB2MNX(sharedIndex))  %ans = 0
numel(find(~cellfun(@isequal,ihuman.rxnKEGG2MNX(sharedIndex),ihuman.rxnBiGGDB2MNX(sharedIndex))))  
% ans = 276, there are 276 conflicting pairs
% Then conduct a comparison between rxnReactome2MNX and rxnBiGGDB2MNX
sharedIndex=intersect(find(~cellfun(@isempty,ihuman.rxnReactome2MNX)),find(~cellfun(@isempty,ihuman.rxnBiGGDB2MNX)));
isequal(ihuman.rxnReactome2MNX(sharedIndex),ihuman.rxnBiGGDB2MNX(sharedIndex))  %ans = 0
numel(find(~cellfun(@isequal,ihuman.rxnReactome2MNX(sharedIndex),ihuman.rxnBiGGDB2MNX(sharedIndex))))
% ans = 32, there are 32 conflicting pairs

% Unify KEGG, Reactome and BiGG (including EHMN and HepatoNet1) associations toward MNX
% If one reaction were associated to several MNX ids, keep them all

rxnMNXID=cell(num,1);
rxnMNXID(:,1)={''};
% Loop through all reactions
count=0;
for i=1:num
		% check out if they are all non-empty
		if ~isempty(ihuman.rxnReactome2MNX{i}) || ~isempty(ihuman.rxnBiGGDB2MNX{i}) || ~isempty(ihuman.rxnKEGG2MNX{i})
				rxnMNXID{i}=strcat(ihuman.rxnReactome2MNX{i},';',ihuman.rxnBiGGDB2MNX{i},';',ihuman.rxnKEGG2MNX{i});
				rxnMNXID{i}=unique(strsplit(rxnMNXID{i},';'));             %convert from string to cell array
				rxnMNXID{i}=rxnMNXID{i}(~cellfun('isempty',rxnMNXID{i}));  %remove empty elements
				if numel(rxnMNXID{i})>1
						count=count+1;
				end
		end	
end
%count=303
ihuman.rxnMNXID=rxnMNXID;
%---There are 303 rxns with multiple MNXref associations (10 rxns with three different assoc)

numel(find(~cellfun(@isempty,ihuman.rxnMNXID)))
% ans = 5588 with single or identical association

save('ihumanRxns2MNX.mat','ihuman');  %2018-05-21
