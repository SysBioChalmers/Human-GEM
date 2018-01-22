%
%   FILE NAME:  rxnAssocByMNX.m  
% 
%   DATE CREATED: 2018-01-06
%
%	
%   PROGRAMMER:   Hao Wang
%                 Department of Biology and Biological Engineering
%                 Chalmers University of Technology
% 
% 
%   PURPOSE: HMR2 reaction association to MNXref based on provided
%            external databasse identifiers
%

% Load HMR model with BiGG association
load('ihuman2BiGG.mat');

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
numel(find(~cellfun(@isempty,ihuman.HMR2BiGG(1:5127))))  % ans = 2531


%===Reaction mapping through MNXref database

% Load NNX reaction references version 3.0
load('MNXrefRxns.mat');

% From KEGG to MNXref
ihuman.KEGG2MNX=cell(num,1);
ihuman.KEGG2MNX(:,1)={''};
[a, b]=ismember(ihuman.rxnKEGGID,MNXrefRxns.KEGGxref);
I=find(a);
ihuman.KEGG2MNX(I)=MNXrefRxns.KEGGMNXid(b(I));
numel(find(~cellfun(@isempty,ihuman.KEGG2MNX)))  % ans = 1730/1745

% From Reactome to MNXref

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
% Map to new Stable IDs and add them to model
ihuman.rxnREACTOMEStableID=cell(num,1);
ihuman.rxnREACTOMEStableID(:,1)={''};
[a, b]=ismember(ihuman.rxnREACTOMEID,ReactomID.oldID);
I=find(a);
ihuman.rxnREACTOMEStableID(I)=ReactomID.stableID(b(I));
numel(find(~cellfun(@isempty,ihuman.rxnREACTOMEStableID)))  % ans = 214/215
% Because cannot associate REACT_22293 to stable Reactom id
% Then associate new Reactome ids to MNXref
ihuman.Reactome2MNX=cell(num,1);
ihuman.Reactome2MNX(:,1)={''};
[a, b]=ismember(ihuman.rxnREACTOMEStableID,MNXrefRxns.Reactomexref);
I=find(a);
ihuman.Reactome2MNX(I)=MNXrefRxns.ReactomeMNXid(b(I));
numel(find(~cellfun(@isempty,ihuman.Reactome2MNX)))  % ans = 208/214

% From BiGG to MNXref
ihuman.BiGG2MNX=cell(num,1);
ihuman.BiGG2MNX(:,1)={''};
[a, b]=ismember(ihuman.rxnBiGGID,MNXrefRxns.BiGGxref);
I=find(a);
ihuman.BiGG2MNX(I)=MNXrefRxns.BiGGMNXid(b(I));
numel(find(~cellfun(@isempty,ihuman.BiGG2MNX)))  % ans = 2281/3061
% Compare BiGG2MNX and HMR2BiGG2MNX, should be identical
sharedIndex=intersect(find(~cellfun(@isempty,ihuman.BiGG2MNX)),find(~cellfun(@isempty,ihuman.HMR2BiGG2MNX)));
isequal(ihuman.BiGG2MNX(sharedIndex),ihuman.HMR2BiGG2MNX(sharedIndex))  %ans = 1
% All BiGG2MNX elements are identical to the corresponding ones
% in HMR2BiGG2MNX, so only HMR2BiGG2MNX will be used for combining

save('ihuman2MNX.mat','ihuman');  %2018-01-19

% Some statistics for unmatched MNX associations
% First conduct a comparison between KEGG2MNX and Reactome2MNX
sharedIndex=intersect(find(~cellfun(@isempty,ihuman.Reactome2MNX)),find(~cellfun(@isempty,ihuman.KEGG2MNX)));
isequal(ihuman.Reactome2MNX(sharedIndex),ihuman.KEGG2MNX(sharedIndex))  %ans = 0
numel(find(~cellfun(@isequal,ihuman.KEGG2MNX(sharedIndex),ihuman.Reactome2MNX(sharedIndex))))  
% ans = 24, there are 24 conflicting pairs
% Then conduct a comparison between KEGG2MNX and HMR2BiGG2MNX
sharedIndex=intersect(find(~cellfun(@isempty,ihuman.KEGG2MNX)),find(~cellfun(@isempty,ihuman.HMR2BiGG2MNX)));
isequal(ihuman.KEGG2MNX(sharedIndex),ihuman.HMR2BiGG2MNX(sharedIndex))  %ans = 0
numel(find(~cellfun(@isequal,ihuman.KEGG2MNX(sharedIndex),ihuman.HMR2BiGG2MNX(sharedIndex))))  
% ans = 252, there are 252 conflicting pairs
% Then conduct a comparison between Reactome2MNX and HMR2BiGG2MNX
sharedIndex=intersect(find(~cellfun(@isempty,ihuman.Reactome2MNX)),find(~cellfun(@isempty,ihuman.HMR2BiGG2MNX)));
isequal(ihuman.Reactome2MNX(sharedIndex),ihuman.HMR2BiGG2MNX(sharedIndex))  %ans = 0
numel(find(~cellfun(@isequal,ihuman.Reactome2MNX(sharedIndex),ihuman.HMR2BiGG2MNX(sharedIndex))))  
% ans = 22, there are 22 conflicting pairs

% Unify KEGG, Reactome and BiGG (including EHMN and HepatoNet1) associations toward MNX
% 1. The only or identical MNX ids will be directly applied
% 2. If associated MNX ids are different, leave them for manual curation
% 3. Add 'conflictMNXAssoc' field to indicate these conflicting match
ihuman.rxnMNXID=cell(num,1);
ihuman.rxnMNXID(:,1)={''};
ihuman.conflictMNXAssoc=zeros(num,1);

% Define output format
outThree='%d, %s: (Reactome)%s-(KEGG)%s-(BiGG)%s\n';
outNoKegg='%d, %s: (Reactome)%s-(BiGG)%s\n';
outNoBigg='%d, %s: (Reactome)%s-(KEGG)%s\n';
outNoReactome='%d, %s: (KEGG)%s-(BiGG)%s\n';
% Loop through all reactions
for i=1:num
	% check out if they are all non-empty
	if ~isempty(ihuman.Reactome2MNX{i}) && ~isempty(ihuman.HMR2BiGG2MNX{i}) && ~isempty(ihuman.KEGG2MNX{i})
		if isequal(ihuman.Reactome2MNX{i},ihuman.HMR2BiGG2MNX{i},ihuman.KEGG2MNX{i})
			ihuman.rxnMNXID{i}=ihuman.Reactome2MNX{i};
		else
			ihuman.conflictMNXAssoc(i,1)=1;
			%fprintf(outThree,i,ihuman.rxns{i},ihuman.Reactome2MNX{i},ihuman.KEGG2MNX{i},ihuman.HMR2BiGG2MNX{i});
		end
	elseif ~isempty(ihuman.Reactome2MNX{i}) && ~isempty(ihuman.HMR2BiGG2MNX{i}) && isempty(ihuman.KEGG2MNX{i})
		if isequal(ihuman.Reactome2MNX{i},ihuman.HMR2BiGG2MNX{i})
			ihuman.rxnMNXID{i}=ihuman.Reactome2MNX{i};
		else
			ihuman.conflictMNXAssoc(i,1)=1;
			%fprintf(outNoKegg,i,ihuman.rxns{i},ihuman.Reactome2MNX{i},ihuman.HMR2BiGG2MNX{i});
		end
	elseif ~isempty(ihuman.Reactome2MNX{i}) && isempty(ihuman.HMR2BiGG2MNX{i}) && ~isempty(ihuman.KEGG2MNX{i})
		if isequal(ihuman.Reactome2MNX{i},ihuman.KEGG2MNX{i})
			ihuman.rxnMNXID{i}=ihuman.Reactome2MNX{i};
		else
			ihuman.conflictMNXAssoc(i,1)=1;
			%fprintf(outNoBigg,i,ihuman.rxns{i},ihuman.Reactome2MNX{i},ihuman.KEGG2MNX{i});
		end
	elseif isempty(ihuman.Reactome2MNX{i}) && ~isempty(ihuman.HMR2BiGG2MNX{i}) && ~isempty(ihuman.KEGG2MNX{i})
		if isequal(ihuman.HMR2BiGG2MNX{i},ihuman.KEGG2MNX{i})
			ihuman.rxnMNXID{i}=ihuman.HMR2BiGG2MNX{i};
		else
			ihuman.conflictMNXAssoc(i,1)=1;
			%fprintf(outNoReactome,i,ihuman.rxns{i},ihuman.KEGG2MNX{i},ihuman.HMR2BiGG2MNX{i});
		end
	elseif isempty(ihuman.Reactome2MNX{i}) && isempty(ihuman.HMR2BiGG2MNX{i}) && ~isempty(ihuman.KEGG2MNX{i})
		ihuman.rxnMNXID{i}=ihuman.KEGG2MNX{i};
	elseif isempty(ihuman.Reactome2MNX{i}) && ~isempty(ihuman.HMR2BiGG2MNX{i}) && isempty(ihuman.KEGG2MNX{i})
		ihuman.rxnMNXID{i}=ihuman.HMR2BiGG2MNX{i};
	elseif ~isempty(ihuman.Reactome2MNX{i}) && isempty(ihuman.HMR2BiGG2MNX{i}) && isempty(ihuman.KEGG2MNX{i})
		ihuman.rxnMNXID{i}=ihuman.Reactome2MNX{i};
	end
end
%---There are 281 rxns with conflicting MNXref associations

numel(find(~cellfun(@isempty,ihuman.rxnMNXID)))  
% ans = 5220 with single or identical association

numel(find(ihuman.conflictMNXAssoc))  
% ans = 281 with conflicting MNX associations

save('ihuman2MNX.mat','ihuman');  %2018-01-22
