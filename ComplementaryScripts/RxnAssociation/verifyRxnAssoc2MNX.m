%
%   FILE NAME:    verifyRxnAssoc2MNX.m
% 
%   DATE CREATED: 2018-03-08
%        
%	
%   PROGRAMMER:   Hao Wang
%                 Department of Biology and Biological Engineering
%                 Chalmers University of Technology
% 
% 
%   PURPOSE: Prepare uniquely associated MNX reactions to HMR for
%            manual curation
%

% Detect duplicated reactions that occurred in multiple comparments
% by using the updated mergeCompartments function in RAVEN.
[mergedModel, deletedRxns, duplicateRxns]=mergeCompartments(metCheck,0,0,0);
mergedModel.duplicateRxns=duplicateRxns;

% Detecting addtional duplications by checking identical columns
% Found two with different reaction direction
%3934: HMR_0688-HMR_5257
%3948: HMR_5257-HMR_0688
%Remove addtional one
mergedModel=removeReactions(mergedModel,{'HMR_5257'},true,true);
mergedModel.duplicateRxns(3948)=[];
save('mergedModel.mat','mergedModel');


%How many uniuqe reactions are from the above 5127
numel(intersect(metCheck.rxns(1:5127),mergedModel.rxns))  %ans= 3906


%Go through these 3906 unique reactions (groups) in mergedModel
%and check reaction association from ihuman2MNX for statistic purpose
%Need to prepare another script to elucidate the detail

% Load the HMR2MNX and mergeModel
load('ihuman2MNX.mat');

singleMNXAssoc=0;
noMNXAssoc=0;
numRxnAssoc=0;
for i=1:3906
		%Determine unique reaction
		if isempty(mergedModel.duplicateRxns{i})
				index=find(strcmp(mergedModel.rxns{i},ihuman.rxns));

				%Check external association
				if ~(ihuman.rxnAssocNum(index)==0)
						numRxnAssoc=numRxnAssoc+1;
				end

				%check MNX association
				if ~isempty(ihuman.rxnMNXID{index})
						singleMNXAssoc=singleMNXAssoc+1;
				elseif ihuman.conflictMNXAssoc(index)==0
						noMNXAssoc=noMNXAssoc+1;
				end

		%Found with duplicate reactions
		else
				%Merge the duplicate reactions into a cell array and loop through
				rxns=[mergedModel.rxns{i}, strsplit(mergedModel.duplicateRxns{i},';')];
				
				%Check external association
				rxnAssoc=false;
				for j=1:numel(rxns)
						index=find(strcmp(rxns{j},ihuman.rxns));
						if ~(ihuman.rxnAssocNum(index)==0)
								rxnAssoc=true;
						end						
				end
				if rxnAssoc
						numRxnAssoc=numRxnAssoc+1;
				end

				%check MNX association
				singleAssoc=false;
				noAssoc=true;
				for j=1:numel(rxns)
						index=find(strcmp(rxns{j},ihuman.rxns));
						if ~isempty(ihuman.rxnMNXID{index})
								singleAssoc=true;
								noAssoc=false;
						elseif ihuman.conflictMNXAssoc(index)==1
								noAssoc=false;
						end
				end
				if singleAssoc
						singleMNXAssoc=singleMNXAssoc+1;
				elseif noAssoc
						noMNXAssoc=noMNXAssoc+1;
				end
		end
end

% numRxnAssoc=3409;
% singleMNXAssoc=2718;
% noMNXAssoc=962;

% Get the reaction equations
equationStrings=constructEquations(ihuman,ihuman.rxns,1,1,1);
ihuman.constructedEquations=equationStrings;
ihuman.constructedEquations=regexprep(equationStrings,'\[\w\]','');  % Clear up the compartment id

% Fetch the corresponding MNX equations
load('MNXRxns.mat');  % Load MNX reactions
ihuman.MNXequations=cell(num,1);
ihuman.MNXequations(:)={''};
[a, b]=ismember(ihuman.rxnMNXID,MNXRxns.MNX_ID);
I=find(a);
ihuman.MNXequations(I)=MNXRxns.Description(b(I));

% Treat the simple cases with single unique MNX association
count=0;
fid=fopen('verifySimpleMNXAssociation.txt','w');
for i=1:3906
		if isempty(mergedModel.duplicateRxns{i})
				[~, index]=ismember(mergedModel.rxns{i},ihuman.rxns);
				if ~isempty(ihuman.rxnMNXID{index})
						count=count+1;
						fprintf(fid,'#%s.\n----------------------------------------\n....
						%s(%s):\t%s\n%s:\t%s\n....
						----------------------------------------n\n',....
						num2str(count),mergedModel.rxns{i},num2str(index),....
						ihuman.constructedEquations{index},ihuman.rxnMNXID{index},....
						ihuman.MNXequations{index});
				end
		end
end
fclose(fid);
% count=2144


% Treat the cases of unique rxn associate to multiple MNX rxns
fid=fopen('verifyUnique2MultipleMNXAssociation.txt','w');
count=0;
for i=1:3906
		if isempty(mergedModel.duplicateRxns{i}) % Unique rxns
				index=find(strcmp(mergedModel.rxns{i},ihuman.rxns));  % Find index in HMR
				if ihuman.conflictMNXAssoc(index)  % With multiple assoc (196)
						count=count+1;
						fprintf(fid,'#%s.\n----------------------------------------\n%s(%s):\t%s\n'....
						,num2str(count),mergedModel.rxns{i},num2str(index),....
						ihuman.constructedEquations{index});
						
						MNXAssoc=strsplit(ihuman.rxnMNXID{index},';');
						for j=1:numel(MNXAssoc)
								MNXID=strsplit(MNXAssoc{j},':');
								I=find(strcmp(MNXID{2},MNXRxns.MNX_ID));
								fprintf(fid,'%s:\t%s\n',MNXAssoc{j},MNXRxns.Description{I});
						end
						fprintf(fid,'----------------------------------------\n\n');								
											
				%elseif isempty(ihuman.rxnMNXID{index})  % No assoc (748)
						
				end
		end
end
fclose(fid);
% count=196
