%
%   FILE NAME:    verifyRxnAssoc2MNX.m
% 
%   PURPOSE: 1. Generate a reduced model by removing exchange reactions and
%               cross-compartment transport reactions, as well as merging
%               duplicate reactions from multiple compartments into single;
%            2. Use the merged model for rxn association to MetaNetX and
%               further manual curation;
%


% Detect duplicated reactions that occurred in multiple comparments
% by using the updated mergeCompartments function in RAVEN.
load('ihumanRxns2MNX.mat');   % load MNX association by external ids
[mergedModel, deletedRxns, duplicateRxns]=mergeCompartments(ihuman,0,0,0);
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


% Get the reaction equations
equationStrings=constructEquations(ihuman,ihuman.rxns,1,1,1);
ihuman.constructedEquations=equationStrings;
%ihuman.constructedEquations=regexprep(equationStrings,'\[\w\]','');  % Clear up the compartment id

% Fetch the corresponding MNX equations
load('MNXRxns.mat');  % Load MNX reactions
ihuman.MNXequations=cell(num,1);
ihuman.MNXequations(:)={''};
[a, b]=ismember(ihuman.rxnMNXID,MNXRxns.MNX_ID);
I=find(a);
ihuman.MNXequations(I)=MNXRxns.Description(b(I));


% load new MNX association through the updated BiGG DB % 2018-05-22
load('ihumanRxns2MNX.mat');

% Regenerate the nested array of MNX association
load('mergedModel.mat');
mergedModel.oldAssocMNXID=mergedModel.rxnAssocMNXID;
mergedModel.rxnAssocMNXID(:)={''};
for i=1:numel(mergedModel.rxns)
		if isempty(mergedModel.duplicateRxns{i})
				hit=find(strcmp(mergedModel.rxns{i},ihuman.rxns));
				mergedModel.rxnAssocMNXID{i}=ihuman.rxnMNXID{hit};
		else
				rxns=[mergedModel.rxns{i};transpose(strsplit(mergedModel.duplicateRxns{i},';'))];
				for j=1:numel(rxns)
						hit=find(strcmp(rxns{j},ihuman.rxns));
						if ~isempty(ihuman.rxnMNXID{hit})
								if isempty(mergedModel.rxnAssocMNXID{i})
										mergedModel.rxnAssocMNXID{i}=ihuman.rxnMNXID{hit};
								else
										mergedModel.rxnAssocMNXID{i}=[mergedModel.rxnAssocMNXID{i},ihuman.rxnMNXID{hit}];
								end
						end
				end
				
		end
		
		if ~isempty(mergedModel.rxnAssocMNXID{i})
				mergedModel.rxnAssocMNXID{i}=unique(mergedModel.rxnAssocMNXID{i});
		end
end

save('mergedModel.mat','mergedModel');   % 2018-05-22
