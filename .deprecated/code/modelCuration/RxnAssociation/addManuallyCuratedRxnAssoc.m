%
%   FILE NAME:    addManuallyCuratedRxnAssoc.m
% 
%   PURPOSE: Incorporate manual curation results of rxn association 
%            with MNXref ids to the mergedModel
%


% Move to the target folder of manual curation files
cd('/Users/haowa/Box Sync/HMR3/Curation Files/manualRxnCuration');

% Load the manual curation results
[~, textData1]=xlsread('conflictAssocByIDsAndViaMets.xlsx','ManualCuration_2018May');
[~, textData2]=xlsread('filledRxnAssocViaMets.xlsx','ManualCuration_2018May');
[~, textData3]=xlsread('rxnAssocOnlyByRxnIDs.xlsx','ManualCuration_2018May');

manuallyCuratedRxnAssoc.rxnHMRID=[textData1(2:end,1);textData2(2:end,1);textData3(2:end,1)];
manuallyCuratedRxnAssoc.rxnMNXID=[textData1(2:end,2);textData2(2:end,2);textData3(2:end,2)];


% Move back to rxnAssoc folder, and save this data structure
%save('manuallyCuratedRxnAssoc.mat','manuallyCuratedRxnAssoc');  % 2018-05-25

% Load the HMR2 model after merging compartments
load('mergedModel.mat');    % 2018-05-25

numel(find(~cellfun(@isempty,mergedModel.confirmedMNXID)))  %ans = 2506
% Add manually curated rxn associations to confirmedMNXID field, some of
% which are mostly not assigned (with empty values)
for i=1:numel(manuallyCuratedRxnAssoc.rxnHMRID)
		% No need to check the duplicateRxns field, because all reactions
		% for curation were selected from the rxns field
		[~, index]=ismember(manuallyCuratedRxnAssoc.rxnHMRID{i},mergedModel.rxns);
		if isempty(mergedModel.confirmedMNXID{index})
				mergedModel.confirmedMNXID{index}{1}=manuallyCuratedRxnAssoc.rxnMNXID{i};
		else
				if ~ismember(manuallyCuratedRxnAssoc.rxnMNXID{i},mergedModel.confirmedMNXID{index})
						mergedModel.confirmedMNXID{index}=[mergedModel.confirmedMNXID{index};manuallyCuratedRxnAssoc.rxnMNXID{i}];
				end
		end
end
numel(find(~cellfun(@isempty,mergedModel.confirmedMNXID)))  %ans = 2988
% A total of 492 (10 already have rxnAssoc) reactions in mergedModel were newly associated to MNX

% Re-processing the filed 'confirmedFilteredMNXID'
mergedModel.confirmedFilteredMNXID=filterBalancedRxns(mergedModel.confirmedMNXID);

% Save the mergedModel
%save('mergedModel.mat','mergedModel');  % 2018-05-25

%% sub functions
%Filter MNX rxns according to their balance status
function filteredMNXrxns = filterBalancedRxns(rxnList)

load('MNXRxns.mat');          % load MNX reactions
filteredMNXrxns=cell(numel(rxnList),1);
filteredMNXrxns(:)={''};

for i=1:numel(rxnList)
		if ~isempty(rxnList{i})
				countNum=numel(rxnList{i});
				if countNum==1  % there is only one MNX id, skip balance check
						filteredMNXrxns{i}=rxnList{i};
				elseif countNum>1   % there are several MNX ids
						[a, b]=ismember(rxnList{i},MNXRxns.MNX_ID);
						% The priority order of Balance status: true > NA > ambiguous
						if find(strcmp('true',MNXRxns.Balance(b)));
								trueHits=find(strcmp('true',MNXRxns.Balance(b)));
								filteredMNXrxns{i}=rxnList{i}(trueHits);
						elseif find(strcmp('NA',MNXRxns.Balance(b)))
								NAHits=find(strcmp('NA',MNXRxns.Balance(b)));
								filteredMNXrxns{i}=rxnList{i}(NAHits);
						elseif find(strcmp('ambiguous',MNXRxns.Balance(b)))
								amHits=find(strcmp('ambiguous',MNXRxns.Balance(b)));
								filteredMNXrxns{i}=rxnList{i}(amHits);
						else
								filteredMNXrxns{i}=rxnList{i};
						end
				end
		end
end

end
