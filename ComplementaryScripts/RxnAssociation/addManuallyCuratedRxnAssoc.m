%
%   FILE NAME:    addManuallyCuratedRxnAssoc.m
% 
%   DATE CREATED: 2018-05-15
%        
%   PROGRAMMER:   Hao Wang
%                 Department of Biology and Biological Engineering
%                 Chalmers University of Technology
% 
% 
%   PURPOSE: Incorporate manual curation results of rxn association 
%            with MNXref ids to the mergedModel
%

% Move to the target folder of manual curation files
cd('/Users/haowa/Box Sync/HMR3/Curation Files/manualRxnCuration');

% Load the manual curation results
[~, textData1]=xlsread('conflictAssocByIDsAndViaMets.xlsx','ManualCuration_2018May');
[~, textData2]=xlsread('filledRxnAssocViaMets.xlsx','ManualCuration_2018May');

manuallyCuratedRxnAssoc.rxnHMRID=[textData1(2:end,1);textData2(2:end,1)];
manuallyCuratedRxnAssoc.rxnMNXID=[textData1(2:end,2);textData2(2:end,2)];


% Move back to rxnAssoc folder, and save this data structure
save('manuallyCuratedRxnAssoc.mat','manuallyCuratedRxnAssoc');  % 2018-05-15

% Load the HMR2 model after merging compartments
load('mergedModel.mat');

numel(find(~cellfun(@isempty,mergedModel.confirmedMNXID)))  %ans = 2499
% Add manually curated rxn associations to confirmedMNXID field, which
% previously are not assigned (with empty values)
for i=1:numel(manuallyCuratedRxnAssoc.rxnHMRID)
		% No need to check the duplicateRxns field, because all reactions
		% for curation were selected from the rxns field
		[~, index]=ismember(manuallyCuratedRxnAssoc.rxnHMRID{i},mergedModel.rxns);
		if isempty(mergedModel.confirmedMNXID{index})
				mergedModel.confirmedMNXID{index}{1}=manuallyCuratedRxnAssoc.rxnMNXID{i};
		else
				mergedModel.confirmedMNXID{index}=[mergedModel.confirmedMNXID{index};manuallyCuratedRxnAssoc.rxnMNXID{i}];
		end
end
numel(find(~cellfun(@isempty,mergedModel.confirmedMNXID)))  %ans = 2945
% A total of 446 reactions in mergedModel were newly associated to MNX

% Save the mergedModel
save('mergedModel.mat','mergedModel');  % 2018-05-15
