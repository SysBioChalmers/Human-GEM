%
%   FILE NAME:    integrateRxnAssocViaIDsAndMets.m
% 
%   PURPOSE: Integrate rxn associations from mets and external rxn ids
%


load('MNXRxns.mat');          % load MNX reactions
load('mergedModel.mat');      % load mergedModel
load('mapRxnResults.mat');    % load automatic rxn association via mets

% Treat the MNX assoc overlapped by rxn-ids and via mets as confirmed
mergedModel.confirmedMNXID=cell(numel(mergedModel.rxns),1);
mergedModel.confirmedMNXID(:)={''};
mergedModel.confirmedFilteredMNXID=cell(numel(mergedModel.rxns),1);
mergedModel.confirmedFilteredMNXID(:)={''};
mergedModel.unResolved=zeros(numel(mergedModel.rxns),1);

numUnresloved=0;
numConfirmed=0;
numToCheck=0;
for i=1:length(mergedModel.rxns)
		if ~isempty(mergedModel.rxnAssocMNXID{i}) || ~isempty(results.rxnMNXID{i})
				overlap=intersect(mergedModel.rxnAssocMNXID{i},results.rxnMNXID{i});

				if ~isempty(overlap)
						numConfirmed=numConfirmed+1;
						mergedModel.confirmedMNXID{i}=overlap;
						mergedModel.confirmedFilteredMNXID{i}=filterBalancedRxns(overlap,MNXRxns);
				else
						numToCheck=numToCheck+1;
				end
		
		else  % These rxns have no association both by rxn-ids and mets
				numUnresloved=numUnresloved+1;
				mergedModel.unResolved(i)=1;
		end
end
%numConfirmed=2506;
%numUnresloved=726;
%numToCheck=739;

% Treat the cases need to be checked
filledRxnAssocViaMets.rxnHMRID=cell(1000,1);
filledRxnAssocViaMets.rxnMNXID=cell(1000,1);
filledRxnAssocViaMets.rxnFilteredMNXID=cell(1000,1);
filledRxnAssocViaMets.notes=cell(1000,1);
needToCheck.rxnHMRID=cell(1000,1);
needToCheck.rxnAssocMNXID=cell(1000,1);
needToCheck.rxnViaMetsMNXID=cell(1000,1);
needToCheck.rxnViaMetsFilteredMNXID=cell(1000,1);
needToCheck.notes=cell(1000,1);

numToCheck=0;
numFilledViaMets=0;
for i=1:numel(mergedModel.rxns)
		if ~isempty(mergedModel.rxnAssocMNXID{i}) || ~isempty(results.rxnMNXID{i})
				overlap=intersect(mergedModel.rxnAssocMNXID{i},results.rxnMNXID{i});
				if isempty(overlap)
						if ~isempty(mergedModel.rxnAssocMNXID{i})
								numToCheck=numToCheck+1;
								needToCheck.rxnHMRID{numToCheck}=mergedModel.rxns{i};
								needToCheck.rxnAssocMNXID{numToCheck}=mergedModel.rxnAssocMNXID{i};
								needToCheck.notes{numToCheck}=results.notes{i};

								% Both have assoc, but non-overlap (202) could be caused by wrong assoc by rxn-ids
								if ~isempty(results.rxnMNXID{i})
										needToCheck.rxnViaMetsMNXID{numToCheck}=results.rxnMNXID{i};
										needToCheck.rxnViaMetsFilteredMNXID{numToCheck}=filterBalancedRxns(results.rxnMNXID{i},MNXRxns);

								% No assoc via mets, keep origianl ones (248) Check the manual curation 
								else
										needToCheck.rxnViaMetsMNXID{numToCheck}='';
										needToCheck.rxnViaMetsFilteredMNXID{numToCheck}='';						
								end
								
						% No rxn assoc, fill by assoc via mets (prefer to true type) (296)
						else
								numFilledViaMets=numFilledViaMets+1;
								filledRxnAssocViaMets.rxnHMRID{numFilledViaMets}=mergedModel.rxns{i};
								filledRxnAssocViaMets.rxnMNXID{numFilledViaMets}=results.rxnMNXID{i};
								filledRxnAssocViaMets.rxnFilteredMNXID{numFilledViaMets}=filterBalancedRxns(results.rxnMNXID{i},MNXRxns);
								filledRxnAssocViaMets.notes{numFilledViaMets}=results.notes{i};
						end
				end
		end
end
%numToCheck=444;
%numFilledViaMets=295;

needToCheck.rxnHMRID=needToCheck.rxnHMRID(1:numToCheck);
needToCheck.rxnAssocMNXID=needToCheck.rxnAssocMNXID(1:numToCheck);
needToCheck.rxnViaMetsMNXID=needToCheck.rxnViaMetsMNXID(1:numToCheck);
needToCheck.rxnViaMetsFilteredMNXID=needToCheck.rxnViaMetsFilteredMNXID(1:numToCheck);
needToCheck.notes=needToCheck.notes(1:numToCheck);

filledRxnAssocViaMets.rxnHMRID=filledRxnAssocViaMets.rxnHMRID(1:numFilledViaMets);
filledRxnAssocViaMets.rxnMNXID=filledRxnAssocViaMets.rxnMNXID(1:numFilledViaMets);
filledRxnAssocViaMets.rxnFilteredMNXID=filledRxnAssocViaMets.rxnFilteredMNXID(1:numFilledViaMets);
filledRxnAssocViaMets.notes=filledRxnAssocViaMets.notes(1:numFilledViaMets);

%save('mergedModel.mat','mergedModel');
%save('filledRxnAssocViaMets_20180522.mat','filledRxnAssocViaMets');
%save('needToCheck_20180522.mat','needToCheck');

%% sub functions
%Filter MNX rxns according to their balance status
function filteredMNXrxns = filterBalancedRxns(MNXrxns,mnx)

		countNum=numel(MNXrxns);
		if countNum==1  % there is only one MNX id, skip balance check
				filteredMNXrxns=MNXrxns;
		elseif countNum>1   % there are several MNX ids
				[a, b]=ismember(MNXrxns,mnx.MNX_ID);
				% The priority order of Balance status: true > NA > ambiguous
				if find(strcmp('true',mnx.Balance(b)));
						trueHits=find(strcmp('true',mnx.Balance(b)));
						filteredMNXrxns=MNXrxns(trueHits);
				elseif find(strcmp('NA',mnx.Balance(b)))
						NAHits=find(strcmp('NA',mnx.Balance(b)));
						filteredMNXrxns=MNXrxns(NAHits);
				elseif find(strcmp('ambiguous',mnx.Balance(b)))
						amHits=find(strcmp('ambiguous',mnx.Balance(b)));
						filteredMNXrxns=MNXrxns(amHits);
				else
						filteredMNXrxns=MNXrxns;
				end
		end
end

