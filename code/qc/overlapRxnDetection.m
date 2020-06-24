function overlapRxns=overlapRxnDetection(modelA, modelB)
% overlapRxnDetection
%
%   Detect overlap reactions between two models
%
%   overlapRxns   cell array of reaction pairs between two query models
%
%   Usage: overlapRxns=overlapRxnDetection(modelA, modelB)
%


% Handle the input
if nargin<2
    disp('Missing input model');
    return;
end

% Check the correctness of model structures

% Remove duplicate reactions according to identifiers that can be easily
% resloved elsewhere, and thus removed from here for simplicity
[a, b]=ismember(modelB.rxns,modelA.rxns);
rxnToRemove=find(a);
reducedB=removeReactions(modelB,rxnToRemove,1,1,1);

% Make sure there is no duplicate reactions within the enqury models
test=detectDuplicateRxns(modelA,0);
if ~isempty(test.group)
    disp('Found duplicate reactions in modelA');
    return;
end
test=detectDuplicateRxns(reducedB,0);
if ~isempty(test.group)
    disp('Found duplicate reactions in modelB');
    return;
end

% Merge two models and detect overlap reactions
if ~isfield(reducedB, 'id')
		reducedB.id='reducedB';
end
if ~isfield(modelA, 'id')
		reducedB.id='modelA';
end
mergedModel=mergeModels({modelA reducedB});

% Ignore the reaction direction, other options may be considered later
overlap=detectDuplicateRxns(mergedModel,0);

% Output overlap reaction pairs
index=find(cellfun(@numel, overlap.group)==2);
overlapRxns.rxnModelA={};
overlapRxns.rxnModelB={};
for i=1:length(index)
		m=index(i);
		if ismember(overlap.group{m}{1},modelA.rxns) && ismember(overlap.group{m}{2},reducedB.rxns)
				overlapRxns.rxnModelA=[overlapRxns.rxnModelA;overlap.group{m}{1}];
				overlapRxns.rxnModelB=[overlapRxns.rxnModelB;overlap.group{m}{2}];
		end
end

end
