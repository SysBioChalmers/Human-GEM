function [draftModel, removedRxns] = getModelFromOrthology(templateModel,orthologPairs)
%getModelFromOrthology  
%   Constructs a draft model based on a template model and provided gene
%   orthology information between the query and template organisms
%
% Usage:
%
%   [draftModel, removedRxns] = getModelFromOrthology(templateModel,orthologPairs)
%
% Inputs:
%
%   templateModel    a model structures to be used as a template
%
%
%   orthologPairs    an Nx2 cell array of the ortholog pairs, where the
%                    first column contains gene IDs from the reference
%                    organism, and the second includes gene IDs of the
%                    query organism
%
% Outputs:
%
%   draftModel       a model structure for the new organism
%
%   removedRxns      an array of removed rxns due to the missing of
%                    orthology information
%
%

% handle input arguments
if nargin < 2
    error('Missing input.');
end


% pre-process the template model

% remove non-standard fields, if any
fieldsToRemove = intersect({'rxnFrom','metFrom'}, fieldnames(templateModel));
templateModel = rmfield(templateModel, fieldsToRemove);

% clean metadata fields
templateModel.id = '';
templateModel.description = '';
templateModel.version = '';
templateModel.annotation = structfun(@(x) '',templateModel.annotation,'UniformOutput',0);
templateModel.annotation.defaultLB = -1000;
templateModel.annotation.defaultUB = 1000;

% find the index of non-empty grRules before replacing genes
preNonEmptyRuleInd = find(~cellfun(@isempty, templateModel.grRules));


% Replace genes according to the mapped ortholog pairs, which should be in
% the defined format (an Nx2 cell array)
draftModel = templateModel;
[grRules,genes,rxnGeneMat] = replaceGrRules(draftModel.grRules,orthologPairs);


% Update with modified gene fields
draftModel.grRules    = grRules;
draftModel.genes      = genes;
draftModel.rxnGeneMat = rxnGeneMat;


% find the index of non-empty grRules after replacing genes
postNonEmptyRuleInd = find(~cellfun(@isempty, draftModel.grRules));


% remove the rxns whose grRules become empty after replacement of orthologs
removedRxns = setdiff(preNonEmptyRuleInd, postNonEmptyRuleInd);
if any(removedRxns)
    draftModel = removeReactions(draftModel, removedRxns, true, true, true);
end


