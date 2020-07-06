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
%   orthologPairs    an Nx2 cell array of the ortholog pairs, , where the
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

% clean template model - to be refined
% now it's a list non-standard fields. Later this part could be replace
% with read in a list core fields, and remove the others not in this list
fieldsToRemove = {'rxnFrom','metFrom'};
templateModel = rmfield(templateModel, fieldsToRemove);

% clear generic fields
templateModel.id = '';
templateModel.description = '';
templateModel.version = '';
templateModel.annotation = '';


% a possible section to remove unnecessary subSystems

% pre-check of grRules
preNonEmptyRuleInd = find(~cellfun(@isempty, templateModel.grRules));


% Replace genes according to the mapped orthologs, that should be in
% designed format
draftModel = templateModel;
[grRules,genes,rxnGeneMat] = replaceGrRules(draftModel.grRules,orthologPairs);


% Update with modified fields
draftModel.grRules    = grRules;
draftModel.genes      = genes;
draftModel.rxnGeneMat = rxnGeneMat;


% post-check of grRules
postNonEmptyRuleInd = find(~cellfun(@isempty, draftModel.grRules));


% find and remove the rxns with empty grRules after ortholog replacement
if ~isequal(preNonEmptyRuleInd, postNonEmptyRuleInd) &&...
    all(ismember(postNonEmptyRuleInd, preNonEmptyRuleInd))

    removedRxns = setdiff(preNonEmptyRuleInd, postNonEmptyRuleInd);
    draftModel = removeReactions(draftModel, removedRxns, true, true, true);
end


