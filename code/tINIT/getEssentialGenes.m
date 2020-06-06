function [essentialGenes, essentialGenesIndexes] = getEssentialGenes(model,ignoreGenes,fluxTol)
% getEssentialGenes
%   Calculate the essential genes for a model to be solvable
%
%   model                   a model structure
%   ignoreGenes             cell array of gene IDs which should not be
%                           checked (opt, default {})
%   fluxTol                 absolute flux value below which a reaction is
%                           considered "off" (opt, default 10^-8)
%
%   essentialGenes          cell array with the IDs of the essential genes
%   essentialGenesIndexes   vector with the indexes of the essential genes
%
%   Essential genes are those which, when deleted, result in an infeasible
%   problem. This function is derived from the similar "getEssentialRxns"
%   function.
%
%   NOTE: To use this function, the model must have a "rules" field. This
%         field can be added using the "generateRules" function from the
%         COBRA toolbox.
%
%   Usage: [essentialGenes, essentialGenesIndexes] = getEssentialGenes(model,ignoreGenes)
%


if nargin < 2
    ignoreGenes = {};
end
if nargin < 3
    fluxTol = 10^-8;
end

% verify that the model has a "rules" field
if ~isfield(model,'rules')
    error('The model must have a "rules" field. Use e.g. the COBRA "generateRules" function to add this field.');
end

% ensure we don't try to optimize for something
model.c = zeros(numel(model.rxns),1);

% first check that the problem is solvable
[sol, hsSolOut] = solveLP(model,1);
if sol.stat == -1 || isempty(sol.x)
    EM = 'No feasible solution to the full model';
    dispEM(EM);
end

% Check which reactions carry flux. Only genes associated with at least one 
% of these reactions can be essential.
activeRxnInds = abs(sol.x) > fluxTol;
genesToCheck = setdiff(model.genes(any(model.rxnGeneMat(activeRxnInds,:) == 1, 1)), ignoreGenes);

% delete each gene and determine if the problem becomes unsolvable
essentialGenes = {};
model.rules(cellfun(@isempty,model.rules)) = {'true'};  % empty rules are always true
[~,geneIndsToCheck] = ismember(genesToCheck,model.genes);
for i = 1:numel(genesToCheck)
    x = true(size(model.genes));  % generate boolean vector representing presence/absence of all genes
    x(geneIndsToCheck(i)) = false;  % indicate deletion of the current gene
    delRxns = model.rxns(~cellfun(@eval,model.rules));  % determine which reactions require the deleted gene
    
    if ~isempty(delRxns)
        sol = solveLP(setParam(model,'eq',delRxns,0),0,[],hsSolOut);
        if sol.stat == -1 || isempty(sol.x)
            essentialGenes = [essentialGenes; genesToCheck(i)];
        end
    end
end

[~, essentialGenesIndexes] = ismember(essentialGenes,model.genes);
end
