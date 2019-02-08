function [newModel,keptGenes] = removeLowScoreGenes(model,geneScores,complexScoring)
%removeLowScoreGenes  Remove low-scoring genes from model.
%
%   This function removes genes from a model based on their scores, a step
%   used by the tINIT package. The function recognizes and differentiates
%   between isozymes and subunits of an enzyme complex. Genes are removed
%   from each grRule, subject to the following conditions:
%       1) At least one gene must remain associated with the reaction
%       2) Genes involved in a complex (joined by ANDs) are not removed
%
%
% Usage:
%
%   [newModel,keptGenes] = removeLowScoreGenes(model,geneScores);
%
% Inputs:
%
%   model           Model structure from which genes are to be removed.
%
%   geneScores      A vector of scores associated with the model genes.
%                   Genes with a positive score will remain in the model,
%                   whereas genes with a negative score will try to be
%                   removed.
%                   
%                   If all genes associated with a reaction have a negative
%                   score, then the least-negative gene will remain; if 
%                   there is a tie, one will be selected at random.
%
%                   If a negative-scoring gene is a subunit in a complex, 
%                   it will not be removed; however, the entire complex may
%                   be removed. See the following example cases:
%
%                    Original: G1 or (G2 and G3 and G4)
%                    Negative scores: G1, G2
%                    New: (G2 and G3 and G4)
%
%                    Original: G1 or (G2 and G3) or (G4 and G5)
%                    Negative 
%
%                       (G2 and G3 and G4)
%
%                   However, if 
%
%
%   complexScoring  Method used to calculate the score of an enzyme complex
%                   from its subunit gene scores: 'min', 'max', 'median',
%                   or 'average'.
%                   (Opt, default 'min').
%
%   keepComplexes   If true (default), remGenes that are part of an enzyme 
%                   complex will not be removed from the model.
%                   For example, if G1 and G2 are to be removed from the
%                   following grRule:
%                   
%                       G1 or (G2 and G3 and G4)
%
%                   the resulting rule would be:
%
%                       G2 and G3 and G4
%
%                   If false, remGenes will be removed from grRules
%                   regardless of their participation in an enzyme complex;
%                   however, the remaining members of the complex will not
%                   be affected. For example, if G1 and G2 are to be
%                   removed from the following grRule:
%
%                       G1 or (G2 and G3 and G4)
%
%                   the resulting rule would be:
%
%                       G3 and G4
%
% Outputs:
%
%   newModel        Model with updated gene, grRules, and rxnGeneMat fields
%                   after attempted gene removals.
%
%   keptGenes       A list of remGenes that the function attempted to
%                   remove from the model, but were kept in at least one
%                   grRule due to its participation in an enzyme complex.
%
%
% Jonathan Robinson, 2019-02-07
%



% NEW PLAN: USE THE ALGORITHM TO BREAK RULES INTO CHUNKS. IF A CHUNK
% CONTAINS AND RELATIONSHIPS, DON'T REMOVE ANYTHING. IF A CHUNK CONTAINS OR
% RELATIONSHIPS, REMOVE ANY NEGATIVE-SCORING GENES, BUT KEEP AT LEAST ONE
% GENE. WILL NEED TO SCORE CHUNKS AS IN THE NEW SCORING FUNCTION.








if nargin < 3
    keepComplexes = true;
end

if ~iscell(remGenes)
    remGenes = {remGenes};
end

% get list of all genes in grRules
genes = getGenesFromGrRules(model.grRules);

% check if any genes to be removed do not exist in grRules
nonExistGenes = setdiff(remGenes,genes);
if ~isempty(nonExistGenes)
    fprintf('\nWARNING: The following genes to be removed were not found in grRules:\n');
    fprintf('\t%s\n',nonExistGenes{1:min([10,numel(nonExistGenes)])});
    if numel(nonExistGenes) > 10
        fprintf('\t... plus %u additional genes.\n',numel(nonExistGenes) - 10);
    end
    fprintf('\n');
end

if ( keepComplexes )
    
    % convert logical operators to symbols
    grRules = regexprep(model.grRules,' and ',' & ');
    grRules = regexprep(grRules,' or ',' | ');
    
    
    
    
else
    % generate new list of genes, where those to be removed are replaced
    % with an empty string {''}
    newGenes = genes;
    newGenes(ismember(newGenes,remGenes)) = {''};
    
    % remove genes from grRules
    model.grRules = translateGeneRules(model.grRules,[genes,newGenes]);
end

% regenerate other gene fields
[genes,rxnGeneMat] = getGenesFromGrRules(model.grRules);
model.genes = genes;
model.rxnGeneMat = rxnGeneMat;

% assign outputs
newModel = model;
keptGenes = remGenes(ismember(remGenes,newModel.genes));











