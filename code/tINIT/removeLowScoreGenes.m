function [newModel,remGenes] = removeLowScoreGenes(model,geneScores,isozymeScoring,complexScoring)
%removeLowScoreGenes  Remove low-scoring genes from model.
%
%   This function removes genes from a model based on their scores, a step
%   used by the tINIT package. The function recognizes and differentiates
%   between isozymes and subunits of an enzyme complex. Genes are removed
%   from each grRule, subject to the following conditions:
%       1) At least one gene must remain associated with the reaction
%       2) Genes involved in a complex (joined by ANDs) are not removed
%
% Usage:
%
%   [newModel,remGenes] = removeLowScoreGenes(model,geneScores,complexScoring);
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
%                    Negative: G1, G2
%                         New: G2 and G3 and G4
%
%                    Original: G1 or (G2 and G3) or (G4 and G5)
%                    Negative: G1, G2
%                         New: G4 and G5     [using default complexScoring]
%
%                    Original: (G1 and (G2 or G3) and G4)
%                    Negative: G2
%                         New: (G1 and G3 and G4)
%
%   isozymeScoring  Method used to combine the scores of multiple isozymes;
%                   'min', 'max', 'median', or 'average'. 
%                   (Opt, default 'max')
%
%   complexScoring  Method used to combine the scores of enzyme complex
%                   subunits: 'min', 'max', 'median', or 'average'. 
%                   (Opt, default 'min')
%
% Outputs:
%
%   newModel        Model with updated gene, grRules, and rxnGeneMat fields
%                   after attemping to remove negative-score genes.
%
%   remGenes        A list of negative-score genes that were fully removed
%                   from the model. Negative-score genes that were removed
%                   from some grRules but not all will not be included in
%                   this list.
%


if nargin < 3 || isempty(isozymeScoring)
    isozymeScoring = 'max';
end
if nargin < 4
    complexScoring = 'min';
end

if ~isequal(size(model.genes),size(geneScores))
    error('The dimensions of model genes and geneScores do not match.');
end

% convert logical operators to symbols
grRules = model.grRules;
if any(contains(grRules,{'&','|'}))
    symbolicRules = true;
else
    symbolicRules = false;
end
grRules = regexprep(grRules,' and ',' & ');
grRules = regexprep(grRules,' or ',' | ');

% get unique list of grRules
[uRules,~,rule_ind] = unique(grRules);

% iterate through each rule
newRules = uRules;  %initialize newRules
for i = 1:numel(uRules)
    if isempty(uRules{i}) || ~contains(uRules{i},'|')
        continue
    elseif contains(uRules{i},'&')
        newRules{i} = processComplexRule(uRules{i},model.genes,geneScores,isozymeScoring,complexScoring);
    else
        newRules{i} = processSimpleRule(uRules{i},model.genes,geneScores,isozymeScoring,complexScoring);
    end
end

% re-map unique rules to model
newModel = model;
newModel.grRules = newRules(rule_ind);

% restore original logical operator formatting if it was changed
if ~symbolicRules
    newModel.grRules = regexprep(newModel.grRules, ' \| ', ' or ');
    newModel.grRules = regexprep(newModel.grRules, ' & ', ' and ');
end

% regenerate "genes" and "rxnGeneMat" model fields
[genes,rxnGeneMat] = getGenesFromGrRules(newModel.grRules);
newModel.genes = genes;
newModel.rxnGeneMat = rxnGeneMat;

% update other gene-related fields
remInd = ~ismember(model.genes,newModel.genes);
remGenes = model.genes(remInd);

if isfield(newModel,'geneShortNames')
    newModel.geneShortNames(remInd) = [];
end
if isfield(newModel,'geneMiriams')
    newModel.geneMiriams(remInd) = [];
end
if isfield(newModel,'geneFrom')
    newModel.geneFrom(remInd) = [];
end
if isfield(newModel,'geneComps')
    newModel.geneComps(remInd) = [];
end


end



function [updatedRule,rScore] = processSimpleRule(rule,genes,gScores,isozymeScoring,complexScoring)
% Either score or modify a reaction gene rule containig only ANDs or ORs.
%
% If the rule contains an enzyme complex (all ANDs), the complex will be
% scored based on the score of its subunits. Subunits without a score (NaN)
% will be excluded from the score calculation.
%
% If the rule contains only isozymes (all ORs), the negative-score genes
% will be removed from the rule. Isozymes without a score (NaN) will not be
% removed from the rule. The resuling rule will then be scored.


% get IDs and indices of genes involved in rule
ruleGenes = unique(regexp(rule,'[^&|\(\) ]+','match'));
[~,geneInd] = ismember(ruleGenes,genes);

% rules with one or no genes remain unchanged
if numel(ruleGenes) < 2
    rScore = gScores(geneInd);
    updatedRule = rule;
    return
end

if ~contains(rule,'&')  % rule contains isozymes
    
    scoreMethod = isozymeScoring;
    negInd = gScores(geneInd) < 0;  % NaNs will return false here
    if all(negInd)
        % get the least negative gene, adding a small random value to avoid a tie
        [~,maxInd] = max(gScores(geneInd) + rand(size(geneInd))*(1e-8));
        updatedRule = ruleGenes{maxInd};
    elseif sum(~negInd) == 1
        updatedRule = ruleGenes{~negInd};
    else
        updatedRule = strjoin(ruleGenes(~negInd),' | ');
        if startsWith(rule,'(')
            updatedRule = ['(',updatedRule,')'];
        end
    end
    
    % update ruleGenes and their indices
    ruleGenes = unique(regexp(updatedRule,'[^&|\(\) ]+','match'));
    [~,geneInd] = ismember(ruleGenes,genes);
    
elseif ~contains(rule,'|')  % rule contains enzyme complex
    scoreMethod = complexScoring;
    updatedRule = rule;
else
    error('This function cannot handle rules with both "OR" and "AND" expressions.');
end

% score rule
switch lower(scoreMethod)
    case 'min'
        rScore = min(gScores(geneInd),[],'omitnan');
    case 'max'
        rScore = max(gScores(geneInd),[],'omitnan');
    case 'median'
        rScore = median(gScores(geneInd),'omitnan');
    case 'average'
        rScore = mean(gScores(geneInd),'omitnan');
end

end



function updatedRule = processComplexRule(rule,genes,gScores,isozymeScoring,complexScoring)
% Update reactions containing both AND and OR expressions.
%
% Negative-score genes will be removed if they are isozymic, whereas they
% will not be removed if they are part of an enzyme complex. However, if
% the enzyme complex has a negative score, the entire complex will be
% removed, as long as it is not the only remaining element in the rule.


% Specify phrases to search for in the grRule. These phrases will find
% genes grouped by all ANDs (first phrase) or all ORs (second phrase).
search_phrases = {'\([^&|\(\) ]+( & [^&|\(\) ]+)+\)', '\([^&|\(\) ]+( \| [^&|\(\) ]+)+\)'};

% initialize some variables
subsets = {};  % subsets are groups of genes grouped by all ANDs or all ORs
c = 1;  % counter to keep track of the group (subset) number
r_orig = rule;  % record original rule to determine when it stops changing
for k = 1:100  % iterate some arbitrarily high number of times
    for j = 1:length(search_phrases)
        new_subset = regexp(rule,search_phrases{j},'match')';  % extract subsets
        if ~isempty(new_subset)
            subsets = [subsets; new_subset];  % append to list of subsets
            subset_nums = arrayfun(@num2str,(c:length(subsets))','UniformOutput',false);  % get group numbers to be assigned to the new subsets, and convert to strings
            rule = regexprep(rule,search_phrases{j},strcat('#',subset_nums,'#'),'once');  % replace the subsets in the expression with their group numbers (enclosed by "#"s)
            c = c + length(new_subset);
        end
    end
    if isequal(rule,r_orig)
        break;  % stop iterating when rule stops changing
    else
        r_orig = rule;
    end
end
subsets{end+1} = rule;  % add final state of rule as the last subset

% score and update each subset, and append to gene list and gene scores
for i = 1:numel(subsets)
    [subsets{i},subset_score] = processSimpleRule(subsets{i},genes,gScores,isozymeScoring,complexScoring);
    gScores = [gScores; subset_score];
    genes = [genes; {strcat('#',num2str(i),'#')}];
end

% reconstruct the rule from its updated subsets
updatedRule = subsets{end};
for i = c-1:-1:1
    updatedRule = regexprep(updatedRule,strcat('#',num2str(i),'#'),subsets{i});
end

end



