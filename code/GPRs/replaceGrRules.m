function [grRules_new,genes,rxnGeneMat] = replaceGrRules(grRules,idMapping)
%replaceGrRules  Replace grRules with another set of gene IDs.
%
% NOTE: This function is adapted from "translateGrRules" and specifically
%       designed for converting grRules from current gene ID to another
%       gene ID type
%
% NOTE: The function also "cleans" the grRules, meaning duplicate genes,
%       extra parentheses, etc. will be removed from each grRule.
%
% Usage:
%
%   [grRules_new,genes,rxnGeneMat] = replaceGrRules(grRules,idMapping)
%
%
% Inputs:
%
%   grRules      A cell array of gene-reaction rules (e.g., model.grRules).
%
%
%   idMapping    must be an Nx2 cell array, where the first column contains
%                gene IDs corresponding to the current grRules, and the
%                second column contains new gene IDs to which the rules
%                will be translated.
%
% Outputs:
%
%   grRules_new  grRules with the gene IDs converted to the new IDs.
%
%   genes        A list of all genes found in the converted grRules.
%
%   rxnGeneMat   A matrix associating reactions to genes, built from the
%                converted grRules and genes list.
%


% handle input arguments
if nargin < 2
    error('Missing input.');
end


% check input data

% get original list of genes from the grRules
genes_orig = getGenesFromGrRules(grRules);
if ismember('and',genes_orig) || ismember('or',genes_orig)
    error('Problem reading grRules. Verify that all "and" and "or" elements are lowercase and surrounded by spaces.');
end

% remove duplicate genes
genes_orig  = unique(genes_orig,'stable');

% idMapping should be a NX2 cell array
if ~(size(idMapping,2) == 2)
    error('The idMapping data structure must be a NX2 cell array.');
% input genes should be consistent with those retrieved from grRules
elseif all(~ismember(idMapping(:,1), genes_orig))
    error('The input genes are Not consistent with those retrieved from grRules.');
end


% conduct replacement

% begin by "cleaning" the original grRules
% this is not necessary but can speed up the translation if the original
% grRules are not in a simplified format.
rules_orig = cleanGrRules(grRules);

% determine logical operator type and change to "&/|" for easy manipulation
textBooleanType = false;
if any(contains(grRules,{' and ',' or '}))
    textBooleanType = true;
    grRules_new = rules_orig;
    grRules_new = regexprep(grRules_new, ' or ', ' | ');
    grRules_new = regexprep(grRules_new, ' and ', ' & ');
end

% define function convertGeneId to convert gene IDs
% use strjoin to combine ids with 'or' operator if multiple hits are found
convertGeneId = @(g) strjoin(idMapping(ismember(idMapping(:,1),g),2), ' | ');

% replace grRules with another ID type
% The next line identifies gene IDs as collections of characters that
% are not spaces, parentheses, or the symbols '&' or '|'. It then replaces
% those gene IDs with the new gene ID type, which calls the convertGeneId
% inline function to retrieve.
grRules_new = regexprep(grRules_new, '[^&|\(\) ]+', '(${convertGeneId($0)})');


% prepare output

% clean up rules (removes extra parentheses, repeated genes, etc.)
grRules_new = cleanGrRules(grRules_new);

% restore "&" as "and" and "|" as "or"
if textBooleanType
    grRules_new = regexprep(grRules_new, ' \| ', ' or ');
    grRules_new = regexprep(grRules_new, ' & ', ' and ');
end

% generate new rxnGeneMat and gene list based on converted grRules
[genes,rxnGeneMat] = getGenesFromGrRules(grRules_new);


