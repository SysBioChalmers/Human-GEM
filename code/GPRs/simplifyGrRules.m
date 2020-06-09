function [simple_rules,skipped] = simplifyGrRules(grRules,expanded)
%simplifyGrRules  Simplify and condense the logic of model grRules.
%
% simplifyGrRules uses Matlab's symbolic toolbox to mathematically simplify
% gene rules. Each rule is converted into a symbolic boolean expression,
% and these expressions are simplified yielding rules that are as reduced/
% simplified as possible, while still maintaining functional equivalency to
% the original grRule logic.
%
% The function will also remove unnecessary parentheses, trailing ANDs/ORs,
% and duplicate genes in model grRules. The function requires that the gene
% names/IDs do not contain spaces, parentheses, "&", "|", or "#".
% Furthermore, grRules must separate gene names/IDs from ANDs and ORs by at
% least one space. For example,
%
%       OK: "GENE1 and GENE2", or "GENE1 & GENE2"
%   NOT OK: "GENE1andGENE2", or "GENE1&GENE2"
% 
%
% *!WARNING!* Simplifying grRules with this function may eliminate some
%             gene-reaction associations. 
%             For example: 'G1 or (G1 and G2)' simplifies to 'G1'
%
%
% USAGE:
%
%   [simple_rules,skipped] = simplifyGrRules(grRules,expanded);
%
%
% INPUT:
%
%   grRules        The grRules field from a genome-scale metabolic model.
%                  Note that grRules can use written ("and" "or") or
%                  symbolic ("&" "|") boolean operators.
%
%   expanded       (Optional, default FALSE). If TRUE, the grRule will be
%                  returned in its expanded form, where groups of ANDs
%                  (enzyme complexes) are separated by ORs, rather than
%                  using a shorter nested format.
%
%                  For example:
%
%                   input: G1 and G2 and (G3 or G4) and (G4 or G3 or G3) and G1
%
%              simplified: G1 and G2 and (G3 or G4)
%
%       simplify + expand: (G1 and G2 and G3) or (G1 and G2 and G4)
%
%                  NOTE: Expanding rules can make them very long, and in
%                        some cases will exceed Matlab's limit on the
%                        length of a symbolic expression, requiring them to
%                        be skipped.
%
%
% OUTPUT:
%
%   simple_rules   Updated/simplified grRules.
%
%   skipped        Logical vector indicating which grRules were skipped
%                  because they were too long to be converted into a
%                  symbolic equation. For these cases, the grRule will be
%                  copied into simple_rules without any simplification.
%

if nargin < 2
    expanded = false;
end

% perform a preliminary "clean" of the gene rules
grRules = cleanGrRules(grRules);


% check if the grRules use written or symbolic boolean operators
if any(contains(grRules,{'&','|'}))
    Btype = 'symbol';  % record rule type, to revert back at the end
else
    % convert "and" to "&" and "or" to "|"
    grRules = regexprep(grRules, ' or ', ' | ');
    grRules = regexprep(grRules, ' and ', ' & ');
    Btype = 'text';  % record rule type, to revert back at the end
end

% get list of all gene IDs present in grRules
rxnGenes = cellfun(@(r) unique(regexp(r,'[^&|\(\) ]+','match')),grRules,'UniformOutput',false);
nonEmpty = ~cellfun(@isempty,rxnGenes);
allGenes = unique([rxnGenes{nonEmpty}]');

% Replace all gene IDs with G1, G2, G3, etc. to avoid potential problems
% with certain gene ID formats (e.g., Entrez IDs are numbers, which will
% cause problems).
tempIDs = strcat('G',arrayfun(@num2str,(1:length(allGenes))','UniformOutput',false));  % construct new ID vector
convertGene = @(g) tempIDs{ismember(allGenes,g)};  % define conversion function
grRules = regexprep(grRules, '[^&|\(\) ]+', '${convertGene($0)}');  % find and convert all gene IDs to new temporary type


% identify all rules containing both ANDs and ORs
AndOrInd = find(contains(grRules,'&') & contains(grRules,'|'));

% initialize output specifying which rules were skipped due to complexity
skipped = false(size(grRules));

% iterate through complex rules (those containing both ANDs and ORs)
for i = 1:length(AndOrInd)
    
    r = grRules{AndOrInd(i)};
    
    % attempt to convert string to symbolic equation
    try
        reqn = str2sym(r);
    catch
        fprintf('grRule #%u was too long (%u characters), and therefore skipped.\n',AndOrInd(i),numel(r));
        skipped(AndOrInd(i)) = true;
        continue
    end
    
    % simplify eqn, allowing for a max of 10 steps
    reqn = simplify(reqn,10);
    
    % attempt to expand equation if requested
    reqn_noexpand = reqn;
    if ( expanded )
        try
            reqn = expand(reqn);
        catch
            fprintf('Failed to expand grRule #%u due to excess length. Rule will remain in simplified form.\n',AndOrInd(i));
        end
    end
    
    % Add parentheses around groups of genes. This may add some
    % unnecessary parentheses, but they will be cleaned up later.
    if is_or(reqn) || is_and(reqn)
        % the input and output of this function are of the type "char"
        try
            % sometimes an expanded equation will fail here
            r = add_parentheses(char(reqn));
        catch
            fprintf('Failed to expand grRule #%u due to excess length. Rule will remain in simplified form.\n',AndOrInd(i));
            reqn = reqn_noexpand;
            r = add_parentheses(char(reqn));
        end
        r = regexprep(r,'^\(|\)$','');  % remove the extra set of parentheses that always enclose the expression
    end
    
    % update rule in grRule vector
    grRules{AndOrInd(i)} = r;
    
end

% restore original gene IDs
restoreGene = @(g) allGenes{ismember(tempIDs,g)};  % define conversion function
grRules = regexprep(grRules, '[^&|\(\) ]+', '${restoreGene($0)}');  % find and restore all gene IDs to original ID

% change boolean operators back to original type
if strcmpi(Btype,'text')
    grRules = regexprep(grRules, ' \| ', ' or ');
    grRules = regexprep(grRules, ' & ', ' and ');
end

% assign output
simple_rules = grRules;

end



%% Add parentheses around groups of equation terms
function eqn_str_new = add_parentheses(eqn_str)
% This is a recursive function, which groups terms in an equation with
% parentheses in order to clearly indicate the order of operations. This is
% necessary to avoid ambiguity when the equations are converted back into
% gene rules.
%
% NOTE: both eqn_str_new and eqn_str are strings (char).

% obtain symbolic form of equation
eqn = str2sym(eqn_str);

% separate equation into terms/pieces ("children")
eqn_pieces = arrayfun(@char,children(eqn),'UniformOutput',false);

% find pieces that are further separable, and recursively separate until
% all pieces are single genes
ind = find(contains(eqn_pieces,{'|','&'}));
for i = 1:length(ind)
    eqn_pieces{ind(i)} = add_parentheses(eqn_pieces{ind(i)});
end

% re-assemble equation pieces, adding | or &, depending on the input eqn
if is_or(eqn)
    eqn_str_new = char(strcat('(',join(eqn_pieces,' | '),')'));
elseif is_and(eqn)
    eqn_str_new = char(strcat('(',join(eqn_pieces,' & '),')'));
else
    % if the input eqn was simply one gene, just return the input
    eqn_str_new = eqn_str;
end

end



%% Functions to test whether an equation is a collection of ORs or ANDs
function res = is_or(eqn)
% check if outermost operations contain ORs
res = length(regexp(char(eqn),'\|')) > length(regexp(char(children(eqn)),'\|'));
end

function res = is_and(eqn)
% check if outermost operations contain ANDs
res = length(regexp(char(eqn),'&')) > length(regexp(char(children(eqn)),'&'));
end




