function [simple_rules,skipped] = simplifyGeneRules(grRules,expanded)
%simplifyGeneRules  Simplify and condense the logic of model grRules.
%
% simplifyGeneRules recasts grRules as mathematical equations, where ORs
% are replaced with addition (+), and ANDs are replaced with multiplication
% (*). These equations are simplified using Matlab's symbolic math toolbox,
% yielding rules that are as reduced/simplified as possible, while still
% maintaining functional equivalency to the original grRule logic.
%
% NOTE: The function performs an extensive simplification process, and can
%       therefore be quite slow if there are many long grRules containing
%       both ANDs and ORs.
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
% Note that grRules can use written ("and" "or") or symbolic ("&" "|")
% boolean operators.
%
%
% USAGE:
%
%   simple_rules = simplifyGeneRules(grRules,expand);
%
%
% INPUT:
%
%   grRules        The grRules field from a genome-scale metabolic model.
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
%
% OUTPUT:
%
%   simple_rules   Updated/simplified grRules.
%
%   skipped        Logical vector indicating which grRules were skipped
%                  because they were too long/complex to be converted into
%                  a symbolic equation. For these cases, the grRule will
%                  be copied into simple_rules without any simplification.
%
%
% Jonathan Robinson, 2018-07-31

if nargin < 2
    expanded = false;
end

if ( expanded )
    warning('The EXPANDED option is not yet functional');
end

% check if the grRules use written or symbolic boolean operators
if any(contains(grRules,{'&','|'}))
    % fix some potential missing spaces between parentheses and &/|
    grRules = regexprep(grRules,'\)&',') &');   % ")&"  ->  ") &"
    grRules = regexprep(grRules,'&\(','& (');   % "&("  ->  "& ("
    grRules = regexprep(grRules,'\)\|',') |');  % ")|"  ->  ") |"
    grRules = regexprep(grRules,'\|\(','| (');  % "|("  ->  "| ("
    Btype = 'symbol';  % record rule type, to revert back at the end
else
    % fix some potential missing spaces between parentheses and AND/OR
    grRules = regexprep(grRules,'\)and',') and');  % ")and" ->  ") and"
    grRules = regexprep(grRules,'and\(','and (');  % "and(" ->  "and ("
    grRules = regexprep(grRules,'\)or',') or');    % ")or"  ->  ") or"
    grRules = regexprep(grRules,'or\(','or (');    % "or("  ->  "or ("
    
    % convert "and" to "&" and "or" to "|" (easier to work with symbols)
    grRules = regexprep(grRules, ' or ', ' | ');
    grRules = regexprep(grRules, ' and ', ' & ');
    Btype = 'text';  % record rule type, to revert back at the end
end

% remove any repeated ORs and ANDs
grRules = regexprep(grRules,'\|(\s+\|)+','|');  % remove repeated ORs
grRules = regexprep(grRules,'&(\s+&)+','&');  % remove repeated ANDs

% remove any trailing ORs and ANDs
grRules = regexprep(grRules,'\s*\|\s*(\))|(\()\s*\|\s*','$1');  % remove orphan ORs
grRules = regexprep(grRules,'^\s*\|\s*|\s*\|\s*$','');  % remove ORs at the beginning or end of rule
grRules = regexprep(grRules,'\s*&\s*(\))|(\()\s*&\s*','$1');  % remove orphan ANDs
grRules = regexprep(grRules,'^\s*&\s*|\s*&\s*$','');  % remove ANDs at the beginning or end of rule


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

% Simplify gene rules by converting them into mathematical equations:
% ANDs are functionally equivalent to multiplication (*), whereas ORs
% are functionally equivalent to addition (+).
h = waitbar(0,'Processing grRules...');
for i = 1:length(AndOrInd)
    waitbar(i/length(AndOrInd),h);
    
    % retrieve rule, and convert '&' to '*', and '|' to '+'
    r = grRules{AndOrInd(i)};
    r = regexprep(r,'&','*');
    r = regexprep(r,'\|','+');
    
    % Iterate through a series of simplifications and cleaning until rule
    % no longer changes
    r_orig = r;
    for k = 1:10  % perform an arbitrarily high number of iterations
        
        try
            % convert string to symbolic equation
            reqn = str2sym(r);
        catch
            fprintf('grRule #%u was too complex, and therefore skipped.\n',AndOrInd(i));
            skipped(AndOrInd(i)) = true;
            break
        end
        
        % Iterate through a series of equation expansions and simplifications, each
        % time removing constant coefficients and exponents, until the equation is
        % no longer changing.
        reqn_orig = reqn;
        for kk = 1:10  % perform an arbitrarily high number of iterations
            
            % for each gene in equation, collect terms for that gene
            g = symvar(reqn);
            for j = 1:length(g)
                
                % collect equation in terms of gene
                reqn = collect(reqn,g(j));
                r = char(reqn);
                
                % remove coefficients and exponents
                r = regexprep(r,'(^| |\()\d+\*','$1');  % remove numerical coeffs
                r = regexprep(r,'\^\d+','');  % remove exponents
                
                if is_sum(reqn) || is_prod(reqn)
                    % group terms, but also remove any terms with a numeric
                    % constant, because they will always be TRUE
                    r = group_terms(r);
                end
                
                % convert back to symbolic eqn
                reqn = str2sym(r);
            end
            
            % expand equation
            r = char(expand(reqn));
            
            % remove coefficients and exponents
            r = regexprep(r,'(^| |\()\d+\*','$1');  % remove numerical coeffs
            r = regexprep(r,'\^\d+','');  % remove exponents
            
            % simplify equation
            try
                reqn = simplify(str2sym(r));
            catch
                fprintf('grRule #%u was too complex, and therefore skipped.\n',AndOrInd(i));
                skipped(AndOrInd(i)) = true;
                break
            end
            
            % check if the equation has changed
            if isequal(char(reqn),char(reqn_orig))  % need to compare eqns in string form, otherwise it compares mathematically
                % if no changes, we can stop iterating
                break;
            else
                % go through the process again if something has changed
                reqn_orig = reqn;
            end
        end
        
        % break outer for-loop if equation was skipped due to complexity
        if skipped(AndOrInd(i))
            break
        end
        
        % Add parentheses around groups of genes. This may add some
        % unnecessary parentheses, but they will be cleaned up later.
        % NOTE: This also removes terms in the equation containing a
        % numeric constant: e.g., "(GENE2 + 1)". In this case, the
        % term will always be TRUE, so it is removed from the rule eqn.
        if is_sum(reqn) || is_prod(reqn)
            % the input and output of this function are of the type "char"
            r = group_terms(char(reqn));
            r = regexprep(r,'^\(|\)$','');  % remove the extra set of parentheses that always enclose the expression
        end
        
        % check if the rule has changed
        if isequal(r,r_orig)  % need to compare eqns in string form, otherwise it compares mathematically
            % if no changes, we can stop iterating
            break;
        else
            % go through the process again if something has changed
            r_orig = r;
        end
        
    end

    if ~skipped(AndOrInd(i))
        % reformat to gene rule logic ("*" becomes "&", and "+" becomes "|")
        r = regexprep(r,'\*',' & ');  % need to add spaces around the "*"s
        r = regexprep(r,'\+','|');
        
        % update rule in grRule vector
        grRules{AndOrInd(i)} = r;
    end
    
end
close(h);


% iterate through the cleaning process until the grRules no longer change
grRules_orig = grRules;
for i = 1:50  % perform a max of 50 iterations (it shouldn't require nearly that many)
    
    % remove extra parentheses
    has_AND_OR = contains(grRules,'&') & contains(grRules,'|');
    grRules(~has_AND_OR) = regexprep(grRules(~has_AND_OR),'\(|\)','');  % remove all parentheses from rules that don't have both ANDs and ORs
    grRules = regexprep(grRules,'\((\S*)\)','$1');  % remove extra parentheses: orphan '()' or those surrounding only one gene '(GENE1)'
    grRules = regexprep(grRules,'\((\([^\(\)]*\))\)','$1');  % remove double parentheses around expressions: '((GENE1 & GENE2))' or '((GENE1 | GENE2))'
    
    % determine if the grRules changed (in length), and if not, exit the loop
    if isequal(grRules,grRules_orig)
        break
    else
        grRules_orig = grRules;
        if i == 50
            warning('After 50 iterations, grRules still changed. There is potentially a problem with grRules or this function.');
        end
    end
    
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



%% Group equation terms, add parentheses, and remove unnecessary terms
function eqn_str_new = group_terms(eqn_str)
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
ind = find(contains(eqn_pieces,{'+','*'}));
for i = 1:length(ind)
    eqn_pieces{ind(i)} = group_terms(eqn_pieces{ind(i)});
end

% remove any empty equation pieces
eqn_pieces(cellfun(@isempty,eqn_pieces)) = [];

% re-assemble equation pieces, adding + or *, depending on the input eqn
if is_sum(eqn)
    if ~all(isnan(str2double(eqn_pieces)))
        % This is a special case where one of the terms of the equation is
        % a numeric constant (e.g., 'GENE1 + 1'). In this case, the
        % expression is always TRUE, so it can be removed.
        eqn_str_new = '';
    else
        eqn_str_new = char(strcat('(',join(eqn_pieces,' + '),')'));
    end
elseif is_prod(eqn)
    eqn_str_new = char(strcat('(',join(eqn_pieces,'*'),')'));
else
    % if the input eqn was simply one gene or empty, just return the input
    eqn_str_new = eqn_str;
end

end



%% Functions to test whether an equation is a sum or product of terms
function res = is_sum(eqn)
res = length(regexp(char(eqn),'\+')) > length(regexp(char(children(eqn)),'\+'));
end

function res = is_prod(eqn)
res = length(regexp(char(eqn),'\*')) > length(regexp(char(children(eqn)),'\*'));
end




