function stdFormulas = standardizeMetFormulas(metFormulas)
%standardizeMetFormulas
%   Standardize metabolite chemical formulas by reordering elements
%   according to the Hill system.
%
%
%   metFormulas    Cell array vector of metabolite atomic formulas.
%                  NOTE: Formulas containing characters other than letters
%                        or digits (A-Z,a-z,0-9) will NOT be reordered.
%
%
%   stdFormulas    Cell array vector of metabolite formulas rearranged in
%                  according to the Hill system: 
%                      C and H first, followed by all other elements in 
%                      alphabetical order. If no C, then all elements 
%                      (including H) are arranged in alphabetical order.
%                      Ex: C17H27N5O4SR2 -> C17H27N5O4R2S
%                      Ex: HBr -> BrH
%
%
%    Usage: stdFormulas = standardizeMetFormulas(metFormulas)
%


% ignore metabolite formulas containing parentheses
ignoreMets = ~cellfun(@isempty,regexp(metFormulas,'[^A-Za-z0-9]','match'));

% initialize alphabetized met formulas with original met formulas
stdFormulas = metFormulas;

% reorder metabolite formulas
stdFormulas(~ignoreMets) = cellfun(@(F) strjoin(sort(regexp(F,'[A-Z][a-z]*\d*','match')),''),stdFormulas(~ignoreMets),'UniformOutput',false);

% place C and H at start of formula if C is present
hasC = ~cellfun(@isempty,regexp(stdFormulas,'C[^a-z]','once')) & ~ignoreMets;
stdFormulas(hasC) = cellfun(@(F) strjoin([sort(regexp(F,'[CH]\d*(?![a-z])','match')), sort(regexp(F,'[A-BD-GI-Z][a-z]*\d*|[A-Z][a-z]+\d*','match'))],''),stdFormulas(hasC),'UniformOutput',false);




