function alphaMets = alphabetizeMetFormulas(metFormulas)
%alphabetizeMetFormulas  Reorder metabolite formulas alphabetically.
%
% USAGE:
%
%   alphaMets = alphabetizeMetFormulas(metFormulas);
%
% INPUTS:
%
%   metFormulas     Cell array vector of metabolite atomic formulas.
%                   NOTE: formulas containing characters other than letters
%                         or digits (A-Z,a-z,0-9), or those containing the
%                         phrase "FULLR", will NOT be reordered.
%
% OUTPUTS:
%
%   alphaMets       Cell array vector of metabolite formulas rearranged in
%                   alphabetical order of the element symbols.
%                   Ex: C17H27N5O4SR2 -> C17H27N5O4R2S
%


% Ignore metabolite formulas containing special characters, or "FULLR", 
% which appears in the Recon3D model.
ignoreMets = ~cellfun(@isempty,regexp(metFormulas,'[^A-Za-z0-9]','match')) | ...
             contains(metFormulas,'FULLR');

% initialize alphabetized met formulas with original met formulas
alphaMets = metFormulas;

% reorder metabolite formulas
alphaMets(~ignoreMets) = cellfun(@(F) strjoin(sort(regexp(F,'[A-Z][a-z]*\d*','match')),''),alphaMets(~ignoreMets),'UniformOutput',false);







