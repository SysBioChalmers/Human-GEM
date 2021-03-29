function complex_mat = getComplexesFromGeneRules(grRules)
%getComplexesFromGeneRules  Extract enzyme complex info from model grRules.
%
% getComplexesFromGeneRules uses the logic in grRules to identify which
% genes are treated as belonging to an enzyme complex, and returns a binary
% matrix describing all identified complexes and their gene constituents.
%
% USAGE:
%
%   complex_mat = getComplexesFromGeneRules(grRules);
%
% INPUT:
%
%   grRules      grRules from a genome-scale metabolic model structure.
%
% OUTPUT:
%
%   complex_mat  A (GxC) binary matrix, where G is the number of genes
%                appearing in grRules, and C is the number of complexes
%                idenditied. Each column of the matrix represents a unique
%                enzyme complex, with 1s indicating which genes belong to
%                the complex. Although the complexes are unique, they can
%                be subsets of one another.
%
