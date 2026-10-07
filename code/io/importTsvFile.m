function structure = importTsvFile(filename, numeric_cols)
% importTsvFile
%
%   Loads content from a tab-separated value (tsv) file into a structure.
%
%   Every column is interpreted as a string, whether or not its entries are
%   quoted (""). Columns that should be numeric (double) are named with the
%   numeric_cols argument.
%
% Input:
%
%   filename      Name of the .tsv annotation file to be loaded.
%
%   numeric_cols  (Optional) Index (or indices) of the columns that should
%                 be interpreted as numeric (double) instead of as a
%                 string.
%
% Output:
%
%   structure     A structure containing the tsv file contents, where field
%                 names of the structure will correspond to column names
%                 from the first line of the tsv file.
%
% Usage:
%
%   structure = importTsvFile(filename, numeric_cols);
%

if nargin < 2
    numeric_cols = [];
end

% detectImportOptions is used to discover the columns without knowing in
% advance how many there are; its guessed types are then replaced, so that
% an identifier column of digits (such as geneEntrezID) is text rather than
% a number, and so that the types do not depend on how the file is quoted.

opt = detectImportOptions(filename, 'FileType', 'text', 'Delimiter', '\t');
opt.VariableTypes(:) = {'char'};
opt.DataLines = [2 Inf];  % data starts from line 2 (readtable sometimes guesses this incorrectly)

if ~isempty(numeric_cols)
    opt.VariableTypes(numeric_cols) = {'double'};
end

% import the file as a table and convert to structure
tab = readtable(filename, opt);
structure = table2struct(tab, 'ToScalar', true);



