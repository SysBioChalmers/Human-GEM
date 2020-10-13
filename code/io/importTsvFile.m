function structure = importTsvFile(filename, numeric_cols)
% importTsvFile
%
%   Loads content from a tab-separated value (tsv) file into a structure.
%
%   importTsvFile will interpret columns as strings if their entries are
%   quoted (""), and will interpret columns without quotes as numeric.
%   If the tsv file does not contain quotes, all columns will be
%   interpreted as strings. This can be overridden with the numeric_cols
%   argument.
%
% Input:
%
%   filename      Name of the .tsv annotation file to be loaded.
%
%   numeric_cols  (Optional) Index (or indices) of the columns that should
%                 be interpreted as numeric (double). This will override
%                 the automatic interpretation of column types.
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

% This loading process is written in this seemingly overcomplicated manner
% because we want to import the data without knowing in advance how many
% columns there are, and to automatically interpret the column type based
% on the presence/absence of quotes "".

% read the second line of the file to determine the column formats
fid = fopen(filename);
fgetl(fid);  % skip first line (column names)
L = strsplit(fgetl(fid), '\t', 'CollapseDelimiters', false);
fclose(fid);

% load file import options and modify with expected column formats
opt = detectImportOptions(filename, 'FileType', 'text', 'Delimiter', '\t');
opt.VariableTypes(contains(L, '"')) = {'char'};
opt.VariableTypes(~contains(L, '"')) = {'double'};
opt.DataLines = [2 Inf];  % data starts from line 2 (readtable sometimes guesses this incorrectly)

% update column types
if all(ismember(opt.VariableTypes, 'double'))
    opt.VariableTypes(:) = {'char'};
end
if ~isempty(numeric_cols)
    opt.VariableTypes(numeric_cols) = {'double'};
end

% import the file as a table and convert to structure
tab = readtable(filename, opt);
structure = table2struct(tab, 'ToScalar', true);



