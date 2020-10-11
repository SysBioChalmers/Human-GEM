function annot = importAnnotations(filename, out_format)
% importAnnotations
%
%   Load Human-GEM reaction, metabolite, or gene annotations from their
%   respective .tsv annotation file.
%
%   NOTE: Only two column types are recognized: 'char' or 'double', based
%         on whether the column entries are quoted ("") or not,
%         respectively.
%
% Input:
%
%   filename    Name of the .tsv annotation file to be loaded.
%
%   out_format  Format of the loaded annotation data.
%               Options are 'struct' (Default) or 'table'
%
% Output:
%
%   annot       A structure or table (depending on OUT_FORMAT) containing
%               the annotation information associated with the GEM
%               reactions, metabolites, or genes.
%
% Usage:
%
%   annot = importAnnotations(filename, out_format);
%

if nargin < 2
    out_format = 'struct';
end

% This loading process was written in this seemingly overcomplicated manner
% because we wanted to import the data without knowing in advance how many
% columns there were, and to automatically interpret the column type based
% on the presence/absence of quotes "".

% read the second line of the file to determine the column formats
fid = fopen(filename);
fgetl(fid);  % skip first line
L = strsplit(fgetl(fid), '\t');
fclose(fid);

% load file import options and modify with expected column formats
opt = detectImportOptions(filename, 'FileType', 'text', 'Delimiter', '\t');
opt.VariableTypes(contains(L, '"')) = {'char'};
opt.VariableTypes(~contains(L, '"')) = {'double'};
opt.DataLines = [2 Inf];  % data starts from second line

% import the file as a table
annot = readtable(filename, opt);

% convert to structure if requested
if startsWith(lower(out_format), 'struct')
    annot = table2struct(annot, 'ToScalar', true);
end



