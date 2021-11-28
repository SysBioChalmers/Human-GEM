function exportTsvFile(data, filename)
% exportTsvFile
%
%   Export structure or table to a tsv file. String (char) entries will be
%   enclosed in double quotes (""), whereas numeric (double) entries will
%   not.
%
% Input:
%
%   data       A structure or table containing the annotation information
%              associated with the GEM reactions, metabolites, or genes.
%
%   filename   Name of the file to which the data will be written,
%              delimited by tabs.
%              If filename does not end in ".tsv", the extension will
%              automatically be appended.
%
% Usage:
%
%   exportTsvFile(data, filename);
%

if isstruct(data)
    data = struct2table(data);
end
if ~istable(data)
    error('Invalid DATA format. Must be a STRUCT or a TABLE.');
end
if ~endsWith(lower(filename), '.tsv') && ~endsWith(lower(filename), '.txt')
    filename = strcat(filename, '.tsv');
end

writetable(data, filename, 'Delimiter', '\t', 'QuoteStrings', true, 'FileType', 'text');

