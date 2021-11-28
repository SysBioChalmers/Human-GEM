function exportTsvFile(data, filename, withQuotes)
% exportTsvFile
%
%   Export structure or table to a tsv (or txt) file. String (char) entries
%   will be enclosed in double quotes ("") by default, whereas numeric
%   (double) entries will not.
%
% Input:
%
%   data       A structure or table containing the annotation information
%              associated with the GEM reactions, metabolites, or genes.
%
%   filename   Name of the file to which the data will be written,
%              delimited by tabs.
%              If filename does not end in ".tsv" or ".txt", the extension
%              will automatically be appended as ".tsv".
%
%   withQuotes Enclose string elements with double quotes (opt, default is
%              TRUE)
%
%
% Usage:
%
%   exportTsvFile(data, filename, withQuotes);
%

if nargin < 3
    withQuotes = true;
end

if isstruct(data)
    data = struct2table(data);
end
if ~istable(data)
    error('Invalid DATA format. Must be a STRUCT or a TABLE.');
end
if ~endsWith(lower(filename), '.tsv') && ~endsWith(lower(filename), '.txt')
    filename = strcat(filename, '.tsv');
end

writetable(data, filename, 'Delimiter', '\t', 'QuoteStrings', withQuotes, 'FileType', 'text');


