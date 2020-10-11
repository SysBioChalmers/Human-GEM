function exportAnnotations(annot, filename)
% exportAnnotations
%
%   Export annotation structure or table to a .tsv file.
%
% Input:
%
%   annot    A structure or table containing the annotation information
%            associated with the GEM reactions, metabolites, or genes.
%
% Usage:
%
%   exportAnnotations(annot, filename);
%

if isstruct(annot)
    annot = struct2table(annot);
end
if ~istable(annot)
    error('Invalid annotation format. Must be a STRUCT or a TABLE.');
end

writetable(annot, filename, 'Delimiter', '\t', 'QuoteStrings', true, 'FileType', 'text');

