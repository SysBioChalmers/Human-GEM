function compModel = compressModelField(model,field,delim)
%compressModelField  Compress multi-column field into single-column field.
%
% For a model with specified field containing multiple columns, the columns
% will be combined into a single column, with the multiple entries for each
% row separated by a delimiter.
%
% USAGE:
%
%   compModel = compressModelField(model,field,delim);
%
% INPUTS:
%
%   model      Model structure.
%
%   field      Model structure field to be compressed into a single column.
%
%   delim      Delimiter by which multiple entries should be separated.
%              (Default = '; ')
%
% OUTPUTS:
%
%   compModel  Model structure with specified field compressed into a
%              single column format.
%


% handle inputs
if nargin < 3
    delim = '; ';
end

% extract data from field
F = model.(field);

% check if model field is of correct format
if size(F,2) == 1
    fprintf('Model field "%s" is already a single column. No changes will be made.\n',field);
    return
elseif ~iscell(F)
    error('Specified model field must be a cell array.');
end

% determine non-empty entries in field
non_empty = ~cellfun(@isempty,F);

% % ensure that empty entries in the field are empty strings
% F(cellfun(@isempty,F)) = {''};

% join entries for each row, separating by delimiter
F = arrayfun(@(i) strjoin(F(i,non_empty(i,:)),delim),1:size(F,1),'UniformOutput',false)';

% add compressed field to output model
compModel = model;
compModel.(field) = F;



