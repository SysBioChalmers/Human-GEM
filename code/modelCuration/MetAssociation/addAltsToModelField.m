function newModel = addAltsToModelField(model,field,newEntries,singleCol)
%addAltsToModelField  Add alternative entries to a model field.
%
% Compares and appends a new set of field values to an existing model
% field. If an existing field entry is empty, it will be overwritten by the
% new entry. If the existing field entry and new entry conflict, both will
% be saved by adding new columns to the model field.
%
% USAGE:
%
%   newModel = addAltsToModelField(model,field,newEntries,singleCol);
%
% INPUT:
%
%   model       model structure
% 
%   field       model structure field to which new entries are to be added
%
%   newEntires  a column or matrix of new entries to add to the specified
%               model field
%
%   singleCol   (optional, default FALSE) If TRUE, the updated field will
%               be compressed to a single column, where multiple entries 
%               per row will be separated by a semicolon and space '; '.
%               If the field contains only one column, this input will 
%               have no effect.
%
% OUTPUT:
%
%   newModel    model structure with the new entries added to the specified
%               field. Note that the updated field may have multiple
%               columns.
%


if nargin < 4
    singleCol = false;
end

if ~isfield(model,field)
    % if the field doesn't yet exist in the model, just add the new entries
    model.(field) = newEntries;
    newModel = model;
    return
elseif isequal(model.(field),newEntries)
    % no changes needed if new entries are identical to existing entries
    newModel = model;
    return
end

if (size(model.(field),2) == 1) && (size(newEntries,2) == 1)
    % if existing and new entries are both column vectors, simply add the
    % new (mismatching) entries as a second column
    mismatch_ind = ~strcmp(model.(field),newEntries);
    newEntries(~mismatch_ind) = {''};
    model.(field) = [model.(field),newEntries];
else
    for i = 1:size(newEntries,1)
        vals = [model.(field)(i,:),newEntries(i,:)];  % combine rows
        vals(cellfun(@isempty,vals)) = [];  % remove empty entries
        vals = unique(vals);  % get all unique entries for row
        model.(field)(i,1:length(vals)) = vals;  % update row in model
    end
end

% replace empty matrices with empty strings
model.(field)(cellfun(@isempty,model.(field))) = {''};

% compress multiple columns into single column, if specified
if ( singleCol )
    model = compressModelField(model,field,'; ');
end

newModel = model;  % assign output


