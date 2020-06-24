function Cnest = nestCell(C,remEmpty)
%NESTCELL  Compress a matrix cell array into a vector of nested cells.
%
%   CNEST = NESTCELL(C) take cell array C and combine all of the columns
%   into nested cells for each row.
%
%   CNEST = NESTCELL(C,REMEMPTY) if TRUE, all empty elements will be
%   removed when combining cells. Default = FALSE.
%


if nargin < 2
    remEmpty = false;
end

% set up indexing variable according to nesting dimension
x = transpose(1:size(C,1));

if ( remEmpty )
    % remove empty indices as cells are nested
    nonempty = ~cellfun(@isempty,C);
    Cnest = arrayfun(@(a) [{C{a,nonempty(a,:)}}],x,'UniformOutput',false);
else
    % keep empty indices as cells are nested
    Cnest = arrayfun(@(a) [{C{a,:}}],x,'UniformOutput',false);
end



