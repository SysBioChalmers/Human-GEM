function list=getNonEmptyList(queryList,type)
% getNonEmptyList
%   Return the index of empty or non-empty elements in a cell array
%
% Input
%   queryList    input cell
%   type         if true (default), get the index of non-empty elements
%                from the input cell, or get the index of empty elements
%                if false
%
% Output:
%   list         the index of empty (or non-empty) elements
%
% Usage: list=getNonEmptyList(queryList,type)
%


list=[];

if nargin<1
	error('Missing input argument');
end

if nargin<2
	type=true;
end

if type
	list=find(~cellfun(@isempty,queryList));
else
	list=find(cellfun(@isempty,queryList));
end

end


