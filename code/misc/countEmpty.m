function countNum = countEmpty(queryList,nonEmpty)
% countEmpty
%   Count number of empty (or non-empty) elements in a cell array
%
% Input:
%   queryList    input cell array
%   nonEmpty     true is for counting non-empty elements, while false (default)
%                for counting empty elements 
% Output:
%   countNum     the number of empty (or non-empty) elements
%
% Usage: countNum = countEmpty(queryList,nonEmpty)
%


if nargin<1
	error('Missing input arguments!');
end

if nargin<2
	nonEmpty = false;
end

if nonEmpty
    countNum = length(find(~cellfun(@isempty,queryList)));
else
	countNum = length(find(cellfun(@isempty,queryList)));
end
