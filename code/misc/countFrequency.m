function report=countFrequency(queryList)
% countFrequency
%   Count occurrences of unique elements in a cell array
%
% Input:
%   queryList     the input cell array
%
% Output:
%   report        the output structure of frequency count
%     uniqueList  the cell array of unique elements
%     frequency   occurrence frequency of unique elements in the input cell array
%                  
% Usage: report=countFrequency(queryList)
%


report={};

if nargin<1
	error('Missing input arguments');
end

if ~iscell(queryList)
	error('Wrong input argument');
end

report.uniqueList=unique(queryList,'stable');

frequency=cellfun(@(x) sum(ismember(queryList,x)),report.uniqueList,'un',0);
report.frequency=cell2mat(frequency);

end
