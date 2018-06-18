function report=countFrequency(queryList)
% 
%   count occurrences of string elements in a cell array:
%
%   queryList    the input cell array
%
%   Usage: report=countFrequency(queryList)
%
%   Hao Wang, 2018-06-02
%

report={};

if nargin<1
		EM='Missing input arguments';
		disp(EM);
end

if ~iscell(queryList)
		EM='Wrong input arguments';
		disp(EM);
end

report.uniqueList=unique(queryList,'stable');
report.frequency=cellfun(@(x) sum(ismember(queryList,x)),report.uniqueList,'un',0);

end