function newCell=reformatElements(inputCell,type,delimiter)
%reformatElements  reformat elements of cell array to desired format
%
%   convert cell array element format between string and nested cell
%
%   inputCell    the input cell array
%   type         two conversion approaches: cell2str and str2cell
%                (default: str2cell)
%   delimiter    specify the delimiter separating values within the
%                element (default: semicolon)
%
%   newCell=reformatElements(inputCell,type,delimiter)
%
%   Hao Wang, 2018-06-03

newCell={};

if ~iscell(inputCell)
		EM='Wrong input arguments';
		disp(EM);
end

if nargin < 2
    type = 'str2cell';
end

if nargin < 3
    delimiter = ';';
end

% get the index of non-empty elements
index=find(~cellfun(@isempty,inputCell));

if length(index)>0
		if (iscell(inputCell{index(1)}) && strcmpi(type,'str2cell')) || (~iscell(inputCell{index(1)})  &&  strcmpi(type,'cell2str'))
		fprintf('The conversion parameter is conflic with input data!\n');
		return;
elseif length(index)==0
		fprintf('The input cell is empty!\n');
		return;		
end

% initilize output
newCell=cell(numel(inputCell),1);
newCell(:)={''};

if strcmp('str2cell',type)
    % from string to nested cell
		newCell(index)=cellfun(@(s) strsplit(s,delimiter),inputCell(index),'UniformOutput', false);
elseif isequal('cell2str',type)
    % combine elements of each cell to string
    newCell(index)=cellfun(@(s) strjoin(s,delimiter),inputCell(index),'UniformOutput', false);
end

end
