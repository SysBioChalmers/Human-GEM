function newCell=reformatElements(inputCell,type,delimiter)
%reformatElements  reformat elements of cell array to desired format
% reformatElements
%   convert cell array element format between string and nested cell
%
% Input:
%   inputCell    the input cell array
%   type         two conversion approaches: cell2str and str2cell
%                (default: str2cell)
%   delimiter    specify the delimiter separating values within the
%                element (default: semicolon)
%
% Output:
%   newCell      the output of cell array with refromatted elements   
%
% Usage: newCell=reformatElements(inputCell,type,delimiter)
%


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
    inputCell = regexprep(inputCell, '\s', '');   % remove space from input
	newCell(index)=cellfun(@(s) strsplit(s,delimiter),inputCell(index),'UniformOutput', false);
elseif isequal('cell2str',type)
    % combine elements of each cell to string
    delimiter = [delimiter(~isspace(delimiter)) ' '];  % append one space only
    newCell(index)=cellfun(@(s) strjoin(s,delimiter),inputCell(index),'UniformOutput', false);
end

end
