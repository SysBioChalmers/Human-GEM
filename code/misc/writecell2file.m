function [] = writecell2file(C,fileName,header,delim,formatSpec,date)
%writecell2file  Write cell array to file.
%
% USAGE:
%
%   writecell2file(C,fileName,header,delim,formatSpec, date);
%
% INPUT:
%
%   C            Cell array to be written to file.
%
%   fileName     Name of file to which cell C will be written.
%
%   header       (Optional, default FALSE) If TRUE, the first row in C will
%                be treated as a header, and processed separately from the
%                remaining rows of C.
%
%                NOTE: if C contains only strings, this argument
%                essentially has no effect.
% 
%   delim        (Optional, default '\t') Delimiter by which to separate
%                columns of C when writing to file.
%
%                NOTE: if formatSpec is included as an input, then DELIM
%                will only be used on the header (if HEADER = TRUE).
%
%   formatSpec   (Optional) A string specifying the formatting rules
%                applied to each row in C. For example, if C contains two
%                text columns and one column of unsigned integers, and
%                columns are to be separated by tabs, the corresponding
%                FORMATSPEC should be: '%s\t%s\t%u\n'
%
%   date         (Optional, default FALSE) If TRUE, print out the date
%                in the first line of output file 
%


% handle input arguments
if nargin < 6
    date = false;
end
if nargin < 5
    formatSpec = [];
end
if nargin < 4 || isempty(delim)
    delim = '\t';
end
if nargin < 3 || isempty(header)
    header = false;
end

[nrows,ncols] = size(C);
if isempty(formatSpec)
    % if formatSpec not provided, assume the cell contains only strings,
    % and will be separated by tabs
    formatSpec = strcat(repmat(strcat('%s',delim),1,ncols-1),'%s','\n');
end

% separate header if specified
if ( header )
    header_data = C(1,:);
    C(1,:) = [];
    nrows = nrows - 1;
    formatSpecHead = strcat(repmat(strcat('%s',delim),1,ncols-1),'%s','\n');
end

% write to file
f = fopen(fileName,'w');
if ( date )
    fprintf(f,'# Date: %s\n',datestr(clock,'yyyy-mm-dd'));
end
if ( header )
    fprintf(f,formatSpecHead,header_data{:});
end
for i = 1:nrows
    fprintf(f,formatSpec,C{i,:});
end
fclose(f);


