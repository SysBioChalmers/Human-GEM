function Padj = adjust_pvalues(P,method,dim)
%adjust_pvalues  Adjust p-values for multiple comparisons/testing.
%
% Usage:
%
%   Padj = adjust_pvalues(P,method,dim);
%
%
% Input:
%
%   P         A vector or matrix of p-values. Does not need to be sorted, 
%             and will be returned in the order provided.
%
%   method    'bonferroni' - Use the Bonferroni correction.
%             'holm' - Use the Holm-Bonferroni correction.
%             'benjamini' - Use the Benjamini-Hochberg correction.
%
%   dim       Dimension along which to group p-values. By default, p-values
%             will be grouped by column, unless P is a row-vector, in which
%             case the row vector will be treated as a single group. The
%             number of p-values in a group corresponds to the number of
%             "tests" for which they will be corrected.
%                 1 = by columns (default)
%                 2 = by rows
%             'all' = treat all p-values as single group
%
%
% Output:
%
%   Padj      A vector or matrix of adjusted p-values. Note that adjusted
%             p-values are restricted to a maximum value of 1.
%
%
% Jonathan Robinson, 2019-04-14


% if DIM is not specified, act along 1st dimension (columns), unless P is a
% vector, then act along the longest dimension.
if nargin < 3
    if all(size(P) > 1)
        dim = 1;
    else
        [~,dim] = max(size(P));
    end
end

% handle DIM input
if ~isnumeric(dim) && strcmpi(dim,'all')
    [nrows,ncols] = size(P);
    P = P(:);
    n = 1;
else
    if dim == 2
        P = P';
    end
    n = size(P,2);
end

Padj = ones(size(P));
for i = 1:n
    
    p = P(:,i);
    m = length(p);
    padj = ones(size(p));
    
    switch lower(method)
        
        case 'bonferroni'
            
            padj = p*m;
            
        case 'holm'
            
            [p,sort_ind] = sort(p);  % sort p-values in ascending order
            padj(1) = p(1)*m;  % adjust first p-value
            
            for j = 2:m  % adjust all other p-values
                padj(j) = max([p(j)*(m-j+1),padj(j-1)]);
            end
            
            % unsort p-values
            [~,unsort_ind] = sort(sort_ind);
            padj = padj(unsort_ind);
            
        case 'benjamini'
            
            [p,sort_ind] = sort(p);  % sort p-values in ascending order
            padj(m) = p(m);  % last p-value is unadjusted
            
            for j = m-1:-1:1  % adjust all other p-values
                padj(j) = min([p(j)*(m/j),padj(j+1)]);
            end
            
            % unsort p-values
            [~,unsort_ind] = sort(sort_ind);
            padj = padj(unsort_ind);
            
        otherwise
            
            error('"%s" is not a valid adjustment method.\n Valid options are: "bonferroni", "holm", or "benjamini".',method);
            
    end
    
    padj(padj > 1) = 1;  % cap adjusted p-values at 1
    Padj(:,i) = padj;
    
end

% return Padj to original dimensions or orientation of input P.
if ~isnumeric(dim) && strcmpi(dim,'all')
    Padj = reshape(Padj,nrows,ncols);
elseif dim == 2
    Padj = Padj';
end

