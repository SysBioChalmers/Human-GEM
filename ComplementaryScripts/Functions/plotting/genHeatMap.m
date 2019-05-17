function h = genHeatMap(data,colnames,rownames,clust_dim,clust_dist,col_map,col_bounds,grid_color)
%genHeatMap  Generate a heatmap for a given matrix of data.
%
% Usage:
%
%   genHeatMap(data,colnames,rownames,,clust_dim,clust_dist,col_map,col_bounds);
%
% Input:
%
%   data        Numerical matrix.
%
%   colnames    Cell array of data column names.
% 
%   rownames    Cell array of data row names.
%
%   clust_dim   'none'  the data will be plotted as provided (DEFAULT)
%               'rows'  cluster/rearrange the rows based on distance
%               'cols'  cluster/rearrange the columns based on distance
%               'both'  cluster/rearrange rows and columns based on distance
%
%   clust_dist  Distance metric to be used for clustering, ignored if
%               clust_dim is 'none'. Options are the same as those for
%               distance in, e.g., PDIST ('euclidean', 'hamming', etc.).
%               (DEFAULT = 'euclidean')
%
%   col_map     Colormap, provided as string (e.g., 'parula', 'hot', etc.)
%               or an Nx3 RGB matrix of N colors.
%               (DEFAULT = 'magma')
%
%   col_bounds  A 2-element vector with min and max values, to manually set
%               the bounds of the colormap.
%               (DEFAULT = min/max of data)
%
%   grid_color  String or 1x3 RGB vector indicating the color of grid lines.
%               (DEFAULT = 'none')
%
% Output:
%
%   h           (Optional) Handle to a pcolor plotting object. The heatmap
%               will be plotted automatically.
%
%
% Jonathan Robinson, 2019-02-27


% handle input arguments
if nargin < 4 || isempty(clust_dim)
    clust_dim = 'none';
elseif ~ismember(clust_dim,{'none','rows','cols','both'})
    error('%s is not a valid CLUST_DIM option. Choose "none", "rows", "cols", or "both".',clust_dim);
end
if nargin < 5 || isempty(clust_dist)
    clust_dist = 'euclidean';
end
if nargin < 6 || isempty(col_map)
    col_map = 'magma';
end
if nargin < 7 || isempty(col_bounds)
    col_bounds = [min(data(:)),max(data(:))];
end
if nargin < 8 || isempty(grid_color)
    grid_color = 'none';
end

% linkage algorithm for hierarchical clustering
% (see "linkage" function for more options)
linkage_method = 'average';

% perform hierarchical clustering to sort rows (if specified)
if ismember(clust_dim,{'rows','both'})
    L = linkage(data,linkage_method,clust_dist);
    row_ind = optimalleaforder(L,pdist(data,clust_dist));
else
    row_ind = 1:size(data,1);
end
% perform hierarchical clustering to sort columns (if specified)
if ismember(clust_dim,{'cols','both'})
    L = linkage(data',linkage_method,clust_dist);
    col_ind = optimalleaforder(L,pdist(data',clust_dist));
else
    col_ind = 1:size(data,2);
end

% reorder data matrix according to clustering results
sortdata = data(row_ind,col_ind);
sortrows = rownames(row_ind);
sortcols = colnames(col_ind);

% check if data is square matrix with identical row and column names
if (length(colnames) == length(rownames)) && all(strcmp(colnames,rownames))
    % flip data so the diagonal is from upper left to lower right
    sortdata = fliplr(sortdata);
    sortcols = flipud(sortcols);
end

% pad data matrix with zeros (pcolor cuts off last row and column)
sortdata(end+1,end+1) = 0;

% generate pcolor plot
a = axes;
set(a,'YAxisLocation','Right','XTick',[],'YTick', (1:size(sortdata,1))+0.5,'YTickLabels',sortrows);
set(a,'TickLength',[0 0],'XLim',[1 size(sortdata,2)],'YLim',[1 size(sortdata,1)]);
hold on

h = pcolor(sortdata);
set(h,'EdgeColor',grid_color);
set(gca,'XTick', (1:size(sortdata,2))+0.5);
set(gca,'YTick', (1:size(sortdata,1))+0.5);
set(gca,'XTickLabels',sortcols,'YTickLabels',sortrows);
set(gca,'XTickLabelRotation',90);

% bad form to have try/catch, but it works for now
try
    colormap(col_map);
catch
    colormap(custom_cmap(col_map));
end

if ~isempty(col_bounds)
    caxis(col_bounds);
end


