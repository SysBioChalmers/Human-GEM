function cmap = custom_cmap(colors,n,fractions)
%custom_cmap  Generate a custom or pre-defined RGB colormap.
%
% Usage:
%
%   cmap = custom_cmap(colors,n,fractions);
%
% Input:
%
%   colors      A Cx3 matrix of C RGB values to include in the colormap, or
%               a string specifying a pre-built custom colormap:
%               'redblue', 'redbluebright', 'green', 'blue', 'darkblue',
%               'magma'
%                   
%               ** custom_cmap('preview'); shows a preview of the pre-built
%                  colormap options
%
%   n           Total number of colors desired in the final colormap.
%               Default = 100.
%
%   fractions   Numeric vector representing the fractions of the colormap
%               dedicated to each color (should sum to 1). If empty or not
%               provided, the colors will be distributed evenly.
%
% Output:
%
%   cmap        An Nx3 matrix of RGB values comprising a colormap.
%
%
% Jonathan Robinson, 2019-02-27


if nargin < 2 || isempty(n)
    n = 100;
end

if ~isnumeric(colors)
    
    % define custom colormaps
    switch lower(colors)
        case 'redblue'
            hotcmap = hot(round(60*n/100));
            cmap = flipud([hotcmap(round(10*n/100+1):end,:); custom_cmap([0.1 0.1 0.8;0 0.7 0.9;1 1 1],round(n/2),[0.25;0.4;0.35])]);
        case 'redbluebright'
            hotcmap = hot(round(80*n/100));
            cmap = flipud([hotcmap(round(30*n/100+1):end,:); custom_cmap([0.2 0.2 1;0 1 1;1 1 1],round(n/2),[0.25;0.4;0.35])]);
        case 'green'
            cmap = custom_cmap([0.0118 0.2392 0.0157;0.2392 0.7490 0.2471;1 1 1],n,[0.25;0.4;0.35]);
        case 'blue'
            cmap = custom_cmap([8,58,106;8,81,156;49,130,189;107,174,214;189,215,231;255,255,255]./255,n);
        case 'darkblue'
            cmap = custom_cmap([0.015,0.106,0.247; 0.031,0.227,0.416; 0.055,0.361,0.580; 0.400,0.592,0.690; 0.706,0.804,0.875],n);
        case 'magma'
            cmap = custom_cmap([15,15,19; 42,0,92; 117,15,110; 211,51,86; 253,144,91; 251,255,178]./255,n);
        case 'redbluedark'
            cmap = custom_cmap([0.573,0.216,0.212; 1,1,1; 0.039,0.345,0.529],n);
        case 'preview'
            % preview all pre-defined colormap options
            cmap_list = {'redblue';'redbluebright';'redbluedark';'green';'blue';'darkblue';'magma'};
            Ncmap = length(cmap_list);
            cmap = [];
            for i = 1:length(cmap_list)
                cmap = [cmap; custom_cmap(cmap_list{i},200)];
            end
            plotmat = flipud(reshape(0:(200*Ncmap-1),200,length(cmap_list))');
            plotmat(end+1,end+1) = 0;
            h = pcolor(plotmat); hold on
            set(h,'EdgeColor','none');
            set(gca,'YTick',(1:Ncmap)+0.5,'YTickLabels',flipud(cmap_list),'FontSize',14);
            plot([zeros(1,Ncmap-1);202*ones(1,Ncmap-1)],repmat(2:Ncmap,2,1),'k-');
            colormap(cmap);
            return
        otherwise
            error('"%s" is not a valid colormap option.',colors);
    end
    
else
    
    [numcolors,z] = size(colors);
    if z ~= 3
        error('COLORS must be an N x 3 matrix.');
    end
    
    if nargin < 3 || isempty(fractions) || length(unique(fractions)) == 1
        cmap = flipud(interp1(1:numcolors, colors, linspace(1,numcolors,n)));
    else
        if n < numcolors
            warning('Too many colors. Not all will be included in CMAP.');
        elseif length(fractions) ~= numcolors
            error('COLORS and FRACTIONS must contain equal number of elements');
        end
        
        cmap = zeros(n,3);
        cx = linspace(0,1,n+1)' + 1/(2*n);
        cx(end) = [];
        
        x = sort([0; cumsum(fractions); cumsum(fractions(1:end-2)) + fractions(2:end-1)/2 ]);
        y = zeros(numcolors*2 - 1,3);
        for i = 1:numcolors
            y(i*2-1,:) = colors(i,:);
            if i < numcolors
                y(i*2,:) = mean([colors(i,:);colors(i+1,:)]);
            end
        end
        for i = 1:3
            cmap(:,i) = interp1(x,y(:,i),cx);
        end
        cmap = flipud(cmap); 
    end
end


    
    
    
    
    
    