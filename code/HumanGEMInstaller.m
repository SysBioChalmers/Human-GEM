classdef HumanGEMInstaller
% HumanGEMInstaller
%   Support for installing and uninstalling
%   Run HumanGEMInstaller.install to install (will set up the path in MATLAB)
%   Run HumanGEMInstaller.uninstall to clear the path from MATLAB
%   To install, you first need to cd to the HumanGEM root.

    methods (Static)
        function install
            sourceDir = fileparts(which(mfilename));
            paths = HumanGEMInstaller.GetFilteredSubPaths(sourceDir);
            addpath(paths);
            savepath;
        end
        function uninstall
            sourceDir = fileparts(which(mfilename));
            paths = HumanGEMInstaller.GetFilteredSubPaths(sourceDir);
            rmpath(paths);
            savepath;
        end
        function path = getHumanGEMMainPath()
			path = fileparts(which(mfilename));
			path = strrep(path, '\', '/'); %get rid of backslashes in Windows
			if ~endsWith(path, '/')
				path = strcat(path,'/');
			end
		end

        function newPaths = GetFilteredSubPaths(path_)
            pathSep = pathsep();
			%Check that there are no separators in the path - that will cause 
            %problems since the separator is used to separate paths in a string
			if contains(path_, pathSep)
				error('The path in which Human-GEM resides may not contain path separator chars such as semicolon for this installation to work!');
			end
            subpaths = arrayfun(@(subf) sprintf('%s', genpath([path_ filesep subf{:}])), HumanGEMInstaller.SUBFOLDERS, 'UniformOutput', false);
            newPaths = strjoin(subpaths, '');
        end
    end
    
    properties (Constant)
      SUBFOLDERS = {'../code', '../data', '../model'};
   end
end
