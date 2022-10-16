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
            paths = HumanGEMInstaller.GetFilteredSubPaths(sourceDir, HumanGEMInstaller.FILE_FILTER);
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

        function newPaths = GetFilteredSubPaths(path_, filter_)
            pathSep = pathsep();
			%Check that there are no separators in the path - that will cause 
            %problems since the separator is used to separate paths in a string
			if contains(path_, pathSep)
				error('The path in which Human-GEM resides may not contain path separator chars such as semicolon for this installation to work!');
			end
            paths = genpath(path_);
            splitPaths = strsplit(paths, pathSep);
            %remove the last, it is empty
            splitPaths = splitPaths(1,1:end-1);
            matches = regexp(splitPaths, filter_, 'match');
            okPaths = cellfun(@isempty, matches);
            pathsLeft = splitPaths(1,okPaths);
            newPaths = strcat(char(join(pathsLeft, pathSep)), pathSep);
        end
    end
    
    properties (Constant)
      FILE_FILTER = '(.*\.git.*)|(.*\.deprecated.*)';
   end
end
