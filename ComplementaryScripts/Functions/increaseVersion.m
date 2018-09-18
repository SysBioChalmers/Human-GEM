function increaseVersion(bumpType)
% increaseVersion
%
%   Increase version for GEM respositories
%
%   Usage: function increaseVersion(bumpType)
%
%   Benjamín J. Sánchez, 2018-07
%   Hao Wang, 2018-09-18
%

%Check if in master:
currentBranch = git('rev-parse --abbrev-ref HEAD');
if ~strcmp(currentBranch,'master')
    error('ERROR: not in master')
end

%Bump version number:
load('../../ModelFiles/mat/humanGEM.mat');
oldVersion = ihuman.version;
oldVersion = str2double(strsplit(oldVersion,'.'));
newVersion = oldVersion;
switch bumpType
    case 'major'
        newVersion(1) = newVersion(1) + 1;
        newVersion(2) = 0;
        newVersion(3) = 0;
    case 'minor'
        newVersion(2) = newVersion(2) + 1;
        newVersion(3) = 0;
    case 'patch'
        newVersion(3) = newVersion(3) + 1;
    otherwise
        error('ERROR: invalid input. Use "major", "minor" or "patch"')
end
newVersion = num2str(newVersion,'%d.%d.%d');

%Check if history has been updated:
%fid     = fopen('../../history.md','r');
%history = fscanf(fid,'%s');
%fclose(fid);
%if ~contains(history,['human' newVersion ':'])
%    error('ERROR: update history.md first')
%end
%To be included

%Load model:
%ihuman = importModel('../../ModelFiles/xml/humanGEM.xml');

%Include tag and save model:
ihuman.version = newVersion;

%Store model as .mat:
save('../../ModelFiles/mat/humanGEM.mat','ihuman');
%deal with other formats later
%exportForGit(ihuman,'humanGEM','..',{'mat', 'txt', 'xlsx', 'xml', 'yml'});

%Update version file:
fid = fopen('../../version.txt','wt');
fprintf(fid,newVersion);
fclose(fid);

end
