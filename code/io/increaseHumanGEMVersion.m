function increaseHumanGEMVersion(bumpType)
% increaseHumanGEMVersion
%   Increase version for the humanGEM respositories
%
% Input:
%   bumpType      the value has to be one selected among 'major', 'minor' or 'patch'
%
% Usage: increaseHumanGEMVersion(bumpType)
%


%Check if in master:
currentBranch = git('rev-parse --abbrev-ref HEAD');
if ~strcmp(currentBranch,'master')
    error('ERROR: not in master')
end

%Get model path
[ST, I]=dbstack('-completenames');
modelPath=fileparts(fileparts(fileparts(ST(I).file)));

%Bump version number:
versionFile=fullfile(modelPath,'version.txt');
fid = fopen(versionFile,'r');
    oldVersion = fscanf(fid, '%s');
fclose(fid);
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
        error('ERROR: invalid input. Use either "major", "minor" or "patch"')
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
ymlFile=fullfile(modelPath,'model','Human-GEM.yml');
ihuman = importHumanYaml(ymlFile);

%Include tag and save model:
ihuman.version = newVersion;

%Export model to multiple formats
exportHumanGEM(ihuman,'Human-GEM',modelPath,{'mat', 'txt', 'xml', 'yml', 'xlsx'});

%Update version file:
fid = fopen(versionFile,'wt');
fprintf(fid,newVersion);
fclose(fid);

%Update readme file:
readmeFile=fullfile(modelPath,'README.md');
content = fileread(readmeFile);
content = strrep(content,'{{DATE}}',datestr(now,29));
content = strrep(content,'{{nRXN}}',num2str(length(ihuman.rxns)));
content = strrep(content,'{{nMET}}',num2str(length(ihuman.mets)));
content = strrep(content,'{{nGENE}}',num2str(length(ihuman.genes)));
fid = fopen(readmeFile,'wt');
fwrite(fid,content);
fclose(fid);

end
