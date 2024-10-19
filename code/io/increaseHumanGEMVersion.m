function increaseHumanGEMVersion(bumpType)
% increaseHumanGEMVersion
%   Increase version for the humanGEM respositories
%
% Input:
%   bumpType      the value has to be one selected among 'major', 'minor' or 'patch'
%
% Usage: increaseHumanGEMVersion(bumpType)
%

%Check if in main:
currentBranch = git('rev-parse --abbrev-ref HEAD');
if ~strcmp(currentBranch,'main')
    error('ERROR: not in main')
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
fid     = fopen(fullfile(modelPath,'history.md'),'r');
history = fscanf(fid,'%s');
fclose(fid);
if ~contains(history,['human' newVersion ':'])
    error('ERROR: update history.md first')
end

%Load model:
ihuman = readYAMLmodel(fullfile(modelPath,'model','Human-GEM.yml'));

%Include tag and save model:
ihuman.version = newVersion;

%Check if it matches reactions.tsv, metabolites.tsv and genes.tsv
fields = {'rxns','reactions';'mets','metabolites';'genes','genes.tsv'};
for i=1:size(fields,1)
    tsvList = importTsvFile(fullfile(modelPath,'model',[fields{i,2} '.tsv']));
    Lia     = ismember(ihuman.(fields{i,1}), tsvList.(fields{i,1}));
    dispEM(['The following ' fields{i,2} ' are in model.' fields{i,1} ...
            ' but not in model/' fields{i,2} '.tsv:'],true,ihuman.(fields{i,1})(~Lia));
    Lia     = ismember(tsvList.(fields{i,1}), ihuman.(fields{i,1}));
    dispEM(['The following ' fields{i,2} ' are in model/' fields{i,2} ...
            '.tsv but not in model.' fields{i,1} ':'],true,tsvList.(fields{i,1})(~Lia));
end

%Export model to multiple formats, without annotation
writeYAMLmodel(ihuman,fullfile(modelPath,'model','Human-GEM.yml'),false,false);
save(ihuman,fullfile(modelPath,'model','Human-GEM.mat'));

ihuman = annotateGEM(ihuman);  % Add annotation data to structure
exportForGit(ihuman,'Human-GEM',modelPath,{'xml', 'xlsx', 'txt'});

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
