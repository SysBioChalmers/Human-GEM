function increaseHumanGEMVersion(bumpType,test)
% increaseHumanGEMVersion
%   Increase version for the humanGEM respositories
%
% Input:
%   bumpType    either 'major', 'minor' or 'patch'
%   test        true if a test run should be done that can take in a
%               development branch (it will not check if current branch is
%               main, nor update the model version, history and readme files).
%
% Usage: increaseHumanGEMVersion(bumpType)
%
if nargin<2
    test=false;
end

%Get model path
[ST, I]=dbstack('-completenames');
modelPath=fileparts(fileparts(fileparts(ST(I).file)));

%Check RAVEN version
currVer = checkInstallation('versionOnly');
if strcmp(currVer,'develop')
    printOrange('WARNING: Cannot determine your RAVEN version as it is in a development branch.\n');
else
    currVerNum = str2double(strsplit(currVer,'.'));
    minmVer = '2.10.3';
    minmVerNum = str2double(strsplit(minmVer,'.'));
    if currVerNum(1) ~= minmVerNum(1)
        wrongVersion = currVerNum(1) < minmVerNum(1);
    elseif currVerNum(2) ~= minmVerNum(2)
        wrongVersion = currVerNum(2) < minmVerNum(2);
    else
        wrongVersion = currVerNum(3) < minmVerNum(3);
    end
end
if wrongVersion
    error('Minimum required RAVEN version is %s.',minmVer);
end

%Check if in main:
if ~test
    currentBranch = git('rev-parse --abbrev-ref HEAD');
    if ~strcmp(currentBranch,'main')
        error('ERROR: not in main')
    end

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
end

%Load model:
humanGEM = readYAMLmodel(fullfile(modelPath,'model','Human-GEM.yml'));

%Include tag and save model:
if ~test
    humanGEM.version = newVersion;
end

%Check if it matches reactions.tsv, metabolites.tsv and genes.tsv
fields = {'rxns','reactions';'mets','metabolites';'genes','genes'};
for i=1:size(fields,1)
    tsvList = importTsvFile(fullfile(modelPath,'model',[fields{i,2} '.tsv']));
    Lia     = ismember(humanGEM.(fields{i,1}), tsvList.(fields{i,1}));
    dispEM(['The following ' fields{i,2} ' are in model.' fields{i,1} ...
        ' but not in model/' fields{i,2} '.tsv:'],true,humanGEM.(fields{i,1})(~Lia),false);
    Lia     = ismember(tsvList.(fields{i,1}), humanGEM.(fields{i,1}));
    dispEM(['The following ' fields{i,2} ' are in model/' fields{i,2} ...
        '.tsv but not in model.' fields{i,1} ':'],true,tsvList.(fields{i,1})(~Lia),false);
end

%Export model to multiple formats, without annotation
writeYAMLmodel(humanGEM,fullfile(modelPath,'model','Human-GEM.yml'),true,false);
save(fullfile(modelPath,'model','Human-GEM.mat'),'humanGEM');
humanGEM = annotateGEM(humanGEM);  % Add annotation data to structure
exportForGit(humanGEM,'Human-GEM',fullfile(modelPath,'model'),{'xml', 'xlsx', 'txt'},'',false);

if ~test
    %Update version file:
    fid = fopen(versionFile,'wt');
    fprintf(fid,newVersion);
    fclose(fid);

    %Update readme file:
    readmeFile=fullfile(modelPath,'README.md');
    content = fileread(readmeFile);
    content = strrep(content,'{{DATE}}',datestr(now,29));
    content = strrep(content,'{{nRXN}}',num2str(length(humanGEM.rxns)));
    content = strrep(content,'{{nMET}}',num2str(length(humanGEM.mets)));
    content = strrep(content,'{{nGENE}}',num2str(length(humanGEM.genes)));
    fid = fopen(readmeFile,'wt');
    fwrite(fid,content);
    fclose(fid);
end
end
