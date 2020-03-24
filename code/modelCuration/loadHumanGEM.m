function model = loadHumanGEM
% loadHumanGEM
%   Load humanGEM and prepare for a simulation-ready model by contraining 
%   problematic reactions archived in incactivationRxns.tsv
%
% Usage: model = loadHumanGEM
%
% Hao Wang, 2018-11-14
%

% get model path
[ST, I]=dbstack('-completenames');
modelPath=fileparts(fileparts(fileparts(ST(I).file)));

% load model
matFile=fullfile(modelPath,'ModelFiles','mat','HumanGEM.mat');
load(matFile);

% get reactions need to be constrained
rxnsToConstrain = '';
inactivateRxnsFile=fullfile(modelPath,'ComplementaryData','modelCuration','inactivationRxns.tsv');
if exist(inactivateRxnsFile, 'file') == 2
    fid = fopen(inactivateRxnsFile,'r');
    input = textscan(fid,'%s %s','Delimiter','\t','Headerlines',1);
    fclose(fid);
    rxnsToConstrain = input{1};
end

% conduct constraining, if any
if ~isempty(rxnsToConstrain)
    model = setParam(ihuman, 'eq', rxnsToConstrain, 0);
else
    model = ihuman;
end

end

