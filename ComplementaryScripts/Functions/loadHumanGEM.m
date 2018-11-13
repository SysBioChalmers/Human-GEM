function model = loadHumanGEM
% loadHumanGEM
%   Load humanGEM and prepare for a simulation-ready model by contraining 
%   problematic reactions archived in incactivationRxns.tsv
%
% Usage: model = loadHumanGEM
%
% Hao Wang, 2018-11-12
%

% get model path
[ST, I]=dbstack('-completenames');
modelPath=fileparts(fileparts(fileparts(ST(I).file)));

% load reactions need to be constrained
inactivateRxnsFile=fullfile(modelPath,'ComplementaryData','modelCuration','inactivationRxns.tsv');
fid = fopen(inactivateRxnsFile,'r');
input = textscan(fid,'%s %s','Delimiter','\t','Headerlines',1);
fclose(fid);
rxnsToConstrain = input{1};

% load model and conduct constraining
matFile=fullfile(modelPath,'ModelFiles','mat','humanGEM.mat');
load(matFile);
model = setParam(ihuman, 'eq', rxnsToConstrain, 0);

end

