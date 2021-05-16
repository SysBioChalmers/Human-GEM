%
%   FILE NAME:    getRxnsFromMNX.m
% 
%   PURPOSE: Generate the data structure for MetaNetX reactions
%


% Move to the target MNX folder

% Load the MNX reactions
T=readtable('reac_prop.xlsx','ReadVariableNames',1);
MNXRxns=table2struct(T,'ToScalar',true);
% Refine the equations by removing compartment suffix
MNXRxns.equations=regexprep(MNXRxns.Equation,'\@MNXD\d+','');
MNXRxns.equations=regexprep(MNXRxns.equations,'=','<=>');
MNXRxns.Description=regexprep(MNXRxns.Description,'1 `','');
MNXRxns.Description=regexprep(MNXRxns.Description,'`','');

% Construct the S matrix and list of metabolites
[S, mets, badRxns, reversible]=constructS(MNXRxns.equations);  %time-consuming
MNXRxns.S=S;
MNXRxns.mets=mets;
MNXRxns.badRxns=badRxns;
% Count metabolite number for reactions
num=numel(MNXRxns.EC);
MNXRxns.metNum=zeros(num,1);
for i=1:num
	MNXRxns.metNum(i,1)=numel(find(MNXRxns.S(:,i)));
end

% Construct the full S matrix and list of metabolites
MNXRxns.fullequations=regexprep(MNXRxns.Equation,'=','<=>');
[fullS, fullMets, badRxns, reversible]=constructS(MNXRxns.fullequations);  %time-consuming
MNXRxns.fullS=fullS;
MNXRxns.fullMets=fullMets;
MNXRxns.fullSbadRxns=badRxns;
% Count metabolite number for reactions
MNXRxns.fullMetNum=zeros(num,1);
for i=1:num
	MNXRxns.fullMetNum(i,1)=numel(find(MNXRxns.fullS(:,i)));
end

save('MNXRxns.mat','MNXRxns');
