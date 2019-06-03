%
% FILE NAME:    getHuman1RxnAssoc.m
%
% DATE CREATED: 2019-05-23
%     MODIFIED: 2019-06-03
%
% PROGRAMMER:   Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
%
% PURPOSE: Provide Human1 with extensively associated exteranl reaction
%          identifieres in JSON format file (#75), as discussed in #107.
%

%% load data

% load reaction association/annotation that were prepared mainly in #8
% and fromatted into JSON in #79
ihumanRxns2MNX = jsondecode(fileread('ihumanRxns2MNX.JSON'));
Recon3Rxns2HMR = jsondecode(fileread('Recon3Rxns2HMR.JSON'));

load('humanGEM.mat'); % HumanGEM v1.0.3


%% generate data structure with comprehensive reaction associations

% align with the rxn index in HumanGEM
clear r;   % clean and use a temp variable
r.rxns = ihuman.rxns;

% prepare associations to a list of external sources
r.rxnKEGGID       = repmat({''},size(r.rxns));  % KEGG
r.rxnBiGGID       = repmat({''},size(r.rxns));  % BiGG
r.rxnEHMNID       = repmat({''},size(r.rxns));  % EHMN
r.rxnHepatoNET1ID = repmat({''},size(r.rxns));  % HepatoNET1
r.rxnREACTOMEID   = repmat({''},size(r.rxns));  % REACTOME Stable ID
r.rxnMNXID        = repmat({''},size(r.rxns));  % MetaNetX

% retrieve information from ihumanRxns2MNX
[a, b] = ismember(r.rxns, ihumanRxns2MNX.rxns);
ind = find(a);

r.rxnKEGGID(ind)       = ihumanRxns2MNX.rxnKEGGID(b(ind));
r.rxnBiGGID(ind)       = ihumanRxns2MNX.rxnBiGGID(b(ind));
r.rxnEHMNID(ind)       = ihumanRxns2MNX.rxnEHMNID(b(ind));
r.rxnHepatoNET1ID(ind) = ihumanRxns2MNX.rxnHepatoNET1ID(b(ind));
r.rxnREACTOMEID(ind)   = ihumanRxns2MNX.rxnREACTOMEStableID(b(ind));
r.rxnMNXID(ind)        = ihumanRxns2MNX.rxnMNXID(b(ind));

% retrieve additional associations from Recon3Rxns2HMR
additionalRxns = setdiff(r.rxns, ihumanRxns2MNX.rxns);
uniqueR3DRxns = intersect(additionalRxns, Recon3Rxns2HMR.rxns);

[~, ind_R3D] = ismember(uniqueR3DRxns, Recon3Rxns2HMR.rxns);
[~, ind_humanGEM] = ismember(uniqueR3DRxns, r.rxns);

r.rxnKEGGID(ind_humanGEM) = Recon3Rxns2HMR.rxnKEGGID(ind_R3D);
r.rxnBiGGID(ind_humanGEM) = Recon3Rxns2HMR.rxnBiGGID(ind_R3D);
r.rxnMNXID(ind_humanGEM)  = Recon3Rxns2HMR.rxnMNXID(ind_R3D);


%% Output rxn association in JSON format

jsonStr = jsonencode(r);
fid = fopen('humanGemRxnAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);

% check the content of JSON file
check = jsondecode(fileread('humanGemRxnAssoc.JSON'));
if isequal(r, check)
    fprintf('\nThe reaction association file is sucessfully exported!\n\n');
end

movefile('humanGemRxnAssoc.JSON','../../ComplementaryData/annotation');
clear;

