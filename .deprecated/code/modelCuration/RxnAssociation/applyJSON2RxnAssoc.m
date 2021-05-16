%
% FILE NAME:    applyJSON2RxnAssoc.m
%
% PURPOSE: Sort out and clean up the previous .mat files generated during
%          reaction association/curation of human-GEM:
%          1. identify redundant rxn-association .mat files and remove them;
%          2. convert the remained .mat files to JSON format as discussed
%          in #75, while retaining the content unchanged.
%


%% sort out and clean the reaction association files to HMR2

% There are two Matlab model structures (ihumanRxns2BiGG.mat and 
% ihumanRxns2MNX.mat) that were generated from the HMR2 model during
% reaction-association curation. Since all the fields in ihumanRxns2BiGG
% are shared in ihumanRxns2MNX and with the same content, and the latter
% includes addtional information. The ihumanRxns2BiGG.mat file is moved to
% `deprecated` subfolder for reducing repetition.

% load the reaction association file and assign its info to a temp variable
load('ihumanRxns2MNX.mat');

r.rxns                = ihuman.rxns;                    % HMR id
r.rxnKEGGID           = ihuman.rxnKEGGID;               % KEGG
r.rxnEHMNID           = ihuman.rxnEHMNID;               % EHMN
r.rxnBiGGID           = ihuman.rxnBiGGID;               % BiGG
r.rxnHepatoNET1ID     = ihuman.rxnHepatoNET1ID;         % HepatoNET1
r.rxnREACTOMEID       = ihuman.rxnREACTOMEID;           % REACTOME
r.BiGG2BiGG           = ihuman.BiGG2BiGG;               % map to BiGG db id
r.HepatoNet12BiGG     = ihuman.HepatoNet12BiGG;         % HepatoNet1 map to BiGG id
r.EHMN2BiGG           = ihuman.EHMN2BiGG;               % EHMN map to BiGG id
r.HMR2BiGG            = ihuman.HMR2BiGG;                % combined BiGG ids from above
r.rxnBiGGDB2MNX       = ihuman.rxnBiGGDB2MNX;           % BiGG id to MNX id
r.rxnKEGG2MNX         = ihuman.rxnKEGG2MNX;             % KEGG id to MNX id
r.rxnREACTOMEStableID = ihuman.rxnREACTOMEStableID;     % REACTOME Stable ID
r.rxnReactome2MNX     = ihuman.rxnReactome2MNX;         % Reactome id 2 MNX id
r.rxnMNXID            = reformatElements(ihuman.rxnMNXID,'cell2str');  % MetaNetX, non-unique association
r.rxnCompIdx          = ihuman.rxnCompIdx;              % index of involved compartments


% encode to JSON format and pretty the layout before saving as plaintext file
jsonStr = jsonencode(r);
fid = fopen('ihumanRxns2MNX.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);

% check if the content of JSON file is identical to that of .mat file
check = jsondecode(fileread('ihumanRxns2MNX.JSON'));
if isequal(r, check)
    fprintf('\nThe file ihumanRxns2MNX.JSON is confirmed with the same content to ihumanRxns2MNX.mat, which therefore can be removed!\n\n');
end


%% sort and clean up reaction association files to Recon3D

% There are two Matlab model structures (Recon3Rxns2MNX.mat and Recon3Rxns2HMR.mat)
% that were generated for Recon3D reaction-association curation. Since
% all the fields of Recon3Rxns2MNX.mat have already been included in
% Recon3Rxns2HMR.mat, which also includes some addtional fields. Hence the
% Recon3Rxns2MNX.mat is thus moved to `deprecated` subfolder for avoidance
% of repetition.

% load the reaction association file and assign its info to a temp variable
load('Recon3Rxns2HMR.mat');
clear r;   % clean and reuse the temp variable

r.rxns         = Recon3D.rxns;          % Recon3D id
r.rxnKEGGID    = Recon3D.rxnKEGGID;     % KEGG
r.rxnBiGGID    = Recon3D.rxnBiGGID;     % BiGG
r.rxnMNXID     = Recon3D.rxnMNXID;      % MetaNetX
r.rxnHMRID     = Recon3D.rxnHMRID;      % HMR id
r.BiGG2HMR     = Recon3D.BiGG2HMR;      % BiGG map to HMR id
r.rxnCompIdx   = Recon3D.rxnCompIdx;    % index of involved compartments
%r.comps       = Recon3D.comps;         % identical field is in Recon3DRAVEN
%r.compNames   = Recon3D.compNames;     % identical field is in Recon3DRAVEN
r.withMNXnoHMR = Recon3D.withMNXnoHMR;  % can be associated to MNX but not HMR

% Note that BiGG2HMR and rxnHMRID are different, the latter is used for the
% generation of rxnAssoc.mat

% encode to JSON format and pretty the layout before saving as plaintext file
jsonStr = jsonencode(r);
fid = fopen('Recon3Rxns2HMR.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);

% check the content of JSON file, make sure they are identical to that of .mat file
check = jsondecode(fileread('Recon3Rxns2HMR.JSON'));
if isequal(r, check)
    fprintf('\nThe file Recon3Rxns2HMR.JSON is confirmed with the same content to Recon3Rxns2HMR.mat, which therefore can be removed!\n\n');
end

