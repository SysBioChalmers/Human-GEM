%
% FILE NAME:    integrateRecon3DToHMR.m
% 
% DATE CREATED: 2018-04-25
% DATE CREATED: 2018-05-30
% 
% PROGRAMMER:   Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% 
% PURPOSE: The master scrit for integrate Recon3D into HMR2 based on
%          reaction/metabolite association to MetaNetX/BiGG databases
%

% 1. Load HMR2 and Recon3D model in RAVEN format
load('Recon3DRaven.mat');     %Recon3D
load('HMRdatabase2_00.mat');  %HMR2


% 2. Get the reduced Recon3D model by removing transport/exchange
% and overlapped rxns

% Load rxn assoc between Recon3D and HMR2
load('Recon3Rxns2HMR.mat');
mappedRecon3DRxns=find(~cellfun(@isempty,Recon3D.rxnHMRID));

% Get the index with transport/exchange rxns independent of subSystems
transportRxnIdx=find(getTransportRxns(Recon3DRaven));          % Transport rxns by getTransprot 4773
subSystems=cellfun(@char,Recon3D.subSystems,'un',0);;
exchangeRxnIdx=find(contains(subSystems, 'Exchange'));         % Exchange rxns 1972
% An addtional step of rxn filtering could be included here
addtionalRxnIdx=find(~cellfun(@isempty,Recon3D.withMNXnoHMR)); % Rxns extended to other comps 259

% Get the reduced Recon3D model by removing reacions
rxnToRemove=unique([mappedRecon3DRxns;transportRxnIdx;exchangeRxnIdx;addtionalRxnIdx]);
uniqueRecon3D=removeReactions(Recon3DRaven,rxnToRemove,1,1,1);
save('uniqueRecon3D.mat','uniqueRecon3D');   % 2018-05-30


% 3. Prepare reduced Recon3D model for combining with HMR2
% Load met association and replace mapped Recon3D mets
load('Recon3Mets2HMR.mat');
[a, b]=ismember(uniqueRecon3D.mets,Recon3Mets2HMR.mets);
uniqueRecon3D.mets(find(a))=Recon3Mets2HMR.metHMRID(b(find(a)));

% An addtional step of met filtering could be included here


% 4. Merge models
uniqueRecon3D.id='uniqueRecon3D';
HMR3=mergeModels({ihuman uniqueRecon3D});


% 5. Retrieve/add transport/exchange reactions
% To be completed


% 6. Detect and remove repetitive reactions and metabolites
% To be completed


% 7. Rebalance reactions under pH 7.3
% To be completed


% 8. The grRules should be considered/resolved between steps 2 to 6
% To be completed
