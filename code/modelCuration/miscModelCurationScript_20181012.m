%
% FILE NAME:    miscModelCurationScript_20181012.m
% 
% PURPOSE: Script to restore two reactions that were previously deleted,
%          but after further investigation were decided that they should be
%          returned to the model. These reactions were deleted as part of
%          the 'miscModelCurationScript_20180921.m' script, so all of their
%          information was retrieved from the model prior to applying that
%          script, and is re-inserted into the model with this script.
%
%          The reactions are added back to the model in such a way as to
%          try and maintain their original indexing in the model, just for
%          convenience.
%


%% Load model

% load latest version of humanGEM
load('humanGEM.mat');  % v0.4.1


%% Restore transport reactions that were removed earlier

% There was insufficient evidence to remove the reaction that allowed 
% transport of HMG-CoA through the peroxisomal membrane (HMGCOAtx). 
% Furthermore, the HMG-CoA transport through the mitochondrial membrane 
% should proceed via the carnitine shuttle, but for now the simplified 
% direct transport is a suitable placeholder reaction (HMR_1572).
%   HMR_1572:  HMG-CoA[c] <=> HMG-CoA[m]
%   HMGCOAtx:  HMG-CoA[c] <=> HMG-CoA[p]

% try to restore them to their original index location in the model (before
% reactions HMR_1284 and HMGCOARr)
[~,r_ind] = ismember({'HMR_1284';'HMGCOARr'},ihuman.rxns);
nMet = length(ihuman.mets);
nGene = length(ihuman.genes);
nProt = length(ihuman.proteins);
met_ind = getIndexes(ihuman,{'HMG-CoA[c]';'HMG-CoA[m]';'HMG-CoA[p]'},'metscomps');

% generate columns for stoich matrix
s1 = zeros(nMet,1);
s1(met_ind(1:2)) = [-1;1];
s2 = zeros(nMet,1);
s2(met_ind([1,3])) = [-1;1];

% restore information in all model fields
ihuman.rxns = [ihuman.rxns(1:r_ind(1)); {'HMR_1572'}; ihuman.rxns((r_ind(1)+1):r_ind(2)); {'HMGCOAtx'}; ihuman.rxns((r_ind(2)+1):end)];
ihuman.S = [ihuman.S(:,1:r_ind(1)), s1, ihuman.S(:,(r_ind(1)+1):r_ind(2)), s2, ihuman.S(:,(r_ind(2)+1):end)];
ihuman.lb = [ihuman.lb(1:r_ind(1)); -1000; ihuman.lb((r_ind(1)+1):r_ind(2)); -1000; ihuman.lb((r_ind(2)+1):end)];
ihuman.ub = [ihuman.ub(1:r_ind(1)); 1000; ihuman.ub((r_ind(1)+1):r_ind(2)); 1000; ihuman.ub((r_ind(2)+1):end)];
ihuman.rev = double((ihuman.lb < 0) & (ihuman.ub > 0));
ihuman.c = [ihuman.c(1:r_ind(1)); 0; ihuman.c((r_ind(1)+1):r_ind(2)); 0; ihuman.c((r_ind(2)+1):end)];
ihuman.rxnNames = [ihuman.rxnNames(1:r_ind(1)); {''}; ihuman.rxnNames((r_ind(1)+1):r_ind(2)); {'Hydroxymethylglutaryl Coenzyme A Reversible Peroxisomal Transport'}; ihuman.rxnNames((r_ind(2)+1):end)];
ihuman.rxnComps = [ihuman.rxnComps(1:r_ind(1)); 3; ihuman.rxnComps((r_ind(1)+1):r_ind(2)); 1; ihuman.rxnComps((r_ind(2)+1):end)];
ihuman.grRules = [ihuman.grRules(1:r_ind(1)); {''}; ihuman.grRules((r_ind(1)+1):r_ind(2)); {''}; ihuman.grRules((r_ind(2)+1):end)];
ihuman.rxnGeneMat = [ihuman.rxnGeneMat(1:r_ind(1),:); zeros(1,nGene); ihuman.rxnGeneMat((r_ind(1)+1):r_ind(2),:); zeros(1,nGene); ihuman.rxnGeneMat((r_ind(2)+1):end,:)];
ihuman.subSystems = [ihuman.subSystems(1:r_ind(1)); {'Transport, mitochondrial'}; ihuman.subSystems((r_ind(1)+1):r_ind(2)); {'Transport, peroxisomal'}; ihuman.subSystems((r_ind(2)+1):end)];
ihuman.eccodes = [ihuman.eccodes(1:r_ind(1)); {''}; ihuman.eccodes((r_ind(1)+1):r_ind(2)); {''}; ihuman.eccodes((r_ind(2)+1):end)];
ihuman.rxnKEGGID = [ihuman.rxnKEGGID(1:r_ind(1)); {''}; ihuman.rxnKEGGID((r_ind(1)+1):r_ind(2)); {''}; ihuman.rxnKEGGID((r_ind(2)+1):end)];
ihuman.rxnEHMNID = [ihuman.rxnEHMNID(1:r_ind(1)); {''}; ihuman.rxnEHMNID((r_ind(1)+1):r_ind(2)); {''}; ihuman.rxnEHMNID((r_ind(2)+1):end)];
ihuman.rxnBiGGID = [ihuman.rxnBiGGID(1:r_ind(1)); {'HMGCOAtm'}; ihuman.rxnBiGGID((r_ind(1)+1):r_ind(2)); {''}; ihuman.rxnBiGGID((r_ind(2)+1):end)];
ihuman.rxnHepatoNET1ID = [ihuman.rxnHepatoNET1ID(1:r_ind(1)); {''}; ihuman.rxnHepatoNET1ID((r_ind(1)+1):r_ind(2)); {''}; ihuman.rxnHepatoNET1ID((r_ind(2)+1):end)];
ihuman.rxnREACTOMEID = [ihuman.rxnREACTOMEID(1:r_ind(1)); {''}; ihuman.rxnREACTOMEID((r_ind(1)+1):r_ind(2)); {''}; ihuman.rxnREACTOMEID((r_ind(2)+1):end)];
ihuman.rxnReferences = [ihuman.rxnReferences(1:r_ind(1)); {''}; ihuman.rxnReferences((r_ind(1)+1):r_ind(2)); {'PMID:14713247'}; ihuman.rxnReferences((r_ind(2)+1):end)];
ihuman.rxnFrom = [ihuman.rxnFrom(1:r_ind(1)); {'HMRdatabase'}; ihuman.rxnFrom((r_ind(1)+1):r_ind(2)); {'Recon3D'}; ihuman.rxnFrom((r_ind(2)+1):end)];
ihuman.rxnMiriams = [ihuman.rxnMiriams(1:r_ind(1)); {''}; ihuman.rxnMiriams((r_ind(1)+1):r_ind(2)); struct('name',{{'pmid'}},'value',{{'PMID:14713247'}}); ihuman.rxnMiriams((r_ind(2)+1):end)];
ihuman.rxnConfidenceScores = [ihuman.rxnConfidenceScores(1:r_ind(1)); NaN; ihuman.rxnConfidenceScores((r_ind(1)+1):r_ind(2)); 0; ihuman.rxnConfidenceScores((r_ind(2)+1):end)];
ihuman.rxnRecon3DID = [ihuman.rxnRecon3DID(1:r_ind(1)); {'HMGCOAtm'}; ihuman.rxnRecon3DID((r_ind(1)+1):r_ind(2)); {''}; ihuman.rxnRecon3DID((r_ind(2)+1):end)];
ihuman.prRules = [ihuman.prRules(1:r_ind(1)); {''}; ihuman.prRules((r_ind(1)+1):r_ind(2)); {''}; ihuman.prRules((r_ind(2)+1):end)];
ihuman.rxnProtMat = [ihuman.rxnProtMat(1:r_ind(1),:); zeros(1,nProt); ihuman.rxnProtMat((r_ind(1)+1):r_ind(2),:); zeros(1,nProt); ihuman.rxnProtMat((r_ind(2)+1):end,:)];
ihuman.priorCombiningGrRules = [ihuman.priorCombiningGrRules(1:r_ind(1)); {''}; ihuman.priorCombiningGrRules((r_ind(1)+1):r_ind(2)); {''}; ihuman.priorCombiningGrRules((r_ind(2)+1):end)];


% clear intermediate variables
clear met_ind nGene nMet nProt r_ind s1 s2



%% Remove non-standard fields (e.g. rxnKEGGID) to comply with defined RAVEN structure
% The following five fields are non-standard. To assist convenient model 
% manipulation, they are removed here given that these information have
% been intactly stored elsewhere.
f = {'rxnKEGGID';'rxnEHMNID';'rxnBiGGID';'rxnHepatoNET1ID';'rxnREACTOMEID'};
ihuman = rmfield(ihuman, f);


%% Save model
save('../../model/Human-GEM.mat','ihuman');

