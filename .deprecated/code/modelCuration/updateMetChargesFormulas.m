%
% FILE NAME:    updateMetChargesFormulas.m
%
% PURPOSE:Update metabolites with curated charges and formulas for humanGEM
%


%% Load models and curated metabolite info
load('humanGEM.mat');
load('Recon3D_301.mat');
load('metAssocHMR2Recon3.mat');  % load curated metabolite info

% rename 'temp006x' to 'm01451x' (both correspond to 'cholesterol-ester pool')
ihuman.mets{ismember(ihuman.mets,'temp006x')} = 'm01451x';

% remove compartment abbrevs
metsNoComp = regexprep(ihuman.mets,'\_\w$','');
metsNoComp = regexprep(metsNoComp,'^(m\d+)\w$','$1');
metsNoComp = regexprep(metsNoComp,'^(temp\d+)\w$','$1');
Recon3D.metsNoComp = regexprep(Recon3D.mets,'\[\w\]$','');


%% Update metabolites with curated charges and formulas

% get index to met curation array and Recon3D
[hit2HMR, indHMRID]=ismember(metsNoComp,metAssocHMR2Recon3.metHMRID);
[hit2R3D, indR3DID]=ismember(metsNoComp,Recon3D.metsNoComp);
IHMR=find(hit2HMR);
IR3D=find(hit2R3D);

% Confirm the met ids are correctly and full mapped
%isequal(numel(metsNoComp),numel([IHMR;IR3D])) % Should be true
%isempty(setdiff(transpose(1:numel(ihuman.mets)),[IHMR;IR3D]))  % Should be true

% update metCharges and metFormulas
metCharges=repmat({''},size(ihuman.mets));
metCharges(IHMR)=metAssocHMR2Recon3.metCuratedCharges(indHMRID(IHMR));
metCharges(IR3D)=num2cell(Recon3D.metCharges(indR3DID(IR3D)));

metFormulas=repmat({''},size(ihuman.mets));
metFormulas(IHMR)=metAssocHMR2Recon3.metCuratedFormulas(indHMRID(IHMR));
metFormulas(IR3D)=regexprep(Recon3D.metFormulas(indR3DID(IR3D)),'FULLR','R');

ihuman.metCharges=cellfun(@int64, metCharges);
ihuman.metFormulas=metFormulas;

%% Save updated model file
save('humanGEM.mat','ihuman');

