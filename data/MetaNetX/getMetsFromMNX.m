%
%   FILE NAME:    getMetsFromMNX.m
% 
%   PURPOSE: Generate the data structure for MetaNetX metabolites
%


% Move to the target MNX folder

% Load the MNX metabolites
T=readtable('chem_prop.xlsx','ReadVariableNames',1);
MNXMets=table2struct(T,'ToScalar',true);

% Rename the fields according to RAVEN specification
MNXMets.mets=MNXMets.MNX_ID;
MNXMets.metNames=MNXMets.Description;
MNXMets.metFormulas=MNXMets.Formula;
MNXMets.metCharges=str2double(MNXMets.Charge);
MNXMets.metMass=MNXMets.Mass;
MNXMets.inchis=MNXMets.InChI;
MNXMets.metSmiles=MNXMets.SMILES;
MNXMets.metInChIKey=MNXMets.InChIKey;
MNXMets.metSource=MNXMets.Source;
MNXMets=rmfield(MNXMets,{'MNX_ID','Description','Mass','Formula',....
'Charge','InChI','SMILES','Source','InChIKey'});

save('MNXMets.mat','MNXMets');
