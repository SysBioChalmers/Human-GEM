%
%   FILE NAME:    HMR2Curation.m
% 
%   PURPOSE: Curating HMR2 database model toward HMR3
%


% Load the original Matlab file of HMR2 database
load('HMRdatabase2_00.mat');

% Import HMR2 reactions from the Excel file 'RXNS' sheet that leave out
% the following five rows that represent reaction classes
% 7709: 'Exchange reactions';
% 8171: 'Fake reactions';
% 8176: 'Biomass reactions';
% 8177: 'HMR_biomass_Renalcancer';
% 8183: 'Included for connectivity for INIT'
T=readtable('HMRdatabase2_00.xlsx','Sheet','RXNS','ReadVariableNames',1);
HMR2=table2struct(T,'ToScalar',true);

% Adding reaction identifiers from external databases
% additional empty spaces were also removed druing this process
if isequal(ihuman.rxns,HMR2.RXNID)  
	ihuman.rxnKEGGID=HMR2.KEGGID;              %KEGG
	ihuman.rxnEHMNID=HMR2.EHMNID;              %EHMN
	ihuman.rxnBiGGID=HMR2.BIGGDATABASEID;      %Recon
	ihuman.rxnHepatoNET1ID=HMR2.HEPATONET1ID;  %Hepatonet1
	ihuman.rxnREACTOMEID=HMR2.REACTOMEID;      %Reactome
	ihuman.rxnReferences=HMR2.REFERENCES;      %References
end

% Some errors were spotted and fixed:
% The reaction HMR_2190 was assocated to four REACTOM reactions.
% The first one (REACT_22270) is correect, the other three
% (REACT_22097; REACT_22133; REACT_22219) should be wrong because
% they are empty ids by searching reactom.org, and thus removed
index=find(strcmp('HMR_2190',ihuman.rxns));
ihuman.rxnREACTOMEID{index}='REACT_22270';

% The KEGG id of HMR_7709 was associated to 'R0302' that should be typo
% It was manually checked and corrected to R03027
index=find(strcmp('HMR_7709',ihuman.rxns));
ihuman.rxnKEGGID{index}='R03027';

% Save as version 2.0.1
save('HMRdatabase2_01.mat','ihuman'); %===2018-01-16
