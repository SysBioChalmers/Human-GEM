% FILE NAME:    getHuman1MetAssoc.m
%
% PURPOSE: To address #107, the previous .mat files generated during mets
%          association/curation were sorted and refined for extensively
%          retrieving exteranl identifiers and saving in JSON format (#75).
%

%% Update metAssocHMR2Recon3 with manual curation results

% Load met association information prepared in #23
load('metAssocHMR2Recon3.mat');


% change variable name and extract field names
m = metAssocHMR2Recon3;
fields = fieldnames(m);

% Convert elements of two fields (metChEBIID, metR3DID) from cell to string
for i=1:numel(fields)
    if iscell(m.(fields{i}){1})
        m.(fields{i}) = reformatElements(m.(fields{i}),'cell2str');
    end
end
m.metChEBIID = regexprep(m.metChEBIID, 'CHEBi:', '');  % small refinements


% some fields need to be reformated for either
% adding space after delimiter (metMNXID, metCuratedMNXID) or
% removing additional spaces (metR3DFormulas,metR3DKEGGID,metCuratedFormulas)
reformatFields = {'metMNXID';'metCuratedMNXID';'metR3DFormulas';'metR3DKEGGID';'metCuratedFormulas'};
for i=1:numel(reformatFields)
    convert2Cell = reformatElements(m.(reformatFields{i}),'str2cell');
    m.(reformatFields{i})  = reformatElements(convert2Cell,'cell2str');
end

% The metLIPIDMAPSID field has three elements that need to manually fixed
% m.metLIPIDMAPSID{336} ='LMFA01050113 LMFA01050349 LMFA01050359 LMFA02000035';
% m.metLIPIDMAPSID{1216}='LMFA01070018 LMFA02000037';
% m.metLIPIDMAPSID{606} ='LMST01010086;LMST01010144';
m.metLIPIDMAPSID{336} ='LMFA01050113; LMFA01050349; LMFA01050359; LMFA02000035';
m.metLIPIDMAPSID{1216}='LMFA01070018; LMFA02000037';
m.metLIPIDMAPSID{606} ='LMST01010086; LMST01010144';


% get the index PAPs and update it formula, which was fixed in #81
metsInd = find(strcmp(m.metHMRID, 'm02682'));
m.metCuratedFormulas(metsInd) = {'C10H11N5O13P2S'};


% save updated information back to metAssocHMR2Recon3
metAssocHMR2Recon3 = m;
save('metAssocHMR2Recon3.mat','metAssocHMR2Recon3');


%% generate new data structure with comprehensive metabolite associations

% load model
load('humanGEM.mat'); % HumanGEM v1.0.3

% align with the met index in HumanGEM
clear metAssoc;   % clean and use a temp variable
metAssoc.mets = ihuman.mets;

% add a field for met id without compartment id
metAssoc.metsNoComp = regexprep(metAssoc.mets, '.$', '', 'lineanchors');
metAssoc.metsNoComp = regexprep(metAssoc.metsNoComp,'\_$',''); % Recon3D

% prepare associations to a list of external sources
metAssoc.metBiGGID       = repmat({''},size(metAssoc.mets));   % BiGG
metAssoc.metKEGGID       = repmat({''},size(metAssoc.mets));   % KEGG
metAssoc.metHMDBID       = repmat({''},size(metAssoc.mets));   % HMDB
metAssoc.metChEBIID      = repmat({''},size(metAssoc.mets));   % ChEBI
metAssoc.metPubChemID    = repmat({''},size(metAssoc.mets));   % PubChem
metAssoc.metLipidMapsID  = repmat({''},size(metAssoc.mets));   % LIPIDMAPS
metAssoc.metEHMNID       = repmat({''},size(metAssoc.mets));   % EHMN
metAssoc.metHepatoNET1ID = repmat({''},size(metAssoc.mets));   % HepatoNET1
metAssoc.metRecon3DID    = repmat({''},size(metAssoc.mets));   % Recon3D
metAssoc.metMNXID        = repmat({''},size(metAssoc.mets));   % MetaNetX


%% extract met association from metAssocHMR2Recon3

[a, b] = ismember(metAssoc.metsNoComp, m.metHMRID);
ind = find(a);   % index in Human1
index = b(ind);  % index in metAssocHMR2Recon3

% mainly use the information prepared for HMR2 here
metAssoc.metBiGGID(ind)       = m.metBiGGID(index);           % BiGG
metAssoc.metKEGGID(ind)       = m.metKEGGID(index);           % KEGG
metAssoc.metHMDBID(ind)       = m.metHMDBID(index);           % HMDB
metAssoc.metChEBIID(ind)      = m.metChEBIID(index);          % ChEBI
metAssoc.metPubChemID(ind)    = m.metR3DPubChemID(index);     % PubChem - Recon3D
metAssoc.metLipidMapsID(ind)  = m.metLIPIDMAPSID(index);      % LIPIDMAPS
metAssoc.metEHMNID(ind)       = m.metEHMNID(index);           % EHMN
metAssoc.metHepatoNET1ID(ind) = m.metHepatoNET1ID(index);     % HepatoNET1
metAssoc.metRecon3DID(ind)    = m.metR3DID(index);            % Recon3D
metAssoc.metMNXID(ind)        = m.metCuratedMNXID(index);     % MetaNetX

% combine the KEGG, HMDB, ChEBI and MNX associations provided by Recon3D

% compare associtions between HMR2 and Recon3D
matchKEGG  = cellfun(@strcmp, m.metKEGGID(index), m.metR3DKEGGID(index));
matchHMDB  = cellfun(@strcmp, m.metHMDBID(index), m.metR3DHMDBID(index));
matchChEBI = cellfun(@strcmp, m.metChEBIID(index), m.metR3DCHEBIID(index));

% get the indexes of newly curated associations in Recon3D
curatedKEGGInd  = intersect(getNonEmptyList(m.metR3DKEGGID(index)), find(~matchKEGG));   % 163
curatedHMDBInd  = intersect(getNonEmptyList(m.metR3DHMDBID(index)), find(~matchHMDB));   % 769
curatedChEBIInd = intersect(getNonEmptyList(m.metR3DCHEBIID(index)), find(~matchChEBI)); % 943
curatedMNXInd   = find(strcmp(m.metCuratedMNXID(index), 'toBeChecked'));       % there are 164 conflict ones

% generate intermediate results for manual inspection
%compareKEGG  = [m.metKEGGID(index(curatedKEGGInd)), m.metR3DKEGGID(index(curatedKEGGInd))];
%compareHMDB  = [m.metHMDBID(index(curatedHMDBInd)), m.metR3DHMDBID(index(curatedHMDBInd))];
%compareChEBI = [m.metChEBIID(index(curatedChEBIInd)), m.metR3DCHEBIID(index(curatedChEBIInd))];
%compareMNX   = [m.metMNXID(index(curatedMNXInd)), m.metR3DMNXID(index(curatedMNXInd)), m.metCuratedMNXID(index(curatedMNXInd))];

% update with the new associations from Recon3D
metAssoc.metKEGGID(ind(curatedKEGGInd))   = m.metR3DKEGGID(index(curatedKEGGInd));
metAssoc.metHMDBID(ind(curatedHMDBInd))   = m.metR3DHMDBID(index(curatedHMDBInd));
metAssoc.metChEBIID(ind(curatedChEBIInd)) = m.metR3DCHEBIID(index(curatedChEBIInd));
metAssoc.metMNXID(ind(curatedMNXInd))     = m.metR3DMNXID(index(curatedMNXInd));

% confirm the changes are correctly made
%isequal(metAssoc.metKEGGID(ind(curatedKEGGInd)), m.metR3DKEGGID(index(curatedKEGGInd)))
%isequal(metAssoc.metHMDBID(ind(curatedHMDBInd)), m.metR3DHMDBID(index(curatedHMDBInd)))
%isequal(metAssoc.metChEBIID(ind(curatedChEBIInd)), m.metR3DCHEBIID(index(curatedChEBIInd)))
%isequal(metAssoc.metMNXID(ind(curatedMNXInd)), m.metR3DMNXID(index(curatedMNXInd)))


%% retrieve additional associations from Recon3Mets2MNX prepared in #6, #8

load('Recon3Mets2MNX.mat');

% resolve the conflicts caused by mismatch of comp ids between HMR and Recon
% first modify Recon compartment ids according to HMR2: e,x -> s,p
Recon3Mets2MNX.metsNew = Recon3Mets2MNX.mets;
Recon3Mets2MNX.metsNew = regexprep(Recon3Mets2MNX.metsNew,'\_e$','\_s');
Recon3Mets2MNX.metsNew = regexprep(Recon3Mets2MNX.metsNew,'\_x$','\_p');
% then convert boundary mets to their extracellular counterparts: x -> s
metsNoBoundary = metAssoc.mets;
metsNoBoundary = regexprep(metsNoBoundary,'\_x$','\_s');

% get the index of unique R3D mets
ind_noMatch = find(~ismember(metAssoc.metsNoComp, m.metHMRID));

% find out the complete association
[c, d] = ismember(metsNoBoundary(ind_noMatch), Recon3Mets2MNX.metsNew);
new_ind = d(find(c));

% extract met associations for unique Recon3D mets
metAssoc.metBiGGID(ind_noMatch(c))       = Recon3Mets2MNX.metBiGGDB2BiGG(new_ind);   % BiGG DB ids
metAssoc.metKEGGID(ind_noMatch(c))       = Recon3Mets2MNX.metKEGGID(new_ind);        % KEGG
metAssoc.metHMDBID(ind_noMatch(c))       = Recon3Mets2MNX.metHMDBID(new_ind);        % HMDB
metAssoc.metChEBIID(ind_noMatch(c))      = Recon3Mets2MNX.metChEBIID(new_ind);       % ChEBI
metAssoc.metPubChemID(ind_noMatch(c))    = Recon3Mets2MNX.metPubChemID(new_ind);     % PubChem
metAssoc.metEHMNID(ind_noMatch(c))       = Recon3Mets2MNX.metEHMNID(new_ind);        % EHMN
metAssoc.metHepatoNET1ID(ind_noMatch(c)) = Recon3Mets2MNX.metHepatoNET1ID(new_ind);  % HepatoNET1
metAssoc.metRecon3DID(ind_noMatch(c))    = Recon3Mets2MNX.mets(new_ind);             % Recon3D
metAssoc.metMNXID(ind_noMatch(c))        = Recon3Mets2MNX.metBiGGDB2MNX(new_ind);    % MNX id through BiGG DB


%% Output rxn association in JSON format

jsonStr = jsonencode(metAssoc);
fid = fopen('humanGEMMetAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);

% check the content of JSON file
check = jsondecode(fileread('humanGEMMetAssoc.JSON'));
if isequal(metAssoc, check)
    fprintf('\nThe metabolite association file is sucessfully exported!\n\n');
end

movefile('humanGEMMetAssoc.JSON','../../ComplementaryData/annotation');
clear;

