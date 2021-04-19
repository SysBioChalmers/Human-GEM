%
% FILE NAME:    curateHMR2Mets.m
%
% PURPOSE: Curate metabolite information in humanGEM for protonation state
%          at physiological pH and associating with external identifiers
%

%............ Obtain information for mets originating from HMR ............

load('metAssocHMR2Recon3.mat');
m=metAssocHMR2Recon3;   % assign a new name

% load HMR model for metabolite information and association
load('ihumanMets2MNX_v2.mat');  % loads as variable "ihuman"
ihuman.metsNoComp = regexprep(ihuman.mets,'\w$','');
[~, ind]=ismember(m.metHMRID,ihuman.metsNoComp);
m.metFormulas = ihuman.metFormulas(ind);           % formulas

% Include exteranl metabolite identifiers 
m.metLIPIDMAPSID  = ihuman.metLIPIDMAPSID(ind);             % LipidMap
m.metEHMNID       = ihuman.metEHMNID(ind);                  % EHMN
m.metKEGGID       = ihuman.metKEGGID(ind);                  % KEGG
m.metHMDBID       = ihuman.metHMDBID(ind);                  % HMDB
m.metHepatoNET1ID = ihuman.metHepatoNET1ID(ind);            % HepatoNet1
m.metChEBIID      = ihuman.metChEBIID(ind);                 % ChEBI
m.metInChI        = repmat({''},size(m.metHMRID));          % InChI
m.metMNXID        = reformatElements(m.metMNXID,'cell2str');% MetaNetX


%.......... Obtain information for mets originating from Recon3D ..........

% load Recon3D model and associations
load('Recon3Mets2MNX.mat');  % loads as variable "Recon3D"
Recon3D.metsNoComp = regexprep(Recon3D.mets,'\_\w$','');  % remove compartment abbrevs from Recon3D met IDs
Recon3D.metFormulas = regexprep(Recon3D.metFormulas,'FULLR','R');  % also replace FULLR with R in met formulas

% retrieve metabolite information from Recon3D
m.metR3DID        = m.metRecon3DID;                 % mets
m.metR3DNames     = repmat({''},size(m.metHMRID));  % metNames
m.metR3DFormulas  = repmat({''},size(m.metHMRID));  % formulas
m.metR3DCharges   = repmat({''},size(m.metHMRID));  % charges
m.metR3DSmiles    = repmat({''},size(m.metHMRID));  % Smiles
m.metR3DHMDBID    = repmat({''},size(m.metHMRID));  % HMDB
m.metR3DInChI     = repmat({''},size(m.metHMRID));  % InChI
m.metR3DKEGGID    = repmat({''},size(m.metHMRID));  % KEGG
m.metR3DPubChemID = repmat({''},size(m.metHMRID));  % PubChem
m.metR3DCHEBIID   = repmat({''},size(m.metHMRID));  % ChEBI
m.metR3DMNXID     = repmat({''},size(m.metHMRID));  % MetaNetX
m=rmfield(m, 'metRecon3DID');


% Resolving associations
tmp=reformatElements(m.metR3DID,'cell2str');   % parepare the Recon3D IDs

% A. Deal with uniquely mapped ids at first
uniqueInd = find(cellfun(@numel,m.metR3DID)==1);
[a, b]=ismember(tmp(uniqueInd),Recon3D.metsNoComp);
I=find(a);
m.metR3DNames(uniqueInd(I)) = Recon3D.metNames(b(I));
m.metR3DFormulas(uniqueInd(I)) = Recon3D.metFormulas(b(I));
m.metR3DCharges(uniqueInd(I)) = num2cell(Recon3D.metCharges(b(I)));
% make metCharges as a cell, so we can have empty entries for unknown charges

% Include exteranl metabolite identifiers
m.metR3DSmiles(uniqueInd(I)) = Recon3D.metSMILES(b(I));
m.metR3DHMDBID(uniqueInd(I)) = Recon3D.metHMDBID(b(I));
m.metR3DInChI(uniqueInd(I)) = Recon3D.metInChI(b(I));
m.metR3DKEGGID(uniqueInd(I)) = Recon3D.metKEGGID(b(I));
m.metR3DPubChemID(uniqueInd(I)) = Recon3D.metPubChemID(b(I));
m.metR3DCHEBIID(uniqueInd(I)) = Recon3D.metChEBIID(b(I));

% Use the metMNXIDs obtained from BiGG DB mapping, unless the association is
% missing, in which case the MNXID(s) obtained via metName and external IDs
% will be used.
BiGGDB2MNX=Recon3D.metBiGGDB2MNX(b(I));
MNXID=Recon3D.metMNXID(b(I));
empty_ind = cellfun(@isempty,BiGGDB2MNX);
BiGGDB2MNX(empty_ind) = MNXID(empty_ind);
m.metR3DMNXID(uniqueInd(I)) = BiGGDB2MNX;

% There are 2 associated Recon3D met ids that aren't found in Recon3D
noHitInd=find(a==0);
m.metHMRID{uniqueInd(noHitInd)}  % HMR met id
% m00077
% m01422
tmp{uniqueInd(noHitInd)}   % Associasted Recon3D met id that are missing from Recon3D!
% CE2416   % This is a EHMN and BiGG met id
% cbtnCCP  % This is a BiGG met id
% Manual correction: empty the Recon3D ids and update to BiGG ids
m.metR3DID{uniqueInd(noHitInd(1))}='';
m.metR3DID{uniqueInd(noHitInd(2))}='';
m.metBiGGID{uniqueInd(noHitInd(1))}='CE2416';
m.metBiGGID{uniqueInd(noHitInd(2))}='cbtnCCP';


% B. Resolve multiplely mapped ids (a lot of exteranl ids, diffcult)
multiInd = find(cellfun(@numel,m.metR3DID)>1);

%fid = fopen('metCuration_HMR2MultiRecon3D_20180910.tsv','w');
%fprintf(fid,['HMRID\tRecon3DID\tmetName\tFormulas\tCharges\tHMDB\tKEGG\tPubChem\tChEBI\tMNX\tBiGG2MNX\tSmiles\tInChI\n']);
%for i=1:length(multiInd)
%		o=multiInd(i);
%		multiRecon3DID=split(tmp{o},';');
%		[c, d]=ismember(multiRecon3DID,Recon3D.metsNoComp);
%		if all(c)
%				% Output for manual check
%				for j=1:length(c)
%						fprintf(fid,'%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',m.metHMRID{o},Recon3D.metsNoComp{d(j)},Recon3D.metNames{d(j)},Recon3D.metFormulas{d(j)},Recon3D.metCharges(d(j)),Recon3D.metHMDBID{d(j)},Recon3D.metKEGGID{d(j)},Recon3D.metPubChemID{d(j)},Recon3D.metChEBIID{d(j)},Recon3D.metMNXID{d(j)},Recon3D.metBiGGDB2MNX{d(j)},Recon3D.metSMILES{d(j)},Recon3D.metInChI{d(j)});
%				end
%		else
%				% fprintf('These cases are not found\n');
%		end
%end
%fclose(fid);

% Update the manual curation results into the array structure
curatedResults=....
{'m00095','c226coa';
'm00099','CE0695'
'm00118','CE0784'
'm00196','M00196'
'm00554','CE5101'
'm00618','cholcoas'
'm00886','CE0782'
'm00894','CE0853'
'm01026','oretn'
'm01448','xoltri27'
'm02839','HC02187'
'm03023','CE2594'};

for i=1:numel(curatedResults(:,1))
		indHMR=find(strcmp(m.metHMRID,curatedResults{i,1}));
		[~, indR3D]=ismember(curatedResults{i,2},Recon3D.metsNoComp);
		m.metR3DNames{indHMR}=Recon3D.metNames{indR3D};         % metNames
		m.metR3DFormulas{indHMR}=Recon3D.metFormulas{indR3D};   % formulas
		m.metR3DCharges(indHMR)=num2cell(Recon3D.metCharges(indR3D));    % charges
		m.metR3DSmiles{indHMR}=Recon3D.metSMILES{indR3D};       % Smiles
		m.metR3DHMDBID{indHMR}=Recon3D.metHMDBID{indR3D};       % HMDB
		m.metR3DInChI{indHMR}=Recon3D.metInChI{indR3D};         % InChI
		m.metR3DKEGGID{indHMR}=Recon3D.metKEGGID{indR3D};       % KEGG
		m.metR3DPubChemID{indHMR}=Recon3D.metPubChemID{indR3D}; % PubChem
		m.metR3DCHEBIID{indHMR}=Recon3D.metChEBIID{indR3D};     % ChEBI
		m.metR3DMNXID{indHMR}=Recon3D.metMNXID{indR3D};         % MetaNetX
end

% Fix some MetaNetX associations based on above manual curation results
m.metR3DMNXID{find(strcmp(m.metHMRID,'m00095'))}='MNXM3234; MNXM91778';        % m00095
m.metR3DMNXID{find(strcmp(m.metHMRID,'m02839'))}='MNXM162627; MNXM690';        % m02839


% C. Deal with the mets without Recon3D association
nullInd = find(cellfun(@numel,m.metR3DID)==0);   %num = 19

% Output HMR mets without Recon3D association for manual curation
%HMRChEBIID=reformatElements(m.metChEBIID,'cell2str','; ');  % parepare ChEBI IDs
%HMRMNXID=reformatElements(m.metMNXID,'cell2str','; ');      % parepare MetaNetX IDs

%fid = fopen('metCuration_NoRecon3DAssoc_20180911.tsv','w');
%fprintf(fid,['HMRID\tmetName\tFormulas\tLipidMap\tEHMN\tBiGG\tHMDB\tKEGG\tHepatoNet1\tChEBI\tMNX\n']);
%for i=1:numel(nullInd)
%		indHMR=nullInd(i);
%		fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',m.metHMRID{indHMR},m.metNames{indHMR},m.metFormulas{indHMR},m.metLIPIDMAPSID{indHMR},m.metEHMNID{indHMR},m.metBiGGID{indHMR},m.metHMDBID{indHMR},m.metKEGGID{indHMR},m.metHepatoNET1ID{indHMR},HMRChEBIID{indHMR},HMRMNXID{indHMR});
%end
%fclose(fid);

% Update above manual curations into the cell arrays
curatedCharges={-1;-1;-1;-5;-1;0;0;-2;-2;0;0;0;0;-1;0;0;0;0;0};
curatedFormulas={'C19H37O2';'C26H35O8';'C20H29O4';'C11H14O19P3R2';....
    'C20H31O3';'C18H26N5O6R2S';'';'C6H11O9P';'C66H111N4O38RCO';'HO';....
    'C25H47NO4';'C6H10N2O2S2R4';'H3N';'C25H45NO11SR';'X';'';'';'';''};
curatedMNXIDs={'MNXM165274';'';'MNXM33400';'MNXM170';'MNXM6760';....
    'MNXM5655';'MNXM9270';'MNXM336';'MNXM11644';'MNXM56888';....
    'MNXM65475';'MNXM96993';'';'MNXM1234';'MNXM165176';'';'MNXM7010';'';'';};

% Update metCharges and metFormulas based on met associations to
% Recon3D assuming mass/charge balance has been resolved there
m.metCuratedCharges=m.metR3DCharges;            % add curatedCharges field
m.metCuratedCharges(nullInd)=curatedCharges;    % update curated charges
m.metCuratedFormulas=m.metR3DFormulas;          % add curatedFormulas field
m.metCuratedFormulas(nullInd)=curatedFormulas;  % update curated formulas

% Add field for curating associated MNXIDs
m.metCuratedMNXID=repmat({''},size(m.metHMRID));% add curatedMNXID field
m.metCuratedMNXID(nullInd)=curatedMNXIDs;       % update MNX ids for non-associated mets
indOthers=setdiff(transpose(1:numel(m.metHMRID)), nullInd);   % Index of the rest MNXIDs
% Get the index of already matched MNXIDs and update as curated MNXIDs
matchedInd=find(strcmp(m.metMNXID(indOthers),m.metR3DMNXID(indOthers)));
m.metCuratedMNXID(indOthers(matchedInd))=m.metMNXID(indOthers(matchedInd));

% Deal with the unmatched MNXIDs by checking MetaNetX database
%load('MNXMets.mat');   % load MNX met infomation
unmatchedInd=indOthers(setdiff(transpose(1:length(indOthers)),matchedInd));   % Unmatched index
HMRMNXID=reformatElements(m.metMNXID(unmatchedInd),'str2cell');
R3DMNXID=reformatElements(m.metR3DMNXID(unmatchedInd),'str2cell','; ');
newMNXID=repmat({''},size(HMRMNXID));
for k = 1:length(unmatchedInd)
		p = unmatchedInd(k);
		overlap=intersect(HMRMNXID{k},R3DMNXID{k});       % intersection of MNX ids
		aggregate=unique([HMRMNXID{k},R3DMNXID{k}]);      % aggregate of MNX ids
		indMNX=find(ismember(MNXMets.mets, aggregate));   % index to MNX database
		chargeValues=num2cell(MNXMets.metCharges(indMNX));% metCharges in cell
		if isempty(m.metMNXID{p}) && ~isempty(m.metR3DMNXID{p})
				newMNXID{k}=R3DMNXID{k};
		elseif ~isempty(m.metMNXID{p}) && isempty(m.metR3DMNXID{p})
				newMNXID{k}=HMRMNXID{k};
		elseif ~isempty(overlap)
        if isequal(MNXMets.metFormulas{indMNX}) && isequal(chargeValues{:})
        		newMNXID{k}=aggregate;
        else
        		newMNXID{k}=overlap;
    		end
		elseif isempty(overlap)
        if isequal(MNXMets.metFormulas{indMNX}) && isequal(chargeValues{:})
        		newMNXID{k}=aggregate;
        else
        		newMNXID{k}{1}='toBeChecked';   % 69 cases
    		end
		end
end
m.metCuratedMNXID(unmatchedInd)=reformatElements(newMNXID,'cell2str');

% Add two additional met association to Recon3D (duplicate mets in HMR2)
m.metR3DID{find(strcmp('m00555',m.metHMRID))}{1}='pail35p_hs';
m.metR3DID{find(strcmp('m02487',m.metHMRID))}{1}='trdrd';

% Save back to metAssocHMR2Recon3.mat
metAssocHMR2Recon3=m;
save('metAssocHMR2Recon3.mat','metAssocHMR2Recon3');
%..........................................................................

% Also sync the curation info with ihumanMets2MNX_v2.mat
% A. Deal with uniquely mapped ids at first
ihuman.metBiGGID{find(strcmp('m00077p',ihuman.mets))}='CE2416';
ihuman.metBiGGID{find(strcmp('m01422c',ihuman.mets))}='cbtnCCP';
ihuman.metRecon3DID{find(strcmp('m00077p',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m01422c',ihuman.mets))}='';
% C. Deal with the mets without Recon3D association
ihuman.metRecon3DID{find(strcmp('m00555c',ihuman.mets))}{1}='pail35p_hs';
ihuman.metRecon3DID{find(strcmp('m00555g',ihuman.mets))}{1}='pail35p_hs';
ihuman.metRecon3DID{find(strcmp('m00555r',ihuman.mets))}{1}='pail35p_hs';
ihuman.metRecon3DID{find(strcmp('m02487c',ihuman.mets))}{1}='trdrd';
ihuman.metRecon3DID{find(strcmp('m02487m',ihuman.mets))}{1}='trdrd';

save('ihumanMets2MNX_v2.mat','ihuman');  % 2018-09-20

