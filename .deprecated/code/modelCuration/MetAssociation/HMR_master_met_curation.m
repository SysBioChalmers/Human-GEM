% Master script for curating and mapping HMR and Recon metabolites
%


%% Load necessary items

% load MNX metabolite data (takes a few minutes)
mnx_met = buildMNXmodel('met');


%% Load HMR2 model and perform initial updates and ID mapping

% load HMR model (will be loaded as variable "ihuman")
load('ModelFiles/mat/HMRdatabase2_02.mat');

% Load metabolite information from spreadsheet and add to model, and
% perform a few minor corrections to various model attributes.
ihuman = HMR_update_met_attributes(ihuman);

% map metabolites to MNXIDs using names and external IDs (takes 1 min)
ihuman = mapModelMets(ihuman,mnx_met);

% distribute mapped met IDs to identical mets in different compartments
ihuman = spreadInfoAcrossComps(ihuman,'met',true);


%% Load Recon3D model and perform minor changes and ID mapping

% load Recon3D model (loaded as variable "Recon3D")
load('ComplementaryScripts/Recon3D_301.mat');

% rename some of the model fields for consistency
Recon3D = renameStructField(Recon3D,'metInChIString','metInChI');
Recon3D = renameStructField(Recon3D,'metSmiles','metSMILES');
Recon3D = renameStructField(Recon3D,'metCHEBIID','metChEBIID');

% remove the "CHEBI:" preceding some of the ChEBI IDs
Recon3D.metChEBIID = regexprep(Recon3D.metChEBIID,'^CHEBI:','','ignorecase');

% Some of the entries in the "mets" field are actually external IDs. 
% Extract these, and add to the corresponding ID field.
mets = regexprep(Recon3D.mets,'\[.\]$','');  % remove trailing compartment info from "mets" entries

% add missing fields
if ~isfield(Recon3D,'metEHMNID')
    Recon3D.metEHMNID = repmat({''},size(mets));
end
if ~isfield(Recon3D,'metHepatoNET1ID')
    Recon3D.metHepatoNET1ID = repmat({''},size(mets));
end

% check for EHMN IDs (format: CE#### or CN####)
ind = ~cellfun(@isempty,regexp(mets,'^C[EN]\d{4}$','match')) & cellfun(@isempty,Recon3D.metEHMNID);
Recon3D.metEHMNID(ind) = mets(ind);

% check for KEGG IDs (format: C#####, D#####, or G#####)
ind = ~cellfun(@isempty,regexp(mets,'^[CDG]\d{5}$','match')) & cellfun(@isempty,Recon3D.metKEGGID);
Recon3D.metKEGGID(ind) = mets(ind);

% check for HepatoNET1 IDs (format: HC#####)
ind = ~cellfun(@isempty,regexp(mets,'^HC\d{5}$','match')) & cellfun(@isempty,Recon3D.metHepatoNET1ID);
Recon3D.metHepatoNET1ID(ind) = mets(ind);

% all of the "mets" entries are technically BiGG IDs
Recon3D.metBiGGID = mets;


% some metabolites have a weird "return" character in the metName field
% that needs to be removed
ind = ismember(mets,'M01870');
Recon3D.metNames(ind) = {'(GlcNAc)7 (Man)3 (Asn)1'};
ind = ismember(mets,'M00304');
Recon3D.metNames(ind) = {'11-Trans-Leukotriene E4'};

% rename "All-Trans-Retinal" to "9-Cis-Retinal" (appears to be naming error)
Recon3D.metNames(ismember(mets,'retinal_cis_9')) = {'9-Cis-Retinal'};

% split metabolite names by semicolon, and add additional names to a new
% metNamesAlt field
metNames = cellfun(@(n) strsplit(n,'; '),Recon3D.metNames,'UniformOutput',false);
metNames = flattenCell(metNames,true);  % flatten nested cell array
Recon3D.metNames = metNames(:,1);
Recon3D.metNamesAlt = metNames(:,2:end);

% add compartment-related fields to model
metComps = regexprep(Recon3D.mets,'^.*\[|\]$','');
Recon3D.comps = unique(metComps);
[~,Recon3D.metComps] = ismember(metComps,Recon3D.comps);

% remove spaces from model metabolite formulas
Recon3D.metFormulas = regexprep(Recon3D.metFormulas,'\s','');

% reorder metabolite formulas so elements are in alphabetical order
% (this is to be consistent with format of MNX database formulas)
% NOTE: this does not re-order formulas containing special characters (i.e., those other than A-Z,0-9).
Recon3D.metFormulas = alphabetizeMetFormulas(Recon3D.metFormulas);

% map metabolites to MNXIDs using names and external IDs
Recon3D = mapModelMets(Recon3D,mnx_met);

% distribute mapped met IDs to identical mets in different compartments
Recon3D = spreadInfoAcrossComps(Recon3D,'met',true);

% clear intermediate variables that are no longer needed
clear ind metComps metNames mets


%% Perform mapping between HMR and Recon3D metabolites

% Specify fields by which metabolites should be mapped, as well as the
% corresponding weight of each field. Metabolites mapped via a field with
% greater weight will be chosen over a field with lower weight (unless the
% "all" option is specified as an input to the mapMetsToAltModel function,
% in which case all matches along all fields will be included).
mapFields = {'mets', 100
             'metNames', 4
             'metBiGGID', 2
             'metChEBIID', 2
             'metHMDBID', 2
             'metEHMNID', 2
             'metHepatoNET1ID', 2
             'metKEGGID', 2
             'metMNXID', 1};

% map HMR mets to Recon3D mets
ihuman = mapMetsToAltModel(ihuman,Recon3D,mapFields,'score');
ihuman = renameStructField(ihuman,'metAltModelID','metRecon3DID');  % rename field

% map Recon3D mets to HMR mets
Recon3D = mapMetsToAltModel(Recon3D,ihuman,mapFields,'score');
Recon3D = renameStructField(Recon3D,'metAltModelID','metHMRID');  % rename field

% clear intermediate variables
clear mapFields


%% Additional MNX matching

% map HMR mets to additional MNXIDs via their Recon3D met associations
R3IDs = nestCell(ihuman.metRecon3DID,true);
R3mets = regexprep(Recon3D.mets,'\[.\]$','');
R32MNX = cellfun(@(m) unique(Recon3D.metMNXID(ismember(R3mets,m),:)),R3IDs,'UniformOutput',false);
ihuman.metRecon3DID2MNX = flattenCell(nestCell(flattenCell(R32MNX,true),true),true);  % very ugly way to remove empty entries

% combine newly mapped HMR MNXIDs with the HMR metMNXID field
non_empty_orig = ~cellfun(@isempty,ihuman.metMNXID);
non_empty_add = ~cellfun(@isempty,ihuman.metRecon3DID2MNX);
ihuman.metMNXID = arrayfun(@(i) unique([ihuman.metMNXID(i,non_empty_orig(i,:)),ihuman.metRecon3DID2MNX(i,non_empty_add(i,:))]),[1:size(ihuman.metMNXID,1)]','UniformOutput',false);
ihuman.metMNXID = flattenCell(ihuman.metMNXID,true);


% map Recon3D mets to additional MNXIDs via their HMR met associations
HMRIDs = nestCell(Recon3D.metHMRID,true);
HMRmets = regexprep(ihuman.mets,'.$','');
HMR2MNX = cellfun(@(m) unique(ihuman.metMNXID(ismember(HMRmets,m),:)),HMRIDs,'UniformOutput',false);
Recon3D.metHMRID2MNX = flattenCell(nestCell(flattenCell(HMR2MNX,true),true),true);  % very ugly way to remove empty entries

% combine newly mapped HMR MNXIDs with the HMR metMNXID field
non_empty_orig = ~cellfun(@isempty,Recon3D.metMNXID);
non_empty_add = ~cellfun(@isempty,Recon3D.metHMRID2MNX);
Recon3D.metMNXID = arrayfun(@(i) unique([Recon3D.metMNXID(i,non_empty_orig(i,:)),Recon3D.metHMRID2MNX(i,non_empty_add(i,:))]),[1:size(Recon3D.metMNXID,1)]','UniformOutput',false);
Recon3D.metMNXID = flattenCell(Recon3D.metMNXID,true);


% distribute mapped met IDs to identical mets in different compartments
ihuman = spreadInfoAcrossComps(ihuman,'met',true);
Recon3D = spreadInfoAcrossComps(Recon3D,'met',true);

% clear intermediate variables
clear HMR2MNX HMRIDs HMRmets R32MNX R3IDs R3mets non_empty_add non_empty_orig

% *** Note (2018-05-02): at this point, ihuman and Recon3D were saved as
% "ihumanMets2MNX.mat" and "Recon3Mets2MNX.mat", respectively.


%% NOTE: For mapping reactions to MNX via metabolite IDs
% (2018-05-02)
% The following steps were taken, after those above:
%   1) Load "mergedModel.mat"
%   2) Remove all reactions in ihuman that are not present in mergedModel.
%   3) Input the resulting filtered ihuman model into the "mapRxnsViaMets"
%      function.


%% Identify new metMNXIDs through rxnMNXID associations

% first load mergedModel.mat
load('ComplementaryScripts/RxnAssociation/mergedModel.mat');

% extract rxnMNXID associations from mergedModel, and add them to ihuman
[~,ind] = ismember(mergedModel.rxns,ihuman.rxns);
ihuman.rxnMNXID = repmat({''},size(ihuman.rxns));
ihuman.rxnMNXID(ind) = mergedModel.confirmedFilteredMNXID;

% ensure that all rxnMNXIDs are row vectors
ihuman.rxnMNXID = cellfun(@(x) x(:)',ihuman.rxnMNXID,'UniformOutput',false);

% consolidate ihuman.metMNXID field into column of nested cells
ihuman.metMNXID = nestCell(ihuman.metMNXID,true);

% run script to identify potential new metMNXID associations to add
mnx_rxn = buildMNXmodel('rxn');  % load MNX database rxn info
results = compareRxnMNXIDsWithMets(ihuman,mnx_rxn,true);

% *** MANUALLY CURATE metMNXIDs BASED ON RESULTS STRUCTURE ***
% Note: this will need to be updated with any changes to the previous code
addMNXIDinfo = {'m00037', {'MNXM165175', 'MNXM560', 'MNXM3395'}
                'm00075', {'MNXM10022', 'MNXM35866', 'MNXM663', 'MNXM108213'}
                'm00084', {'MNXM1251', 'MNXM163785'}
                'm00134', {'MNXM1203', 'MNXM1513', 'MNXM47387', 'MNXM5306', 'MNXM169215'}
                'm00200', {'MNXM73306', 'MNXM6404'}
                'm00681', {'MNXM163157', 'MNXM145987'}
                'm00989', {'MNXM894', 'MNXM97048'}
                'm01633', {'MNXM786', 'MNXM162802'}
                'm01657', {'MNXM148196', 'MNXM48918', 'MNXM507601', 'MNXM7283', 'MNXM148197'}
                'm01733', {'MNXM167079', 'MNXM278', 'MNXM51501', 'MNXM145597'}
                'm01734', {'MNXM296', 'MNXM7322', 'MNXM145684'}
                'm01995', {'MNXM12747'}
                'm02998', {'MNXM114114', 'MNXM162591'}
                'm03052', {'MNXM1013', 'MNXM162560', 'MNXM162627', 'MNXM690'} };

% add manually-curated updated metMNXID associations to the model
mets_noComp = regexprep(ihuman.mets,'.$','');
for i = 1:size(addMNXIDinfo,1)
    ind = ismember(mets_noComp,addMNXIDinfo(i,1));
    ihuman.metMNXID(ind) = addMNXIDinfo(i,2);
end


%% Identify metMNXIDs that should be removed
% Find metMNXIDs that don't appear in any of the MNX reactions currently
% associated with the model, and remove them from the model.

% run filter function
[fmodel,removed] = filterMetMNXIDsViaRxns(ihuman,mnx_rxn,true);

% remove rxnMNXID field from fmodel (no longer needed)
fmodel = rmfield(fmodel,'rxnMNXID');
  
% consolidate model fields that contain multiple columns into single column
fmodel.metNamesAlt = nestCell(fmodel.metNamesAlt,true);
fmodel.metChEBIID = nestCell(fmodel.metChEBIID,true);
fmodel.metName2MNX = nestCell(fmodel.metName2MNX,true);
fmodel.metChEBIID2MNX = nestCell(fmodel.metChEBIID2MNX,true);
fmodel.metRecon3DID = nestCell(fmodel.metRecon3DID,true);
fmodel.metRecon3DID2MNX = nestCell(fmodel.metRecon3DID2MNX,true);
  
% remove and re-add metMNXID field so it is listed as the last field
metMNXID = fmodel.metMNXID;
fmodel = rmfield(fmodel,'metMNXID');
fmodel.metMNXID = metMNXID;

% rename and save original (pre-filter) model
ihuman_orig = ihuman;
ihuman = fmodel;


%% Update mapping of HMR mets to Recon3D mets

% remove original recon3D met ID field from HMR model
ihuman = rmfield(ihuman,'metRecon3DID');
ihuman = rmfield(ihuman,'metRecon3DID2MNX');

% flatten nested Recon3D met fields
metFields = fields(Recon3D);
metFields = metFields(startsWith(metFields,'met'));
for i = 1:length(metFields)
    if iscell(Recon3D.(metFields{i})) && any(any(cellfun(@iscell,Recon3D.(metFields{i}))))
        Recon3D.(metFields{i}) = flattenCell(Recon3D.(metFields{i}),true);
    end
end
% flatten nested ihuman met fields
metFields = fields(ihuman);
metFields = metFields(startsWith(metFields,'met'));
for i = 1:length(metFields)
    if iscell(ihuman.(metFields{i})) && any(any(cellfun(@iscell,ihuman.(metFields{i}))))
        ihuman.(metFields{i}) = flattenCell(ihuman.(metFields{i}),true);
    end
end

% regenerate mapping to Recon3D mets
mapFields = {'mets', 100
             'metNames', 4
             'metBiGGID', 2
             'metChEBIID', 2
             'metHMDBID', 2
             'metEHMNID', 2
             'metHepatoNET1ID', 2
             'metKEGGID', 2
             'metMNXID', 1};

% map HMR mets to Recon3D mets
ihuman = mapMetsToAltModel(ihuman,Recon3D,mapFields,'score');
ihuman = renameStructField(ihuman,'metAltModelID','metRecon3DID');

% re-consolidate Recon3D met fields into single columns
metFields = fields(Recon3D);
metFields = metFields(startsWith(metFields,'met'));
for i = 1:length(metFields)
    if iscell(Recon3D.(metFields{i})) && size(Recon3D.(metFields{i}),2) > 1
        Recon3D.(metFields{i}) = nestCell(Recon3D.(metFields{i}),true);
    end
end
% re-consolidate ihuman met fields into single columns
metFields = fields(ihuman);
metFields = metFields(startsWith(metFields,'met'));
for i = 1:length(metFields)
    if iscell(ihuman.(metFields{i})) && size(ihuman.(metFields{i}),2) > 1
        ihuman.(metFields{i}) = nestCell(ihuman.(metFields{i}),true);
    end
end

% save latest ihuman model file
save('ihumanMets2MNX_v2.mat','ihuman');



%% Update mapping of Recon3D mets to HMR mets
% (2018-05-30)

% load ihuman model with latest metMNXID mapping (if not yet loaded)
load('ihumanMets2MNX_v2.mat');

% rename Recon3D structure generated from above commands so it is not
% overwritten when loading the Recon3D model to which new fields have been
% added (metBiGGDB2BiGG and metBiGGDB2MMNX)
Recon3D_orig = Recon3D;
load('Recon3Mets2MNX.mat');  % load model

% remove metHMRID2MNX field from Recon3D_orig model, to avoid any confusion
Recon3D_orig = rmfield(Recon3D_orig,'metHMRID2MNX');

% remove the existing metHMRID field from Recon3D_orig, since it will be
% regenerated with the updated Recon3D model (with added MNX associations)
Recon3D_orig = rmfield(Recon3D_orig,'metHMRID');

% consolidate met fields in Recon3D_orig into single columns
Recon3D_orig.metNamesAlt = nestCell(Recon3D_orig.metNamesAlt,true);
Recon3D_orig.metName2MNX = nestCell(Recon3D_orig.metName2MNX,true);
Recon3D_orig.metMNXID = nestCell(Recon3D_orig.metMNXID,true);

% merge associations by importing MNXIDs from Recon3D_orig.metMNXID into 
% ONLY the EMPTY entries in Recon3D.metBiGGDB2MNX
mergedMNXID = Recon3D.metBiGGDB2MNX;
empty_inds = cellfun(@isempty,mergedMNXID);
mergedMNXID(empty_inds) = Recon3D_orig.metMNXID(empty_inds);

% flatten and re-nest array so that all cells are in same format
mergedMNXID = nestCell(flattenCell(mergedMNXID,true),true);

% replace Recon3D_orig.metMNXID field with mergedMNXID, and rename
% structure to Recon3D.
Recon3D_orig.metMNXID = mergedMNXID;
Recon3D = Recon3D_orig;

% flatten nested Recon3D met fields
metFields = fields(Recon3D);
metFields = metFields(startsWith(metFields,'met'));
for i = 1:length(metFields)
    if iscell(Recon3D.(metFields{i})) && any(any(cellfun(@iscell,Recon3D.(metFields{i}))))
        Recon3D.(metFields{i}) = flattenCell(Recon3D.(metFields{i}),true);
    end
end
% flatten nested ihuman met fields
metFields = fields(ihuman);
metFields = metFields(startsWith(metFields,'met'));
for i = 1:length(metFields)
    if iscell(ihuman.(metFields{i})) && any(any(cellfun(@iscell,ihuman.(metFields{i}))))
        ihuman.(metFields{i}) = flattenCell(ihuman.(metFields{i}),true);
    end
end

% Re-associate Recon3D mets to HMR mets
mapFields = {'mets', 100
             'metNames', 4
             'metBiGGID', 2
             'metChEBIID', 2
             'metHMDBID', 2
             'metEHMNID', 2
             'metHepatoNET1ID', 2
             'metKEGGID', 2
             'metMNXID', 1};
         
Recon3D = mapMetsToAltModel(Recon3D,ihuman,mapFields,'score');
Recon3D = renameStructField(Recon3D,'metAltModelID','metHMRID');  % rename field
Recon3D = spreadInfoAcrossComps(Recon3D,'met',true);  % distribute information across compartments

% To determine which HMRIDs are associated to multiple Recon3 mets,
% generate the inverse of the Recon3D-HMR met map.

% get unique metabolite list for Recon3D
mets_noComp = regexprep(Recon3D.mets,'\[\w\]','');
[uniq_R3mets,uniq_ind] = unique(mets_noComp);
uniq_R3hmrID = Recon3D.metHMRID(uniq_ind,:);

% obtain list of all HMR metabolites that have been mapped to Recon3D mets
hmr_mets = unique(Recon3D.metHMRID(:));
hmr_mets(1) = [];  % remove the empty string

% now determine the Recon3D mets that have been mapped to each of these HMR mets
hmr2r3 = {};  % intialize variable
for i = 1:length(hmr_mets)
    r3mets = uniq_R3mets(any(ismember(uniq_R3hmrID,hmr_mets(i)),2));
    hmr2r3(i,1:length(r3mets)) = r3mets;
end
hmr2r3(cellfun(@isempty,hmr2r3)) = {''};


% re-consolidate Recon3D met fields into single columns
metFields = fields(Recon3D);
metFields = metFields(startsWith(metFields,'met'));
for i = 1:length(metFields)
    if iscell(Recon3D.(metFields{i})) && size(Recon3D.(metFields{i}),2) > 1
        Recon3D.(metFields{i}) = nestCell(Recon3D.(metFields{i}),true);
    end
end


% NOTE: at this point, this version of the Recon3D model has not yet been
% used in any further analyses/curation.



%% Further model processing (single compartment, writing to file, etc.)
% NOTE: This section was NOT used for any analyses/curation, only for 
% conventient output of model metabolite information.

%***** SELECT MODEL *****
model = ihuman;
% model = Recon3D;
%************************

% distribute mapped met IDs to identical mets in different compartments
model_spread = spreadInfoAcrossComps(model,'met',true);

% extract unique metabolite information from model (merge compartments)
% WARNING: THIS MODEL WILL LOSE PROPER CONNECTION BETWEEN METABOLITES AND
% REACTIONS - USE ONLY THE METABOLITE INFORMATION.
model_nocomp = model_spread;  % intialize
if strcmp(model_nocomp.mets{1}(end),']')
    % remove compartment abbreviation in brackets from metID
    model_nocomp.mets = regexprep(model_nocomp.mets,'\[.\]$','');
else
    % assume compartment info is single character at the end of metID
    model_nocomp.mets = regexprep(model_nocomp.mets,'.$','');
end
[~,uniq_ind] = unique(model_nocomp.mets);  % get unique metabolite indices
f = fields(model_nocomp);  % get model fields
for i = 1:length(f)
    % find fields corresponding to met information, and only keep those
    % corresponding to unique metabolic species
    if size(model_nocomp.(f{i}),1) == length(model.mets)
        model_nocomp.(f{i}) = model_nocomp.(f{i})(uniq_ind,:);
    end
end

%............... write metabolite data to text file ...............
m = {};
metFields = fields(model_nocomp);
metFields(~startsWith(metFields,'met')) = [];
% combine metIDs into single column, separating multiple entries by semicolons
for i = 1:length(metFields)
    if size(model_nocomp.(metFields{i}),2) == 1
        % skip fields that are already a single column
        m.(metFields{i}) = model_nocomp.(metFields{i});
        continue
    end
    ids = model_nocomp.(metFields{i});
    empty_inds = cellfun(@isempty,ids);
    m.(metFields{i}) = arrayfun(@(x) strjoin(ids(x,~empty_inds(x,:)),'; '),1:size(ids,1),'UniformOutput',false)';
end
f = {'mets','metNames','metNamesAlt','metFormulas','metLIPIDMAPSID','metBiGGID',...
    'metEHMNID','metKEGGID','metPubChemID','metHMDBID','metHepatoNET1ID','metChEBIID',...
    'metSMILES','metInChI','metHMRID','metRecon3DID','metMNXID'};
f(~isfield(m,f)) = [];  % remove fields that don't exist in the model
z = {};
for i = 1:length(f)
    z = [z,m.(f{i})];
end
z = [f;z];
writecell2file(z,'metdata.txt',true,'\t');
%..........................................................




