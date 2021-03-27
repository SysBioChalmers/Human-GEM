function model_new = HMR_update_met_attributes(model)
%HMR_update_met_attributes  Add and correct HMR metabolite information.
%
% HMR_update_met_attributes loads metabolite-related information from the 
% HMRdatabase2_00.xls spreadsheet, and performs a variety of miscellaneous
% corrections and formatting changes to the model.
%
% USAGE:
%
%   model_new = HMR_update_met_attributes(model);
%
% INPUT:
%
%   model    HMR2 model structure.
%
% OUTPUT:
%
%   model_new   HMR2 model structure with added metabolite information and
%               minor corrections applied.
%


% reorder metabolite formulas so elements are in alphabetical order
% NOTE: this does not re-order formulas containing special characters 
%       (i.e., those other than A-Z,0-9).
model.metFormulas = alphabetizeMetFormulas(model.metFormulas);

% load data from spreadsheet
warning('off','MATLAB:table:ModifiedAndSavedVarnames');  % disable warning
mdata = readtable('ComplementaryScripts/HMRdatabase2_00.xlsx','Sheet','METS');
warning('on','MATLAB:table:ModifiedAndSavedVarnames');  % re-enable warning

% remove table rows with "#" in first column
mdata(ismember(mdata.x_,'#'),:) = [];

% add new met fields to model, or merge ID data with existing fields
model = addAltsToModelField(model,'metLIPIDMAPSID',mdata.LM_ID);
model = addAltsToModelField(model,'metBiGGID',mdata.BIGGID);
model = addAltsToModelField(model,'metEHMNID',mdata.EHMNID);
model = addAltsToModelField(model,'metKEGGID',mdata.KEGG_ID);
model = addAltsToModelField(model,'metHMDBID',mdata.HMDB_ID);
model = addAltsToModelField(model,'metHepatoNET1ID',mdata.HepatoNETID);

% merge "systematic name" and "synonyms" in new "metNamesAlt" field
model = addAltsToModelField(model,'metNamesAlt',mdata.SYSTEMATIC_NAME);
if ismember('Sedoheptulose 1-phosphate;',model.metNamesAlt)
    % there is one metabolite in this field with a trailing semi-colon
    model.metNamesAlt(ismember(model.metNamesAlt,'Sedoheptulose 1-phosphate;')) = {'Sedoheptulose 1-phosphate'};
end

% The synonyms field contains multiple entries for some mets, separated by
% a semicolon. These need to be split into separate columns before
% appending to the metNamesAlt field.
metSynon = cellfun(@(x) strsplit(x,'; '),mdata.SYNONYMS,'UniformOutput',false);
metSynon = flattenCell(metSynon,true);  % flatten cell array
model = addAltsToModelField(model,'metNamesAlt',metSynon);

% Some of the glycans have problems matching based on their names, so
% add their formulas as an alternative metabolite name
glycan_formulas = repmat({''},size(model.mets));
glycan_ind = startsWith(model.metFormulas,{'(Gal)','(GalNAc)','(Glc)','(GlcNAc)'});
glycan_formulas(glycan_ind) = model.metFormulas(glycan_ind);
model = addAltsToModelField(model,'metNamesAlt',glycan_formulas);

% note that there are two ChEBI ID columns in the spreadsheet
% remove preceing "CHEBI:" string from ChEBI IDs
mdata.CHEBI_ID = regexprep(mdata.CHEBI_ID,'CHEBI:','');
mdata.CHEBI_ID_1 = regexprep(mdata.CHEBI_ID_1,'CHEBI:','');
% add to model
model = addAltsToModelField(model,'metChEBIID',mdata.CHEBI_ID);
model = addAltsToModelField(model,'metChEBIID',mdata.CHEBI_ID_1);


% the metabolite formula of lepidimoide (m02357c and m02357s) contains a
% period "." at the end, which needs to be removed
ind = ismember(model.mets,{'m02357c','m02357s'});
model.metFormulas(ind) = regexprep(model.metFormulas(ind),'\.','');

% update protonation of a few metabolites
ind = contains(model.mets,'m02751');
model.metFormulas(ind) = {'HO4P'};  % originally listed as "H3PO4"
ind = contains(model.mets,'m02949');
model.metFormulas(ind) = {'O3S'};  % originally listed as "H2SO3"

% rename 'temp006x' to 'm01451x' (both correspond to 'cholesterol-ester pool')
model.mets(ismember(model.mets,'temp006x')) = {'m01451x'};

% rename some BiGG IDs to be consistent with the current BiGG database
model.metBiGGID = regexprep(model.metBiGGID,',|-|/','_');  % replace commas, dashes, and slashes with underscores
model.metBiGGID(ismember(model.metBiGGID,'Xyl_L_Ser_(protein)')) = {'xser'};

% assign output
model_new = model;

