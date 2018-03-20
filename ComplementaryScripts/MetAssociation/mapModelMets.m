function mappedModel = mapModelMets(model,mnx)
%mapModelMets  Retrieve and assign standard IDs to model metabolites.
%
%
%
%
%
% Jonathan Robinson, 2018-05-19

% handle input arguments
if nargin < 2
    mnx = [];
end


% HMR2-specific section: import data from original model excel sheet
if isfield(model,'id') && strcmp(model.id,'HMRdatabase')
    
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');  % disable warning
    mdata = readtable('HMRdatabase2_00.xlsx','Sheet','METS');
    warning('on','MATLAB:table:ModifiedAndSavedVarnames');  % re-enable warning
    
    % remove table rows with "#" in first column
    mdata(ismember(mdata.x_,'#'),:) = [];
    
    
    % add new met fields to model, or merge ID data with existing fields
%     model = addAltsToField(model,'metFormulas',mdata.COMPOSITION);
    model = addAltsToField(model,'metLIPIDMAPSID',mdata.LM_ID);
    model = addAltsToField(model,'metBiGGID',mdata.BIGGID);
    model = addAltsToField(model,'metEHMNID',mdata.EHMNID);
    model = addAltsToField(model,'metKEGGID',mdata.KEGG_ID);
    model = addAltsToField(model,'metHMDBID',mdata.HMDB_ID);
    model = addAltsToField(model,'metHepatoNETID',mdata.HepatoNETID);
    
    % merge "systematic name" and "synonyms" in new "metNamesAlt" field
    model = addAltsToField(model,'metNamesAlt',mdata.SYSTEMATIC_NAME);
    if ismember('Sedoheptulose 1-phosphate;',model.metNamesAlt)
        % there is one metabolite in this field with a trailing semi-colon
        model.metNamesAlt(ismember(model.metNamesAlt,'Sedoheptulose 1-phosphate;')) = {'Sedoheptulose 1-phosphate'};
    end
    % The synonyms field contains multiple entries for some mets, separated by
    % a semicolon. These need to be split into separate columns before
    % appending to the metNamesAlt field.
    metSynon = cellfun(@(x) strsplit(x,'; '),mdata.SYNONYMS,'UniformOutput',false);
    metSynon = flattenCell(metSynon,true);  % flatten cell array
    model = addAltsToField(model,'metNamesAlt',metSynon);
    
    % note that there are two ChEBI ID columns in the spreadsheet
    % remove preceing "CHEBI:" string from ChEBI IDs
    mdata.CHEBI_ID = regexprep(mdata.CHEBI_ID,'CHEBI:','');
    mdata.CHEBI_ID_1 = regexprep(mdata.CHEBI_ID_1,'CHEBI:','');
    % add to model
    model = addAltsToField(model,'metChEBIID',mdata.CHEBI_ID);
    model = addAltsToField(model,'metChEBIID',mdata.CHEBI_ID_1);
    
end

% get list of metID fields
ignoreFields = {'metFormulas','metMiriams','metComps','mets', ...
                'metNamesAlt','metCharges','metSMILES','metPdMap'};
metIDfields = fields(model);
metIDfields(~startsWith(metIDfields,'met') | ismember(lower(metIDfields),lower(ignoreFields))) = [];

% load metabolite information from MNX database file
if isempty(mnx)
    mnx = buildMNXmodel('met');
end

% associate each set of IDs to MNX IDs
fprintf('Mapping metabolite external IDs to MNX IDs:\n');
for i = 1:length(metIDfields)
    
    % skip field if it isn't present in the MNX database model structure
    if ~isfield(mnx,metIDfields{i})
        continue
    end
    
    fprintf('\t%s\n',metIDfields{i});
    if strcmp(metIDfields{i},'metNames')  % the 'metNames' field is handled differently than others
        
        % combine names and alternative names into single cell array
        if isfield(model,'metNamesAlt')
            metNames = lower([model.metNames,model.metNamesAlt]);
        else
            metNames = lower(model.metNames);
        end
        
        % extract subset of MNXID-name pairs containing matching names (for faster processing)
        mnx.mnxID2name(:,2) = lower(mnx.mnxID2name(:,2));  % make names lowercase
        keep_ind = ismember(mnx.mnxID2name(:,2),metNames);
        mnx_ids = mnx.mnxID2name(keep_ind,1);
        mnx_names = mnx.mnxID2name(keep_ind,2);
        
        % convert metNames from cell array to column vector of nested cells (for next processing step)
        metNames = nestCell(metNames,true);
        empty_inds = cellfun(@isempty,metNames);
        
        % retrieve all matching IDs for each met
        model.metName2MNX = repmat({''},size(model.mets));
        model.metName2MNX(~empty_inds) = cellfun(@(x) mnx_ids(ismember(mnx_names,x)),metNames(~empty_inds),'UniformOutput',false);
        
        % remove entries that have matched to too many IDs (>100)
        model.metName2MNX(cellfun(@numel,model.metName2MNX) > 100) = {''};
        
        % flatten cell array
        model.metName2MNX = flattenCell(model.metName2MNX,true);
        
    else
        % get model met IDs, and compress each row into nested cells
        % (this is to deal with fields that have multiple columns)
        model_ids = nestCell(model.(metIDfields{i}),true);
        
        % extract only subset of MNX model containing matching IDs (for faster processing)
        keep_ind = ismember(mnx.mnxID2extID(:,2),metIDfields(i)) & ismember(mnx.mnxID2extID(:,3),model.(metIDfields{i}));
        ext_ids = mnx.mnxID2extID(keep_ind,3);
        mnx_ids = mnx.mnxID2extID(keep_ind,1);
        
        % get empty indices for current field
        empty_inds = cellfun(@isempty,model_ids);
        
        % retrieve all matching IDs for each met
        newField = strcat(metIDfields{i},'2MNX');
        model.(newField) = repmat({''},size(model.mets));
        model.(newField)(~empty_inds) = cellfun(@(x) mnx_ids(ismember(ext_ids,x)),model_ids(~empty_inds),'UniformOutput',false);
        
        % flatten cell array
        model.(newField) = flattenCell(model.(newField),true);
    end
    
end
fprintf('Done.\n');


% now combine all the MNX IDs for each metabolite
mnxIDfields = fields(model);
mnxIDfields(~(startsWith(mnxIDfields,'met') & endsWith(mnxIDfields,'2MNX'))) = [];
metMNXIDs = {};  % intialize cell array of met MNX IDs
for i = 1:length(mnxIDfields)
    % append each field as new column(s)
    metMNXIDs = [metMNXIDs,model.(mnxIDfields{i})];
end
empty_inds = cellfun(@isempty,metMNXIDs);
% metMNXIDs(empty_inds) = {''};


% obtain unique set of MNX IDs for each metabolite
met_index = transpose(1:length(model.mets));
model.metMNXID = arrayfun(@(i) unique(metMNXIDs(i,~empty_inds(i,:))),met_index,'UniformOutput',false);
model.metMNXID = flattenCell(model.metMNXID,true);

% assign output
mappedModel = model;

end





function newModel = addAltsToField(model,field,newEntries)
% Compares and appends a new set of field values to an existing model
% field. If an existing field entry is empty, it will be overwritten by the
% new entry. If the existing field entry and new entry conflict, both will
% be saved by adding new columns to the model field.

    if ~isfield(model,field)
        % if the field doesn't yet exist in the model, just add the new entries
        model.(field) = newEntries;
        newModel = model;
        return
    elseif isequal(model.(field),newEntries)
        % no changes needed if new entries are identical to existing entries
        newModel = model;
        return
    end

    if (size(model.(field),2) == 1) && (size(newEntries,2) == 1)
        % if existing and new entries are both column vectors, simply add the
        % new (mismatching) entries as a second column
        mismatch_ind = ~strcmp(model.(field),newEntries);
        newEntries(~mismatch_ind) = {''};
        model.(field) = [model.(field),newEntries];
    else
        for i = 1:size(newEntries,1)
            vals = [model.(field)(i,:),newEntries(i,:)];  % combine rows
            vals(cellfun(@isempty,vals)) = [];  % remove empty entries
            vals = unique(vals);  % get all unique entries for row
            model.(field)(i,1:length(vals)) = vals;  % update row in model
        end
    end
    
    % replace empty matrices with empty strings
    model.(field)(cellfun(@isempty,model.(field))) = {''};
        
    newModel = model;  % assign output

end






