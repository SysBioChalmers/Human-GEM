function mappedModel = mapModelMets(model,mnx)
%mapModelMets  Retrieve and assign standard IDs to model metabolites.
%
% USAGE:
%
%   mappedModel = mapModelMets(model,mnx);
%
% INPUTS:
%
%   model   A genome scale model structure, containing metabolite related
%           fields (e.g., mets, metNames, etc.).
%
%   mnx     An MNX metabolite database structure generated using the
%           buildMNXmodel('met') function.
%
% OUTPUS:
%
%   mappedModel   The model returned with additional metabolite-related
%                 fields, associating various metabolite IDs to MNX IDs, as
%                 well as a metMNXID field which combines the MNX IDs
%                 obtained by mapping metabolites along each of its
%                 available ID fields.
%                 *NOTE: Some of the met-related fields in mappedModel may
%                        contain multiple columns, which contain the
%                        multiple IDs that matched to one or more of the
%                        metabolites.
%

% handle input arguments
if nargin < 2
    mnx = [];
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
            metNames = [model.metNames,model.metNamesAlt];
        else
            metNames = model.metNames;
        end
        
        % ignore case
        metNames = lower(metNames);
        mnx.mnxID2name(:,2) = lower(mnx.mnxID2name(:,2));
        
        % Perform the name-matching process twice. The first pass will
        % search for exact name matches (ignoring case), whereas the second
        % pass will loosen the criteria by removing all special characters
        % (e.g., hyphens, parentheses, spaces, etc.) from the met names,
        % and search again for any mets that were not matched during the
        % first pass.
        for ii = 1:2
           
            if ii == 2
                % remove special characters from metabolite names
                metNames = regexprep(metNames,'[^a-zA-Z0-9]','');
                mnx.mnxID2name(:,2) = regexprep(mnx.mnxID2name(:,2),'[^a-zA-Z0-9]','');
            else
                model.metName2MNX = repmat({''},size(model.mets));
            end
            
            % extract subset of MNXID-name pairs containing matching names (for faster processing)
            keep_ind = ismember(mnx.mnxID2name(:,2),metNames);
            mnx_ids = mnx.mnxID2name(keep_ind,1);
            mnx_names = mnx.mnxID2name(keep_ind,2);
            
            % convert metNames from cell array to column vector of nested cells (for next processing step)
            metNamesNest = nestCell(metNames,true);
            ignore_inds = ~cellfun(@isempty,model.metName2MNX) | cellfun(@isempty,metNamesNest);
            
            % retrieve all matching IDs for each met
            model.metName2MNX(~ignore_inds) = cellfun(@(x) mnx_ids(ismember(mnx_names,x)),metNamesNest(~ignore_inds),'UniformOutput',false);
            
        end
        
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









