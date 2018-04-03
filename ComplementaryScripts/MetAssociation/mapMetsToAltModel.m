function model = mapMetsToAltModel(refModel,mapModel,metFields,prioritize)
%mapMetsToAltModel  Map metabolites from one model to another.
%
% USAGE:
%
%   model = mapMetsToAltModel(refModel,mapModel,metFields,allmatches);
%
%
% INPUTS:
%
%   refModel   Model structure containing metabolites that are to be mapped
%              to mapModel metabolite identifiers.
%
%   mapModel   Model structure containing metabolites to which refModel 
%              will be mapped.
%
%   metFields  One or more metabolite-related fields that will be compared 
%              between the two models to map metabolites from one model to 
%              the other. If "metNames" is included in metFields,
%              metabolite names will be processed in three steps:
%               1) Map refModel metNames to mapModel metNames and
%                  metNamesAlt
%               2) Map refModel metNames and metNamesAlt to mapModel
%                  metNames and metNamesAlt
%               3) Make specific modifications to met names based on known
%                  differences (unique to situations where Recon3D and HMR
%                  are the input models)
%
%   prioritize (Optional, default FALSE) If TRUE, then the function will
%              prioritize fields based on the order in which they occur in
%              metFields, and will not try to map metabolites via a field
%              that occurs later on in the metFields list if they have
%              already been mapped via a field earlier in the list.
%              If FALSE, the function will try to map metabolites via all
%              listed fields in metFields, regardless of the order in which
%              they are processed, and will combine all matches in the
%              resulting metAltModelID field.
%              NOTE: This prioritization is enforced within the metName
%                    mapping stages, where mets that have been mapped in
%                    stage 1 are excluded from mapping in stage 2, and
%                    those mapped in stage 1 or 2 are excluded from mapping
%                    during stage 3.
%
% OUTPUTS:
%
%   model     Model structure with added field "metAltModelID", which 
%             contains the mapModel metIDs that were matched to the 
%             refModel metabolites.
%
%
% Jonathan Robinson 2018-04-03


% handle input arguments
if nargin < 4
    prioritize = false;
end

if any(~isfield(refModel,metFields)) || any(~isfield(mapModel,metFields))
    error('One or more of the specified altFields is not present in one or both models.');
elseif ischar(metFields)
    metFields = {metFields};
else
    metFields = metFields(:);  % force into column vector
end

% remove metNamesAlt if present in metFields - this field is automatically
% included when metNames is specified
if ismember(metFields,'metNamesAlt')
    fprintf('Note: the "metNamesAlt" entry is ignored, because it is included in the "metNames" processing step.\n');
    metFields(ismember(metFields,'metNamesAlt')) = [];
    if ~ismember('metNames',metFields)
        % warn user if they include metNamesAlt, but not metNames, in metFields
        fprintf('WARNING: comparison by "metNames" has not been specified - mets will not be compared by names.\n');
    end
end

% Add empty metNamesAlt field to either of the models if they are missing
% the field. This just makes things easier later on.
if ~isfield(refModel,'metNamesAlt')
    refModel.metNamesAlt = repmat({''},size(refModel.mets));
end
if ~isfield(mapModel,'metNamesAlt')
    mapModel.metNamesAlt = repmat({''},size(mapModel.mets));
end

% initialize output model
model = refModel;

% strip compartments from mapModel met IDs to obtain compartment-free
% met-to-name (or met-to-ID) pairs
if strcmp(mapModel.mets{1}(end),']')
    % compartment information is formatted in brackets at end of ID
    map_mets = regexprep(mapModel.mets,'\[.\]$','');
elseif length(unique(mapModel.mets)) == length(mapModel.mets)
    % If the met IDs are not all unique, assume the compartment has already
    % been removed; otherwise, assume that the last character of the ID is
    % the compartment abbreviation, and remove it.
    map_mets = regexprep(mapModel.mets,'.$','');
end

% obtain unique list of compartment-free met IDs, and only keep field
% entries associated with this unique set of metabolites
[~,uniq_met_ind] = unique(map_mets);
map_mets = map_mets(uniq_met_ind);
mapModel.metNamesAlt = mapModel.metNamesAlt(uniq_met_ind,:);
for i = 1:length(metFields)
    % also remove rows correponding to non-unique mets from metFields
    mapModel.(metFields{i}) = mapModel.(metFields{i})(uniq_met_ind,:);
end

% if "metNames" is one of the metFields, convert it to three entries
% for the three stages of metabolite name comparisons
[~,ind] = ismember('metNames',metFields);
if ind > 0
    metFields = [metFields(1:ind-1);{'metNames (stage 1)';'metNames (stage 2)';'metNames (stage 3)'};metFields(ind+1:end)];
end

% map metabolites along each field in metFields
metAltModelID = {};  % initialize
for f = 1:length(metFields)
    
    fprintf('Mapping metabolites via %s... ',metFields{f});
    ids = repmat({''},size(refModel.mets));  % initialize
    
    if  prioritize && ~isempty(metAltModelID)
        % only map metabolites that have not yet been mapped
        unmapped = all(cellfun(@isempty,metAltModelID),2);
    else
        % map all metabolites, regardless if they have already been mapped
        % previously via some other field
        unmapped = true(length(refModel.mets),1);
    end
    
    % extract IDs from model field
    switch metFields{f}
        
        case 'mets'
            
            % If comparing the "mets" field, remove the compartment
            is_name = false;
            if strcmp(refModel.mets{1}(end),']')
                ref_mets = regexprep(refModel.mets,'\[.\]$','');
            elseif length(unique(refModel.mets)) == length(refModel.mets)
                ref_mets = regexprep(refModel.mets,'.$','');
            end
            ref_ids = lower(ref_mets);
            map_ids = lower(map_mets);
            
        case 'metNames (stage 1)'
            
            % METNAMES STAGE 1: Map mets via metNames field
            is_name = true;
            ref_ids = refModel.metNames;
            map_ids = [mapModel.metNames,mapModel.metNamesAlt];
            
            % initialize variable to keep track of mapping via met names
            unmapped_viaNames = true(size(unmapped));
            
        case 'metNames (stage 2)'
            
            % METNAMES STAGE 2: Map remaining mets via metNamesAlt field
            is_name = true;
            ref_ids = [refModel.metNames,refModel.metNamesAlt];
            map_ids = [mapModel.metNames,mapModel.metNamesAlt];
            
            % do not try to map metabolites that have already been mapped
            % during the first metNames mapping stage
            unmapped = unmapped & unmapped_viaNames;
            
        case 'metNames (stage 3)'
            
            % METNAMES STAGE 3: Make manual changes to metNames, and map.
            % These changes are specific to the mapping between HMR2 and 
            % Recon3D, and are based on observations in naming differences.
            is_name = true;
            ref_ids = [refModel.metNames,refModel.metNamesAlt];
            map_ids = [mapModel.metNames,mapModel.metNamesAlt];
            
            % make name adjustments to Recon3D model if present
            if isfield(refModel,'modelID') && strcmpi(refModel.modelID,'Recon3D')
                ref_ids = applyMetNameChanges(ref_ids);
            elseif isfield(mapModel,'modelID') && strcmpi(mapModel.modelID,'Recon3D')
                map_ids = applyMetNameChanges(map_ids);
            else
                % neither model is recognized as Recon3D
                fprintf('No Recon3D recognized, skipped.\n');
                continue
            end
            
            % do not try to map metabolites that have already been mapped
            % during the first or second metNames mapping stages
            unmapped = unmapped & unmapped_viaNames;
            
        otherwise
            
            is_name = false;
            ref_ids = refModel.(metFields{f});
            map_ids = mapModel.(metFields{f});
            
    end
    
    % if map_ids contains multiple columns, flatten into a single
    % column, and repeat entries of mapModel mets to maintain alignment
    non_empty = ~cellfun(@isempty,map_ids);
    if size(map_ids,2) > 1
        repMets = arrayfun(@(i) repmat(map_mets(i),sum(non_empty(i,:),2),1),[1:length(map_mets)]','UniformOutput',false);
        repMets = vertcat(repMets{:});
        map_ids = map_ids';
        map_ids = map_ids(non_empty');
    else
        map_ids = map_ids(non_empty);
        repMets = map_mets(non_empty);
    end
    
    % assemble mapping array
    met2id = [repMets,map_ids];
    
    % map metabolites
    ids(unmapped) = matchIDs(ref_ids(unmapped,:),met2id,is_name);
    
    % to keep track of name matching process
    if startsWith(metFields{f},'metNames (stage')
        unmapped_viaNames(~cellfun(@isempty,ids)) = false;
    end
    
    % flatten and append newly mapped IDs
    ids = flattenCell(ids,true);
    metAltModelID = [metAltModelID,ids];
    
    fprintf('Done.\n');
end

% remove duplicate mapped IDs
empty_inds = cellfun(@isempty,metAltModelID);
metAltModelID = arrayfun(@(i) unique(metAltModelID(i,~empty_inds(i,:))),[1:length(refModel.mets)]','UniformOutput',false);
metAltModelID = flattenCell(metAltModelID,true);

% add mapped IDs to model structure
model.metAltModelID = metAltModelID;

end


function chgNames = applyMetNameChanges(names)

% met name changes specific to Recon3D-HMR2 metabolite mapping
chgNames = regexprep(names,'Coenzyme A','CoA');
chgNames = regexprep(chgNames,'\(.*-Density Lipoprotein\)','');
chgNames = regexprep(chgNames,'13-Cis-Retinoyl Glucuronide','13-Cis-Retinoyl-beta-D-Glucuronide');
chgNames = regexprep(chgNames,"Cytidine-5'-Diphosphate",'CDP');
chgNames = regexprep(chgNames,'Human Liver Homolog','');
chgNames = regexprep(chgNames,'Glutathionyl-Leuc4','glutathionyl-leukotriene C4');
chgNames = regexprep(chgNames,'Omega-Cooh-Tetranor-Leukotriene E3','Omega-Cooh-Tetranor-LTE3');

end


function metID = matchIDs(ids,met2id,is_name)

if is_name
    % make all met names lowercase, and remove special characters (non-word or
    % digit characters, e.g., dashes, parentheses, spaces, etc.)
    met2id(:,2) = lower(regexprep(met2id(:,2),'[^a-zA-Z0-9]',''));
    ids = lower(regexprep(ids,'[^a-zA-Z0-9]',''));
end

% compress refModel met IDs into single column of nested entries
if size(ids,2) > 1
    ids = nestCell(ids,true);
end

% find empty met ID indices to ignore
ignore_ind = cellfun(@isempty,ids);

% map refModel met IDs to mapModel met IDs
metID = repmat({''},size(ids,1),1);
metID(~ignore_ind) = cellfun(@(x) unique(met2id(ismember(met2id(:,2),x),1)),ids(~ignore_ind),'UniformOutput',false);

end



