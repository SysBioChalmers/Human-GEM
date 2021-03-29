function model = mapMetsToAltModel(refModel,mapModel,metFields,mapMethod)
%mapMetsToAltModel  Map metabolites from one model to another.
%
% USAGE:
%
%   model = mapMetsToAltModel(refModel,mapModel,metFields,mapMethod);
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
%   mapMethod  (Optional) Specify the method used to handle metabolites in
%              refModel that map to multiple IDs in mapModel across 
%              multiple metFields.
%                 'all'  (Default) For each metabolite, met IDs will be
%                        returned for all mets in mapModel that were mapped
%                        along any of the specified metFields.
%               'order'  The order in which metFields are listed will
%                        determine their priority in mapping met IDs (where
%                        first is most important/confident, and last is the
%                        least). Once a met is mapped via one of the fields
%                        listed in metFields, there will be no further
%                        attempts to map that metabolite via all subsequent
%                        fields in metFields.
%               'score'  Each met will be mapped along all fields listed in
%                        metFields. The mapped met ID(s) with the highest 
%                        score for each metabolite will be kept, whereas
%                        all other mapped IDs with lower scores will be
%                        removed as potential matches. The score is
%                        calculated based on the number of fields in
%                        metFields by which the refModel met is mapped to
%                        the mapModel met. If metFields includes a second
%                        column of field weights, then the score will be
%                        weighted by these values, where mets mapped along
%                        fields with greater weights will receive higher
%                        scores than those mapped along fields with lower
%                        weights.
%
% OUTPUTS:
%
%   model     Model structure with added field "metAltModelID", which 
%             contains the mapModel metIDs that were matched to the 
%             refModel metabolites.
%


% handle input arguments
if nargin < 4
    mapMethod = 'all';
end

if ischar(metFields)
    % convert char to cell
    metFields = {metFields};
    metFieldWeights = 1;  % weight is irrelevant if only one metField is provided
elseif isvector(metFields)
    % ensure that metFields is a column vector (probably unnecessary)
    metFields = metFields(:);
    metFieldWeights = ones(size(metFields));  % weight all metFields equally
else
    % if metFields contains a second column, it is treated as the vector of
    % weights corresponding to the list of met fields, and is extracted
    metFieldWeights = cell2mat(metFields(:,2));
    metFields = metFields(:,1);
end

if strcmpi(mapMethod,'order')
    % if mapMethod is "order", then weight fields based on their order of appearance in metFields
    metFieldWeights = (length(metFields):-1:1)';
end

% verify that all metFields are present in both models
if any(~isfield(refModel,metFields)) || any(~isfield(mapModel,metFields))
    error('One or more of the specified metFields is not present in one or both models.');
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

% initialize variables
matchScores = [];
name_stage = 1;
f = 1;

% map metabolites along each field in metFields
while f <= length(metFields)
    
    fprintf('Mapping metabolites via %s... ',metFields{f});
    
    % extract IDs from model field
    switch metFields{f}
        
        case 'mets'
            
            % If comparing the "mets" field, remove the compartment
            if strcmp(refModel.mets{1}(end),']')
                ref_mets = regexprep(refModel.mets,'\[.\]$','');
            elseif length(unique(refModel.mets)) == length(refModel.mets)
                ref_mets = regexprep(refModel.mets,'.$','');
            end
            ref_ids = lower(ref_mets);
            map_ids = lower(map_mets);
            is_name = false;
            ids = repmat({''},size(refModel.mets));  % initialize matches
            
        case 'metNames'
            
            if name_stage == 1
                
                % METNAMES STAGE 1: Map mets via metNames field
                fprintf('(stage 1) ');
                is_name = true;
                ref_ids = refModel.metNames;
                map_ids = [mapModel.metNames,mapModel.metNamesAlt];
                ids = repmat({''},size(refModel.mets));  % initialize matches
                
            elseif name_stage == 2
                
                % METNAMES STAGE 2: Map remaining mets via metNamesAlt field
                fprintf('(stage 2) ');
                is_name = true;
                ref_ids = [refModel.metNames,refModel.metNamesAlt];
                map_ids = [mapModel.metNames,mapModel.metNamesAlt];
                
            elseif name_stage == 3
                
                % METNAMES STAGE 3: Make manual changes to metNames, and map.
                % These changes are specific to the mapping between HMR2 and
                % Recon3D, and are based on observations in naming differences.
                fprintf('(stage 3) ');
                is_name = true;
                ref_ids = [refModel.metNames,refModel.metNamesAlt];
                map_ids = [mapModel.metNames,mapModel.metNamesAlt];
                
                % make name adjustments to Recon3D model if present
                if isfield(refModel,'modelID') && strcmpi(refModel.modelID,'Recon3D')
                    ref_ids = applyMetNameChanges(ref_ids);
                elseif isfield(mapModel,'modelID') && strcmpi(mapModel.modelID,'Recon3D')
                    map_ids = applyMetNameChanges(map_ids);
                end
                
            end
            
        otherwise
            
            % comparing by field other than "mets" or "metNames"
            is_name = false;  % not a met name
            ref_ids = refModel.(metFields{f});
            map_ids = mapModel.(metFields{f});
            ids = repmat({''},size(refModel.mets));  % initialize matches
            
    end
    
    % if map_ids contains multiple columns, flatten into a single
    % column, and repeat entries of mapModel mets to maintain alignment
    non_empty = ~cellfun(@isempty,map_ids);
    if size(map_ids,2) > 1
        repMets = arrayfun(@(i) repmat(map_mets(i),sum(non_empty(i,:),2),1),(1:length(map_mets))','UniformOutput',false);
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
    ignore_mets = ~cellfun(@isempty,ids);  % for the name-matching process, ignore mets that have been mapped in previous name-matching stages
    ids(~ignore_mets) = matchIDs(ref_ids(~ignore_mets,:),met2id,is_name);
    fprintf('Done.\n');
    
    if strcmpi(metFields{f},'metNames') && name_stage < 3
        % if processing metNames, move to next stage
        name_stage = name_stage + 1;
    else
        
        % get indices of mapped mets in each of the respective models
        ref_met_inds = arrayfun(@(i) repmat(i,numel(ids{i}),1),(1:length(ids))','UniformOutput',false);
        ref_met_inds = vertcat(ref_met_inds{:});
        [~,map_met_inds] = ismember(vertcat(ids{:}),map_mets);
        
        % append matches to match score matrix
        matchScores(:,:,f) = accumarray(unique([ref_met_inds,map_met_inds],'rows'),metFieldWeights(f),[length(refModel.mets),length(map_mets)]);
        f = f + 1;  % proceed to next metField
        
    end
    
end


if strcmpi(mapMethod,'order')
    % only use the match from the top-scoring (earliest-listed) metField
    matchScores(matchScores < max(matchScores,[],3)) = 0;
end

% determine the total score for each match by adding scores among all fields
matchScores = sum(matchScores,3);

if ~strcmpi(mapMethod,'all')
    % For each metabolite, only keep the mapped IDs with the highest score. 
    % If there is a tie, keep all IDs that are tied.
    matchScores(matchScores < max(matchScores,[],2)) = 0;
end

% retrieve mapModel met IDs and add to output model structure
metAltModelID = arrayfun(@(i) map_mets(matchScores(i,:) > 0),(1:length(refModel.mets))','UniformOutput',false);
model.metAltModelID = flattenCell(metAltModelID,true);

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



