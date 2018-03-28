function model = mapMetsToAltModel(refModel,mapModel,altFields)
%mapMetsToAltModel  Map metabolites from one model to another via names.
%
% USAGE:
%
%   model = mapMetsToAltModel(refModel,mapModel,altFields);
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
%   altFields  (Optional) A list of one or more metabolite-related field(s)
%              that, in addition to metabolite names, will be compared 
%              between the two models to try and map metabolites from one
%              model to the other.
%
% OUTPUTS:
%
%   model     Model structure with added field "metAltModelID", which 
%             contains the mapModel metIDs that were matched to the 
%             refModel metabolites.
%
%
% Jonathan Robinson 2018-03-22


% handle input arguments
if nargin < 3
    altFields = [];
elseif any(~isfield(refModel,altFields)) || any(~isfield(mapModel,altFields))
    error('One or more of the specified altFields is not present in one or both models.');
elseif ischar(altFields)
    altFields = {altFields};
end


% initialize output
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

% do the same for refModel mets
if strcmp(refModel.mets{1}(end),']')
    % compartment information is formatted in brackets at end of ID
    ref_mets = regexprep(refModel.mets,'\[.\]$','');
elseif length(unique(refModel.mets)) == length(refModel.mets)
    % If the met IDs are not all unique, assume the compartment has already
    % been removed; otherwise, assume that the last character of the ID is
    % the compartment abbreviation, and remove it.
    ref_mets = regexprep(refModel.mets,'.$','');
end

% obtain unique list of compartment-free met IDs and met names
[~,uniq_met_ind] = unique(map_mets);
map_mets = map_mets(uniq_met_ind);
if isfield(mapModel,'metNamesAlt')
    metNames = [mapModel.metNames(uniq_met_ind,:),mapModel.metNamesAlt(uniq_met_ind,:)];
else
    metNames = mapModel.metNames(uniq_met_ind,:);
end
for i = 1:length(altFields)
    % also remove rows correponding to non-unique mets from altFields, if specified
    mapModel.(altFields{i}) = mapModel.(altFields{i})(uniq_met_ind,:);
end

if size(metNames,2) > 1
    % if metNames contains multiple columns, flatten into a single column
    % of names, and repeat entries of mets to keep alignment
    non_empty = ~cellfun(@isempty,metNames);
    repMets = arrayfun(@(i) repmat(map_mets(i),sum(non_empty(i,:),2),1),[1:length(map_mets)]','UniformOutput',false);
    repMets = vertcat(repMets{:});
    metNames = metNames';
    metNames = metNames(non_empty');
elseif any(contains(mapModel.metNames,'; '))
    % if metNames is a single column, search for semicolons and split names
    % by semicolons, and repeat entries of mets to keep alignment
    metNames = cellfun(@(s) strsplit(s,'; ')',metNames,'UniformOutput',false);
    repMets = arrayfun(@(i) repmat(map_mets(i),numel(metNames{i}),1),[1:length(map_mets)]','UniformOutput',false);
    repMets = vertcat(repMets{:});
    metNames = vertcat(metNames{:});
end
met2name = [repMets,metNames];


% ***** STAGE 1: Map mets via metNames field *****
fprintf('Mapping metabolites via metNames... ');
metAltModelID = matchIDs(refModel.metNames,met2name);
unmapped = cellfun(@isempty,metAltModelID);  % find unmapped metabolites
fprintf('Done.\n');


% ***** STAGE 2: Map remaining mets via metNamesAlt field *****
if isfield(refModel,'metNamesAlt')
    fprintf('Mapping metabolites via metNamesAlt... ');
    refModel.metNames = [refModel.metNames,refModel.metNamesAlt];
    metAltModelID(unmapped) = matchIDs(refModel.metNames(unmapped,:),met2name);
    unmapped = cellfun(@isempty,metAltModelID);  % find unmapped metabolites
    fprintf('Done.\n');
end


% ***** STAGE 3: Make manual changes to metNames, and repeat *****
% These changes are specific to the mapping between HMR2 and Recon3D, and
% are based on observations in naming differences.
if isfield(refModel,'modelID') && strcmpi(refModel.modelID,'Recon3D')
    fprintf('Mapping metabolites via metNames after name adjustments... ');
    refModel.metNames = applyMetNameChanges(refModel.metNames);
    metAltModelID(unmapped) = matchIDs(refModel.metNames(unmapped,:),met2name);
    fprintf('Done.\n');
elseif isfield(mapModel,'modelID') && strcmpi(mapModel.modelID,'Recon3D')
    fprintf('Mapping metabolites via metNames after name adjustments... ');
    met2name(:,2) = applyMetNameChanges(met2name(:,2));
    metAltModelID(unmapped) = matchIDs(refModel.metNames(unmapped,:),met2name);
    fprintf('Done.\n');
end

% flatten cell column into 2D cell array
metAltModelID = flattenCell(metAltModelID,true);


% ***** OPTIONAL STAGE: Map metabolites through other fields *****
% Map metabolites using additional model fields, if provided.
if ~isempty(altFields)
    newMappedIDs = repmat({''},size(refModel.mets));
    for f = 1:length(altFields)
        fprintf('Mapping metabolites via %s... ',altFields{f});
        
        % extract IDs from model field
        if strcmp(altFields{f},'mets')
            % If comparing the "mets" field, use the compartment-free
            % format that was generated earlier. Also ignore case.
            ref_ids = lower(ref_mets);
            map_ids = lower(map_mets);
        else
            ref_ids = refModel.(altFields{f});
            map_ids = mapModel.(altFields{f});
        end
        
        % if the field contains multiple columns, flatten into a single 
        % column, and repeat entries of mets to maintain alignment
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
        ids = matchIDs(ref_ids,met2id,true);
        
        % flatten and append newly mapped IDs
        ids = flattenCell(ids,true);
        newMappedIDs = [newMappedIDs,ids];
        
        fprintf('Done.\n');
    end
    
    % append newly mapped IDs to name-matched IDs, and remove duplicates
    metAltModelID = [metAltModelID,newMappedIDs];
    empty_inds = cellfun(@isempty,metAltModelID);
    metAltModelID = arrayfun(@(i) unique(metAltModelID(i,~empty_inds(i,:))),[1:length(refModel.mets)]','UniformOutput',false);
    metAltModelID = flattenCell(metAltModelID,true);
    
end


% assign output
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


function metID = matchIDs(names,met2name,is_id)

if nargin < 3
    is_id = false;
end

if ~is_id
    % make all met names lowercase, and remove special characters (non-word or
    % digit characters, e.g., dashes, parentheses, spaces, etc.)
    met2name(:,2) = lower(regexprep(met2name(:,2),'[^a-zA-Z0-9]',''));
    names = lower(regexprep(names,'[^a-zA-Z0-9]',''));
end

% compress refModel met names into single column of nested entries
if size(names,2) > 1
    names = nestCell(names,true);
end

% find empty met name indices to ignore
ignore_ind = cellfun(@isempty,names);

% map refModel met names to mapModel met names
metID = repmat({''},size(names,1),1);
metID(~ignore_ind) = cellfun(@(x) unique(met2name(ismember(met2name(:,2),x),1)),names(~ignore_ind),'UniformOutput',false);

end



