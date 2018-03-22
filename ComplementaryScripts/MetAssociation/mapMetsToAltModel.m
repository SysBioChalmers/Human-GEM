function model = mapMetsToAltModel(refModel,mapModel)
%mapMetsToAltModel  Map metabolites from one model to another via names.
%
% USAGE:
%
%   model = mapMetsToAltModel(refModel,mapModel);
%
%
% INPUTS:
%
%   refModel  Model structure containing metabolites that are to be mapped
%             to mapModel metabolite identifiers.
%
%   mapModel  Model structure containing metabolites to which refModel will
%             be mapped.
%
%
% OUTPUTS:
%
%   model     Model structure with added field "metAltModelID", which 
%             contains the mapModel metIDs that were matched to the 
%             refModel metabolites.
%
%
% Jonathan Robinson 2018-03-21


% initialize output
model = refModel;

% strip compartments from mapModel met IDs to obtain compartment-free
% met-to-metName pairs
if strcmp(mapModel.mets{1}(end),']')
    % compartment information is formatted in brackets at end of ID
    mets = regexprep(mapModel.mets,'\[.\]$','');
elseif length(unique(mapModel.mets)) == length(mapModel.mets)
    % If the met IDs are not all unique, assume the compartment has already
    % been removed; otherwise, assume that the last character of the ID is
    % the compartment abbreviation, and remove it.
    mets = regexprep(mapModel.mets,'.$','');
end

% obtain unique list of compartment-free met IDs and met names
[~,ind] = unique(mets);
mets = mets(ind);
if isfield(mapModel,'metNamesAlt')
    metNames = [mapModel.metNames(ind,:),mapModel.metNamesAlt(ind,:)];
else
    metNames = mapModel.metNames(ind,:);
end

% met2Name = [mapModel.mets(ind),mapModel.metNames(ind,:)];


if size(metNames,2) > 1
    % if metNames contains multiple columns, flatten into a single column
    % of names, and repeat entries of mets to keep alignment
    non_empty = ~cellfun(@isempty,metNames);
    mets = arrayfun(@(i) repmat(mets(i),sum(non_empty(i,:),2),1),[1:length(mets)]','UniformOutput',false);
    mets = vertcat(mets{:});
    metNames = metNames';
    metNames = metNames(non_empty');
elseif any(contains(mapModel.metNames,'; '))
    % if metNames is a single column, search for semicolons and split names
    % by semicolons, and repeat entries of mets to keep alignment
    metNames = cellfun(@(s) strsplit(s,'; ')',metNames,'UniformOutput',false);
    mets = arrayfun(@(i) repmat(mets(i),numel(metNames{i}),1),[1:length(mets)]','UniformOutput',false);
    mets = vertcat(mets{:});
    metNames = vertcat(metNames{:});
end
met2name = [mets,metNames];


% ***** STAGE 1: Map mets via metNames field *****
metAltModelID = matchIDs(refModel.metNames,met2name);
unmapped = cellfun(@isempty,metAltModelID);  % find unmapped metabolites


% ***** STAGE 2: Map remaining mets via metNamesAlt field *****
if isfield(refModel,'metNamesAlt')
    refModel.metNames = [refModel.metNames,refModel.metNamesAlt];
    metAltModelID(unmapped) = matchIDs(refModel.metNames(unmapped,:),met2name);
    unmapped = cellfun(@isempty,metAltModelID);  % find unmapped metabolites
end


% ***** STAGE 3: Make manual changes to metNames, and repeat *****
% These changes are specific to the mapping between HMR2 and Recon3D, and
% are based on observations in naming differences.
if isfield(refModel,'modelID') && strcmpi(refModel.modelID,'Recon3D')
    refModel.metNames = applyMetNameChanges(refModel.metNames);
elseif isfield(mapModel,'modelID') && strcmpi(mapModel.modelID,'Recon3D')
    met2name(:,2) = applyMetNameChanges(met2name(:,2));
end
metAltModelID(unmapped) = matchIDs(refModel.metNames(unmapped,:),met2name);


% flatten cell column into 2D cell array
metAltModelID = flattenCell(metAltModelID,true);

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


function metID = matchIDs(names,met2name)

% make all met names lowercase, and remove special characters (non-word or 
% digit characters, e.g., dashes, parentheses, spaces, etc.)
met2name(:,2) = lower(regexprep(met2name(:,2),'[^a-zA-Z0-9]',''));
names = lower(regexprep(names,'[^a-zA-Z0-9]',''));

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



