function annModel = annotateModel(model,annType,addMiriams,addFields,overwrite)
% Add reaction, metabolite, and/or gene annotation to a model.
%
% Input:
%
%   model        Model structure.
%
%   annType      String or cell array of strings specifying the type(s) of 
%                annotation data to add: 'rxn', 'met', and/or 'gene'. To
%                add all annotation types, use 'all'.
%                For example, annType={'rxn','gene'} will add reaction and
%                gene annotation information, but not metabolite.
%                (Opt, default 'all')
%
%   addMiriams   If TRUE, annotation information will be added to the 
%                metMiriams field. If the metMiriams field does not exist,
%                it will be created.
%                (Opt, default TRUE)   
%
%   addFields    If TRUE, separate model fields will be added for each 
%                added ID type, e.g., reaction KEGG IDs will be added to a
%                rxnKEGGID field.
%                (Opt, default TRUE)
%
%   overwrite    If TRUE, any existing metMiriam or model ID field entries
%                will be overwritten. If FALSE, they will be appended with
%                the new data.
%                (Opt, default TRUE)
%
%
% Output:
%
%   annModel     Model structure with added annotation information.
%
%
% Usage:
%
%   annModel = annotateMets(model,annType,addMiriams,addFields);
%
%
% Jonathan Robinson, 2019-07-18
%


%% Inputs and setup

if nargin < 2 || isempty(annType) || strcmpi(annType,'all')
    annType = {'rxn','met','gene'};
elseif ~all(ismember(annType,{'rxn','met','gene','reaction','metabolite'}))
    error('annType input(s) not recognized. Valid options are "rxn", "met", and/or "gene", or "all"');
end

if nargin < 3 || isempty(addMiriams)
    addMiriams = true;
end

if nargin < 4 || isempty(addFields)
    addFields = true;
end

if nargin < 5
    overwrite = true;
end

if (~addMiriams) && (~addFields)
    error('addMiriams and addFields are both FALSE - no annotation will be added.')
end

% map field names to those on identifiers.org (used for building Miriams)
id2miriam = {%reactions
             'rxnKEGGID'        'kegg.reaction'
             'rxnBiGGID'        'bigg.reaction'
             'rxnREACTOMEID'    'reactome'
             'rxnMetaNetXID'    'metanetx.reaction'
             % metabolites
             'metBiGGID'        'bigg.metabolite'
             'metKEGGID'        'kegg.compound'
             'metHMDBID'        'hmdb'
             'metChEBIID'       'chebi'
             'metPubChemID'     'pubchem.compound'
             'metLipidMapsID'   'lipidmaps'
             'metMetaNetXID'    'metanetx.chemical'
             % genes
             'geneNames'        'hgnc.symbol'
             'geneEnsemblID'    'ensembl'
             'geneEntrezID'     'ncbigene'
             'geneUniProtID'    'uniprot'};


%% Load and organize annotation data

% load reaction annotation data
if any(ismember({'rxn','reaction'},lower(annType)))
    [ST, I] = dbstack('-completenames');
    path = fileparts(ST(I).file);
    tmpfile = fullfile(path,'../../ComplementaryData/annotation','humanGEMRxnAssoc.JSON');
    rxnAssoc = jsondecode(fileread(tmpfile));
    
    rxnAssocArray = struct2cell(rxnAssoc);
    rxnAssocArray = horzcat(rxnAssocArray{:});
else
    rxnAssoc = [];
end

% load metabolite annotation data
if any(ismember({'met','metabolite'},lower(annType)))
    [ST, I] = dbstack('-completenames');
    path = fileparts(ST(I).file);
    tmpfile = fullfile(path,'../../ComplementaryData/annotation','humanGEMMetAssoc.JSON');
    metAssoc = jsondecode(fileread(tmpfile));
    
    metAssocArray = struct2cell(metAssoc);
    metAssocArray = horzcat(metAssocArray{:});
else
    metAssoc = [];
end

% load and organize gene annotation data
if ismember('gene',lower(annType))
    % there is not a .JSON annotation file for genes because it is
    % unnecessary - we will instead use the grRule translation function,
    % which uses gene annotations retrieved from the Ensembl database.
    geneAssoc.genes = model.genes;
    geneAssoc.geneNames = regexprep(translateGrRules(geneAssoc.genes,'Name'),' or ','; ');  % HGNC Symbols
    geneAssoc.geneEnsemblID = model.genes;
    geneAssoc.geneEntrezID = regexprep(translateGrRules(geneAssoc.genes,'Entrez'),' or ','; ');
    geneAssoc.geneUniProtID = regexprep(translateGrRules(geneAssoc.genes,'UniProt'),' or ','; ');
    
    geneAssocArray = struct2cell(geneAssoc);
    geneAssocArray = horzcat(geneAssocArray{:});
else
    geneAssoc = [];
end


%% Add annotation data to model MIRIAM fields

if ( addMiriams )
    
    % Reactions
    if ~isempty(rxnAssoc)
        % map model rxns to those in the association structures
        [hasRxn,rxnInd] = ismember(model.rxns,rxnAssoc.rxns);
        
        % add rxnMiriams field if it does not exist or if overwrite is TRUE
        if ~isfield(model,'rxnMiriams') || overwrite
            model.rxnMiriams = repmat({''},size(model.rxns));
        end
        
        % map rxn ID fields to the proper identifiers.org terms
        [hasRxnField,rxnFieldInd] = ismember(fieldnames(rxnAssoc),id2miriam(:,1));
        miriamNames = id2miriam(rxnFieldInd(hasRxnField),2);
        miriamValues = rxnAssocArray(:,hasRxnField);
        
        % add new data to the rxnMiriams field
        for i = 1:numel(model.rxns)
            if hasRxn(i)
                model.rxnMiriams{i} = appendMiriamData(model.rxnMiriams{i}, miriamNames, miriamValues(rxnInd(i),:)');
            end
        end
    end
    
    % Metabolites
    if ~isempty(metAssoc)
        % map model mets to those in the association structures
        [hasMet,metInd] = ismember(model.mets,metAssoc.mets);
        
        % add metMiriams field if it does not exist or if overwrite is TRUE
        if ~isfield(model,'metMiriams') || overwrite
            model.metMiriams = repmat({''},size(model.mets));
        end
        
        % map met ID fields to the proper identifiers.org terms
        [hasMetField,metFieldInd] = ismember(fieldnames(metAssoc),id2miriam(:,1));
        miriamNames = id2miriam(metFieldInd(hasMetField),2);
        miriamValues = metAssocArray(:,hasMetField);
        
        % add new data to the metMiriams field
        for i = 1:numel(model.mets)
            if hasMet(i)
                model.metMiriams{i} = appendMiriamData(model.metMiriams{i}, miriamNames, miriamValues(metInd(i),:)');
            end
        end
    end
    
    % Genes
    if ~isempty(geneAssoc)
        % Note: because geneAssoc was built using geneIDs directly from the
        %       model, we do not need to map model.genes to the geneAssoc
        %       structure - the indexing should be identical.
        
        % add geneMiriams field if it does not exist or if overwrite is TRUE
        if ~isfield(model,'geneMiriams') || overwrite
            model.geneMiriams = repmat({''},size(model.genes));
        end
        
        % map gene ID fields to the proper identifiers.org terms
        [hasGeneField,geneFieldInd] = ismember(fieldnames(geneAssoc),id2miriam(:,1));
        miriamNames = id2miriam(geneFieldInd(hasGeneField),2);
        miriamValues = geneAssocArray(:,hasGeneField);
        
        % add new data to the geneMiriams field
        for i = 1:numel(model.genes)
            model.geneMiriams{i} = appendMiriamData(model.geneMiriams{i}, miriamNames, miriamValues(i,:)');
        end
    end
    
end


%% Add annotation data to individual model fields

if ( addFields )
    
    % merge association structures
    allAssoc = mergeStructures({rxnAssoc,metAssoc,geneAssoc});
    
    % some fields should not be added to the model
    remFields = {'rxns','mets','metsNoComp','genes'};
    allAssoc = rmfield(allAssoc, intersect(fieldnames(allAssoc),remFields));
    
    % add individual ID fields to the model
    f = fieldnames(allAssoc);
    for i = 1:numel(f)
        
        if startsWith(f{i},'rxn')
            [hasRxn,rxnInd] = ismember(model.rxns,rxnAssoc.rxns);
            ids = repmat({''},size(model.rxns));
            ids(hasRxn) = allAssoc.(f{i})(rxnInd(hasRxn));
        elseif startsWith(f{i},'met')
            [hasMet,metInd] = ismember(model.mets,metAssoc.mets);
            ids = repmat({''},size(model.mets));
            ids(hasMet) = allAssoc.(f{i})(metInd(hasMet));
        elseif startsWith(f{i},'gene')
            ids = allAssoc.(f{i});
        else
            continue
        end
        
        if ~isfield(model,f{i}) || overwrite
            model.(f{i}) = ids;
        else
            % merge new IDs with existing IDs
            merged_ids = join([model.(f{i}), ids],'; ');
            merged_ids = arrayfun(@(id) strjoin(unique(split(id,'; ')),'; '), merged_ids, 'UniformOutput', false);
            model.(f{i}) = regexprep(merged_ids,'^;\s*','');  % remove any preceding semicolons
        end
        
    end
    
end

%%

% assign output
annModel = model;

end



%% Additional functions

function miriam_new = appendMiriamData(miriam,names,values)
% Append new entries to an existing miriam structure.

% ignore empty entries
remove = cellfun(@isempty,values);
names(remove) = [];
values(remove) = [];
if isempty(names)
    miriam_new = miriam;
    return
end

% convert existing miriam entry to cell array
if ~isempty(miriam)
    miriam = [miriam.name(:), miriam.value(:)];
end

% append miriam entry with new data
for i = 1:numel(names)
    if contains(values{i},';')
        % handle multiple values separated by semicolon (;)
        val = strtrim(split(values{i},';'));
        name = repmat(names(i),numel(val),1);
        miriam = [miriam; [name, val]];
    else
        miriam = [miriam; [names(i), values(i)]];
    end
end

% remove any repeated entries
[~,miriam_num] = ismember(miriam,miriam);
[~,uniq_ind] = unique(miriam_num,'rows');
miriam = miriam(uniq_ind,:);

% convert back to structure
miriam_new.name = miriam(:,1);
miriam_new.value = miriam(:,2);

end


function mergedStruct = mergeStructures(structs)
% Merge multiple structures into a single structure.

mergedStruct = {};
for i = 1:numel(structs)
    
    s = structs{i};
    
    if isempty(s)
        continue
    end
    
    f = fieldnames(s);
    for j = 1:numel(f)
        mergedStruct.(f{j}) = s.(f{j});
    end
    
end

end



