function writeYAMLmodel(model,fileName,preserveQuotes,sortIds)
% writeYAMLmodel
%   Writes a yaml file matching (roughly) the cobrapy yaml structure
%
%   model           a model structure
%   fileName        name that the file will have.  A dialog window will 
%                   open if no file name is specified.
%   preserveQuotes  if quotes should be preserved for strings
%                   (logical, default=true)
%   sortIds         if metabolites, reactions, genes and compartments
%                   should be sorted alphabetically by their identifier,
%                   otherwise they are kept in their original order
%                   (logical, default=false)   
%
% Usage: writeYAMLmodel(model,fileName,preserveQuotes,sortIds)
if nargin<2|| isempty(fileName)
    [fileName, pathName] = uiputfile({'*.yml;*.yaml'}, 'Select file for model export',[model.id '.yml']);
    if fileName == 0
        error('You should provide a file location')
    else
        fileName = fullfile(pathName,fileName);
    end
end
fileName=char(fileName);

if nargin < 3
    preserveQuotes = true;
end
if nargin < 4
    sortIds = false;
end
if ~endsWith(fileName,{'.yml','.yaml'})
    fileName = strcat(fileName,'.yml');
end

%Check that model is in RAVEN format:
if isfield(model,'rules')
    model = ravenCobraWrapper(model);
end

%Sort identifiers alphabetically
if sortIds == true
    model = sortIdentifiers(model);
end

%Simplify Miriam fields:
if isfield(model,'metMiriams')
    [model.newMetMiriams,model.newMetMiriamNames]   = extractMiriam(model.metMiriams);
end
if isfield(model,'rxnMiriams')
    [model.newRxnMiriams,model.newRxnMiriamNames]   = extractMiriam(model.rxnMiriams);
end
if isfield(model,'geneMiriams')
    [model.newGeneMiriams,model.newGeneMiriamNames] = extractMiriam(model.geneMiriams);
end
if isfield(model,'compMiriams')
    [model.newCompMiriams,model.newCompMiriamNames] = extractMiriam(model.compMiriams);
end

%Open file:
fid = fopen(fileName,'wt');
if fid == -1
    error(['Cannot write to ' fileName ', does the directory exist?'])
end
fprintf(fid,'---\n!!omap\n');

%Insert file header (metadata)
writeMetadata(model,fid);

%Metabolites:
fprintf(fid,'- metabolites:\n');
for i = 1:length(model.mets)
    fprintf(fid,'    - !!omap\n');
    writeField(model, fid, 'mets',        'txt', i, '  - id',          preserveQuotes)
    writeField(model, fid, 'metNames',    'txt', i, '  - name',        preserveQuotes)
    writeField(model, fid, 'metComps',    'txt', i, '  - compartment', preserveQuotes)
    writeField(model, fid, 'metFormulas', 'txt', i, '  - formula',     preserveQuotes)
    writeField(model, fid, 'metCharges',  'num', i, '  - charge',      preserveQuotes)
    writeField(model, fid, 'inchis',      'txt', i, '  - inchis',      preserveQuotes)
    writeField(model, fid, 'metSmiles',   'txt', i, '  - smiles',      preserveQuotes)
    writeField(model, fid, 'metMiriams',  'txt', i, '  - annotation',  preserveQuotes)
    writeField(model, fid, 'metDeltaG',   'num', i, '  - deltaG',      preserveQuotes)
    writeField(model, fid, 'metNotes',    'txt', i, '  - notes',       preserveQuotes)
    writeField(model, fid, 'metFrom',     'txt', i, '  - metFrom',     preserveQuotes)
end

%Reactions:
fprintf(fid,'- reactions:\n');
for i = 1:length(model.rxns)
    fprintf(fid,'    - !!omap\n');
    writeField(model, fid, 'rxns',                 'txt', i, '  - id',                    preserveQuotes)
    writeField(model, fid, 'rxnNames',             'txt', i, '  - name',                  preserveQuotes)
    writeField(model, fid, 'S',                    'txt', i, '  - metabolites',           preserveQuotes)
    writeField(model, fid, 'lb',                   'num', i, '  - lower_bound',           preserveQuotes)
    writeField(model, fid, 'ub',                   'num', i, '  - upper_bound',           preserveQuotes)
    writeField(model, fid, 'grRules',              'txt', i, '  - gene_reaction_rule',    preserveQuotes)
    writeField(model, fid, 'rxnFrom',              'txt', i, '  - rxnFrom',               preserveQuotes)
    if model.c(i)~=0
        writeField(model, fid, 'c',                    'num', i, '  - objective_coefficient', preserveQuotes)    
    end
    writeField(model, fid, 'eccodes',              'txt', i, '  - eccodes',               preserveQuotes)
    writeField(model, fid, 'rxnReferences',        'txt', i, '  - references',            preserveQuotes)
    writeField(model, fid, 'subSystems',           'txt', i, '  - subsystem',             preserveQuotes)
    writeField(model, fid, 'rxnMiriams',           'txt', i, '  - annotation',            preserveQuotes)
    writeField(model, fid, 'rxnDeltaG',            'num', i, '  - deltaG',                preserveQuotes)
    writeField(model, fid, 'rxnConfidenceScores',  'num', i, '  - confidence_score',      preserveQuotes)
    writeField(model, fid, 'rxnNotes',             'txt', i, '  - rxnNotes',              preserveQuotes)
end

%Genes:
fprintf(fid,'- genes:\n');
for i = 1:length(model.genes)
    fprintf(fid,'    - !!omap\n');
    writeField(model, fid, 'genes',          'txt', i, '  - id',         preserveQuotes)
    writeField(model, fid, 'geneShortNames', 'txt', i, '  - name',       preserveQuotes)
    writeField(model, fid, 'proteins',   'txt', i, '  - protein',    preserveQuotes)
    writeField(model, fid, 'geneMiriams',    'txt', i, '  - annotation', preserveQuotes)
end

%Compartments:
fprintf(fid,'- compartments: !!omap\n');
for i = 1:length(model.comps)
    writeField(model, fid, 'compNames',   'txt', i, ['- ' model.comps{i}], preserveQuotes)
    writeField(model, fid, 'compMiriams', 'txt', i, '- annotation',             preserveQuotes)
end


%EC-model:
if isfield(model,'ec')
    fprintf(fid,'- ec-rxns:\n');
    for i = 1:length(model.ec.rxns)
        fprintf(fid,'  - !!omap\n');
        writeField(model.ec, fid, 'rxns',      'txt', i, '- id',      preserveQuotes)
        writeField(model.ec, fid, 'kcat',      'num', i, '- kcat',    preserveQuotes)
        writeField(model.ec, fid, 'source',    'txt', i, '- source',  preserveQuotes)
        writeField(model.ec, fid, 'notes',     'txt', i, '- notes',   preserveQuotes)
        writeField(model.ec, fid, 'eccodes',   'txt', i, '- eccodes', preserveQuotes)
        writeField(model.ec, fid, 'rxnEnzMat', 'txt', i, '- enzymes', preserveQuotes)
    end

    fprintf(fid,'- ec-enzymes:\n');
    for i = 1:length(model.ec.genes)
        fprintf(fid,'  - !!omap\n');
        writeField(model.ec, fid, 'genes',    'txt', i, '- genes',    preserveQuotes)
        writeField(model.ec, fid, 'enzymes',  'txt', i, '- enzymes',  preserveQuotes)
        writeField(model.ec, fid, 'mw',       'num', i, '- mw',       preserveQuotes)
        writeField(model.ec, fid, 'sequence', 'txt', i, '- sequence', preserveQuotes)
        writeField(model.ec, fid, 'concs',    'num', i, '- concs',    preserveQuotes)
    end
end

%Close file:
fclose(fid);

end

function writeField(model,fid,fieldName,type,pos,name,preserveQuotes)
%Writes a new line in the yaml file if the field exists and the field is
%not empty at the correspoinding position. It's recursive for some fields
%(metMiriams, rxnMiriams, and S)

if isfield(model,fieldName)
    if strcmp(fieldName,'metComps')
        %metComps: write full name
        fieldName = 'comps';
        pos       = model.metComps(pos);
    end
    
    field = model.(fieldName);
    
    if strcmp(fieldName,'metMiriams')
        if ~isempty(model.metMiriams{pos})
            fprintf(fid,'    %s: !!omap\n',name);
            for i=1:size(model.newMetMiriams,2)
                %'i' represents the different miriam names, e.g.
                %kegg.compound or chebi
                if ~isempty(model.newMetMiriams{pos,i})
                    %As during the following writeField call the value of
                    %'i' would be lost, it is temporarily concatenated to
                    %'name' parameter, which will be edited later
                    writeField(model, fid, 'newMetMiriams', 'txt', pos, ['      - ' model.newMetMiriamNames{i} '_' sprintf('%d',i)], preserveQuotes)
                end
            end
        end
        
    elseif strcmp(fieldName,'rxnMiriams')
        if ~isempty(model.rxnMiriams{pos})
            fprintf(fid,'    %s: !!omap\n',name);
            for i=1:size(model.newRxnMiriams,2)
                if ~isempty(model.newRxnMiriams{pos,i})
                    writeField(model, fid, 'newRxnMiriams', 'txt', pos, ['      - ' model.newRxnMiriamNames{i} '_' sprintf('%d',i)], preserveQuotes)
                end
            end
        end
        
    elseif strcmp(fieldName,'geneMiriams')
        if ~isempty(model.geneMiriams{pos})
            fprintf(fid,'    %s: !!omap\n',name);
            for i=1:size(model.newGeneMiriams,2)
                if ~isempty(model.newGeneMiriams{pos,i})
                    writeField(model, fid, 'newGeneMiriams', 'txt', pos, ['      - ' model.newGeneMiriamNames{i} '_' sprintf('%d',i)], preserveQuotes)
                end
            end
        end
        
    elseif strcmp(fieldName,'compMiriams')
        if ~isempty(model.compMiriams{pos})
            fprintf(fid,'    %s: !!omap\n',name);
            for i=1:size(model.newCompMiriams,2)
                if ~isempty(model.newCompMiriams{pos,i})
                    writeField(model, fid, 'newCompMiriams', 'txt', pos, ['      - ' model.newCompMiriamNames{i} '_' sprintf('%d',i)], preserveQuotes)
                end
            end
        end
        
    elseif strcmp(fieldName,'S')
        %S: create header & write each metabolite in a new line
        fprintf(fid,'    %s: !!omap\n',name);
        if sum(field(:,pos) ~= 0) > 0
            model.mets   = model.mets(field(:,pos) ~= 0);
            model.coeffs = field(field(:,pos) ~= 0,pos);
            %Sort metabolites:
            [model.mets,order] = sort(model.mets);
            model.coeffs       = model.coeffs(order);
            for i = 1:length(model.mets)
                writeField(model, fid, 'coeffs',  'num', i, ['      - ' model.mets{i}], preserveQuotes)
            end
        end

    elseif strcmp(fieldName,'rxnEnzMat')
        %S: create header & write each enzyme in a new line
        fprintf(fid,'    %s: !!omap\n',name);
        if sum(field(pos,:) ~= 0) > 0
            model.enzymes = model.enzymes(field(pos,:) ~= 0);
            model.coeffs  = field(pos,field(pos,:) ~= 0);
            %Sort metabolites:
            [model.enzymes,order] = sort(model.enzymes);
            model.coeffs          = model.coeffs(order);
            for i = 1:length(model.enzymes)
                writeField(model, fid, 'coeffs',  'num', i, ['    - ' model.enzymes{i}], preserveQuotes)
            end
        end        

    elseif sum(strcmp({'subSystems','newMetMiriams','newRxnMiriams','newGeneMiriams','newCompMiriams','eccodes'},fieldName)) > 0
        %eccodes/rxnNotes: if 1 write in 1 line, if more create header and list
        if strcmp(fieldName,'subSystems')
            list = field{pos};  %subSystems already comes in a cell array
            if isempty(list)
                return
            end
        elseif strcmp(fieldName,'newMetMiriams')
            index = str2double(regexprep(name,'^.+_',''));
            name  = regexprep(name,'_\d+$','');
            list = strsplit(model.newMetMiriams{pos,index},'; ');
        elseif strcmp(fieldName,'newRxnMiriams')
            index = str2double(regexprep(name,'^.+_',''));
            name  = regexprep(name,'_\d+$','');
            list = strsplit(model.newRxnMiriams{pos,index},'; ');
        elseif strcmp(fieldName,'newGeneMiriams')
            index = str2double(regexprep(name,'^.+_',''));
            name  = regexprep(name,'_\d+$','');
            list = strsplit(model.newGeneMiriams{pos,index},'; ');
        elseif strcmp(fieldName,'newCompMiriams')
            index = str2double(regexprep(name,'^.+_',''));
            name  = regexprep(name,'_\d+$','');
            list = strsplit(model.newCompMiriams{pos,index},'; ');
        elseif ~isempty(field{pos})
            list = strrep(field{pos},' ','');
            list = strsplit(list,';');
        else
            return % empty, needs no line in file
        end
        list=strip(list);

        if length(list) == 1 && ~strcmp(list{1},'') && ~strcmp(fieldName,'subSystems')
            if preserveQuotes
                list = ['"' list{1} '"'];
            end
            if iscell(list)
                list=list{1};
            end
            fprintf(fid,'    %s: %s\n',name,list);
        elseif ischar(list) && strcmp(fieldName,'subSystems')
            if preserveQuotes
                list = ['"' list '"'];
            end
            fprintf(fid,'    %s: %s\n',name,list);            
        elseif length(list) > 1 || strcmp(fieldName,'subSystems')
            if preserveQuotes
                for j=1:numel(list)
                    list{j} = ['"' list{j} '"'];
                end
            end
            fprintf(fid,'    %s:\n',name);
            for i = 1:length(list)
                fprintf(fid,'%s        - %s\n',regexprep(name,'(^\s*).*','$1'),list{i});
            end
        end
        
    elseif sum(pos) > 0
        %All other fields:
        if strcmp(type,'txt')
            value = field{pos};
            if preserveQuotes && ~isempty(value)
                value = ['"',value,'"'];
            end
        elseif strcmp(type,'num')
            if isnan(field(pos))
                value = [];
            else
                value = sprintf('%.15g',full(field(pos)));
            end
        end
        if ~isempty(value)
            fprintf(fid,'    %s: %s\n',name,value);
        end
    end
end
end

function writeMetadata(model,fid)
% Writes model metadata to the yaml file. This information will eventually
% be extracted entirely from the model, but for now, many of the entries
% are hard-coded defaults for HumanGEM.

fprintf(fid, '- metaData:\n');
if isfield(model,'id')
    fprintf(fid, '    id: "%s"\n',  model.id);
else
    fprintf(fid, '    id: "blankID"\n');
end
if isfield(model,'name')
    fprintf(fid, '    name: "%s"\n',model.name);
else
    fprintf(fid, '    name: "blankName"\n');
end
if isfield(model,'version')
    fprintf(fid, '    version: "%s"\n',model.version);
end
fprintf(fid, '    date: "%s"\n',datestr(now,29));  % 29=YYYY-MM-DD
if isfield(model,'annotation')
    if isfield(model.annotation,'defaultLB')
        fprintf(fid, '    defaultLB: "%g"\n',   model.annotation.defaultLB);
    end
    if isfield(model.annotation,'defaultUB')
        fprintf(fid, '    defaultUB: "%g"\n',   model.annotation.defaultUB);
    end
    if isfield(model.annotation,'givenName')
        fprintf(fid, '    givenName: "%s"\n',   model.annotation.givenName);
    end
    if isfield(model.annotation,'familyName')
        fprintf(fid, '    familyName: "%s"\n',  model.annotation.familyName);
    end
    if isfield(model.annotation,'authors')
        fprintf(fid, '    authors: "%s"\n',     model.annotation.authors);
    end
    if isfield(model.annotation,'email')
        fprintf(fid, '    email: "%s"\n',       model.annotation.email);
    end
    if isfield(model.annotation,'organization')
        fprintf(fid, '    organization: "%s"\n',model.annotation.organization);
    end
    if isfield(model.annotation,'taxonomy')
        fprintf(fid, '    taxonomy: "%s"\n',    model.annotation.taxonomy);
    end
    if isfield(model.annotation,'note')
        fprintf(fid, '    note: "%s"\n',        model.annotation.note);
    end
    if isfield(model.annotation,'sourceUrl')
        fprintf(fid, '    sourceUrl: "%s"\n',   model.annotation.sourceUrl);
    end
end
if isfield(model,'ec')
    if model.ec.geckoLight
        geckoLight = 'true';
    else
        geckoLight = 'false';
    end
    fprintf(fid,'    geckoLight: "%s"\n',geckoLight);
end
end
