function writeHumanYaml(model,name)
% writeHumanYaml  Write a yaml file for HumanGEM.
%
%   Writes a yaml file that is similar to the cobrapy yaml structure, but
%   contains some minor changes and a header file for compatibility with
%   the Metabolic Atlas. Adapted from RAVEN's "writeYaml" function.
%
% Usage: 
%   writeYaml(model,name);
%
% Input:
%   model       a model structure
%   name        name of the yaml file to write
%

%{
% check for yaml extension
if ~endsWith(name,{'.yml','.yaml'})
    name = strcat(name,'.yaml');
end

% check that model is in RAVEN format
if isfield(model,'rules')
    model = ravenCobraWrapper(model);
end

% simplify Miriam fields
if isfield(model,'metMiriams')
    [model.newMetMiriams,model.newMetMiriamNames]   = extractMiriam(model.metMiriams);
    model.newMetMiriams                             = regexprep(model.newMetMiriams,'^.+/','');
end
if isfield(model,'rxnMiriams')
    [model.newRxnMiriams,model.newRxnMiriamNames]   = extractMiriam(model.rxnMiriams);
    model.newRxnMiriams                             = regexprep(model.newRxnMiriams,'^.+/','');
end
if isfield(model,'geneMiriams')
    [model.newGeneMiriams,model.newGeneMiriamNames] = extractMiriam(model.geneMiriams);
    model.newGeneMiriams                            = regexprep(model.newGeneMiriams,'^.+/','');
end
if isfield(model,'compMiriams')
    [model.newCompMiriams,model.newCompMiriamNames] = extractMiriam(model.compMiriams);
    model.newCompMiriams                            = regexprep(model.newCompMiriams,'^.+/','');
end
%}

% open file
fid = fopen(name,'wt');
fprintf(fid,'---\n!!omap\n');

% insert file header (metadata)
writeMetadata(model,fid); 

% metabolites
fprintf(fid,'- metabolites:\n');
%[~,pos] = sort(model.mets);
for i = 1:length(model.mets)
    fprintf(fid,'    - !!omap\n');
    writeField(model, fid, 'mets',        'txt', i, '  - id')
    writeField(model, fid, 'metNames',    'txt', i, '  - name')
    writeField(model, fid, 'metComps',    'txt', i, '  - compartment')
    writeField(model, fid, 'metFormulas', 'txt', i, '  - formula')
    writeField(model, fid, 'metCharges',  'num', i, '  - charge')
    writeField(model, fid, 'inchis',      'txt', i, '  - inchis')
    writeField(model, fid, 'metFrom',     'txt', i, '  - metFrom')
%     writeField(model, fid, 'metMiriams',  'txt', i, '  - annotation')
end

% reactions
fprintf(fid,'- reactions:\n');
%[~,pos] = sort(model.rxns);
for i = 1:length(model.rxns)
    fprintf(fid,'    - !!omap\n');
    writeField(model, fid, 'rxns',                 'txt', i, '  - id')
    writeField(model, fid, 'rxnNames',             'txt', i, '  - name')
    writeField(model, fid, 'S',                    'txt', i, '  - metabolites')
    writeField(model, fid, 'lb',                   'num', i, '  - lower_bound')
    writeField(model, fid, 'ub',                   'num', i, '  - upper_bound')
    writeField(model, fid, 'grRules',              'txt', i, '  - gene_reaction_rule')
    writeField(model, fid, 'priorCombiningGrRules','txt', i, '  - HMR2_grRule')    
    writeField(model, fid, 'rxnFrom',              'txt', i, '  - rxnFrom')
    if model.c(i)
        writeField(model, fid, 'c',                'num', i, '  - objective_coefficient')
    end
    writeField(model, fid, 'eccodes',              'txt', i, '  - eccodes')
    writeField(model, fid, 'rxnReferences',        'txt', i, '  - references')
    writeField(model, fid, 'subSystems',           'txt', i, '  - subsystem')
%     writeField(model, fid, 'rxnMiriams',          'txt', i, '  - annotation')
    writeField(model, fid, 'rxnConfidenceScores',  'num', i, '  - confidence_score')
end

% genes
fprintf(fid,'- genes:\n');
%[~,pos] = sort(model.genes);
for i = 1:length(model.genes)
    fprintf(fid,'    - !!omap\n');
    writeField(model, fid, 'genes',          'txt', i, '  - id')
    writeField(model, fid, 'geneShortNames', 'txt', i, '  - name')
%     writeField(model, fid, 'geneMiriams',    'txt', i, '  - annotation')
end

% compartments
fprintf(fid,'- compartments: !!omap\n');
%[~,pos] = sort(model.comps);
for i = 1:length(model.comps)
    writeField(model, fid, 'compNames',   'txt', i, ['- ' model.comps{i}])
%     writeField(model, fid, 'compMiriams', 'txt', i, '- annotation')
end

% close file:
fclose(fid);

end

function writeField(model,fid,fieldName,type,pos,name)
% It's recursive for some fields (metMiriams, rxnMiriams, and S)

if isfield(model,fieldName)
    if strcmp(fieldName,'metComps')
        % metComps: write full name
        fieldName = 'comps';
        pos       = model.metComps(pos);
    end
    
    field = eval(['model.' fieldName]);
    
    if strcmp(fieldName,'metMiriams')
        if ~isempty(model.metMiriams{pos})
            fprintf(fid,['    ' name ': !!omap\n']);
            for i=1:size(model.newMetMiriams,2)
                % 'i' represents the different miriam names, e.g.
                % kegg.compound or chebi
                if ~isempty(model.newMetMiriams{pos,i})
                    % As during the following writeField call the value of
                    % 'i' would be lost, it is temporarily concatenated to
                    % 'name' parameter, which will be edited later
                    writeField(model, fid, 'newMetMiriams', 'txt', pos, ['  - ' model.newMetMiriamNames{i} '_' num2str(i)])
                end
            end
        end
        
    elseif strcmp(fieldName,'rxnMiriams')
        if ~isempty(model.eccodes{pos}) || ~isempty(model.rxnMiriams{pos})
            fprintf(fid,['    ' name ': !!omap\n']);
            writeField(model, fid, 'eccodes',  'txt', pos, '  - ec-code')
            for i=1:size(model.newRxnMiriams,2)
                if ~isempty(model.newRxnMiriams{pos,i})
                    writeField(model, fid, 'newRxnMiriams', 'txt', pos, ['  - ' model.newRxnMiriamNames{i} '_' num2str(i)])
                end
            end
        end
        
    elseif strcmp(fieldName,'geneMiriams')
        if ~isempty(model.geneMiriams{pos})
            fprintf(fid,['    ' name ': !!omap\n']);
            for i=1:size(model.newGeneMiriams,2)
                if ~isempty(model.newGeneMiriams{pos,i})
                    writeField(model, fid, 'newGeneMiriams', 'txt', pos, ['  - ' model.newGeneMiriamNames{i} '_' num2str(i)])
                end
            end
        end
        
    elseif strcmp(fieldName,'compMiriams')
        if ~isempty(model.compMiriams{pos})
            fprintf(fid,['    ' name ': !!omap\n']);
            for i=1:size(model.newCompMiriams,2)
                if ~isempty(model.newCompMiriams{pos,i})
                    writeField(model, fid, 'newCompMiriams', 'txt', pos, ['  - ' model.newCompMiriamNames{i} '_' num2str(i)])
                end
            end
        end
        
    elseif strcmp(fieldName,'S')
        % S: create header & write each metabolite in a new line
        fprintf(fid,['    ' name ': !!omap\n']);
        if sum(field(:,pos) ~= 0) > 0
            model.mets   = model.mets(field(:,pos) ~= 0);
            model.coeffs = field(field(:,pos) ~= 0,pos);
            % sort metabolites
            [model.mets,order] = sort(model.mets);
            model.coeffs       = model.coeffs(order);
            for i = 1:length(model.mets)
                writeField(model, fid, 'coeffs',  'num', i, ['      - ' model.mets{i}])
            end
        end
        
    %elseif sum(strcmp({'eccodes','rxnReferences','subSystems','newMetMiriams','newRxnMiriams','newGeneMiriams','newCompMiriams'},fieldName)) > 0
    elseif sum(strcmp({'subSystems','newMetMiriams','newRxnMiriams','newGeneMiriams','newCompMiriams'},fieldName)) > 0
        % eccodes/rxnNotes: if 1 write in 1 line, if more create header and list
        if strcmp(fieldName,'subSystems')
            list = field{pos};  %subSystems already comes in a cell array
        elseif strcmp(fieldName,'newMetMiriams')
            index = str2double(regexprep(name,'^.+_',''));
            name  = regexprep(name,'_\d+$','');
            list = strsplit(model.newMetMiriams{pos,index},';');
        elseif strcmp(fieldName,'newRxnMiriams')
            index = str2double(regexprep(name,'^.+_',''));
            name  = regexprep(name,'_\d+$','');
            list = strsplit(model.newRxnMiriams{pos,index},';');
        elseif strcmp(fieldName,'newGeneMiriams')
            index = str2double(regexprep(name,'^.+_',''));
            name  = regexprep(name,'_\d+$','');
            list = strsplit(model.newGeneMiriams{pos,index},';');
        elseif strcmp(fieldName,'newCompMiriams')
            index = str2double(regexprep(name,'^.+_',''));
            name  = regexprep(name,'_\d+$','');
            list = strsplit(model.newCompMiriams{pos,index},';');
        else
            list = strrep(field{pos},' ','');     % exception for eccodes
            list = strsplit(list,';');
        end
        if length(list) == 1 && ~strcmp(list{1},'') && ~strcmp(fieldName,'subSystems')
            fprintf(fid,['    ' name ': ' list{1} '\n']);
        elseif length(list) > 1 || strcmp(fieldName,'subSystems')
            fprintf(fid,['    ' name ':\n']);
            for i = 1:length(list)
                fprintf(fid,['          - "' list{i} '"\n']);
            end
        end
        
    elseif sum(pos) > 0
        % all other fields
        if strcmp(type,'txt')
            % enclose all string elements with double quotes
            value = ['"',field{pos},'"'];
        elseif strcmp(type,'num')
            if isnan(field(pos))
                value = [];
            else
                value = num2str(field(pos), 12);
            end
        end
        %if ~isempty(value)
            fprintf(fid,['    ' name ': ' value '\n']);
        %end
    end
end

end


function writeMetadata(model,fid)
% Writes model metadata to the yaml file. This information will eventually
% be extracted entirely from the model, but for now, many of the entries
% are hard-coded defaults for HumanGEM.

fprintf(fid, '- metaData:\n');
fprintf(fid, ['    short_name: "', model.id, '"\n']);
fprintf(fid, ['    full_name: "', model.description, '"\n']);
fprintf(fid, ['    version: "', model.version, '"\n']);
fprintf(fid, ['    date: "', datestr(now,29), '"\n']);  % 29=YYYY-MM-DD
fprintf(fid, ['    authors: "', model.annotation.authorList, '"\n']);
fprintf(fid, ['    email: "', model.annotation.email, '"\n']);
fprintf(fid, ['    organization: "', model.annotation.organization, '"\n']);
fprintf(fid, ['    taxonomy: "', model.annotation.taxonomy, '"\n']);
fprintf(fid, ['    github: "', model.annotation.sourceUrl, '"\n']);
fprintf(fid, ['    description: "', model.annotation.note, '"\n']);

end
