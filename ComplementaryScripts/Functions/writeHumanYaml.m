function writeYaml(model,name)
% writeYaml
%   Writes a yaml file matching (roughly) the cobrapy yaml structure
%
%   model       a model structure
%   name        name that the file will have
%
%   Usage: writeYaml(model,name)
%
%   Benjamin Sanchez, 2018-04-11
%   Simonas Marcisauskas, 2018-04-13
%

%Check that model is in RAVEN format:
if isfield(model,'rules')
    model = ravenCobraWrapper(model);
end

%Simplify Miriam fields:
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

%Open file:
fid = fopen(name,'wt');
fprintf(fid,'!!omap\n');

%Metabolites:
fprintf(fid,'- metabolites:\n');
[~,pos] = sort(model.mets);
for i = 1:length(model.mets)
    fprintf(fid,'  - !!omap\n');
    writeField(model, fid, 'mets',        'txt', pos(i), '- id')
    writeField(model, fid, 'metNames',    'txt', pos(i), '- name')
    writeField(model, fid, 'metComps',    'txt', pos(i), '- compartment')
    writeField(model, fid, 'metFormulas', 'txt', pos(i), '- formula')
    writeField(model, fid, 'metCharges',  'num', pos(i), '- charge')
    writeField(model, fid, 'metMiriams',  'txt', pos(i), '- annotation')
end

%Reactions:
fprintf(fid,'- reactions:\n');
[~,pos] = sort(model.rxns);
for i = 1:length(model.rxns)
    fprintf(fid,'  - !!omap\n');
    writeField(model, fid, 'rxns',                'txt', pos(i), '- id')
    writeField(model, fid, 'rxnNames',            'txt', pos(i), '- name')
    writeField(model, fid, 'S',                   'txt', pos(i), '- metabolites')
    writeField(model, fid, 'lb',                  'num', pos(i), '- lower_bound')
    writeField(model, fid, 'ub',                  'num', pos(i), '- upper_bound')
    writeField(model, fid, 'grRules',             'txt', pos(i), '- gene_reaction_rule')
    writeField(model, fid, 'subSystems',          'txt', pos(i), '- subsystem')
    writeField(model, fid, 'rxnMiriams',          'txt', pos(i), '- annotation')
    writeField(model, fid, 'rxnConfidenceScores', 'num', pos(i), '- confidence_score')
end

%Genes:
fprintf(fid,'- genes:\n');
[~,pos] = sort(model.genes);
for i = 1:length(model.genes)
    fprintf(fid,'  - !!omap\n');
    writeField(model, fid, 'genes',          'txt', pos(i), '- id')
    writeField(model, fid, 'geneShortNames', 'txt', pos(i), '- name')
    writeField(model, fid, 'geneMiriams',    'txt', pos(i), '- annotation')
end

%Compartments:
fprintf(fid,'- compartments: !!omap\n');
[~,pos] = sort(model.comps);
for i = 1:length(model.comps)
    writeField(model, fid, 'compNames',   'txt', pos(i), ['- ' model.comps{pos(i)}])
    writeField(model, fid, 'compMiriams', 'txt', pos(i), '- annotation')
end

%TODO: include id, name & version (lost in RAVEN)

%Close file:
fclose(fid);

end

function writeField(model,fid,fieldName,type,pos,name)
%Writes a new line in the yaml file if the field exists and the field is
%not empty at the correspoinding position. It's recursive for some fields
%(metMiriams, rxnMiriams, and S)

if isfield(model,fieldName)
    if strcmp(fieldName,'metComps')
        %metComps: write full name
        fieldName = 'comps';
        pos       = model.metComps(pos);
    end
    
    field = eval(['model.' fieldName]);
    
    if strcmp(fieldName,'metMiriams')
        if ~isempty(model.metMiriams{pos})
            fprintf(fid,['    ' name ': !!omap\n']);
            for i=1:size(model.newMetMiriams,2)
                %'i' represents the different miriam names, e.g.
                %kegg.compound or chebi
                if ~isempty(model.newMetMiriams{pos,i})
                    %As during the following writeField call the value of
                    %'i' would be lost, it is temporarily concatenated to
                    %'name' parameter, which will be edited later
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
        %S: create header & write each metabolite in a new line
        fprintf(fid,['    ' name ': !!omap\n']);
        if sum(field(:,pos) ~= 0) > 0
            model.mets   = model.mets(field(:,pos) ~= 0);
            model.coeffs = field(field(:,pos) ~= 0,pos);
            %Sort metabolites:
            [model.mets,order] = sort(model.mets);
            model.coeffs       = model.coeffs(order);
            for i = 1:length(model.mets)
                writeField(model, fid, 'coeffs',  'num', i, ['  - ' model.mets{i}])
            end
        end
        
    elseif sum(strcmp({'eccodes','subSystems','newMetMiriams','newRxnMiriams','newGeneMiriams','newCompMiriams'},fieldName)) > 0
        %eccodes/rxnNotes/subSystems: if 1 write in 1 line, if more create
        %header and list
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
            list = strrep(field{pos},' ','');     %Exception for eccodes
            list = strsplit(list,';');
        end
        if length(list) == 1 && ~strcmp(list{1},'')
            fprintf(fid,['    ' name ': ' list{1} '\n']);
        elseif length(list) > 1
            fprintf(fid,['    ' name ':\n']);
            for i = 1:length(list)
                fprintf(fid,['        - ' list{i} '\n']);
            end
        end
        
    elseif sum(pos) > 0
        %All other fields:
        if strcmp(type,'txt')
            value = field{pos};
        elseif strcmp(type,'num')
            if isnan(field(pos))
                value = [];
            else
                value = num2str(field(pos));
            end
        end
        if ~isempty(value)
            fprintf(fid,['    ' name ': ' value '\n']);
        end
    end
end


end