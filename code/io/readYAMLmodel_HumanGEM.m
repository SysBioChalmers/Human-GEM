function model=readYAMLmodel(fileName, verbose)
% readYAMLmodel
%   Reads a yaml file matching (roughly) the cobrapy yaml structure
%
%   Input:
%   fileName    a model file in yaml file format. A dialog window will open
%               if no file name is specified.
%   verbose     set as true to monitor progress (optional, default false)
%
%   Output:
%   model       a model structure
%
% Usage: model = readYAMLmodel(fileName, verbose)
if nargin<1 || isempty(fileName)
    [fileName, pathName] = uigetfile({'*.yml;*.yaml'}, 'Please select the model file');
    if fileName == 0
        error('You should select a model file')
    else
        fileName = fullfile(pathName,fileName);
    end
end
if nargin < 2
    verbose = false;
end

if ~isfile(fileName)
    error('Yaml file %s cannot be found', string(fileName));
end

if verLessThan('matlab','9.9') %readlines introduced 2020b
    fid=fopen(fileName);
    line_raw=cell(1000000,1);
    i=1;
    while ~feof(fid)
        line_raw{i}=fgetl(fid);
        i=i+1;
    end
    line_raw(i:end)=[];
    line_raw=string(line_raw);
    fclose(fid);
else
    line_raw=readlines(fileName);
end
% If entry is broken of multiple lines, concatenate. Assumes at least 6
% leading spaces to avoid metaData to be concatenated.
newLine=regexp(line_raw,'^ {6,}([\w\(\)].*)','tokenExtents');
brokenLine=find(~cellfun('isempty',newLine));
for i=1:numel(brokenLine)
    extraLine = char(line_raw(brokenLine(i)));
    extraLine = extraLine(newLine{brokenLine(i)}{1}(1):end);
    line_raw{brokenLine(i)-1} = strjoin({line_raw{brokenLine(i)-1},extraLine},' ');
end
line_raw(brokenLine)=[];

line_key = regexprep(line_raw,'^ *-? ([^:]+)(:)($| .*)','$1');
line_key = regexprep(line_key,'(.*!!omap)|(---)|( {4,}.*)','');

line_value = regexprep(line_raw, '.*:$','');
line_value = regexprep(line_value, '[^":]+: "?(.+)"?$','$1');
line_value = regexprep(line_value, '(")|(^ {4,}- )','');

modelFields =   {'id',char();...
               'name',char();...
        'description',char();...
            'version',char();...
               'date',char();...
         'annotation',struct();...
               'rxns',{};...
           'rxnNames',{};...
               'mets',{};...
           'metNames',{};...
                  'S',sparse([]);...
                 'lb',{};... %Changed to double in the end.
                 'ub',{};... %Changed to double in the end.
                'rev',{};... %Changed to double in the end.
                  'c',[];...
                  'b',cell(0,0);... %Changed to double in the end.
              'genes',cell(0,0);...
            'grRules',cell(0,0);...
         'rxnGeneMat',sparse([]);...
           'rxnComps',cell(0,0);... %Changed to double in the end.
         'subSystems',cell(0,0);...
            'eccodes',cell(0,0);...
         'rxnMiriams',cell(0,0);...
          'rxnDeltaG',{};... %Changed to double in the end.
           'rxnNotes',cell(0,0);...
      'rxnReferences',cell(0,0);...
'rxnConfidenceScores',cell(0,0);...
           'metComps',cell(0,0);... %Changed to double in the end.
             'inchis',cell(0,0);...
          'metSmiles',cell(0,0);...
        'metFormulas',cell(0,0);...
         'metMiriams',cell(0,0);...
          'metDeltaG',{};... %Changed to double in the end.
         'metCharges',cell(0,0);... %Changed to double in the end.
           'metNotes',cell(0,0);...
              'comps',cell(0,0);...
          'compNames',cell(0,0);...
        'compOutside',cell(0,0);...
          'geneComps',cell(0,0);... %Changed to double in the end.
        'geneMiriams',cell(0,0);...
     'geneShortNames',cell(0,0);...
       'proteins',cell(0,0);...
      'unconstrained',cell(0,0);... %Changed to double in the end.
            'metFrom',cell(0,0);...
            'rxnFrom',cell(0,0)};
for i=1:size(modelFields,1)
    model.(modelFields{i,1})=modelFields{i,2};
end

% If GECKO model
if any(contains(line_key,'geckoLight'))
    isGECKO=true;
    ecFields = {'geckoLight', false;...
                      'rxns', {};...
                      'kcat', {};...
                    'source', cell(0,0);...
                     'notes', cell(0,0);...
                   'eccodes', cell(0,0);...
                     'genes', cell(0,0);...
                   'enzymes', cell(0,0);...
                        'mw', cell(0,0);...
                  'sequence', cell(0,0);...
                     'concs', cell(0,0);...
                 'rxnEnzMat', []};
    for i=1:size(ecFields,1)
        model.ec.(ecFields{i,1})=ecFields{i,2};
    end
    ecGecko=cell(25000,2);      ecGeckoNo=1;
    enzStoich=cell(100000,3);   enzStoichNo=1;
else
    isGECKO=false;
end

section = 0;
metMiriams=cell(100000,3);   metMirNo=1;
rxnMiriams=cell(100000,3);   rxnMirNo=1;
geneMiriams=cell(100000,3);  genMirNo=1;
subSystems=cell(100000,2);   subSysNo=1;
eccodes=cell(100000,2);      ecCodeNo=1;
equations=cell(100000,3);   equatiNo=1;

for i=1:numel(line_key)
    tline_raw = line_raw{i};
    tline_key = line_key{i};
    tline_value = line_value{i};
    % import different sections
    switch tline_raw
        case '- metaData:'
            section = 1;
            if verbose
                fprintf('\t%d\n', section);
            end
            continue % Go to next line
        case '- metabolites:'
            section = 2;
            if verbose
                fprintf('\t%d\n', section);
            end
            pos=0;
            continue
        case '- reactions:'
            section = 3;
            if verbose
                fprintf('\t%d\n', section);
            end
            pos=0;
            continue
        case '- genes:'
            section = 4;
            if verbose
                fprintf('\t%d\n', section);
            end
            pos=0;
            continue
        case '- compartments: !!omap'
            section = 5;
            if verbose
                fprintf('\t%d\n', section);
            end
            pos=0;
            continue
        case '- ec-rxns:'
            section = 6;
            if verbose
                fprintf('\t%d\n', section);
            end
            pos=0;
            continue
        case '- ec-enzymes:'
            section = 7;
            if verbose
                fprintf('\t%d\n', section);
            end
            pos=0;
            continue
    end

    % skip over empty keys
    if isempty(tline_raw) || (isempty(tline_key) && contains(tline_raw,'!!omap'))
        continue;
    end
    
    % import metaData
    if section == 1
        switch tline_key
            case {'short_name','id'} %short_name used by human-GEM
                model.id = tline_value;
            case 'name'
                model.name = tline_value;                
            case 'full_name' %used by human-GEM
                model.description = tline_value;
            case 'version'
                model.version = tline_value;
            case 'date'
                model.date = tline_value;                
            case 'taxonomy'
                model.annotation.taxonomy = tline_value;
            case {'description','note'} %description used by human-GEM
                model.annotation.note = tline_value;
            case 'github'
                model.annotation.sourceUrl = tline_value;
            case 'sourceUrl'
                model.annotation.sourceUrl = tline_value;                
            case 'givenName'
                model.annotation.givenName = tline_value;
            case 'familyName'
                model.annotation.familyName = tline_value;
            case 'authors'
                model.annotation.authorList = tline_value;
            case 'email'
                model.annotation.email = tline_value;
            case 'organization'
                model.annotation.organization = tline_value;
            case 'geckoLight'
                if strcmp(tline_value,'true')
                    model.ec.geckoLight = true;
                end
        end; continue
    end

    % import metabolites:
    if section == 2
        switch tline_key
            case 'id'
                pos = pos + 1;
                model = readFieldValue(model, 'mets', tline_value,pos);
                readList=''; miriamKey='';
            case 'name'
                model = readFieldValue(model, 'metNames', tline_value, pos);
                readList=''; miriamKey='';
            case 'compartment'
                model = readFieldValue(model, 'metComps', tline_value, pos);
                readList=''; miriamKey='';
            case 'formula'
                model = readFieldValue(model, 'metFormulas', tline_value, pos);
                readList=''; miriamKey='';
            case 'charge'
                model = readFieldValue(model, 'metCharges', tline_value, pos);
                readList=''; miriamKey='';
            case 'notes'
                model = readFieldValue(model, 'metNotes', tline_value, pos);
                readList=''; miriamKey='';
            case 'inchis'
                model = readFieldValue(model, 'inchis', tline_value, pos);
                readList=''; miriamKey='';
            case 'smiles'
                model = readFieldValue(model, 'metSmiles', tline_value, pos);
                readList=''; miriamKey='';    
            case 'deltaG'
                model = readFieldValue(model, 'metDeltaG', tline_value, pos);
                readList=''; miriamKey='';                                
            case 'metFrom'
                model = readFieldValue(model, 'metFrom', tline_value, pos);
                readList=''; miriamKey='';
            case 'annotation'
                readList = 'annotation';
            otherwise
                switch readList
                    case 'annotation'
                        [metMiriams, miriamKey, metMirNo] = gatherAnnotation(pos,metMiriams,tline_key,tline_value,miriamKey,metMirNo);
                    otherwise
                        error(['Unknown entry in yaml file: ' tline_raw])
                end
        end; continue
    end

    % import reactions:
    if section == 3
        switch tline_key
            case 'id'
                pos = pos + 1;
                model = readFieldValue(model, 'rxns', tline_value,pos);
                readList=''; miriamKey='';
            case 'name'
                model = readFieldValue(model, 'rxnNames', tline_value, pos);
                readList=''; miriamKey='';
            case 'lower_bound'
                model.lb(pos,1) = {tline_value};
                readList=''; miriamKey='';
            case 'upper_bound'
                model.ub(pos,1) = {tline_value};
                readList=''; miriamKey='';
            case 'rev'
                model.rev(pos,1) = {tline_value};
                readList=''; miriamKey='';
            case 'gene_reaction_rule'
                model = readFieldValue(model, 'grRules', tline_value, pos);
                readList=''; miriamKey='';
            case 'rxnNotes'
                model = readFieldValue(model, 'rxnNotes', tline_value, pos);
                readList=''; miriamKey='';
            case 'rxnFrom'
                model = readFieldValue(model, 'rxnFrom', tline_value, pos);
                readList=''; miriamKey='';
            case 'deltaG'
                model = readFieldValue(model, 'rxnDeltaG', tline_value, pos);
                readList=''; miriamKey='';                
            case 'objective_coefficient'
                model.c(pos,1) = 1;
                readList=''; miriamKey='';
            case 'references'
                model = readFieldValue(model, 'rxnReferences', tline_value, pos);
                readList=''; miriamKey='';
            case 'confidence_score'
                model = readFieldValue(model, 'rxnConfidenceScores', tline_value, pos);
                readList=''; miriamKey='';
            case 'eccodes'
                if isempty(tline_value)
                    readList = 'eccodes';
                else
                    eccodes(ecCodeNo,1:2)={pos,tline_value};
                    ecCodeNo=ecCodeNo+1;
                end
            case 'subsystem'
                if isempty(tline_value)
                    readList = 'subsystem';
                else
                    subSystems(subSysNo,1:2)={pos,tline_value};
                    subSysNo=subSysNo+1;
                end
            case 'metabolites'
                readList = 'equation';
            case 'annotation'
                readList = 'annotation';
                
            otherwise
                switch readList
                    case 'eccodes'
                        eccodes(ecCodeNo,1:2)={pos,regexprep(tline_value,'^ +- "?(.*)"?$','$1')};
                        ecCodeNo=ecCodeNo+1;
                    case 'subsystem'
                        subSystems(subSysNo,1:2)={pos,regexprep(tline_value,'^ +- "?(.*)"?$','$1')};
                        subSysNo=subSysNo+1;
                    case 'annotation'
                        [rxnMiriams, miriamKey,rxnMirNo] = gatherAnnotation(pos,rxnMiriams,tline_key,tline_value,miriamKey,rxnMirNo);
                    case 'equation'
                        coeff = sscanf(tline_value,'%f');
                        equations(equatiNo,1:3)={pos,tline_key,coeff};
                        equatiNo=equatiNo+1;
                    otherwise
                        error(['Unknown entry in yaml file: ' tline_raw])
                end                
        end; continue
    end

    % import genes:
    if section == 4
        switch tline_key
            case 'id'
                pos = pos + 1;
                model = readFieldValue(model, 'genes', tline_value, pos);
                readList = '';
                miriamKey = '';
            case 'name'
                model = readFieldValue(model, 'geneShortNames', tline_value, pos);
            case 'protein'
                model = readFieldValue(model, 'proteins', tline_value, pos);
            case 'annotation'
                readList = 'annotation';
            otherwise
                switch readList
                    case 'annotation'
                        [geneMiriams, miriamKey,genMirNo] = gatherAnnotation(pos,geneMiriams,tline_key,tline_value,miriamKey,genMirNo);
                    otherwise
                        error(['Unknown entry in yaml file: ' tline_raw])
                end
        end; continue
    end

    % import compartments:
    if section == 5
        model.comps(end+1,1) = {tline_key};
        model.compNames(end+1,1) = {tline_value};
    end

    % import ec reaction info
    if section == 6
        switch tline_key
            case 'id'
                pos = pos + 1;
                model.ec = readFieldValue(model.ec, 'rxns', tline_value, pos);
                readList='';
            case 'kcat'
                model.ec = readFieldValue(model.ec, 'kcat', tline_value, pos);
                readList='';
            case 'source'
                model.ec = readFieldValue(model.ec, 'source', tline_value, pos);
                readList='';
            case 'notes'
                model.ec = readFieldValue(model.ec, 'notes', tline_value, pos);
                readList='';
            case 'eccodes'
                if isempty(tline_value)
                    readList = 'eccodes';
                else
                    ecGecko(ecGeckoNo,1:2)={pos,tline_value};
                    ecGeckoNo=ecGeckoNo+1;
                end
            case 'enzymes'
                readList = 'enzStoich';
            otherwise
                switch readList
                    case 'eccodes'
                        ecGecko(ecGeckoNo,1:2)={pos,regexprep(tline_value,'^ +- "?(.*)"?$','$1')};
                        ecGeckoNo=ecGeckoNo+1;
                    case 'enzStoich'
                        coeff = sscanf(tline_value,'%f');
                        enzStoich(enzStoichNo,1:3)={pos,tline_key,coeff};
                        enzStoichNo=enzStoichNo+1;
                    otherwise
                        error(['Unknown entry in yaml file: ' tline_raw])
                end
        end; continue
    end

    % import ec enzyme info
    if section == 7
        switch tline_key
            case 'genes'
                pos = pos + 1;
                model.ec = readFieldValue(model.ec, 'genes', tline_value, pos);
            case 'enzymes'
                model.ec = readFieldValue(model.ec, 'enzymes', tline_value, pos);
            case 'mw'
                model.ec = readFieldValue(model.ec, 'mw', tline_value, pos);
            case 'sequence'
                model.ec = readFieldValue(model.ec, 'sequence', tline_value, pos);
            case 'concs'
                model.ec = readFieldValue(model.ec, 'concs', tline_value, pos);
            otherwise
                error(['Unknown entry in yaml file: ' tline_raw])
        end; continue
    end
end

%Parse annotations
if ~isempty(metMiriams)
    locs = cell2mat(metMiriams(:,1));
    for i=unique(locs)'
        model.metMiriams{i,1}.name=metMiriams(locs==i,2);
        model.metMiriams{i,1}.value=metMiriams(locs==i,3);
    end
end
if ~isempty(rxnMiriams)
    locs = cell2mat(rxnMiriams(:,1));
    for i=unique(locs)'
        model.rxnMiriams{i,1}.name=rxnMiriams(locs==i,2);
        model.rxnMiriams{i,1}.value=rxnMiriams(locs==i,3);
    end
end
if ~isempty(geneMiriams)
    locs = cell2mat(geneMiriams(:,1));
    for i=unique(locs)'
        model.geneMiriams{i,1}.name=geneMiriams(locs==i,2);
        model.geneMiriams{i,1}.value=geneMiriams(locs==i,3);
    end
end

%Parse subSystems
if ~isempty(subSystems)
    locs = cell2mat(subSystems(:,1));
    for i=unique(locs)'
        model.subSystems{i,1}=subSystems(locs==i,2);
    end
end

%Parse ec-codes
if ~isempty(eccodes)
    locs = cell2mat(eccodes(:,1));
    for i=unique(locs)'
        eccodesCat=strjoin(eccodes(locs==i,2),';');
        model.eccodes{i,1}=eccodesCat;
    end
    emptyEc=cellfun('isempty',model.eccodes);
    model.eccodes(emptyEc)={''};
end

% follow-up data processing
if verbose
    fprintf('\nimporting completed\nfollow-up processing...');
end
[~, model.metComps] = ismember(model.metComps, model.comps);
[~, model.geneComps] = ismember(model.geneComps, model.comps);
[~, model.rxnComps] = ismember(model.rxnComps, model.comps);

% Fill S-matrix
rxnIdx = cellfun('isempty', equations(:,1));
equations(rxnIdx,:) = '';
rxnIdx = cell2mat(equations(:,1));
[~,metIdx] = ismember(equations(:,2),model.mets);
coeffs = cell2mat(equations(:,3));
model.S=sparse(max(metIdx),max(rxnIdx));
linearIndices = sub2ind([max(metIdx), max(rxnIdx)],metIdx,rxnIdx);
model.S(linearIndices) = coeffs;

% Convert strings to numeric
model.metCharges = str2double(model.metCharges);
model.lb = str2double(model.lb);
model.ub = str2double(model.ub);
model.rxnConfidenceScores = str2double(model.rxnConfidenceScores);
model.b = zeros(length(model.mets),1);
model.metDeltaG = str2double(model.metDeltaG);
model.rxnDeltaG = str2double(model.rxnDeltaG);

% Fill some other fields
model.annotation.defaultLB = min(model.lb);
model.annotation.defaultUB = max(model.ub);
if numel(model.lb)<numel(model.rxns) %No LB reported = min
    model.lb(end+1:numel(model.rxns)-numel(model.lb),1) = double(model.annotation.defaultLB);
end
if numel(model.ub)<numel(model.rxns) %No UB reported = max
    model.ub(end+1:numel(model.rxns)-numel(model.ub),1) = double(model.annotation.defaultUB);
end
if ~all(cellfun('isempty',model.rev))
    model.rev = str2double(model.rev);
else
    model.rev = [];
end
if numel(model.rev)<numel(model.rxns) %No rev reported, assume from LB and UB
    model.rev(end+1:numel(model.rxns)-numel(model.rev),1) = double(model.lb<0 & model.ub>0);
end

% Remove empty fields, otherwise fill to correct length
% Reactions
for i={'rxnNames','grRules','eccodes','rxnNotes','rxnReferences',...
       'rxnFrom','subSystems','rxnMiriams'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'rxns');
end
for i={'c'} % Zeros
   model = emptyOrFill(model,i{1},0,'rxns',true);
end
for i={'rxnConfidenceScores','rxnDeltaG'} % NaNs
   model = emptyOrFill(model,i{1},NaN,'rxns');
end
for i={'rxnComps'} % Ones, assume first compartment
   model = emptyOrFill(model,i{1},1,'rxns');
end
% Metabolites
for i={'metNames','inchis','metFormulas','metMiriams','metFrom','metSmiles','metNotes'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'mets');
end
for i={'metCharges','unconstrained'} % Zeros
   model = emptyOrFill(model,i{1},0,'mets');
end
for i={'metDeltaG'} % % NaNs
    model = emptyOrFill(model,i{1},NaN,'mets');
 end
for i={'metComps'} % Ones, assume first compartment
   model = emptyOrFill(model,i{1},1,'mets');
end
% Genes
for i={'geneMiriams','geneShortNames','proteins'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'genes');
end
for i={'geneComps'} % Ones, assume first compartment
   model = emptyOrFill(model,i{1},1,'genes');
end
% Comps
for i={'compNames'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'comps');
end
for i={'compOutside'} % First comp
   model = emptyOrFill(model,i{1},model.comps{1},'comps');
end
% Single fields are kept, even if empty
% for i={'description','name','version','date','annotation'}
%     if isempty(model.(i{1}))
%         model = rmfield(model,i{1});
%     end
% end

% Make rxnGeneMat fields and map to the existing model.genes field
[genes, rxnGeneMat] = getGenesFromGrRules(model.grRules);
model.rxnGeneMat = sparse(numel(model.rxns),numel(model.genes));
[~,geneOrder] = ismember(genes,model.genes);
if any(geneOrder == 0)
    error(['The grRules includes the following gene(s), that are not in '...
           'the list of model genes: ', genes{~geneOrder}])
end
model.rxnGeneMat(:,geneOrder) = rxnGeneMat;

% Finalize GECKO model
if isGECKO
    % Fill in empty fields and empty entries
    for i={'kcat','source','notes','eccodes'} % Even keep empty
        model.ec = emptyOrFill(model.ec,i{1},{''},'rxns',true);
    end
    for i={'enzymes','mw','sequence'}
        model.ec = emptyOrFill(model.ec,i{1},{''},'genes',true);
    end
    model.ec = emptyOrFill(model.ec,'concs',{'NaN'},'genes',true);
    model.ec = emptyOrFill(model.ec,'kcat',{'0'},'genes',true);
    % Change string to double
    for i={'kcat','mw','concs'}
        if isfield(model.ec,i{1})
            model.ec.(i{1}) = str2double(model.ec.(i{1}));
        end
    end
    % Fill rxnEnzMat
    rxnIdx              = cellfun('isempty', enzStoich(:,1));
    enzStoich(rxnIdx,:) = '';
    rxnIdx              = cell2mat(enzStoich(:,1));
    [~,enzIdx]          = ismember(enzStoich(:,2),model.ec.enzymes);
    coeffs              = cell2mat(enzStoich(:,3));
    model.ec.rxnEnzMat  = zeros(numel(model.ec.rxns), numel(model.ec.genes));
    linearIndices       = sub2ind([numel(model.ec.rxns), numel(model.ec.genes)], rxnIdx, enzIdx);
    model.ec.rxnEnzMat(linearIndices) = coeffs;
    %Parse ec-codes
    if ~isempty(ecGecko)
        locs = cell2mat(ecGecko(:,1));
        for i=unique(locs)'
            ecGeckoCat=strjoin(ecGecko(locs==i,2),';');
            model.ec.eccodes{i,1}=ecGeckoCat;
        end
        emptyEc=cellfun('isempty',model.ec.eccodes);
        model.ec.eccodes(emptyEc)={''};
    end
end

if verbose
    fprintf(' Done!\n');
end
end

function model = emptyOrFill(model,field,emptyEntry,type,keepEmpty)
if nargin<5
    keepEmpty=false;
end
if isnumeric(emptyEntry)
    emptyCells=isempty(model.(field));
else
    emptyCells=cellfun('isempty',model.(field));
end
if all(emptyCells) && ~keepEmpty
    model = rmfield(model, field);
elseif numel(model.(field))<numel(model.(type))
    model.(field)(end+1:numel(model.(type)),1)=emptyEntry;
end
end

function model = readFieldValue(model, fieldName, value, pos)
if numel(model.(fieldName))<pos-1
    model.(fieldName)(end+1:pos,1) = {''};
end
model.(fieldName)(pos,1) = {value};
end

function [miriams,miriamKey,entryNumber] = gatherAnnotation(pos,miriams,key,value,miriamKey,entryNumber)
if isempty(key)
    key=miriamKey;
else
    miriamKey=key;
end
if ~isempty(value)
    miriams(entryNumber,1:3) = {pos, key, strip(value)};
    entryNumber = entryNumber + 1;
end
end
