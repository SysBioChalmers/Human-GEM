function model=importYaml(yamlFilename, silentMode)
% importYaml
%   Imports a yaml file matching (roughly) the cobrapy yaml structure
%
%   Input:
%   yamlFile    a file in yaml model structure. As defined in Human-GEM, the
%               yaml file contains 5 sections: metaData, metabolites,
%               reactions, genes and compartments
%   silentMode  set as true to turn off notificaiton messages (opt, default
%               false)
%
%   Output:
%   model       a model structure
%
%   Usage: model = importYaml(yamlFilename, silentMode)
%
% This function is to reverse engineer the RAVEN function `writeYaml`
%

if nargin < 2
    silentMode = false;
end

if ~(exist(yamlFilename,'file')==2)
    error('Yaml file %s cannot be found', string(yamlFilename));
end

if verLessThan('matlab','9.9') %readlines introduced 2020b
    fid=fopen(yamlFilename);
    line_raw=cell(1000000,1);
    while ~feof(fid)
        line_raw{i}=fgetl(fid);
        i=i+1;
    end
    line_raw(i:end)=[];
    line_raw=string(line_raw);
else
    line_raw=readlines(yamlFile');
end

line_key=regexprep(line_raw,'^ *-? ([^:]+)(:).*','$1');
line_key=regexprep(line_key,'(.*!!omap)|(---)','');

line_value = regexprep(line_raw, '[^":]+: "?(.+)"?$','$1');
line_value = regexprep(line_value, '"','');

% Define the required fields of humanGEM
% There are a total of 37 fields in the model so far, the non-generic ones
% are excluded here
model=[];
model.id=[];
model.description=[];
model.rxns={};
model.mets={};
model.S=[];
model.lb={};
model.ub={};
model.rev=[];
model.c=[];
model.b=[];
model.comps={};
model.compNames={};
model.rxnNames={};
model.grRules={};
model.rxnGeneMat=[];
model.subSystems={};
model.eccodes={};
model.rxnNotes={};
model.genes={};
model.metNames={};
model.metComps={};
model.inchis={};
model.metFormulas={};
%model.unconstrained=[]; %abandoned
model.rxnReferences={};
model.rxnFrom={};
model.metFrom={};
model.rxnConfidenceScores={};
model.metCharges={};
model.version='';
model.annotation=[];
equations={};
leftEqns={};
rightEqns={};
objRxns={};

section = 0;
for i=1:numel(line_key)
    tline_raw = line_raw{i};
    tline_key = line_key{i};
    tline_value = line_value{i};
    % import different sections
    switch tline_raw
        case '- metaData:'
            section = 1;
            if ~silentMode
                fprintf('\t%d\n', section);
            end
            continue % Go to next line
        case '- metabolites:'
            section = 2;
            if ~silentMode
                fprintf('\t%d\n', section);
            end
            continue
        case '- reactions:'
            section = 3;
            readSubsystems = false;
            readEquation = false;
            rxnId = '';
            if ~silentMode
                fprintf('\t%d\n', section);
            end
            continue
        case '- genes:'
            section = 4;
            if ~silentMode
                fprintf('\t%d\n', section);
            end
            continue
        case '- compartments: !!omap'
            section = 5;
            if ~silentMode
                fprintf('\t%d\n', section);
            end
            continue
    end

    % skip over empty keys
    if isempty(tline_key)
        continue;
    end
    
    % import metaData
    if section == 1
        switch tline_key
            case 'short_name'
                model.id = tline_value;
            case 'full_name'
                model.description = tline_value;
            case 'version'
                model.version = tline_value;
            case 'taxonomy'
                model.annotation.taxonomy = tline_value;
            case 'description'
                model.annotation.note = tline_value;
            case 'github'
                model.annotation.sourceUrl = tline_value;
            case 'authors'
                model.annotation.authorList = tline_value;
            case 'email'
                model.annotation.email = tline_value;
            case 'organization'
                model.annotation.organization = tline_value;
        end
    end

    % import metabolites:
    if section == 2
        switch tline_key
            case 'id'
                model = readFieldValue(model, 'mets', tline_value);
            case 'name'
                model = readFieldValue(model, 'metNames', tline_value);
            case 'compartment'
                model = readFieldValue(model, 'metComps', tline_value);
            case 'formula'
                model = readFieldValue(model, 'metFormulas', tline_value);
            case 'charge'
                model = readFieldValue(model, 'metCharges', tline_value);
            case 'inchis'
                model = readFieldValue(model, 'inchis', tline_value);
            case 'metFrom'
                model = readFieldValue(model, 'metFrom', tline_value);
        end
    end

    % import reactions:
    if section == 3
        switch tline_key
            case 'id'
                model = readFieldValue(model, 'rxns', tline_value);
                rxnId = tline_value;
            case 'name'
                model = readFieldValue(model, 'rxnNames', tline_value);

            case 'lower_bound'
                model.lb(end+1,1)  = {tline_value};
                leftEqns(end+1,1)  = {leftEquation};
                rightEqns(end+1,1) = {rightEquation};
                readEquation = false;

            case 'upper_bound'
                model.ub(end+1,1)  = {tline_value};

            case 'gene_reaction_rule'
                model = readFieldValue(model, 'grRules', tline_value);

            case 'rxnNotes'
                model = readFieldValue(model, 'rxnNotes', tline_value);

            case 'rxnFrom'
                model = readFieldValue(model, 'rxnFrom', tline_value);

            case 'objective_coefficient'
                objRxns(end+1,1) = {rxnId};

            case 'eccodes'
                model = readFieldValue(model, 'eccodes', tline_value);

            case 'references'
                model = readFieldValue(model, 'rxnReferences', tline_value);

            case 'subsystem'
                readSubsystems = true;
                subSystems = {};

            case 'confidence_score'
                model = readFieldValue(model, 'rxnConfidenceScores', tline_value);
                model.subSystems(end+1,1) = {subSystems};
                readSubsystems = false;

            case 'metabolites'
                readEquation = true;
                leftEquation  = '';
                rightEquation = '';

            otherwise
                if readSubsystems
                    subSystems(end+1,1) = {regexprep(tline_value, '^ *- (.+)$','$1')};

                % resolve the equation
                elseif readEquation
                    metCoeffi = regexp(regexprep(tline_raw, ' +- ', ''), ': ', 'split');
                    coeffi = metCoeffi{2};
                    if str2double(coeffi) < 0
                        if strcmp(leftEquation, '')
                            leftEquation = [coeffi(2:end),' ',metCoeffi{1}]; %Remove minus sign from coefficient
                        else
                            leftEquation = [leftEquation,' + ',coeffi(2:end),' ',metCoeffi{1}];
                        end
                    else
                        if strcmp(rightEquation, '')
                            rightEquation = [' ',coeffi,' ',metCoeffi{1}];
                        else
                            rightEquation = [rightEquation,' + ',coeffi,' ',metCoeffi{1}];
                        end
                    end
                end
        end
    end

    % import genes:
    if section == 4
        model = readFieldValue(model, 'genes', tline_value);
    end

    % import compartments:
    if section == 5
        model.comps(end+1,1) = {tline_key};
        model.compNames(end+1,1) = {tline_value};
    end

end

% follow-up data processing
if ~silentMode
    fprintf('\nimporting completed\nfollow-up processing...');
end
[~, model.metComps] = ismember(model.metComps, model.comps);
model.metCharges = int64(str2double(model.metCharges));
model.lb = str2double(model.lb);
model.ub = str2double(model.ub);
model.annotation.defaultLB = min(model.lb);
model.annotation.defaultUB = max(model.ub);
model.rev = double(model.lb<0 & model.ub>0);
model.rxnConfidenceScores = str2double(model.rxnConfidenceScores);
model.b = zeros(length(model.mets),1);
model.c = double(ismember(model.rxns, objRxns));

[genes, rxnGeneMat] = getGenesFromGrRules(model.grRules);
if isequal(sort(genes), sort(model.genes))
    model.rxnGeneMat = rxnGeneMat;
    model.genes = genes;
else
    error('The gene list and grRules are inconsistent.');
end
% regenerate equations
equations = cell(length(model.rxns), 1);
revInd = find(model.rev);
irrevInd = setdiff(transpose([1: length(model.rxns)]), revInd);
equations(revInd)   = strcat(leftEqns(revInd), ' <=>', rightEqns(revInd));
equations(irrevInd) = strcat(leftEqns(irrevInd), ' =>', rightEqns(irrevInd));

% regenerate S matrix
[S, newMets, ~, ~] = constructS(equations, model.mets, model.rxns);
[~, metIdx] = ismember(model.mets, newMets);
model.S = S(metIdx, :);

if ~silentMode
    fprintf(' Done!\n');
end
end

function model = readFieldValue(model, fieldName, value)
    model.(fieldName)(end+1,1) = {value};
end
