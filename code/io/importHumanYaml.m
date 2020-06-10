function model=importHumanYaml(yamlFilename, silentMode)
% importHumanYaml
%   Imports a yaml file matching (roughly) the cobrapy yaml structure
%
%   Input:
%   yamlFile    a file in yaml model structure. As defined in HumanGEM, the
%               yaml file contains 5 sections: metaData, metabolites,
%               reactions, genes and compartments
%   silentMode  set as true to turn off notificaiton messages (opt, default
%               false)
%
%   Output:
%   model       a model structure
%
%   Usage: model=importYaml(yamlFilename, silentMode)
%
% This function is to reverse engineer the RAVEN function `writeYaml`
%

if nargin < 2
    silentMode = false;
end

if ~(exist(yamlFilename,'file')==2)
    error('Yaml file %s cannot be found', string(yamlFilename));
end

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
%model.rxnNotes={}; %not sure
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


% Load Yaml format model

fid = fopen(yamlFilename);
if ~silentMode
    fprintf('Start importing...\n');
end

section = 0;
while ~feof(fid)
    tline = fgetl(fid);
    
    % import different sections
    change_to_section = 0;
    switch tline
        case '- metaData:'
            change_to_section = 1;
        case '- metabolites:'
            change_to_section = 2;
        case '- reactions:'
            change_to_section = 3;
            readSubsystems = false;
            readEquation = false;
            rxnId = '';
        case '- genes:'
            change_to_section = 4;
        case '- compartments: !!omap'
            change_to_section = 5;
    end
    if logical(change_to_section)
        section = change_to_section;
        tline = fgetl(fid);
        if ~silentMode
            fprintf('\t%d\n', section);
        end
    end

    % skip over lines containing only omap
    if any(regexp(tline, "- !!omap"))
        tline = fgetl(fid);
    end
    
    % import metaData
    if section == 1
        [tline_key, tline_value] = tokenizeYamlLine(tline);
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
        [tline_key, tline_value] = tokenizeYamlLine(tline);
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
        [tline_key, tline_value] = tokenizeYamlLine(tline);
        switch tline_key
            case 'id'
                model = readFieldValue(model, 'rxns', tline_value);
                rxnId = tline_value;
            case 'name'
                model = readFieldValue(model, 'rxnNames', tline_value);

            case 'lower_bound'
                model.lb = [model.lb; tline_value];
                leftEqns  = [leftEqns; leftEquation];
                rightEqns = [rightEqns; rightEquation];
                readEquation = false;

            case 'upper_bound'
                model.ub = [model.ub; tline_value];

            case 'gene_reaction_rule'
                model = readFieldValue(model, 'grRules', tline_value);

            case 'rxnFrom'
                model = readFieldValue(model, 'rxnFrom', tline_value);

            case 'objective_coefficient'
                objRxns = [objRxns; rxnId];

            case 'eccodes'
                model = readFieldValue(model, 'eccodes', tline_value);

            case 'references'
                model = readFieldValue(model, 'rxnReferences', tline_value);

            case 'subsystem'
                readSubsystems = true;
                subSystems = {};

            case 'confidence_score'
                model = readFieldValue(model, 'rxnConfidenceScores', tline_value);
                model.subSystems = [model.subSystems; {subSystems}];
                readSubsystems = false;

            case 'metabolites'
                readEquation = true;
                leftEquation  = {''};
                rightEquation = {''};

            otherwise
                if readSubsystems
                    subSystems = [subSystems; regexprep(tline_key, '"', '')];
                    
                % resolve the equation
                elseif readEquation
                    metCoeffi = regexp(regexprep(tline, ' +- ', ''), ': ', 'split');
                    coeffi = str2num(metCoeffi{2});
                    if coeffi < 0
                        if strcmp(leftEquation, '')
                            leftEquation = strcat(num2str(abs(coeffi), 12),32,metCoeffi{1});
                        else
                            leftEquation = strcat(leftEquation,' +',32,num2str(abs(coeffi), 12),32,metCoeffi{1});
                        end
                    else
                        if strcmp(rightEquation, '')
                            rightEquation = strcat(32,num2str(coeffi, 12),32,metCoeffi{1});
                        else
                            rightEquation = strcat(rightEquation,' +',32,num2str(coeffi, 12),32,metCoeffi{1});
                        end
                    end
                end
        end
    end

    % import genes:
    if section == 4
        [tline_key, tline_value] = tokenizeYamlLine(tline);
        model = readFieldValue(model, 'genes', tline_value);
    end

    % import compartments:
    if section == 5
        [tline_key, tline_value] = tokenizeYamlLine(tline);
        model.comps = [model.comps; tline_key];
        model.compNames = [model.compNames; tline_value];
    end

end
fclose(fid);


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

% Now the `unconstrained` field is abandoned in HumanGEM
%model.unconstrained = double(endsWith(model.mets, 'x'));

if ~silentMode
    fprintf(' Done!\n');
end
end

function model = readFieldValue(model, fieldName, value)
    model.(fieldName) = [model.(fieldName); {value}];
end

function [line_key, line_value]= tokenizeYamlLine(line)
    line_key = regexp(line, '^ *-? ([^:]+)', 'tokens');
    line_key = char(line_key{1});
    line_value = regexp(line, '^ [^:]+: "?(.+)"?$', 'tokens');
    if isempty(line_value)
        line_value = '';
    else
        line_value = regexprep(line_value{1}, '"', '');
        line_value = char(line_value{1});
    end
end 