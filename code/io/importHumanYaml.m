function model=importHumanYaml(yamlFilename)
% importHumanYaml
%   Imports a yaml file matching (roughly) the cobrapy yaml structure
%
%   Input:
%   yamlFile    a file in yaml model structure. As defined in HumanGEM, the
%               yaml file contains 5 sections: metaData, metabolites,
%               reactions, genes and compartments
%
%   Output:
%   model       a model structure
%
%   Usage: model=importYaml(yamlFilename)
%
%   Hao Wang, 2019-07-26
%
% This function is to reverse engineer the RAVEN function `writeYaml`
%

if ~(exist(yamlFilename,'file')==2)
    error('Yaml file %s cannot be found',string(yamlFilename));
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
model.unconstrained=[];
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
fprintf('Start importing...\n');

section = 0;
while ~feof(fid)
    tline = fgetl(fid);
    % import different sections
    
    % import metaData
    if isequal(tline, '- metaData:')
        fprintf('\tmetaData\n');
        section = 1;
    end

    if section == 1 && numel(tline) > 17
        switch tline(1:17)
            case '    short_name  :'
                model.id = strip(tline(19:end),'"');

            case '    full_name   :'
                model.description = strip(tline(19:end),'"');

            case '    version     :'
                model.version = strip(tline(19:end),'"');

            case '    taxonomy    :'
                model.annotation.taxonomy = strip(tline(19:end),'"');

            case '    description :'
                model.annotation.note = strip(tline(19:end),'"');

            case '    github      :'
                model.annotation.sourceUrl = strip(tline(19:end),'"');

            case '    authors     :'
                model.annotation.authorList = strip(tline(19:end),'"');

            case '    email       :'
                model.annotation.email = strip(tline(19:end),'"');

            case '    organization:'
                model.annotation.organization = strip(tline(19:end),'"');
        end 
    end


    % import metabolites:
    if isequal(tline, '- metabolites:')
        fprintf('\tmetabolites\n');
        section = 2;
    end

    if section == 2
        if startsWith(tline,'    - id: ')
            model = readFieldElement(model, tline, 'mets', '    - id: ');

        elseif startsWith(tline,'    - name: ')
            model = readFieldElement(model, tline, 'metNames', '    - name: ');

        elseif startsWith(tline,'    - compartment: ')
            model = readFieldElement(model, tline, 'metComps', '    - compartment: ');

        elseif startsWith(tline,'    - formula: ')
            model = readFieldElement(model, tline, 'metFormulas','    - formula: ');

        elseif startsWith(tline,'    - charge: ')
            model = readFieldElement(model, tline, 'metCharges','    - charge: ');
 
        elseif startsWith(tline,'    - inchis: ')
            model = readFieldElement(model, tline, 'inchis','    - inchis: ');

        elseif startsWith(tline,'    - metFrom: ')
            model = readFieldElement(model, tline, 'metFrom','    - metFrom: ');

        end
    end
    
    
    % import reactions:
    if isequal(tline, '- reactions:')
        fprintf('\treactions\n');
        section = 3;
        readSubsystems = false;
        readEquation = false;
        rxnId = '';
    end

    if section == 3
        if startsWith(tline,'    - id: ')
            model = readFieldElement(model, tline, 'rxns','    - id: ');
            rxnId = tline(12:end-1);

        elseif startsWith(tline,'    - name: ')
            model = readFieldElement(model, tline, 'rxnNames','    - name: ');
           
        elseif startsWith(tline,'    - lower_bound: ')
            model.lb = [model.lb; tline(20:end)];
            leftEqns  = [leftEqns; leftEquation];
            rightEqns = [rightEqns; rightEquation];
            readEquation = false;
            
        elseif startsWith(tline,'    - upper_bound: ')
            model.ub = [model.ub; tline(20:end)];
            
        elseif startsWith(tline,'    - gene_reaction_rule: ')
            model = readFieldElement(model, tline, 'grRules','    - gene_reaction_rule: ');
            
        elseif startsWith(tline,'    - rxnFrom: ')
            model = readFieldElement(model, tline, 'rxnFrom','    - rxnFrom: ');
            
        elseif startsWith(tline,'    - objective_coefficient: ')
            objRxns = [objRxns; rxnId];

        elseif startsWith(tline,'    - eccodes: ')
            model = readFieldElement(model, tline, 'eccodes','    - eccodes: ');
            
        elseif startsWith(tline,'    - references: ')
            model = readFieldElement(model, tline, 'rxnReferences','    - references: ');
                        
        elseif isequal(tline,'    - subsystem:')
            readSubsystems = true;
            subSystems = {};
            
        elseif startsWith(tline,'    - confidence_score: ')
            model = readFieldElement(model, tline, 'rxnConfidenceScores','    - confidence_score: ');
            model.subSystems = [model.subSystems; {subSystems}];
            readSubsystems = false;
        
        elseif isequal(tline,'    - metabolites: !!omap')
            readEquation = true;
            leftEquation  = '';
            rightEquation = '';
        else
            if readSubsystems
                subSystems = [subSystems; tline(12:end-1)];
                
            % resolve the equation
            elseif readEquation
                metCoeffi = split(tline(9:end), ': ');
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
    if isequal(tline, '- genes:')
        fprintf('\tgenes\n');
        section = 4;
    end
       
    if section == 4 && startsWith(tline,'    - id: ')
        model = readFieldElement(model, tline, 'genes','    - id: ');
    end


    % import compartments:
    if isequal(tline, '- compartments: !!omap')
        fprintf('\tcompartments\n');
        section = 5;
    end

    if section == 5 && numel(tline) > 7 && isequal(tline(1:6),'    - ')
        str = split(tline(7:end-1),': "');
        model.comps = [model.comps; str{1}];
        model.compNames = [model.compNames; str{2}];
    end
    
end
fclose(fid);


% follow-up data processing
fprintf('\nimporting completed\nfollow-up processing...');
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

% Although this works with HumanGEM, but it is NOT a generic solution of
% dealing with the `unconstrained` field for other models!
model.unconstrained = double(endsWith(model.mets, 'x'));

fprintf(' Done!\n');

end

function model = readFieldElement(model, lineStr, fieldName, keyStr)
%disp('reach here');    
keyStrLen = length(keyStr);
if isequal(lineStr, keyStr)
    model.(fieldName) = [model.(fieldName); {''}];
elseif length(lineStr) > keyStrLen && isequal(lineStr(1:keyStrLen), keyStr)
    newElement = strip(lineStr(keyStrLen+1:end),'"');
    model.(fieldName) = [model.(fieldName); {newElement}];
end

end
