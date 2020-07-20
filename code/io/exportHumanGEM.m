function out=exportHumanGEM(ihuman,prefix,path,formats,masterFlag,dependencies)
% exportHumanGEM
%   Generates a directory structure and populates this with model files, ready
%   to be commited to a Git(Hub) maintained model repository. Writes the model
%   as SBML L3V1 FBCv2 (both XML and YAML), COBRA text, Matlab MAT-file
%   orthologies in KEGG
%
% Input:
%   ihuman              humanGEM model structure in RAVEN format
%   prefix              prefix for all filenames (opt, default 'model')
%   path                path where the directory structure should be generated
%                       and populated with all files (opt, default to current
%                       working directory)
%   formats             cell array of strings specifying in what file formats
%                       the model should be exported (opt, default to all
%                       formats as {'mat', 'txt', 'xlsx', 'xml', 'yml'})
%   masterFlag          logical, if true, function will error if RAVEN (and
%                       COBRA if detected) is/are not on the master branch.
%                       (opt, default false)
%   dependencies        logical, if false, will not output the dependency
%                       information (dependencies.txt). (opt, default true)
%
% Usage: exportHumanGEM(ihuman,prefix,path,formats,masterFlag,dependencies)
%


if nargin<6
    dependencies=true;
end
if nargin<5
    masterFlag=false;
end
if nargin<4
    formats={'mat', 'txt', 'xlsx', 'xml', 'yml'};
end
if ischar(formats)
    formats={formats};
end
if any(~ismember(formats, {'mat', 'txt', 'xlsx', 'xml', 'yml'}))
    EM='Unknown file format defined. Only mat, txt, xlsx, xml and yml are allowed file formats.';
    error(EM)
end
if nargin<3
    path='.';
end
if nargin<2
    prefix='model';
end

%Get versions or commits of toolboxes:
RAVENver = getToolboxVersion('RAVEN','ravenCobraWrapper.m',masterFlag);
COBRAver = getToolboxVersion('COBRA','initCobraToolbox.m',masterFlag);

%Retrieve libSBML version:
try % 5.17.0 and newer
    libSBMLver=OutputSBML;
    libSBMLver=libSBMLver.libSBML_version_string;
catch % before 5.17.0
    fid = fopen('tempModelForLibSBMLversion.xml','w+');
    fclose(fid);
    evalc('[~,~,libSBMLver]=TranslateSBML(''tempModelForLibSBMLversion.xml'',0,0)');
    libSBMLver=libSBMLver.libSBML_version_string;
    delete('tempModelForLibSBMLversion.xml');
end

% Make model folder, no warnings if folder already exists
[~,~,~]=mkdir(fullfile(path,'model'));
for i = 1:length(formats)
    [~,~,~]=mkdir(fullfile(path,'model',formats{i}));
end

% Write TXT format
if ismember('txt', formats)
    fid=fopen(fullfile(path,'model',strcat(prefix,'.txt')),'w');
    eqns=constructEquations(ihuman,ihuman.rxns,false,false,false,true);
    eqns=strrep(eqns,' => ','  -> ');
    eqns=strrep(eqns,' <=> ','  <=> ');
    eqns=regexprep(eqns,'> $','>');
    grRules=regexprep(ihuman.grRules,'\((?!\()','( ');
    grRules=regexprep(grRules,'(?<!\))\)',' )');
    fprintf(fid, 'Rxn name\tFormula\tGene-reaction association\tLB\tUB\tObjective\n');
    for i = 1:numel(ihuman.rxns)
        fprintf(fid, '%s\t', ihuman.rxns{i});
        fprintf(fid, '%s \t', eqns{i});
        fprintf(fid, '%s\t', grRules{i});
        fprintf(fid, '%6.2f\t%6.2f\t%6.2f\n', ihuman.lb(i), ihuman.ub(i), ihuman.c(i));
    end
    fclose(fid);
end

% Write YML format
if ismember('yml', formats)
    writeHumanYaml(ihuman,fullfile(path,'model',strcat(prefix,'.yml')));
end

% Write MAT format
if ismember('mat', formats)
    save(fullfile(path,'model',strcat(prefix,'.mat')),'ihuman');
end

% Write XLSX format
if ismember('xlsx', formats)
    exportToExcelFormat(ihuman,fullfile(path,'model',strcat(prefix,'.xlsx')));
end

% Write XML format
if ismember('xml', formats)
    model = simplifyModel(ihuman,false,false,true);  % remove inactivated rxns
    model = annotateModel(model);  % add annotation data to structure
    model = rmfield(model,'inchis');  % temporarily remove inchis until export function is updated
    model.id = regexprep(model.id,'-','');  % remove dash from model ID since it causes problems with SBML I/O
    exportModel(model,fullfile(path,'model',strcat(prefix,'.xml')));
end

%Save file with versions:
if dependencies
    fid = fopen(fullfile(path,'model','dependencies.txt'),'wt');
    fprintf(fid,['MATLAB\t' version '\n']);
    fprintf(fid,['libSBML\t' libSBMLver '\n']);
    fprintf(fid,['RAVEN_toolbox\t' RAVENver '\n']);
    if ~isempty(COBRAver)
        fprintf(fid,['COBRA_toolbox\t' COBRAver '\n']);
    end
    %if isfield(ihuman,'modelVersion')
    %    fields = fieldnames(ihuman.modelVersion);
    %    for i = 1:length(fields)
    %        value = ihuman.modelVersion.(fields{i});
    %        fprintf(fid,[fields{i} '\t' num2str(value) '\n']);
    %    end
    %end
    fclose(fid);
end

end
