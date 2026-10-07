function model = generateAnimalGEM(species, varargin)
% generateAnimalGEM  Regenerate an animal GEM from Human-GEM.
%
% MATLAB counterpart of generateAnimalGEM.py, with the same steps, inputs,
% outputs and options. The species-specific inputs are read from the
% Animal-GEM repository:
%
%     <repo>/data/human2<Species>Orthologs.tsv   Alliance of Genome Resources orthologs
%     <repo>/data/<species>SpecificRxns.tsv       reactions Human-GEM does not have
%     <repo>/data/<species>SpecificMets.tsv       metabolites those reactions need
%     <repo>/model/                               model files and annotation tables (rewritten)
%     <repo>/version.txt                          the model version
%
% Steps, in order:
%
% 1. Rewrite every Human-GEM grRule through the ortholog table (Ensembl id
%    -> human symbol -> species gene) with getModelFromHomology. Reactions
%    whose grRule becomes empty are removed; gene-free reactions stay.
% 2. Add the species-specific metabolites and reactions.
% 3. Gap-fill the essential metabolic tasks from Human-GEM. Gap-filled
%    reactions carry no grRule and a note saying so.
% 4. Stamp the version and date into the model, merge the annotation
%    tables and write the model files.
%
% Parameters
% ----------
% species : char
%     'Mouse', 'Rat', 'Worm', 'Fruitfly' or 'Zebrafish'.
%
% Name-Value Arguments
% --------------------
% repo : char
%     the Animal-GEM repository (default: <Species>-GEM next to this
%     Human-GEM repository).
% humanRepo : char
%     the Human-GEM checkout to build on, normally a release tag (default:
%     this repository). A model without a version (develop) is refused
%     unless allowUnreleased is true.
% allowUnreleased : logical
%     build on a Human-GEM that has no version (default false).
% version : char
%     model version, also written to version.txt. Its major.minor must be
%     the Human-GEM release's; see animalVersion (default: the repository's
%     version if it matches, otherwise <major>.<minor>.0).
% date : char
%     model date, YYYY-MM-DD (default: today).
% formats : cell
%     output formats, any of 'yml', 'mat', 'xml' (default all three).
% keepBiomass : logical
%     keep Human-GEM's biomass reaction (MAR13082) as the objective instead
%     of the generic cell components (MAR00021) (default false).
%
% Returns
% -------
% model : struct
%     the generated model.
%
% Examples
% --------
%     generateAnimalGEM('Mouse', 'repo', '../Mouse-GEM', 'humanRepo', '../Human-GEM-release');
%
% Notes
% -----
% Requires RAVEN 3 with the getModelFromHomology options complexPolicy,
% keepGeneFree and preserveNotes. The gap-filling MILP needs a MILP solver
% (Gurobi or SCIP); set it with setRavenSolver.

taxonomy = struct('Mouse','10090','Rat','10116','Worm','6239','Fruitfly','7227','Zebrafish','7955');
species = char(species);
if ~isfield(taxonomy, species)
    error('Unknown species "%s"; use one of %s.', species, strjoin(fieldnames(taxonomy)', ', '));
end

thisRepo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
p = parseRAVENargs(varargin, {'repo', fullfile(fileparts(thisRepo), [species '-GEM']); ...
    'humanRepo', thisRepo; 'allowUnreleased', false; 'version', ''; 'date', ''; ...
    'formats', {'yml','mat','xml'}; 'keepBiomass', false});
repo = char(p.repo);
human = char(p.humanRepo);
formats = cellstr(p.formats);
modelId = [species '-GEM'];
dataDir = fullfile(repo, 'data');
modelDir = fullfile(repo, 'model');
orthologFile = fullfile(dataDir, ['human2' species 'Orthologs.tsv']);
rxnFile = fullfile(dataDir, [lower(species) 'SpecificRxns.tsv']);
metFile = fullfile(dataDir, [lower(species) 'SpecificMets.tsv']);
for f = {orthologFile, rxnFile, metFile}
    if ~isfile(f{1})
        error('File not found: %s', f{1});
    end
end
versionFile = fullfile(repo, 'version.txt');
current = '';
if isfile(versionFile)
    current = strtrim(fileread(versionFile));
end
modelDate = char(p.date);
if isempty(modelDate)
    modelDate = char(datetime('today', 'Format', 'yyyy-MM-dd'));
end

%% Template
fprintf('Reading template %s\n', fullfile(human, 'model', 'Human-GEM.yml'));
template = readYAMLmodel(fullfile(human, 'model', 'Human-GEM.yml'));
templateVersion = '';
if isfield(template, 'version') && ~isempty(template.version)
    templateVersion = char(template.version);
end
if isempty(templateVersion) && ~p.allowUnreleased
    error(['%s holds an unreleased Human-GEM (no version in the model). The animal GEMs are ' ...
        'built on a stable release: pass humanRepo a checkout of a release tag, or ' ...
        'allowUnreleased to build on this one anyway.'], human);
end
version = animalVersion(templateVersion, current, char(p.version));
if ~p.keepBiomass && ismember(templateVersion, {'2.0.0','2.0.1','2.1.0'})
    [template, changed] = fixLipoylBiomass(template);
    if changed
        fprintf(['Human-GEM %s: MAR00022 uses lipoic acid instead of ' ...
            '[protein]-N6-(lipoyl)lysine (see issue #1140)\n'], templateVersion);
    end
end

%% 1. Ortholog draft
orthologs = readAllianceOrthologs(orthologFile);
geneMap = ensemblToSpeciesGenes(fullfile(human, 'model', 'genes.tsv'), orthologs);
fprintf('%d ortholog pairs, %d Human-GEM genes with an ortholog\n', size(orthologs.pairs, 1), geneMap.Count);
model = buildOrthologDraft(template, geneMap, modelId);
fprintf('Ortholog draft: %d of %d reactions\n', numel(model.rxns), numel(template.rxns));

%% 2. Species-specific network
rxnTable = readTsv(rxnFile);
metTable = readTsv(metFile);
[model, added] = addSpeciesNetwork(model, rxnTable, metTable);
fprintf('Added %d species-specific reactions, %d metabolites\n', numel(added), height(metTable));
warnUnbalanced(model, added, templateVersion);

%% 3. Gap-filling
[model, filled] = fillEssentialTasks(model, template, ...
    parseTaskList(fullfile(human, 'data', 'metabolicTasks', 'metabolicTasks_Essential.txt')), ~p.keepBiomass);
fprintf('Gap-filled %d reactions\n', numel(filled));
model = rmfield(model, intersect({'rxnFrom','metFrom','geneFrom'}, fieldnames(model)));

%% 4. Metadata, annotation tables and output
symbols = [orthologs.pairs(:,2); speciesRuleGeneIds(rxnTable.grRules)];
names = [orthologs.symbols; speciesRuleGeneSymbols(rxnTable.grRules)];
[found, idx] = ismember(model.genes, symbols);
model.geneShortNames = model.genes;
model.geneShortNames(found) = names(idx(found));

model.id = modelId;
model.name = modelId;
model.version = version;
model.date = modelDate;
model.annotation.taxonomy = taxonomy.(species);
model.annotation.sourceUrl = ['https://github.com/SysBioChalmers/' modelId];

rxnAnnotation = mergeAnnotation(readTsv(fullfile(human, 'model', 'reactions.tsv')), rxnTable, 'rxns', model.rxns);
metAnnotation = mergeAnnotation(readTsv(fullfile(human, 'model', 'metabolites.tsv')), metTable, 'mets', model.mets);
writeOutputs(model, modelDir, modelId, species, rxnAnnotation, metAnnotation, formats);
if ~strcmp(current, version)
    fid = fopen(versionFile, 'w');
    fprintf(fid, '%s', version);
    fclose(fid);
end
fprintf('%s %s written to %s\n', modelId, version, modelDir);
end


%% Version ------------------------------------------------------------------

function version = animalVersion(templateVersion, current, requested)
% The animal GEM version: its major.minor follows the Human-GEM release it
% is built on, and the patch number is its own. Same rule as animal_version
% in generateAnimalGEM.py.
if isempty(templateVersion)
    if ~isempty(requested)
        version = requested;
    elseif ~isempty(current)
        version = current;
    else
        error('No version: pass version.');
    end
    return
end
base = majorMinor(templateVersion);
if ~isempty(requested)
    if ~strcmp(majorMinor(requested), base)
        error(['Version %s does not follow Human-GEM %s: the animal GEMs share ' ...
            'Human-GEM''s major.minor (%s.z).'], requested, templateVersion, base);
    end
    version = requested;
elseif ~isempty(current) && strcmp(majorMinor(current), base)
    version = current;
else
    version = [base '.0'];
end
end

function mm = majorMinor(v)
parts = strsplit(v, '.');
if numel(parts) ~= 3 || any(cellfun(@(x) isempty(regexp(x, '^\d+$', 'once')), parts))
    error('Not an x.y.z version: "%s"', v);
end
mm = [parts{1} '.' parts{2}];
end


%% Orthologs ----------------------------------------------------------------

function id = safeGeneId(symbol)
% Gene id for a gene symbol: characters a grRule cannot hold become "_";
% letters (including non-ASCII), digits, ".", ":" and "-" stay. A symbol
% that is a grRule operator gets a "_gene" suffix. Same rule as
% safe_gene_id in generateAnimalGEM.py.
keep = isstrprop(symbol, 'alphanum') | ismember(symbol, '_.:-');
id = symbol;
id(~keep) = char(1);
id = regexprep(id, [char(1) '+'], '_');
id = regexprep(id, '^_+|_+$', '');
if ismember(lower(id), {'and','or','not'})
    id = [id '_gene'];
end
end

function ids = speciesRuleGeneIds(rules)
ids = cellfun(@safeGeneId, speciesRuleGeneSymbols(rules), 'UniformOutput', false);
end

function symbols = speciesRuleGeneSymbols(rules)
symbols = {};
for i = 1:numel(rules)
    tokens = regexp(rules{i}, '[^\s()]+', 'match');
    symbols = [symbols; tokens(~ismember(lower(tokens), {'and','or'}))']; %#ok<AGROW>
end
end

function rule = safeRule(rule)
% rule with each gene symbol replaced by its safeGeneId.
[tokens, separators] = regexp(rule, '[^\s()]+', 'match', 'split');
out = separators{1};
for i = 1:numel(tokens)
    if ~ismember(lower(tokens{i}), {'and','or'})
        tokens{i} = safeGeneId(tokens{i});
    end
    out = [out tokens{i} separators{i+1}]; %#ok<AGROW>
end
rule = out;
end

function orthologs = readAllianceOrthologs(file)
% Reduce an Alliance of Genome Resources ortholog table to human/species
% pairs, as read_alliance_orthologs in generateAnimalGEM.py:
% 1. drop pairs that are neither best forward nor best reverse;
% 2. keep every human gene with a single remaining hit;
% 3. for the others, keep the hits that are both best forward and reverse;
% 4. if none is, keep the hit supported by the most methods (first on a tie).
t = readTsv(file);
t = t(~(strcmp(t.best, 'No') & strcmp(t.bestReverse, 'No')), :);
methodCount = str2double(t.methodCount);
[groups, ~, g] = unique(t.fromGeneId, 'stable');
keep = false(height(t), 1);
for i = 1:numel(groups)
    rows = find(g == i);
    if numel(rows) == 1
        keep(rows) = true;
        continue
    end
    both = rows(strcmp(t.best(rows), 'Yes') & strcmp(t.bestReverse(rows), 'Yes'));
    if ~isempty(both)
        keep(both) = true;
    else
        [~, k] = max(methodCount(rows));
        keep(rows(k)) = true;
    end
end
kept = t(keep, :);
ids = cellfun(@safeGeneId, kept.toSymbol, 'UniformOutput', false);
[uIds, ~, j] = unique(ids);
for i = 1:numel(uIds)
    if numel(unique(kept.toSymbol(j == i))) > 1
        error('Species symbols that map to one gene id: %s', strjoin(unique(kept.toSymbol(j == i))', ', '));
    end
end
orthologs.pairs = [kept.fromSymbol, ids];
orthologs.symbols = kept.toSymbol;
end

function geneMap = ensemblToSpeciesGenes(genesFile, orthologs)
% Map each Human-GEM gene (Ensembl id) to its species orthologs, via the
% gene symbol.
bySymbol = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:size(orthologs.pairs, 1)
    h = orthologs.pairs{i,1};
    if isKey(bySymbol, h)
        bySymbol(h) = unique([bySymbol(h); orthologs.pairs(i,2)], 'stable');
    else
        bySymbol(h) = orthologs.pairs(i,2);
    end
end
genes = readTsv(genesFile);
geneMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:height(genes)
    found = {};
    for s = strtrim(strsplit(genes.geneSymbols{i}, ';'))
        if ~isempty(s{1}) && isKey(bySymbol, s{1})
            found = unique([found; bySymbol(s{1})], 'stable');
        end
    end
    if ~isempty(found)
        geneMap(genes.genes{i}) = found;
    end
end
end

function draft = buildOrthologDraft(template, geneMap, modelId)
% Draft model from template with every grRule rewritten through geneMap,
% using RAVEN's getModelFromHomology with a table of perfect hits.
if ~isfield(template, 'id') || isempty(template.id)
    template.id = 'HumanGEM';   % older releases carry no id
end
keys = geneMap.keys;
pairs = cell(0, 2);
for i = 1:numel(keys)
    v = geneMap(keys{i});
    pairs = [pairs; repmat(keys(i), numel(v), 1), v(:)]; %#ok<AGROW>
end
blast = makeFakeBlastStructure(pairs, template.id, modelId);
evalc(['draft = getModelFromHomology({template}, blast, modelId, ''complexPolicy'', ''keep'', ' ...
    '''keepGeneFree'', true, ''preserveNotes'', true);']);
% Bookkeeping of the template and of the homology call, not properties of
% the animal model
draft = rmfield(draft, intersect({'rxnFrom','metFrom','geneFrom'}, fieldnames(draft)));
end


%% Species-specific network -------------------------------------------------

function [model, added] = addSpeciesNetwork(model, rxns, mets)
% Add the species-specific metabolites and reactions; returns the new
% reaction ids. Equations name their metabolites as name[compartment].
needM = {'mets','metNames','metFormulas','metCharges','compartments'};
needR = {'rxns','equations','subSystems','grRules'};
if ~all(ismember(needM, mets.Properties.VariableNames)) || ~all(ismember(needR, rxns.Properties.VariableNames))
    error('A species-specific table is missing a required column.');
end
clash = [mets.mets(ismember(mets.mets, model.mets)); rxns.rxns(ismember(rxns.rxns, model.rxns))];
if ~isempty(clash)
    error('Already in the model, cannot be added: %s', strjoin(clash', ', '));
end
known = [strcat(model.metNames, '[', model.comps(model.metComps), ']'); ...
    strcat(mets.metNames, '[', mets.compartments, ']')];
missing = {};
for i = 1:height(rxns)
    terms = regexp(rxns.equations{i}, '(?:^|\s\+\s|=>|<=>|-->)\s*(?:\d+(?:\.\d+)?\s+)?(.+?)\[(\w+)\]', 'tokens');
    for k = 1:numel(terms)
        m = [strtrim(terms{k}{1}) '[' terms{k}{2} ']'];
        if ~ismember(m, known)
            missing{end+1} = [m ' in ' rxns.rxns{i}]; %#ok<AGROW>
        end
    end
end
if ~isempty(missing)
    error(['%d metabolite(s) in the species-specific reactions are neither in the model nor ' ...
        'in the species-specific metabolite table (renamed in Human-GEM?): %s'], ...
        numel(missing), strjoin(unique(missing), '; '));
end

metsToAdd.mets = mets.mets;
metsToAdd.metNames = mets.metNames;
metsToAdd.metFormulas = mets.metFormulas;
metsToAdd.metCharges = toNumber(mets.metCharges, 0);
metsToAdd.compartments = mets.compartments;
model = addMets(model, metsToAdd, 'copyInfo', false);

rxnsToAdd.rxns = rxns.rxns;
rxnsToAdd.equations = rxns.equations;
rxnsToAdd.grRules = cellfun(@safeRule, rxns.grRules, 'UniformOutput', false);
rxnsToAdd.subSystems = cellfun(@(s) {s}, rxns.subSystems, 'UniformOutput', false);
rxnsToAdd.subSystems(cellfun(@(s) isempty(s{1}), rxnsToAdd.subSystems)) = {{}};
optional = {'rxnNames','rxnNames'; 'eccodes','eccodes'; 'rxnReferences','rxnReferences'};
for k = 1:size(optional, 1)
    if ismember(optional{k,1}, rxns.Properties.VariableNames)
        rxnsToAdd.(optional{k,2}) = rxns.(optional{k,1});
    end
end
if ismember('lb', rxns.Properties.VariableNames)
    rxnsToAdd.lb = toNumber(rxns.lb, -1000);
end
if ismember('ub', rxns.Properties.VariableNames)
    rxnsToAdd.ub = toNumber(rxns.ub, 1000);
end
if ismember('rxnConfidenceScores', rxns.Properties.VariableNames)
    rxnsToAdd.rxnConfidenceScores = toNumber(rxns.rxnConfidenceScores, NaN);
end
model = addRxns(model, rxnsToAdd, 'eqnType', 3, 'allowNewGenes', true);
added = rxns.rxns;
end

function warnUnbalanced(model, added, templateVersion)
% Report species-specific reactions that are not mass or charge balanced.
% Sinks and demands (one metabolite) are unbalanced by definition.
idx = find(ismember(model.rxns, added));
idx = idx(sum(model.S(:, idx) ~= 0, 1)' > 1);
balance = getElementalBalance(model, 'rxns', idx);
bad = model.rxns(idx(balance.balanceStatus == 0 | balance.chargeStatus == 0));
if ~isempty(bad)
    warning('%d of %d species-specific reactions are not mass or charge balanced against Human-GEM %s: %s', ...
        numel(bad), numel(added), templateVersion, strjoin(bad', ', '));
end
end


%% Gap-filling --------------------------------------------------------------

function [template, changed] = fixLipoylBiomass(template)
% Let the cofactor pool of MAR00021 use lipoic acid instead of
% [protein]-N6-(lipoyl)lysine, which Human-GEM 2.0.0-2.1.0 cannot make (as
% the current human cofactor pool MAR10065 does). The curated fix belongs in
% Human-GEM (issue #1140).
changed = false;
r22 = find(strcmp(template.rxns, 'MAR00022'));
r65 = find(strcmp(template.rxns, 'MAR10065'));
if isempty(r22) || isempty(r65)
    return
end
lipoyl = find(template.S(:, r22) ~= 0 & strcmp(template.metNames, '[protein]-N6-(lipoyl)lysine'));
acid = find(template.S(:, r65) ~= 0 & strcmp(template.metNames, 'lipoic acid'));
if numel(lipoyl) ~= 1 || numel(acid) ~= 1
    return
end
template.S(acid, r22) = template.S(lipoyl, r22);
template.S(lipoyl, r22) = 0;
if ~isfield(template, 'rxnNotes')
    template.rxnNotes = repmat({''}, numel(template.rxns), 1);
end
note = regexprep(template.rxnNotes{r22}, ';$', '');
if ~isempty(note)
    note = [note ';'];
end
template.rxnNotes{r22} = [note '[protein]-N6-(lipoyl)lysine replaced by lipoic acid by the ' ...
    'animal GEM generator (Human-GEM issue 1140)'];
changed = true;
end

function model = resetBiomass(model)
% Block the human biomass reaction and make the generic cell components
% (MAR00021) the objective.
human = strcmp(model.rxns, 'MAR13082');
components = strcmp(model.rxns, 'MAR00021');
if ~any(human) || ~any(components)
    error('MAR13082 or MAR00021 is not in the model; is the template Human-GEM?');
end
model.lb(human) = 0;
model.ub(human) = 0;
model.c(:) = 0;
model.ub(components) = 1000;
model.c(components) = 1;
end

function [model, filled] = fillEssentialTasks(model, template, tasks, reset)
% Gap-fill model from template until the essential tasks pass. Added
% reactions keep no grRule and are marked in their notes. resolveTies pins
% each degenerate minimum-cost fill to the fewest, lowest-id reactions, so
% the result matches generateAnimalGEM.py and does not depend on the solver.
if numel(intersect(model.rxns, template.rxns)) < 0.5 * numel(model.rxns)
    error('The model shares under half of its reactions with the template.');
end
if reset
    model = resetBiomass(model);
    template = resetBiomass(template);
end
[~, addedRxns, failed] = fitTasks(closeModel(model), closeModel(template), [], ...
    'taskStructure', tasks, 'gapFillMode', 'preMerged', 'resolveTies', true, 'printOutput', false);
if any(failed)
    msg = sprintf('Gap-filling could not satisfy tasks: %s.', strjoin(unique({tasks(failed).id}), ', '));
    if reset
        msg = [msg ' The reference cannot run them with MAR00021 as the biomass reaction ' ...
            '(Human-GEM issue #1140); keepBiomass keeps the human biomass instead.'];
    end
    error(msg); %#ok<SPERR>
end
filled = template.rxns(any(addedRxns, 2));
filled = filled(~ismember(filled, model.rxns));
if ~isempty(filled)
    model = addRxnsGenesMets(model, template, filled, 'addGene', false);
    idx = find(ismember(model.rxns, filled));
    for k = idx'
        [~, t] = ismember(model.rxns{k}, template.rxns);
        note = '';
        if isfield(template, 'rxnNotes')
            note = regexprep(template.rxnNotes{t}, ';$', '');
        end
        if ~isempty(note)
            note = [note ';'];
        end
        model.rxnNotes{k} = [note 'reaction added by gap filling'];
        if isfield(template, 'rxnConfidenceScores') && isfield(model, 'rxnConfidenceScores')
            model.rxnConfidenceScores(k) = template.rxnConfidenceScores(t);
        end
        model.grRules{k} = '';
    end
    [model.grRules, model.rxnGeneMat] = standardizeGrRules(model, 'embedded', true);
    model = deleteUnusedGenes(model, 'verbose', 0);
end
report = checkTasks(closeModel(model), [], 'taskStructure', tasks, 'printOutput', false);
if ~all(report.ok)
    error('Tasks still fail after gap-filling: %s', strjoin(report.id(~report.ok)', ', '));
end
end


%% Annotation tables and output ---------------------------------------------

function merged = mergeAnnotation(humanTable, speciesTable, idCol, ids)
% Annotation table of the animal model: Human-GEM's columns, species rows
% appended (species-only columns are not carried over), rows in model order.
cols = humanTable.Properties.VariableNames;
extra = table();
for c = cols
    if ismember(c{1}, speciesTable.Properties.VariableNames)
        extra.(c{1}) = speciesTable.(c{1});
    else
        extra.(c{1}) = repmat({''}, height(speciesTable), 1);
    end
end
combined = [humanTable; extra];
[~, last] = unique(combined.(idCol), 'last');
combined = combined(sort(last), :);
[found, idx] = ismember(ids, combined.(idCol));
if ~all(found)
    error('Model components with no annotation row: %s', strjoin(ids(~found)', ', '));
end
merged = combined(idx, :);
end

function writeOutputs(model, modelDir, modelId, species, rxnAnnotation, metAnnotation, formats)
if ~isfolder(modelDir)
    mkdir(modelDir);
end
writetable(rxnAnnotation, fullfile(modelDir, 'reactions.tsv'), 'FileType', 'text', 'Delimiter', '\t', 'QuoteStrings', 'none');
writetable(metAnnotation, fullfile(modelDir, 'metabolites.tsv'), 'FileType', 'text', 'Delimiter', '\t', 'QuoteStrings', 'none');
if ismember('yml', formats)
    writeYAMLmodel(model, 'fileName', fullfile(modelDir, [modelId '.yml']));
end
if ismember('mat', formats)
    s.([lower(species) 'GEM']) = model; %#ok<STRNU>
    save(fullfile(modelDir, [modelId '.mat']), '-struct', 's');
end
if ismember('xml', formats)
    % SBML ids cannot contain a dash. Cross-references are merged into the
    % SBML copy only, as in Human-GEM.
    annotated = annotateGEM(model, modelDir, {'rxn','met'});
    annotated.id = strrep(modelId, '-', '');
    exportModel(annotated, 'fileName', fullfile(modelDir, [modelId '.xml']));
end
end


%% Helpers ------------------------------------------------------------------

function t = readTsv(file)
% A tab-separated table with every column as text, whitespace kept as in the
% file; empty cells are ''.
opts = detectImportOptions(file, 'FileType', 'text', 'Delimiter', '\t');
opts.VariableNamingRule = 'preserve';
opts = setvartype(opts, 'char');
opts = setvaropts(opts, 'FillValue', '', 'WhitespaceRule', 'preserve');
t = readtable(file, opts);
end

function v = toNumber(c, default)
v = str2double(c);
v(cellfun(@(x) isempty(strtrim(x)), c)) = default;
end
