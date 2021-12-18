% script to change model compartment abbreviations
%   Extracellular:  s -> e
%   Boundary:       x -> b
%   Peroxisome:     p -> x

% specify model name
modelname = 'Human-GEM';

% load model and other data files
model = importYaml(fullfile('..', '..', 'model', [modelname '.yml']));
metAssoc = importTsvFile(fullfile('..', '..', 'model', 'metabolites.tsv'));
metAssoc_table = struct2table(metAssoc);
metsDep = importTsvFile(fullfile('..', '..', 'data', 'deprecatedIdentifiers', 'deprecatedMetabolites.tsv'));
metsDep_table = struct2table(metsDep);

% update comps field
model.comps(ismember(model.compNames, 'Extracellular')) = {'e'};
model.comps(ismember(model.compNames, 'Peroxisome')) = {'x'};

% update mets field
model.mets = regexprep(model.mets, 's$', 'e');
model.mets = regexprep(model.mets, 'p$', 'x');

% copy deprecated metabolites to deprecated metabolite table
dep_ind = endsWith(metAssoc.mets, {'s', 'p'});
metsDep_table = [metsDep_table; metAssoc_table(dep_ind, :)];

% update the metabolite compartments in the annotation file
metAssoc.mets = regexprep(metAssoc.mets, 's$', 'e');
metAssoc.mets = regexprep(metAssoc.mets, 'p$', 'x');

% export the new model and files
exportYaml(model, fullfile('..', '..', 'model', [modelname '.yml']));
exportTsvFile(metAssoc, fullfile('..', '..', 'model', 'metabolites.tsv'));
exportTsvFile(metsDep_table, fullfile('..', '..', 'data', 'deprecatedIdentifiers', 'deprecatedMetabolites.tsv'));


% other tasks to do in addition to running this script:
%   - replace all instances of [s], [x], [p], with [e], [b], [x],
%     respectively, in all metabolic task files
%   - replace the use of "x" with "b" in the addBoundaryMets function



