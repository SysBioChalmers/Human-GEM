%we are using depmap version 2021_Q3

%Prepares the DepMap data, specifically by filtering out RNA-Seq samples for which no CRISPR data exists
%cd C:/Work/MatlabCode/components/human-GEM/Human-GEMDepMapEval/Human-GEM/code/DepMapGeneEss %to be used if not running the whole file, may need to change this
cd(fileparts(which(mfilename))); %to be used if running the whole file

%% Load and prepare DepMap RNA-Seq data (cell lines)

% load RNA-Seq data from txt file
rna_data = readtable('data/DepMap_tpm_gene_symbols.txt');

% load gene essentiality data (Achilles gene effect)
ach_data = readtable('data/Achilles_gene_effect.csv');
samples = ach_data.DepMap_ID;  % extract sample IDs

% filter RNA-Seq data to only include samples for which we have
% essentiality data
cellLineNames = rna_data.Properties.VariableNames;
%now replace '_' with '-'
cellLineNames = strrep(cellLineNames, '_', '-');

%{'Original column heading: 'ACH-001113''}

keep = ismember(cellLineNames, samples);

sum(keep) %891
sum(keep)/length(keep) % 65%, seems reasonable

% add RNA-Seq data to arrayData
arrayDataDepMap.genes = rna_data.gene;
arrayDataDepMap.tissues = cellLineNames(keep)';
arrayDataDepMap.levels = table2array(rna_data(:, keep));
arrayDataDepMap.threshold = 1;


% save tINIT inputs
save('data/arrayDataDepMap.mat','arrayDataDepMap');

%Generate ftINIT prepData - only needs to be done once. Can take up to an hour to run
model = importYaml('../../model/Human-GEM.yml');
[model.grRules, skipped] = simplifyGrRules(model.grRules, true);%takes a few minutes to run
prepData = prepHumanModelForftINIT(model, true, '../../data/metabolicTasks/metabolicTasks_Essential.txt', '../../model/reactions.tsv');
save('data/prepDataGeneSymbols.mat', 'prepData')
