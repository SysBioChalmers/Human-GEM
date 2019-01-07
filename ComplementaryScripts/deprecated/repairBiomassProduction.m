%
% FILE NAME:    curateBiomassReactions.m
% 
% DATE CREATED: 2018-11-06
%     MODIFIED: 2018-11-06
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script to curate biomass-related reactions in the model.
%


%% Run some initial scripts to get the latest version of the model

constrainVariableMassReactions
repairModelLeaks



%% Enable export of fatty acids

model = ihuman;

% load list of reactions to add
tmp = readtable('ComplementaryScripts/modelCuration/fatty_acid_export_rxns.txt','Delimiter','\t','HeaderLines',0);

% extract list of new metabolites that need to be added
% (note: these mets already exist in the model, but not in these comps)
newMets = split(tmp.eqn,' => ');
newMets = split(newMets(:,2),{'[',']'});
newMetNames = newMets(:,1);
newMetComps = newMets(:,2);

% get their existing met IDs
[~,ind] = ismember(newMetNames,model.metNames);
newMetIDs = strcat(regexprep(model.mets(ind),'[a-z]$',''),newMetComps);

% remove any that already exist
rem_ind = ismember(newMetIDs,model.mets);
newMetNames(rem_ind) = [];
newMetComps(rem_ind) = [];
newMetIDs(rem_ind) = [];

% construct structure of metabolite information to add
metsToAdd = {};
metsToAdd.mets = newMetIDs;
metsToAdd.metNames = newMetNames;
metsToAdd.compartments = newMetComps;
metsToAdd.unconstrained = double(ismember(newMetComps,'x'));

% add mets to model
model = addMets(model,metsToAdd,true);


% construct structure
rxnsToAdd = {};
rxnsToAdd.rxns = tmp.rxn;
rxnsToAdd.equations = tmp.eqn;
rxnsToAdd.ub = 1000*ones(size(tmp.rxn));
rxnsToAdd.lb = zeros(size(tmp.rxn));
rxnsToAdd.subSystems = repmat({'Transport reactions'},size(tmp.rxn));
prevNumRxns = length(model.rxns);  % get number of rxns before adding new

% this is a workaround to the addRxns function
grRules = model.grRules;
model.grRules(:) = {''};

% add reactions
model = addRxns(model,rxnsToAdd,3);

% restore grRules
model.grRules = repmat({''},length(model.rxns),1);
model.grRules(1:prevNumRxns) = grRules;
model.rxnGeneMat(length(prevNumRxns)+1:length(model.rxns),:) = 0;  % update rxnGeneMat as well, just in case


% update other non-standard fields
ind = prevNumRxns:length(model.rxns);
% model.rxnKEGGID(ind) = {''};
% model.rxnEHMNID(ind) = {''};
% model.rxnBiGGID(ind) = {''};
% model.rxnHepatoNET1ID(ind) = {''};
% model.rxnREACTOMEID(ind) = {''};
model.rxnRecon3DID(ind) = {''};
model.prRules(ind) = {''};
model.rxnProtMat(ind,:) = 0;
model.priorCombiningGrRules(ind) = {''};


% update ihuman variable
ihuman = model;


%% Analyze biomass production

model = simplifyModel(ihuman);
exch_inds = sum(model.S ~= 0) == 1;
model.S(:,exch_inds) = -abs(model.S(:,exch_inds));  % make all exch rxns the same direction (met --> 0)
model.ub(exch_inds) = 1000;   % specify export bound
model.lb(exch_inds) = 0;  % specify import bound

% allow import of specified mets
model.lb(ismember(model.rxns,{'HMR_9034'})) = -1000;  % glucose
model.lb(ismember(model.rxns,{'HMR_9048'})) = -1000;  % oxygen


% Set the new objective
% [~,obj_ind] = ismember('biomass_reaction',model.rxns);         % 0.50563 alanine[c] + 0.35926 arginine[c] + 0.27942 asparagine[c] + 0.35261 aspartate[c] + 20.7045 ATP[c] + 0.020401 cholesterol[c] + 0.011658 CL pool[c] + 0.039036 CTP[c] + 0.046571 cysteine[c] + 0.013183 dATP[n] + 0.009442 dCTP[n] + 0.009898 dGTP[n] + 0.013091 dTTP[n] + 0.27519 glucose-6-phosphate[c] + 0.38587 glutamate[c] + 0.326 glutamine[c] + 0.53889 glycine[c] + 0.036117 GTP[c] + 20.6508 H2O[c] + 0.12641 histidine[c] + 0.28608 isoleucine[c] + 0.54554 leucine[c] + 0.59211 lysine[c] + 0.15302 methionine[c] + 0.15446 PC-LD pool[c] + 0.055374 PE-LD pool[c] + 0.25947 phenylalanine[c] + 0.023315 PI pool[c] + 0.41248 proline[c] + 0.005829 PS-LD pool[c] + 0.39253 serine[c] + 0.017486 SM pool[c] + 0.31269 threonine[c] + 0.013306 tryptophan[c] + 0.15967 tyrosine[c] + 0.053446 UTP[c] + 0.35261 valine[c] + 0.002914 PG-CL pool[c] => 20.6508 ADP[c] + 20.6508 H+[c] + 20.6508 Pi[c]
% [~,obj_ind] = ismember('biomass_components',model.rxns);       % alanine[c] + arginine[c] + asparagine[c] + aspartate[c] + cholesterol[c] + cholesterol-ester pool[r] + CL pool[c] + cysteine[c] + DNA[n] + DNA-5-methylcytosine[n] + glutamate[c] + glutamine[c] + glycine[c] + histidine[c] + isoleucine[c] + leucine[c] + lipid droplet[c] + lysine[c] + methionine[c] + phenylalanine[c] + phosphatidate-LD-TAG pool[c] + PI pool[c] + proline[c] + RNA[c] + serine[c] + SM pool[c] + threonine[c] + tryptophan[c] + tyrosine[c] + valine[c] + glycogen[c] + cofactors and vitamins[c] => biomass[x]
% [~,obj_ind] = ismember('HMR_biomass_Renalcancer',model.rxns);  % 0.0031466 1,2-diacylglycerol-LD-TAG pool[c] + 0.0038601 1-acylglycerol-3P-LD-TG1 pool[c] + 0.0012642 2-lysolecithin pool[c] + 0.11871 alanine[c] + 0.0009274 arginine[c] + 0.012056 asparagine[c] + 0.69277 aspartate[c] + 0.0082708 cholesterol[c] + 0.0065713 cholesterol-ester pool[r] + 0.00057356 CL pool[c] + 0.0012056 cysteine[c] + 0.007581 DNA[n] + 0.16786 glutamate[c] + 0.19104 glutamine[c] + 0.13911 glycine[c] + 0.028749 histidine[c] + 0.0037096 isoleucine[c] + 0.011129 leucine[c] + 0.009274 lysine[c] + 0.0018548 methionine[c] + 0.017171 PC-LD pool[c] + 0.020051 PE-LD pool[c] + 0.0037096 phenylalanine[c] + 0.016414 phosphatidate-LD-TAG pool[c] + 0.0030972 PI pool[c] + 0.017992 proline[c] + 0.005891 PS-LD pool[c] + 0.030815 RNA[c] + 0.037096 serine[c] + 0.0024143 SM pool[c] + 0.027482 TAG-LD pool[c] + 0.020403 threonine[c] + 0.0009274 tryptophan[c] + 0.0055644 tyrosine[c] + 0.012056 valine[c] + 0.26551 glycogen[c] + 0.001 cofactors and vitamins[c] => biomass[c]
% [~,obj_ind] = ismember('cofactors_vitamins',model.rxns);         % [protein]-N6-(lipoyl)lysine[c] + biotin[c] + CoA[c] + cobamide-coenzyme[m] + cytochrome-C[m] + FADH2[c] + L-carnitine[c] + NADH[c] + NADPH[c] + riboflavin[c] + tetrahydrobiopterin[c] + THF[c] + ubiquinol[m] + 0.1 vitamin A derivatives[c] + 0.1 vitamin D derivatives[c] + 0.1 vitamin E derivatives[c] => cofactors and vitamins[c]
% [~,obj_ind] = ismember('HMR_3964',model.rxns);  % ATP hydrolysis
% [~,obj_ind] = ismember('HMR_6916',model.rxns);  % ATP synthase
% [~,obj_ind] = ismember('HMR_9034',model.rxns);  % glucose export
[~,obj_ind] = ismember('HMR_9058',model.rxns);  % CO2 export
% [~,obj_ind] = ismember('HMR_9048',model.rxns);  % O2 export
% [~,obj_ind] = ismember('HMR_9079',model.rxns);  % H+ export

model.c(:) = 0;
model.c(obj_ind) = 1;







% delete all sink and demand reactions
sink_dm = find(startsWith(ihuman.rxns,{'sink_','DM_'}));
model = removeReactionsFull(ihuman,sink_dm);

% simplify model
model = simplifyModel(model);

% Specify metUptake
% (these are the components of Ham's media)
mets = {'arginine[s]';'histidine[s]';'lysine[s]';'methionine[s]';'phenylalanine[s]';'tryptophan[s]';'tyrosine[s]';'alanine[s]';'glycine[s]';'serine[s]';'threonine[s]';'aspartate[s]';'glutamate[s]';'asparagine[s]';'glutamine[s]';'isoleucine[s]';'leucine[s]';'proline[s]';'valine[s]';'cysteine[s]';'thiamin[s]';'hypoxanthine[s]';'folate[s]';'biotin[s]';'pantothenate[s]';'choline[s]';'inositol[s]';'nicotinamide[s]';'pyridoxine[s]';'riboflavin[s]';'thymidine[s]';'aquacob(III)alamin[s]';'lipoic acid[s]';'glucose[s]';'sulfate[s]';'NEFA blood pool in[s]';'linoleate[s]';'linolenate[s]';'O2[s]';'H2O[s]';'retinoate[s]';'Fe2+[s]';'Pi[s]'};
mets = [mets; {'alpha-tocopherol[s]';'gamma-tocopherol[s]'}];  % add required vitamin E precursors (tocopherols)
model = specifyMetUptake(model,mets);

% also allow consumption of apocytochrome-C[c] and apoA1[c]
rxnsToAdd = {};
rxnsToAdd.equations = {' => apocytochrome-C[c]';' => apoA1[c]'};
rxnsToAdd.rxns = {'gen_apocytochrome_C';'gen_apoA1'};
rxnsToAdd.ub = [1;1];
rxnsToAdd.lb = [0;0];
model.grRules(:) = {''};
model = addRxns(model,rxnsToAdd,3);

% model.b(getIndexes(model,'cofactors and vitamins[c]','metscomps'),:) = [-1 0];
% model.b(getIndexes(model,'lipid droplet[c]','metscomps'),:) = [-1 0];


% Set the new objective
% [~,obj_ind] = ismember('biomass_reaction',model.rxns);         % 0.50563 alanine[c] + 0.35926 arginine[c] + 0.27942 asparagine[c] + 0.35261 aspartate[c] + 20.7045 ATP[c] + 0.020401 cholesterol[c] + 0.011658 CL pool[c] + 0.039036 CTP[c] + 0.046571 cysteine[c] + 0.013183 dATP[n] + 0.009442 dCTP[n] + 0.009898 dGTP[n] + 0.013091 dTTP[n] + 0.27519 glucose-6-phosphate[c] + 0.38587 glutamate[c] + 0.326 glutamine[c] + 0.53889 glycine[c] + 0.036117 GTP[c] + 20.6508 H2O[c] + 0.12641 histidine[c] + 0.28608 isoleucine[c] + 0.54554 leucine[c] + 0.59211 lysine[c] + 0.15302 methionine[c] + 0.15446 PC-LD pool[c] + 0.055374 PE-LD pool[c] + 0.25947 phenylalanine[c] + 0.023315 PI pool[c] + 0.41248 proline[c] + 0.005829 PS-LD pool[c] + 0.39253 serine[c] + 0.017486 SM pool[c] + 0.31269 threonine[c] + 0.013306 tryptophan[c] + 0.15967 tyrosine[c] + 0.053446 UTP[c] + 0.35261 valine[c] + 0.002914 PG-CL pool[c] => 20.6508 ADP[c] + 20.6508 H+[c] + 20.6508 Pi[c]
% [~,obj_ind] = ismember('biomass_components',model.rxns);       % alanine[c] + arginine[c] + asparagine[c] + aspartate[c] + cholesterol[c] + cholesterol-ester pool[r] + CL pool[c] + cysteine[c] + DNA[n] + DNA-5-methylcytosine[n] + glutamate[c] + glutamine[c] + glycine[c] + histidine[c] + isoleucine[c] + leucine[c] + lipid droplet[c] + lysine[c] + methionine[c] + phenylalanine[c] + phosphatidate-LD-TAG pool[c] + PI pool[c] + proline[c] + RNA[c] + serine[c] + SM pool[c] + threonine[c] + tryptophan[c] + tyrosine[c] + valine[c] + glycogen[c] + cofactors and vitamins[c] => biomass[x]
% [~,obj_ind] = ismember('HMR_biomass_Renalcancer',model.rxns);  % 0.0031466 1,2-diacylglycerol-LD-TAG pool[c] + 0.0038601 1-acylglycerol-3P-LD-TG1 pool[c] + 0.0012642 2-lysolecithin pool[c] + 0.11871 alanine[c] + 0.0009274 arginine[c] + 0.012056 asparagine[c] + 0.69277 aspartate[c] + 0.0082708 cholesterol[c] + 0.0065713 cholesterol-ester pool[r] + 0.00057356 CL pool[c] + 0.0012056 cysteine[c] + 0.007581 DNA[n] + 0.16786 glutamate[c] + 0.19104 glutamine[c] + 0.13911 glycine[c] + 0.028749 histidine[c] + 0.0037096 isoleucine[c] + 0.011129 leucine[c] + 0.009274 lysine[c] + 0.0018548 methionine[c] + 0.017171 PC-LD pool[c] + 0.020051 PE-LD pool[c] + 0.0037096 phenylalanine[c] + 0.016414 phosphatidate-LD-TAG pool[c] + 0.0030972 PI pool[c] + 0.017992 proline[c] + 0.005891 PS-LD pool[c] + 0.030815 RNA[c] + 0.037096 serine[c] + 0.0024143 SM pool[c] + 0.027482 TAG-LD pool[c] + 0.020403 threonine[c] + 0.0009274 tryptophan[c] + 0.0055644 tyrosine[c] + 0.012056 valine[c] + 0.26551 glycogen[c] + 0.001 cofactors and vitamins[c] => biomass[c]
% [~,obj_ind] = ismember('cofactors_vitamins',model.rxns);       % [protein]-N6-(lipoyl)lysine[c] + biotin[c] + CoA[c] + cobamide-coenzyme[m] + cytochrome-C[m] + FADH2[c] + L-carnitine[c] + NADH[c] + NADPH[c] + riboflavin[c] + tetrahydrobiopterin[c] + THF[c] + ubiquinol[m] + 0.1 vitamin A derivatives[c] + 0.1 vitamin D derivatives[c] + 0.1 vitamin E derivatives[c] => cofactors and vitamins[c]
% [~,obj_ind] = ismember('vitaminE',model.rxns);                 % 3-carboxy-alpha-chromanol[c] + gamma-CEHC-glucuronide[c] => vitamin E derivatives[c]

% [~,obj_ind] = ismember('HMR_3964',model.rxns);  % ATP hydrolysis
% [~,obj_ind] = ismember('HMR_6916',model.rxns);  % ATP synthase
% [~,obj_ind] = ismember('HMR_9034',model.rxns);  % glucose export
% [~,obj_ind] = ismember('HMR_9058',model.rxns);  % CO2 export
% [~,obj_ind] = ismember('HMR_9048',model.rxns);  % O2 export
% [~,obj_ind] = ismember('HMR_9079',model.rxns);  % H+ export
[~,obj_ind] = ismember('HMR_9036',model.rxns);  % linolenate export

model.c(:) = 0;
model.c(obj_ind) = 1;



% determine if removal of any metabolites makes the biomass reaction work
met_inds = find(model.S(:,obj_ind) ~= 0);
for i = 1:length(met_inds)
    m = model;
    m.S(met_inds(i),obj_ind) = 0;
    sol = solveLP(m);
    fprintf('met: %s, obj = %f\n',model.metNames{met_inds(i)},sol.f);
end


% determine which metabolites can be synthesized
met_inds = find(model.S(:,obj_ind) ~= 0);
rxnsToAdd = {};
for i = 1:length(met_inds)
    rxnsToAdd.equations = {[model.metNames{met_inds(i)},'[',model.comps{model.metComps(met_inds(i))},'] => ']};
    rxnsToAdd.rxns = {'fakerxn'};
    rxnsToAdd.ub = 1000;
    rxnsToAdd.lb = 0;
    m = addRxns(model,rxnsToAdd,3);
    m.c(:) = 0;
    m.c(end) = 1;
    
    m.b = [m.b,m.b];
    
    sol = solveLP(m);
    fprintf('met: %s, obj = %f\n',model.metNames{met_inds(i)},sol.f);
end
    


%% Check metabolic tasks


% delete all sink and demand reactions
sink_dm = find(startsWith(ihuman.rxns,{'sink_','DM_'}));
model = removeReactionsFull(ihuman,sink_dm);

% run checkTasks using large task list
% [taskReport, essentialReactions, taskStructure]=checkTasks(model, ...
%     '/Users/jonrob/Documents/PostDoc/HMR3/CurationFiles/metabolic_tasks/iHuman_simulations.xls');
[taskReport, essentialReactions, taskStructure]=checkTasks(model, ...
    '/Users/jonrob/Documents/PostDoc/HMR3/CurationFiles/metabolic_tasks/INIT_tasks.xls');


% investigating individual tasks
inMetInd = getIndexes(model,{'O2[s]';'glucose[s]'},'metscomps');
outMetInd = getIndexes(model,{'CO2[s]';'H2O[s]';'linolenate[c]'},'metscomps');
reqMetInd = getIndexes(model,{'linolenate[c]'},'metscomps');

m = model;
m.c(:) = 0;  % don't need an objective here
m.b = zeros(length(model.mets),2);
m.b(inMetInd,1) = -1000;
m.b(outMetInd,2) = 1000;
m.b(reqMetInd,1) = 1;

sol = solveLP(m)



