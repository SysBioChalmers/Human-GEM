%
% FILE NAME:    curateBiomassReactions.m
% 
% DATE CREATED: 2018-11-06
%     MODIFIED: 2018-11-08
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script to curate biomass-related reactions in the model. The
%          changes are summarized as follows:
%
%          1. Modification of existing biomass reactions
%
%             - The existing biomass reactions in humanGEM were renamed as
%               follows:
%
%               Original Name
%               'biomass_reaction'            'biomass_Recon3D'
%               'biomass_maintenance'         'biomass_maintenance_Recon3D'
%               'biomass_maintenance_noTrTr'  'biomass_maintenance_noTrTr_Recon3D'
%               'HMR_biomass_Renalcancer'     'biomass_HMR_RenalCancer'
%
%             - In addition, each of these biomass reactions was modified
%               such that they produce the 'biomass[c]' metabolite.
%
%             - "biomass_transport" ([c] to [s]) and "biomass_exchange"
%               ([s] to [x]) reactions were added.
%
%          2. A new biomass reaction for HepG2 cells was added
%
%             - The new reaction is named 'biomass_HepG2', and is
%               modularized to facilitate interpretation/modification of
%               the biomass reaction (i.e., it uses pooling reactions to
%               consolidate protein, DNA, RNA, etc. rather than including
%               all individual metabolites in the biomass reaction).
%
%             - The additional reactions pre-pool metabolites used for the
%               new HepG2 biomass reaction were also added to the model.
%
%          3. Addition of fatty acid export reactions
%
%             - For the biomass_components reaction to carry flux, the 
%               model needs to be able to export some fatty acid pool
%               metabolites (LD-TG1, LD-TG2, LD-TG3, LD-PC, LD-PE, LD-PI,
%               LD-PS, and LD-SM). Therefore, reactions enabling their
%               transport from the cytoplasm to extracellular, and from
%               extracellular to boundary, were added to the model.
%


%% Run some initial scripts to get the latest version of the model

constrainVariableMassReactions
repairModelLeaks


%% Load model

% load current version of humanGEM
if ~exist('ihuman','var')
    load('humanGEM.mat');  % version X.X.X
end
ihuman_orig = ihuman;  % to keep track of changes made


%% Modify existing biomass reactions

% modify each biomass reaction to produce the "biomass" metabolite in the
% cytoplasm compartment (i.e., biomass[c]), and ensure that none are
% producing the boundary biomass metabolite ([x]). This is because some
% RAVEN functions can erroneously classify reactions with boundary mets as 
% exchange reactions, and can be problematic if the biomass reaction is
% classified as such.
biomass_x_ind = getIndexes(ihuman,'biomass[x]','metscomps');
biomass_c_ind = getIndexes(ihuman,'biomass[c]','metscomps');
biomass_rxns = {'biomass_components';'biomass_reaction';'biomass_maintenance';'biomass_maintenance_noTrTr';'HMR_biomass_Renalcancer'};
bm_rxn_ind = getIndexes(ihuman,biomass_rxns,'rxns');
ihuman.S(biomass_c_ind,bm_rxn_ind) = 1;
ihuman.S(biomass_x_ind,bm_rxn_ind) = 0;

% rename biomass reactions
new_rxn_names = {'biomass_components';'biomass_Recon3D';'biomass_maintenance_Recon3D';'biomass_maintenance_noTrTr_Recon3D';'biomass_HMR_RenalCancer'};
[~,ind] = ismember(biomass_rxns,ihuman.rxns);
ihuman.rxns(ind) = new_rxn_names;


%% Add new reactions related to biomass production
% The HepG2 biomass reaction, and it's associated pool rxns, will be added
% to the model. In addition, reactions enabling export of some fatty acids
% will be incorporated to enable the biomass_components reaction to carry
% flux.

% load new reaction information from file
fid = fopen('../../ComplementaryData/modelCuration/new_rxns_for_biomass.tsv');
rxnData = textscan(fid,'%s%s%f%f%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);

% verify that none of the reactions exist in the current model
if any(ismember(rxnData{1},ihuman.rxns))
    error('One or more reactions to be added already exist in the model.');
end


%...................... Add new metabolites to model ......................

% extract all metabolites that are involved in the new reactions
mets = parseRxnEqu(rxnData{6});
tmp = split(mets,{'[',']'});
metNames = tmp(:,1);
metComps = tmp(:,2);
[~,metCompInds] = ismember(metComps,ihuman.comps);

% determine which mets don't yet exist in the model
modelMetsWithComp = strcat(ihuman.metNames,'[',ihuman.comps(ihuman.metComps),']');
new_ind = ~ismember(mets,modelMetsWithComp);
new_mets = mets(new_ind);
new_metNames = metNames(new_ind);
new_metComps = metComps(new_ind);

% get metIDs for mets that exist in the model, but in another compartment
[hasmatch,ind] = ismember(new_metNames,ihuman.metNames);
new_metIDs = repmat({''},size(new_metNames));
new_metIDs(hasmatch) = strcat(regexprep(ihuman.mets(ind(hasmatch)),'[a-z]$',''),new_metComps(hasmatch));

% assign new metIDs to mets that don't yet exist anywhere in the model
%                 metName                metID
nameIDpairs = {'fatty acid pool',   'fattyAcidPool'
               'heparan sulfate',   'heparanSulfate'
               'metabolite pool',   'metabolitePool'
               'phosphatidate',     'phosphatidate'
               'phosphatidyl pool', 'phosphatidylPool'
               'protein pool',      'proteinPool'};
if any(ismember(nameIDpairs(:,1),ihuman.metNames)) || any(ismember(nameIDpairs(:,2),ihuman.mets))
    error('One or more met names and/or IDs to be added already exist in the model!');
end
[hasmatch,ind] = ismember(new_metNames,nameIDpairs(:,1));
new_metIDs(hasmatch) = strcat(nameIDpairs(ind(hasmatch),2),'_',new_metComps(hasmatch));

% construct structure of metabolite information to add
metsToAdd = {};
metsToAdd.mets = new_metIDs;
metsToAdd.metNames = new_metNames;
metsToAdd.compartments = new_metComps;
metsToAdd.unconstrained = double(ismember(new_metComps,'x'));

% add mets to model
ihuman = addMets(ihuman,metsToAdd,true);


%....................... Add new reactions to model .......................

% construct structure
rxnsToAdd = {};
rxnsToAdd.rxns = rxnData{1};
rxnsToAdd.equations = rxnData{6};
rxnsToAdd.lb = rxnData{3};
rxnsToAdd.ub = rxnData{4};
rxnsToAdd.subSystems = rxnData{5};

% add reactions
prevNumRxns = length(ihuman.rxns);  % record original number of rxns
ihuman = addRxns(ihuman,rxnsToAdd,3);

% update other non-standard fields
ind = prevNumRxns+1:length(ihuman.rxns);
ihuman.rxnRecon3DID(ind) = {''};
ihuman.prRules(ind) = {''};
ihuman.rxnProtMat(ind,:) = 0;
ihuman.priorCombiningGrRules(ind) = {''};

% assign the biomass_components reaction as default objective
ihuman.c(:) = 0;
ihuman.c(ismember(ihuman.rxns,'biomass_components')) = 1;


%% Save model and clear intermediate variables

% clear intermediate varaibles
clearvars -except ihuman

% save model file
save('../../ModelFiles/mat/humanGEM.mat','ihuman');


