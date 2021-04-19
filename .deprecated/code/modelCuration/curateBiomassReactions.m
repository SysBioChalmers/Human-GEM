%
% FILE NAME:    curateBiomassReactions.m
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
%          4. Convert the subSystems field to an array of cell arrays, in
%             order to comply with model spec of both RAVEN and COBRA
%


%% Load model

% load current version of humanGEM
if ~exist('ihuman','var')
    load('humanGEM.mat');  % version 0.7.0
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
new_rxn_ids = {'biomass_components';'biomass_Recon3D';'biomass_maintenance_Recon3D';'biomass_maintenance_noTrTr_Recon3D';'biomass_HMR_RenalCancer'};
[~,ind] = ismember(biomass_rxns,ihuman.rxns);
ihuman.rxns(ind) = new_rxn_ids;

% update biomass reaction subsystem to "Artificial reactions"
ihuman.subSystems(ind) = {'Artificial reactions'}; 

% by default, activate only the biomass_components rxn
ihuman = setParam(ihuman,'eq',new_rxn_ids(2:end),0);


%% Add new reactions related to biomass production
% The HepG2 biomass reaction, and it's associated pool rxns, will be added
% to the model. In addition, reactions enabling export of some fatty acids
% will be incorporated to enable the biomass_components reaction to carry
% flux.

% load new reaction information from file
fid = fopen('../../ComplementaryData/modelCuration/rxns4biomass_20181129.tsv');
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
%                 metName            metID
nameIDpairs = {'fatty acid pool',   'm03171'
               'heparan sulfate',   'm03172'
               'metabolite pool',   'm03173'
               'phosphatidate',     'm03174'
               'phosphatidyl pool', 'm03175'
               'protein pool',      'm03176'};
if any(ismember(nameIDpairs(:,1),ihuman.metNames)) || any(ismember(nameIDpairs(:,2),ihuman.mets))
    error('One or more met names and/or IDs to be added already exist in the model!');
end
[hasmatch,ind] = ismember(new_metNames,nameIDpairs(:,1));
new_metIDs(hasmatch) = strcat(nameIDpairs(ind(hasmatch),2),new_metComps(hasmatch));

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


%% Update model subSystems field to an array of cell arrays
% this change is necessary for compatibility with RAVEN and COBRA packages
non_cell = ~cellfun(@iscell, ihuman.subSystems);
ihuman.subSystems(non_cell) = cellfun(@(x) {{x}}, ihuman.subSystems(non_cell));


%% Save model and clear intermediate variables

% clear intermediate varaibles
clearvars -except ihuman

% save model file
save('../../model/Human-GEM.mat','ihuman');


