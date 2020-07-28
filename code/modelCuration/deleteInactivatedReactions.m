%
% FILE NAME:    deleteInactivatedReactions.m
% 
% PURPOSE: This script performs a full removal ("hard deletion") of
%          inactivated (LB = UB = 0) reactions. These reactions were
%          previously identified as invalid.
%
%          The exceptions are the alternative biomass reactions, which will
%          remain inactivated but will not be removed from the model.
%
%              biomass_components
%              biomass_Recon3D
%              biomass_maintenance_Recon3D
%              biomass_maintenance_noTrTr_Recon3D
%              biomass_HMR_RenalCancer
%


%% Process Human-GEM model

% load model
ihuman = importHumanYaml('../../model/Human-GEM.yml');

% find all inactivated reactions
inact_ind = find( (ihuman.ub == 0) & (ihuman.lb == 0) );

% exclude biomass reactions
inact_ind = setdiff(inact_ind, find(startsWith(ihuman.rxns, 'biomass')));

% delete reactions, and remove metabolites or genes that become unused
% after the reaction deletion
model_del = removeReactions(ihuman, inact_ind, true, true, true);

% print some info
n_rxns = numel(ihuman.rxns) - numel(model_del.rxns);
n_mets = numel(ihuman.mets) - numel(model_del.mets);
n_genes = numel(ihuman.genes) - numel(model_del.genes);
n_comps = numel(ihuman.comps) - numel(model_del.comps);
fprintf('Deleted %u inactivated reactions from the model.\n', n_rxns);
if n_mets > 0
    fprintf('Deleted %u now-unused metabolites from the model.\n', n_mets);
end
if n_genes > 0
    fprintf('Deleted %u now-unused genes from the model.\n', n_genes);
end
if n_comps > 0
    fprintf('Deleted %u now-unused compartments from the model.\n', n_comps);
end
% Note: No genes or compartments end up being removed

% write new model to yml
writeHumanYaml(model_del, '../../modelFiles/Human-GEM.yml');


%% Process Human-GEM annotation JSONs

% load annotation JSONs
metAssoc = jsondecode(fileread('../../data/annotation/humanGEMMetAssoc.JSON'));
rxnAssoc = jsondecode(fileread('../../data/annotation/humanGEMRxnAssoc.JSON'));

% remove deleted reactions from the reaction JSON
f = fieldnames(rxnAssoc);
for i = 1:numel(f)
    rxnAssoc.(f{i})(inact_ind) = [];
end

% remove deleted metabolites from the metabolite JSON
remove_met_ind = find(~ismember(ihuman.mets, model_del.mets));
f = fieldnames(metAssoc);
for i = 1:numel(f)
    metAssoc.(f{i})(remove_met_ind) = [];
end

% export reaction and metabolite annotation structures to JSON
jsonStr = jsonencode(rxnAssoc);
fid = fopen('../../data/annotation/humanGEMRxnAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);
jsonStr = jsonencode(metAssoc);
fid = fopen('../../data/annotation/humanGEMMetAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);


