%
% FILE NAME:    curateMitochMembraneComp.m
% 
% PURPOSE: Script for analyzing/curating the inner mitochondrial matrix
%          compartment "[i]". This is done in a few steps:
%
%          1. Remove duplicated electron transport chain reactions
%             - These reactions are duplicated because one version comes
%               from Recon3D, whereas the other from HMR. The HMR version
%               of these rxns will be kept, and the Recon3D version
%               deleted.
%             - Due to proton/compartment differences, these reactions were
%               not identified previously as duplicated, and therefore are
%               removed here.
%             - The rxnAssoc.mat file is updated accordingly.
%
%          2. Two unused/dead-end rxns were constrained to zero, and
%             flagged for future deletion.
%             - Both of these reactions came from Recon3D, and appeared to
%               be related to the proton gradient:
%                   r1330  'H+[m] => H+[c] + Proton-Gradient[m]'
%                   r1331  'H+[c] => H+[s] + Proton-Gradient[c]'
%
%          3. Some rxn bounds were updated to prevent energy-generating
%             proton pumping (M to C)
%
%          4. All electrogenic rxns involving proton transport between
%             mitochondria and cytoplasm were updated to take place between
%             the mitochondria and inner mitochondrial membrane, "i".
%
%          5. Update the bounds on ATP-driven transport reactions to
%             prevent generation of ATP.
%             - These reactions should not be reversible, as they lead to
%               artificial production.
%
%          6. All changes to reactions (removal, bounds, stoichiometry,
%             etc.) are written to a .txt file for documentation purposes.
%             - Written to "curateMitochMembraneComp_rxnChanges.txt"
%


%% Load Model

% load HumanGEM model (if not already loaded)
if ~exist('ihuman','var')
    load('humanGEM.mat');  % version 0.4.2
end
ihuman_orig = ihuman;  % to keep track of changes later


%% Remove Recon3D ATP synthase, complex I,  reaction

% the HMR ATP synthase reaction will be kept, but we need to document the
% reaction association in the model and the rxnAssoc.mat file
%
%  *ATP synthase
%     HMR_6916: ADP[m] + 4 H+[c] + Pi[m] => ATP[m] + 4 H+[m] + H2O[m]
%      ATPS4mi: ADP[m] + 4 H+[i] + Pi[m] => ATP[m] + 3 H+[m] + H2O[m]
%
%  *Complex I
%     HMR_6921: 5 H+[m] + NADH[m] + ubiquinone[m] => NAD+[m] + ubiquinol[m] + 4 H+[c]
%  NADH2_u10mi: 5 H+[m] + NADH[m] + ubiquinone[m] => NAD+[m] + ubiquinol[m] + 4 H+[i]
%
%  *Complex III
%     HMR_6918: 2 ferricytochrome C[m] + 2 H+[m] + ubiquinol[m] => 2 ferrocytochrome C[m] + ubiquinone[m] + 4 H+[c]
%   CYOR_u10mi: 2 ferricytochrome C[m] + 2 H+[m] + ubiquinol[m] => 2 ferrocytochrome C[m] + ubiquinone[m] + 4 H+[i]
%
%  *Complex IV
%     HMR_6914: 4 ferrocytochrome C[m] + 8 H+[m] + O2[m] => 4 ferricytochrome C[m] + 2 H2O[m] + 4 H+[c]
%      CYOOm2i: 4 ferrocytochrome C[m] + 8 H+[m] + O2[m] => 4 ferricytochrome C[m] + 2 H2O[m] + 4 H+[i]
%
hmr_rxns = {'HMR_6916';'HMR_6921';'HMR_6918';'HMR_6914'};
[~,hmr_ind] = ismember(hmr_rxns,ihuman.rxns);
r3_rxns = {'ATPS4mi';'NADH2_u10mi';'CYOR_u10mi';'CYOOm2i'};
[~,r3_ind] = ismember(r3_rxns,ihuman.rxns);

% if some of the Recon3D reactions have already been removed from the
% model, then exclude those from the process
remove_ind = (r3_ind == 0);
hmr_rxns(remove_ind) = [];
hmr_ind(remove_ind) = [];
r3_rxns(remove_ind) = [];
r3_ind(remove_ind) = [];

if ~isempty(r3_ind)
    
    fprintf('\n%u reactions (%s) will be removed from the model\n',length(r3_ind),strjoin(r3_rxns,', '));
    fprintf('because they are duplicates of existing HMR reactions.\n\n');
    
    % update "rxnRecon3DID" field
    for i = 1:length(hmr_rxns)
        if isempty(ihuman.rxnRecon3DID{hmr_ind(i)})
            ihuman.rxnRecon3DID(hmr_ind(i)) = r3_rxns(i);
        else
            ihuman.rxnRecon3DID(hmr_ind(i)) = join(unique([strsplit(ihuman.rxnRecon3DID{hmr_ind(i)},';'), r3_rxns(i)]), ';');
        end
    end
    
    % load rxnAssoc.mat file
    load('ComplementaryScripts/modelIntegration/rxnAssoc.mat');
    
    % ignore any associations that already exist in rxnAssoc.mat
    remove_ind = ismember(join([hmr_rxns,r3_rxns],'***'),join([rxnAssoc.rxnHMRID,rxnAssoc.rxnRecon3DID],'***'));
    add_hmr_rxns = hmr_rxns(~remove_ind);
    add_hmr_ind = hmr_ind(~remove_ind);
    add_r3_rxns = r3_rxns(~remove_ind);
    add_r3_ind = r3_ind(~remove_ind);
    
    % add associations to rxnAssoc structure
    rxnAssoc.rxnHMRID = [rxnAssoc.rxnHMRID; add_hmr_rxns];
    rxnAssoc.rxnRecon3DID = [rxnAssoc.rxnRecon3DID; add_r3_rxns];
    rxnAssoc.lbHMR = [rxnAssoc.lbHMR; ihuman.lb(add_hmr_ind)];
    rxnAssoc.ubHMR = [rxnAssoc.ubHMR; ihuman.ub(add_hmr_ind)];
    rxnAssoc.lbRecon3D = [rxnAssoc.lbRecon3D; ihuman.lb(add_r3_ind)];
    rxnAssoc.ubRecon3D = [rxnAssoc.ubRecon3D; ihuman.ub(add_r3_ind)];

    fprintf('The rxnAssoc.mat file has been updated with rxns related to the electron transport chain.\n\n');
    % update rxnAssoc.mat file in the end
    
    % delete Recon3D reactions from model
    ihuman = removeReactionsFull(ihuman,r3_rxns);
    
    % record changes in rxnNotes array
    rxnNotes = [r3_rxns, join([repmat({'reaction is a duplicate of'},length(r3_rxns),1) ,hmr_rxns])];
    
end


%% Fix other reactions
% Reactions: r1330  'H+[m] => H+[c] + Proton-Gradient[m]'
%            r1331  'H+[c] => H+[s] + Proton-Gradient[c]'
% These reactions originate from Recon3D. They appear to be an attempt to
% fix ATP leakage, but will not be necessary/compatible with the updated
% handling of mitochondrial proton transport.

del_rxns = {'r1330'; 'r1331'};

% first check if reactions have already been removed
del_rxns(~ismember(del_rxns,ihuman.rxns)) = [];

if ~isempty(del_rxns)
    
    % perform soft deletion of reactions
    [~,del_ind] = ismember(del_rxns,ihuman.rxns);
    ihuman.lb(del_ind) = 0;
    ihuman.ub(del_ind) = 0;
    
    % note that these should be scheduled for hard deletion in the future
    rxnNotes = [rxnNotes; [del_rxns,repmat({'reaction is dead-end and uneccessary, should be DELETED'},length(del_rxns),1)]];
    
    % output stats
    fprintf('An additional %u dead-end/unused reactions were constrained to zero in the model:\n',length(del_rxns));
    fprintf('\t%s\n',del_rxns{:});
    fprintf('\n');
end



%% Analysis of reactions transporting protons into/out of the mitochondria
% Identify all reactions that involve proton transport between the
% cytoplasm [c] and mitochondria [m].

% construct rxn equations
eqns = constructEquations(ihuman);

% get indices of protons in relevant compartments
Hc = getIndexes(ihuman,'H+[c]','metscomps');
Hm = getIndexes(ihuman,'H+[m]','metscomps');
Hi = getIndexes(ihuman,'H+[i]','metscomps');

% find all rxns that can transport protons from [m] to [c]
fwd_rxn_inds = find((ihuman.S(Hm,:) < 0) & (ihuman.S(Hc,:) > 0) & (ihuman.ub > 0)')';  % written in forward direction
rev_rxn_inds = find((ihuman.S(Hm,:) > 0) & (ihuman.S(Hc,:) < 0) & (ihuman.lb < 0)')';  % written in reverse direction
rxn_inds = [fwd_rxn_inds; rev_rxn_inds];

fprintf('A total of %u reactions are able to transport protons from [m] to [c].\n\n',length(rxn_inds));

% This resulted in 8 fwd rxns, and 21 rev rxns, for a total of 29 rxns



%% Update rxn bounds to prevent M --> C proton pumping for electrogenic rxns
% Reactions that involve proton transport between the mitochondria [m] and
% the cytoplasm [c] will not be modified as long as they are
% non-electrogenic (i.e., do not involve a net change in either compartment
% charge). If they ARE electrogenic, then proton pumping from [m] to [c]
% should only be permitted for the electron transport chain, but restricted
% otherwise.

% determine which of the M --> C proton pump reactions are electrogenic
charge_diff = crossCompChargeDiff(ihuman,rxn_inds,false);
non_elec = all(charge_diff.mat == 0,2);
fprintf('%u of the %u rxns able to transport protons from [m] to [c] are non-electrogenic,\n',sum(non_elec),length(rxn_inds));
fprintf('and will therefore be permitted.\n\n');
rxn_inds(non_elec) = [];  % don't constrain non-electrogenic transport rxns

% Three reactions are allowed to pump protons from [m] to [c]:
% Complex IV (HMR_6914), Complex III (HMR_6918), Complex I (HMR_6921).
[~,allowed_inds] = ismember({'HMR_6914';'HMR_6918';'HMR_6921'},ihuman.rxns);
non_elec(ismember(rxn_inds,allowed_inds)) = [];
rxn_inds(ismember(rxn_inds,allowed_inds)) = [];  % don't constrain these reactions
fprintf('ETC rxns involving Complex IV (HMR_6914), Complex III (HMR_6918), and Complex I (HMR_6921)\n');
fprintf('are also permitted to transport protons from [m] to [c].\n\n');


% constrain the remaining electrogenic reactions such that protons can only
% be moved from [c] to [m]
fprintf('The remaining %u electrogenic, proton-transporting rxns will be constrained\n',length(rxn_inds));
fprintf('to prevent their proton transport from [m] to [c].\n\n');
ihuman.ub(rxn_inds(ismember(rxn_inds,fwd_rxn_inds))) = 0;
ihuman.lb(rxn_inds(ismember(rxn_inds,rev_rxn_inds))) = 0;
ihuman.rev = double(ihuman.lb < 0 & ihuman.ub > 0);  % update ".rev" field
rxnNotes = [rxnNotes; [ihuman.rxns(rxn_inds), repmat({'constrained to prevent electrogenic proton transport from [m] to [c]'},length(rxn_inds),1)]];

% For reactions that are now irreversible, but only in the direction that
% is backward from the way it is written, turn those reactions around, so
% they are written in the direction they proceed.
flip_ind = ihuman.lb < 0 & ihuman.ub == 0;
ihuman.ub(flip_ind) = -ihuman.lb(flip_ind);
ihuman.lb(flip_ind) = 0;
ihuman.S(:,flip_ind) = -ihuman.S(:,flip_ind);
fprintf('The following %u rxns have been re-written in the reverse direction:\n',sum(flip_ind));
fprintf('\t%s\n',ihuman.rxns{flip_ind});
fprintf('\n');
rxnNotes = [rxnNotes; [ihuman.rxns(flip_ind), repmat({'turned reaction around to reflect new bounds'},sum(flip_ind),1)]];


%% Update electrogenic rxns involving C <--> M proton transport to be I <--> M
% Protons in the mitochondrial compartment [m] will now transport to/from 
% the inner mitochrondrial membrane compartment [i] for reactions that are
% electrogenic.
% NOTE: Reactions that transport protons from [c] to [m] (and vice versa)
% but are NOT electrogenic will remain in the original form (i.e., the
% cytoplasmic proton will NOT be changed to the [i] compartment).

% identify all electrogenic rxns involving C <--> M proton transport (either direction)
charge_diff = crossCompChargeDiff(ihuman,[],false);
electrogenic = any(charge_diff.mat ~= 0, 2);
rxn_inds = find( (ihuman.S(Hm,:) ~= 0) & (sign(ihuman.S(Hm,:)) == -sign(ihuman.S(Hc,:))) & electrogenic' );
num_non_elec = full(sum( (ihuman.S(Hm,:) ~= 0) & (sign(ihuman.S(Hm,:)) == -sign(ihuman.S(Hc,:))) & ~electrogenic' ));

% change the compartment of all cytoplasmic protons (C) to inner mitochondrial membrane (I)
ihuman.S(Hi,rxn_inds) = ihuman.S(Hc,rxn_inds);
ihuman.S(Hc,rxn_inds) = 0;
fprintf('%u electrogenic rxns involving proton transport between [m] and [c] (either direction) were modified\n',length(rxn_inds));
fprintf('to replace cytoplasmic protons (H+[c]) with inner mitochondrial membrane protons (H+[i]).\n\n');
fprintf('The remaining %u non-electrogenic rxns involving proton transport between [m] and [c] will not be\n',num_non_elec);
fprintf('modified, and their cytoplasmic protons (H+[c]) will be left as is.\n\n');

% record changes in notes array
rxnNotes = [rxnNotes; [ihuman.rxns(rxn_inds), repmat({'replaced cytoplasmic protons with inner mitochondrial membrane protons'},length(rxn_inds),1)]];


%% Ensure that all rxns involving ATP-driven transport cannot generate ATP
% For example, met[c] + ATP[c] + H2O[c] --> met[p] + ADP[c] + Pi[c] should
% not be able to operate in the reverse direction.

% identify all rxns involving transport of any met(s) between compartments
reactants = arrayfun(@(i) ihuman.metNames(ihuman.S(:,i) < 0),(1:length(ihuman.rxns))','UniformOutput',false);
products  = arrayfun(@(i) ihuman.metNames(ihuman.S(:,i) > 0),(1:length(ihuman.rxns))','UniformOutput',false);
trans_inds = arrayfun(@(i) any(ismember(reactants{i},products{i})),(1:length(ihuman.rxns))');

% identify all rxns involving ATP <--> ADP conversion
atp_inds = find(ismember(ihuman.metNames,'ATP'));
adp_inds = find(ismember(ihuman.metNames,'ADP'));
atp_consume = (any(ihuman.S(atp_inds,:) < 0) & any(ihuman.S(adp_inds,:) > 0))';
atp_produce = (any(ihuman.S(atp_inds,:) > 0) & any(ihuman.S(adp_inds,:) < 0))';

% identify transport rxns that can produce ATP
atp_trans_inds = find(trans_inds & (atp_consume & (ihuman.lb < 0)) | ...
                      trans_inds & (atp_produce & (ihuman.ub > 0)) );
% Results in 8 reactions.

% some reactions are allowed (ATP synthase, and ADP/ADP translocases)
allowed_rxns = {'HMR_6916';'HMR_6328';'HMR_4908';'HMR_4907';'HMR_4906';'ATPtg'};
atp_trans_inds(ismember(atp_trans_inds,find(ismember(ihuman.rxns,allowed_rxns)))) = [];

% Results in 2 remaining reactions: 
%        HMR_1696: 3alpha,7alpha-dihydroxy-5beta-cholestanate[c] + ATP[c] + H2O[c] <=> 3alpha,7alpha-dihydroxy-5beta-cholestanate[p] + ADP[c] + Pi[c]
%  4GLU56DIHDINDt: ATP[c] + H2O[c] + 4-S-Glutathionyl-5,6-Dihydroxyindoline[c] <=> ADP[c] + H+[c] + Pi[c] + 4-S-Glutathionyl-5,6-Dihydroxyindoline[s]

% restrict the bounds on the remaining reactions to prevent ATP generation
ihuman.lb(atp_trans_inds(ismember(atp_trans_inds,find(atp_consume)))) = 0;
ihuman.ub(atp_trans_inds(ismember(atp_trans_inds,find(atp_produce)))) = 0;

% update ".rev" field in the model, if necessary
ihuman.rev = double(ihuman.lb < 0 & ihuman.ub > 0);

% record changes in notes array
rxnNotes = [rxnNotes; [ihuman.rxns(atp_trans_inds), repmat({'changed bounds to prevent ATP generation'},length(atp_trans_inds),1)]];


%% write reaction change documentation file

rxnChanges = docRxnChanges(ihuman_orig,ihuman,rxnNotes);
writeRxnChanges(rxnChanges,'curateMitochMembraneComp_rxnChanges');


%% clear intermediate variables and save final results

clearvars -except ihuman rxnAssoc rxnChanges
save('../../model/Human-GEM.mat','ihuman');
save('../modelIntegration/rxnAssoc.mat','rxnAssoc');
movefile('curateMitochMembraneComp_rxnChanges.tsv','../../ComplementaryData/modelCuration/');

