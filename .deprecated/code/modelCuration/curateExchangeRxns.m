%
% FILE NAME:    curateExchangeRxns.m
% 
% PURPOSE: Script to curate exchange, sink, and demand reactions.
%          The current model has a mix of exchange reactions that are
%          formulated such that negative flux corresponds to export:
%
%             HMR_9730: albumin[x] <=> albumin[s]
%
%          whereas others have positive flux corresponding to export:
%
%             EX_adrn[e]: adrenic acid[s] <=> adrenic acid[x]
%
%          In order to standardize this, all exchange reactions will be
%          converted to the second format (positive flux = export).
%
%          However, before this can be incorporated, there are 7 exchange
%          reactions in the model that are duplicated but are written in
%          opposite directions:
%
%             HMR_9025: VLDL[x] => VLDL[s]
%             HMR_9049: VLDL[s] => VLDL[x]
% 
%             HMR_9026: HDL[x] => HDL[s]
%             HMR_9050: HDL[s] => HDL[x]
%             
%             HMR_9027: LDL[x] => LDL[s]
%             HMR_9051: LDL[s] => LDL[x]
%             
%             HMR_9028: chylomicron remnant[x] => chylomicron remnant[s]
%             HMR_9052: chylomicron remnant[s] => chylomicron remnant[x]
%             
%             HMR_9029: VLDL remnant[x] => VLDL remnant[s]
%             HMR_9053: VLDL remnant[s] => VLDL remnant[x]
%             
%             HMR_9030: HDL remnant[x] => HDL remnant[s]
%             HMR_9054: HDL remnant[s] => HDL remnant[x]
%             
%             HMR_9031: LDL remnant[x] => LDL remnant[s]
%             HMR_9055: LDL remnant[s] => LDL remnant[x]
%
%          These reactions must first be merged, otherwise switching
%          reaction directionality will result in these rxn pairs being
%          completely identical. Therefore, for each pair, one reaction was
%          deleted, and the second was made reversible. It was confirmed in
%          advance that none of these reactions are associated with any
%          genes.
%
%          To keep a clear and transparent curation process, a new array
%          structure "redundantRxns" is generated for storing such nearly
%          identical reaction pairs (#43) and clarifying which to keep or
%          remove, as well as some explanations. Finally, this structure is
%          saved as a plaintext JSON file, to improve this repo toward a
%          binary-free mode (#27).
%
%          ----------------------------------------------------------------
%
%          The script also inactivates all sink and demand (DM) reactions
%          by constraining their upper and lower bounds to zero.
%          These reactions will be considered for full deletion in future
%          model versions.
%
%          ----------------------------------------------------------------
%
%          Finally, the upper and lower bounds of all exchange reactions
%          are set to +/-1000, respectively. As a result, the model is by
%          default completely "open", allowing free exchange of all
%          metabolites.
%


%% Load model and initialize variables

% load current version of humanGEM
if ~exist('ihuman','var')
    load('humanGEM.mat');  % version 0.6.1
end
ihuman_orig = ihuman;  % to keep track of changes made


%% Remove duplicated exchange reactions
% The current model contains 7 pairs of duplicated exchange reactions
% involving boundary metabolites. Each of these pairs involves irreversible
% reactions, but in opposite directions. Therefore, they will be merged
% into a single reversible reaction. There are no genes associated with any
% of these reactions.

% Define pairs of duplicated (but opposite direction) exchange reactions.
% The second reaction in each pair (those in the second column) will be
% kept, whereas the first column will be removed. This is because the rxns
% in the second column are written in the preferred direction ([s] => [x])
% and have a corresponding Recon3D rxnID.
rxnPairs = {'HMR_9025','HMR_9049';   % 'VLDL[x] => VLDL[s]', 'VLDL[s] => VLDL[x]'
            'HMR_9026','HMR_9050';   % 'HDL[x] => HDL[s]',   'HDL[s] => HDL[x]'
            'HMR_9027','HMR_9051';   % 'LDL[x] => LDL[s]',   'LDL[s] => LDL[x]'
            'HMR_9028','HMR_9052';   % 'chylomicron remnant[x] => chylomicron remnant[s]', 'chylomicron remnant[s] => chylomicron remnant[x]'
            'HMR_9029','HMR_9053';   % 'VLDL remnant[x] => VLDL remnant[s]', 'VLDL remnant[s] => VLDL remnant[x]'
            'HMR_9030','HMR_9054';   % 'HDL remnant[x] => HDL remnant[s]',   'HDL remnant[s] => HDL remnant[x]'
            'HMR_9031','HMR_9055'};  % 'LDL remnant[x] => LDL remnant[s]',   'LDL remnant[s] => LDL remnant[x]'

% get corresponding reaction indices
[~,rxnInds] = ismember(rxnPairs,ihuman.rxns);

% make the reactions reversible
ihuman.lb(rxnInds(:,2)) = -1000;
ihuman.ub(rxnInds(:,2)) = 1000;

% generate annotation notes for changed/deleted reactions
rxnNotes = strcat('reaction is identical to ',rxnPairs(:,2),', but in opposite direction; these rxns were therefore merged into ',rxnPairs(:,2),', which was made reversible.');
rxnNotes = regexprep(rxnNotes,'toHMR','to HMR');  % fix loss of spaces
rxnNotes = [rxnPairs(:,1), rxnNotes];

addNotes = strcat('reaction is identical to ',rxnPairs(:,1),', but in opposite direction; these rxns were therefore merged into ',rxnPairs(:,2),', which was made reversible.');
addNotes = regexprep(addNotes,'toHMR','to HMR');  % fix loss of spaces
addNotes = [rxnPairs(:,2), addNotes];
rxnNotes = [rxnNotes; addNotes];

% save these redundant reaction pairs to a designed array structure
redundantRxns.rxnKeep = rxnPairs(:,2);
redundantRxns.rxnRemove = rxnPairs(:,1);
rxnEqns = constructEquations(ihuman);
redundantRxns.eqnKeep = rxnEqns(rxnInds(:,2));
redundantRxns.eqnRemove = rxnEqns(rxnInds(:,1));
redundantRxns.grRuleKeep = ihuman.grRules(rxnInds(:,2));
redundantRxns.grRuleRemove = ihuman.grRules(rxnInds(:,1));
redundantRxns.notes=repmat({'This pair of exchange reactions are the same but in opposite directions, they will be merged into a single reversible reaction'},size(rxnInds,1),1);

% delete reactions
ihuman = removeReactionsFull(ihuman,rxnPairs(:,1),false,false,false);


%% Identify and flip directionality of backward exchange reactions

% generate simplified model, which has all boundary metabolites removed
smodel = simplifyModel(ihuman);

% Identify all exchange, sink, and demand reactions. These are identified
% as any reaction that contains only one metabolite.
exch_ind = find(sum(smodel.S ~= 0) == 1)';

% Adjust directionality in these reactions so that negative flux
% corresponds to a metabolite entering the system (import), whereas
% positive flux corresponds to metabolite export.
flip_ind = exch_ind(sum(smodel.S(:,exch_ind)) > 0);  % identify rxns that will change direction
ihuman.S(:,flip_ind) = -ihuman.S(:,flip_ind);
rxnNotes = [rxnNotes; [ihuman.rxns(flip_ind), repmat({'changed reaction directionality so positive flux corresponds to metabolite export'},length(flip_ind),1)]];

% swap upper and lower bounds of flipped reactions
lb = ihuman.lb(flip_ind);
ihuman.lb(flip_ind) = -ihuman.ub(flip_ind);
ihuman.ub(flip_ind) = -lb;


%% Inactivate all sink and demand reactions

% identify sink and demand reactions
sd_rxns = startsWith(ihuman.rxns,{'sink_','DM_'});

% constrain reaction lower and upper bounds to zero
ihuman.lb(sd_rxns) = 0;
ihuman.ub(sd_rxns) = 0;

% add notes to rxnNotes
rxnNotes = [rxnNotes; [ihuman.rxns(sd_rxns), repmat({'sink/demand reactions were inactivated, and are scheduled for future DELETION'},sum(sd_rxns),1)]];


%% Set upper and lower bounds of all exchange reactions to +/-1000

% find exchange reactions (but not sink/demand reactions)
exch_rxns = exch_ind(~ismember(exch_ind,find(sd_rxns)));

% determine which reactions will be affected by the change, and annotate
% the change in rxnNotes
affected_ind = (ihuman.lb(exch_rxns) ~= -1000) | (ihuman.ub(exch_rxns) ~= 1000);
affected_rxns = ihuman.rxns(exch_rxns(affected_ind));
rxnNotes = [rxnNotes; [affected_rxns, repmat({'bounds of all exchange reactions by default set to +/-1000'},length(affected_rxns),1)]];

% set upper and lower bounds to 1000 and -1000, respectively
ihuman.lb(exch_rxns) = -1000;
ihuman.ub(exch_rxns) = 1000;



%% Save and export results

% save the redundantRxns structure to a JSON file
fid = fopen('redundantRxns.JSON', 'w');
fwrite(fid, jsonencode(redundantRxns));
fclose(fid);
movefile('redundantRxns.JSON','../../ComplementaryData/modelCuration/');

% generate report on modified/deleted reactions
rxnChanges = docRxnChanges(ihuman_orig,ihuman,rxnNotes);
writeRxnChanges(rxnChanges,'curateExchangeRxns_rxnChanges',true);
movefile('curateExchangeRxns_rxnChanges.tsv','../../ComplementaryData/modelCuration/');


%% Clear intermediate vars and save model file

% clear intermediate varaibles
clearvars -except ihuman

% update rev field
ihuman.rev = double(ihuman.lb < 0 & ihuman.ub > 0);

% save model file
save('../../model/Human-GEM.mat','ihuman');



