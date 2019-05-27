%
% FILE NAME:    curateExchangeReactions2.m
% 
% DATE CREATED: 2019-05-27
%     MODIFIED: 2019-05-27
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script to add exchange reactions for four metabolites that are
%          present in the extracellular compartment, but do not currently
%          have exchange reactions transporting them to/from the boundary.
%          [Addresses Issue #117]
%
%          The four metabolites for which exchange reactions are added:
%          1)  20-hydroxy-arachidonate
%          2)  chenodiol
%          3)  LacCer pool
%          4)  Chylomicron Lipoprotein
%
%
%          In addition, the script addresses the issue of five reactions
%          that involve transport between the boundary compartment and a
%          non-extracellular compartment. The reaction IDs are:
%          1)  HMR_9736
%          2)  xenobiotics
%          3)  arachidonates
%          4)  steroids
%          5)  others
%
%          The HMR_9736 reaction exchanges 'cholesterol-ester pool[l]' with
%          the boundary compartment. Since several other reactions exist in
%          the model that transport this metabolite to other compartments,
%          including from the extracellular compartment to the boundary,
%          this reaction is redundant and should be deleted.
%
%          The other reactions are all pool reactions in which several
%          metabolites combine to form the metabolites 'xenobiotics',
%          'arachidonate derivatives', 'steroids', or 'others' in the
%          boundary compartment. To avoid exchange between
%          non-extracellular compartments and the boundary, these reactions
%          will be modified such that they form the pool metabolites in the
%          extracellular compartment, and additional transport reactions
%          will be added to exchange these metabolites between the
%          extracellular compartment and the boundary.
%
%          Information on new reactions to be added to the model were
%          organized in the newExchangeRxns.tsv file.



%% Load model and initialize variables

% load HumanGEM
load('humanGEM.mat');
ihuman_orig = ihuman;  % to track changes

% initialize reaction change notes
changeNotes = {};



%% Update metabolite IDs and add new metabolites
% the following metabolites have unusual IDs:
%
%   temp002x: others[x]
%   temp003x: steroids[x]
%   temp004x: xenobiotics[x]
%   temp005x: arachidonate derivatives[x]
%
% In addition, the model does not contain an extracellular form of these
% metabolites, which will be necessary to add their exchange reactions
%
% Therefore, these metabolites will be renamed as follows:
%
%   temp002x -> m10000x
%   temp003x -> m10001x
%   temp004x -> m10002x
%   temp005x -> m10003x
%
% And the following corresponding new metabolites will be added:
%
%   m10000s: others[s]
%   m10001s: steroids[s]
%   m10002s: xenobiotics[s]
%   m10003s: arachidonate derivatives[s]
%
% Furthermore, a boundary version of the following metabolites will be
% added to the model:
%
%   m00591x:    20-hydroxy-arachidonate[x]
%   m01435x:    chenodiol[x]
%   m02328x:    LacCer pool[x]
%   chylo_hs_x: Chylomicron Lipoprotein[x]
%

% revise metabolite IDs
mets = {'temp002x';'temp003x';'temp004x';'temp005x'};
mets_ind = getIndexes(ihuman,mets,'mets');
ihuman.mets(mets_ind) = {'m10000x';'m10001x';'m10002x';'m10003x'};
changeNotes = [changeNotes; [[mets;ihuman.mets(mets_ind)], repmat({'Updated metabolite ID to avoid using "temp" IDs.'},numel(mets)*2,1)]];


% add new metabolites
metsToAdd = {};
metsToAdd.mets = {'m10000s';'m10001s';'m10002s';'m10003s';'m00591x';'m01435x';'m02328x';'chylo_hs_x'};
metsToAdd.metNames = {'others';'steroids';'xenobiotics';'arachidonate derivatives';'20-hydroxy-arachidonate';'chenodiol';'LacCer pool';'Chylomicron Lipoprotein'};
metsToAdd.compartments = {'s';'s';'s';'s';'x';'x';'x';'x'};

ihuman = addMets(ihuman,metsToAdd);
changeNotes = [changeNotes; [metsToAdd.mets, repmat({'New compartment version of metabolite added to enable the addition of an exchange reaction.'},numel(metsToAdd.mets),1)]];


%% Add new exchange reactions to the model

% retrieve new reaction information from file
fid = fopen('../../ComplementaryData/modelCuration/newExchangeRxns.tsv');
rxnData = textscan(fid,'%s%s%f%f%s%s','Delimiter','\t','Headerlines',1);
fclose(fid);
nRxns = numel(rxnData{1});

% verify that none of the reactions exist in the current model
if any(ismember(rxnData{1},ihuman.rxns))
    error('One or more reactions to be added already exist in the model.');
end

% construct structure
rxnsToAdd = {};
rxnsToAdd.rxns = rxnData{1};
rxnsToAdd.rxnNames = rxnData{2};
rxnsToAdd.lb = rxnData{3};
rxnsToAdd.ub = rxnData{4};
rxnsToAdd.subSystems = rxnData{5};
rxnsToAdd.equations = rxnData{6};

% add reactions
ihuman = addRxns(ihuman,rxnsToAdd,3);

% update other non-standard fields
ihuman.rxnRecon3DID(end+1:end+nRxns) = {''};
ihuman.prRules(end+1:end+nRxns) = {''};
ihuman.rxnProtMat(end+1:end+nRxns,:) = 0;
ihuman.priorCombiningGrRules(end+1:end+nRxns) = {''};

changeNotes = [changeNotes; [rxnData{1}, repmat({'New exchange reaction to facilitate transport of metabolite between boundadry and extracellular compartment.'},nRxns,1)]];


%% Delete the HMR_9736 reaction
% The following reaction:
%
%   HMR_9736: cholesterol-ester pool[l] <=> cholesterol-ester pool[x]
%
% involves transport between a non-extracellular compartment (the lysosome)
% and the boundary compartment, but only extracellular metabolites should
% be moved to/from the boundary. This reaction is anyway redundant, due to
% the existence of other choleserol-ester pool transport reactions in the
% model:
%
%   HMR_0019:          cholesterol-ester pool[s] => cholesterol-ester pool[l]
%   HMR_3597:          cholesterol-ester pool[l] => cholesterol-ester pool[r]
%   HMR_0020:          cholesterol-ester pool[c] <=> cholesterol-ester pool[r]
%   EX_xolest2_hs[e]:  cholesterol-ester pool[s] <=> cholesterol-ester pool[x]
%   XOLEST2te:         cholesterol-ester pool[s] <=> cholesterol-ester pool[c]
%   XOLEST2HSTDle:     cholesterol-ester pool[c] => cholesterol-ester pool[s]
%
% Therefore, the reaction will be removed from the model.
ihuman = removeReactionsFull(ihuman,'HMR_9736');
changeNotes = [changeNotes; [{'HMR_9736'},{'Reaction is redundant and involves transport between lysosome and boundary, and was therefore DELETED.'}]];


%% Revise pool metabolite compartment in four reactions
% The following pool reactions:
%
%   xenobiotics:    (many metabolites) => xenobiotics[x]
%   arachidonates:  (many metabolites) => arachidonate derivatives[x]
%   steroids:       (many metabolites) => steroids[x]
%   others:         (many metabolites) => others[x]
%
% Involve transport from non-extracellular compartment(s) to the boundary
% compartment. To correct this, the compartment of the pool metabolites
% generated in each of these reactions will be changed from [x] to [s].
rxns = {'xenobiotics';'arachidonates';'steroids';'others'};
rxn_inds = getIndexes(ihuman,rxns,'rxns');

pool_mets = {'xenobiotics';'arachidonate derivatives';'steroids';'others'};
pool_inds_x = getIndexes(ihuman,strcat(pool_mets,'[x]'),'metscomps');
pool_inds_s = getIndexes(ihuman,strcat(pool_mets,'[s]'),'metscomps');

for i = 1:length(rxns)
    ihuman.S(pool_inds_x(i),rxn_inds(i)) = 0;
    ihuman.S(pool_inds_s(i),rxn_inds(i)) = 1;
end

changeNotes = [changeNotes; [rxns, repmat({'updated reaction to generate pool metabolite in extracellular, to avoid transport between non-extracellular compartments and the boundary.'},numel(rxns),1)]];







% %% Document reaction changes
% 
% modelChanges = docModelChanges(ihuman_orig,ihuman,changeNotes);
% writeRxnChanges(rxnChanges,'../../ComplementaryData/modelCuration/curateExchangeReactions2_rxnChanges.tsv');
% 
% 
% %% Export updated HumanGEM
% 
% exportHumanGEM(ihuman,'humanGEM','../../',{'mat','yml'},false,false);







