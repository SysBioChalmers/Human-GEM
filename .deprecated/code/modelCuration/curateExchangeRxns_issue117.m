%
% FILE NAME:    curateExchangeRxns_issue117.m
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
%          Information on new reactions to be added to the model were
%          organized in the newExchangeRxns.tsv file.


%% Load model and initialize variables

% load HumanGEM
load('HumanGEM.mat');
ihuman_orig = ihuman;  % to track changes

% load reaction and metabolite annotation structures
metAssoc = jsondecode(fileread('../../ComplementaryData/annotation/humanGEMMetAssoc.JSON'));
rxnAssoc = jsondecode(fileread('../../ComplementaryData/annotation/humanGEMRxnAssoc.JSON'));
metAssocFields = fieldnames(metAssoc);
rxnAssocFields = fieldnames(rxnAssoc);

% verify that model and annotation structures are aligned
if ~isequal(metAssoc.mets, ihuman.mets) || ~isequal(rxnAssoc.rxns, ihuman.rxns)
    error('HumanGEM mets and rxns are not aligned with those in the annotation structures!');
end

% initialize model change notes
changeNotes = {};



%% Update metabolite IDs and add new metabolites
% the following metabolites have unusual/temporary IDs:
%   temp002x: others[x]
%   temp003x: steroids[x]
%   temp004x: xenobiotics[x]
%   temp005x: arachidonate derivatives[x]
%
% Therefore, these metabolites will be renamed as follows:
%   temp002x -> m10000x
%   temp003x -> m10001x
%   temp004x -> m10002x
%   temp005x -> m10003x

mets_orig = {'temp002x';'temp003x';'temp004x';'temp005x'};
mets_new  = {'m10000x' ;'m10001x' ;'m10002x' ; 'm10003x'};
mets_ind = getIndexes(ihuman,mets_orig,'mets');
ihuman.mets(mets_ind) = mets_new;
metAssoc.mets(mets_ind) = mets_new;

% append change notes
changeNotes = [changeNotes; [[mets_orig;mets_new], repmat({'Updated "temp" metabolite IDs.'},numel([mets_orig;mets_new]),1)]];


% In addition, the model does not contain extracellular forms of these
% metabolites, which will be necessary to add their exchange reactions to
% the model. The following new metabolites will therefore be added:
%   m10000s: others[s]
%   m10001s: steroids[s]
%   m10002s: xenobiotics[s]
%   m10003s: arachidonate derivatives[s]

metsToAdd = {};
metsToAdd.mets     = {'m10000s';'m10001s' ;'m10002s'    ;'m10003s'                 };
metsToAdd.metNames = {'others' ;'steroids';'xenobiotics';'arachidonate derivatives'};
metsToAdd.compartments = 's';

% add new metabolites to model
ihuman = addMets(ihuman,metsToAdd);

% add new metabolites to metAssoc structure
newMetInd = numel(metAssoc.mets) + (1:numel(metsToAdd.mets));
for i = 1:numel(metAssocFields)
    metAssoc.(metAssocFields{i})(newMetInd) = {''};
end
metAssoc.mets(newMetInd) = metsToAdd.mets;
metAssoc.metsNoComp(newMetInd) = regexprep(metsToAdd.mets, '_*\w$', '');

% append change notes
changeNotes = [changeNotes; [metsToAdd.mets, repmat({'Extracellular version of metabolite added to enable the addition of an exchange reaction.'},numel(metsToAdd.mets),1)]];


% Furthermore, a boundary version of the following metabolites do not yet
% exist in HumanGEM, but need to be added in order to allow their exchange:
%   m00591x:    20-hydroxy-arachidonate[x]
%   m01435x:    chenodiol[x]
%   m02328x:    LacCer pool[x]
%   chylo_hs_x: Chylomicron Lipoprotein[x]

metsToAdd = {};
metsToAdd.mets     = {'m00591x'                ;'m01435x'  ;'m02328x'    ;'chylo_hs_x'             };
metsToAdd.metNames = {'20-hydroxy-arachidonate';'chenodiol';'LacCer pool';'Chylomicron Lipoprotein'};
metsToAdd.compartments = 'x';

% add new metabolites to model
ihuman = addMets(ihuman,metsToAdd);

% add new metabolites to metAssoc structure
newMetInd = numel(metAssoc.mets) + (1:numel(metsToAdd.mets));
for i = 1:numel(metAssocFields)
    metAssoc.(metAssocFields{i})(newMetInd) = {''};
end
metAssoc.mets(newMetInd) = metsToAdd.mets;
metAssoc.metsNoComp(newMetInd) = regexprep(metsToAdd.mets, '_*\w$', '');

% append change notes
changeNotes = [changeNotes; [metsToAdd.mets, repmat({'Boundary compartment version of metabolite added to enable the addition of an exchange reaction.'},numel(metsToAdd.mets),1)]];



%% Add new exchange reactions to the model

% retrieve new reaction information from file
fid = fopen('../../ComplementaryData/modelCuration/newExchangeRxns_issue117.tsv');
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
rxnsToAdd.subSystems = cellfun(@(s) {{s}}, rxnData{5});
rxnsToAdd.equations = rxnData{6};

% add new reactions to model
ihuman = addRxns(ihuman,rxnsToAdd,3);
ihuman.prRules(end+1:end+nRxns) = {''};
ihuman.rxnProtMat(end+1:end+nRxns,:) = 0;
ihuman.priorCombiningGrRules(end+1:end+nRxns) = {''};

% add new reactions to rxnAssoc structure
newRxnInd = numel(rxnAssoc.rxns) + (1:numel(rxnsToAdd.rxns));
for i = 1:numel(rxnAssocFields)
    rxnAssoc.(rxnAssocFields{i})(newRxnInd) = {''};
end
rxnAssoc.rxns(newRxnInd) = rxnsToAdd.rxns;

% append change notes
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
rxnInd = getIndexes(ihuman,'HMR_9736','rxns');
ihuman = removeReactionsFull(ihuman,rxnInd);

% remove reaction from rxnAssoc structure
for i = 1:numel(rxnAssocFields)
    rxnAssoc.(rxnAssocFields{i})(rxnInd) = [];
end

% append change notes
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


%% Document reaction changes

modelChanges = docModelChanges(ihuman_orig,ihuman,changeNotes);
writeModelChanges(modelChanges,'../../ComplementaryData/modelCuration/curateExchangeReactions_issue117.tsv');


%% Export reaction and metabolite annotation JSON files

% first verify that annotation structures are still aligned with the model
if ~isequal(metAssoc.mets, ihuman.mets) || ~isequal(rxnAssoc.rxns, ihuman.rxns)
    error('HumanGEM mets and rxns are not aligned with those in the annotation structures!');
end

% write metAssoc to JSON
jsonStr = jsonencode(metAssoc);
fid = fopen('../../ComplementaryData/annotation/humanGEMMetAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);

% write rxnAssoc to JSON
jsonStr = jsonencode(rxnAssoc);
fid = fopen('../../ComplementaryData/annotation/humanGEMRxnAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);


%% Export updated HumanGEM

exportHumanGEM(ihuman,'HumanGEM','../../',{'mat','yml'},false,false);







