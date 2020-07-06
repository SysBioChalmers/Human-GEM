%
% FILE NAME:    removeSinkDMRxns.m
% 
% PURPOSE: HumanGEM currently contains many "sink" and "demand" (DM)
%          reactions that originate from Recon3D. These reactions are
%          artificial, and involve the transport of a metabolite between
%          the boundary compartment [x] and a non-extracellular
%          compartment. This is unlike exchange reactions, which only
%          involve transport between the boundary compartment and the
%          extracellular compartment.
%
%          These sink and demand reactions were previously inactivated
%          (upper and lower bounds fixed to zero), but will now be
%          completely removed from the model. These reactions can be
%          identified by their rxn IDs, which all start with "sink_" or
%          "DM_". Logs with notes from this reaction removal are written to
%          the file "removedSinkDMrxns.tsv".
%
%          Note that the script also removes all unused metabolites and
%          genes after deleting the reactions. No genes were removed, as
%          none of the sink or demand reactions had gene associations.
%          However, 52 metabolites are removed because they appear ONLY in
%          the sink/demand reactions, and nowhere else in the model.
%


% load latest version of humanGEM
load('humanGEM.mat');  % version 1.0.0-beta

% find all sink and demand reactions
remInd = startsWith(ihuman.rxns,{'sink_','DM_'});
remRxns = ihuman.rxns(remInd);

% remove reactions and unused metabolites from humanGEM
reducedModel = removeReactionsFull(ihuman,remRxns,true);

% print changes to user
fprintf('Removed %u sink/demand reactions from humanGEM.\n',numel(remRxns));
fprintf('Subsequently removed %u now-unused metabolites participating only in those sink/demand reactions.\n\n',numel(ihuman.mets)-numel(reducedModel.mets));

% document model changes
rxnNotes = repmat({'Sink/Demand reaction removed because it is artificial and unnecessary.'},numel(remRxns),1);
rxnChanges = docRxnChanges(ihuman,reducedModel,[remRxns,rxnNotes]);
writeRxnChanges(rxnChanges,'../../ComplementaryData/modelCuration/removedSinkDMrxns.tsv');

% save new version of humanGEM
ihuman = reducedModel;
save('../../model/Human-GEM.mat','ihuman');



