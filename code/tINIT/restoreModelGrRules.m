function restoredModel = restoreModelGrRules(model,refModel)
%restoreModelGrRules  Recover grRules (and genes) from reference model.
%
% The current implementation of tINIT removes genes from grRules, and will
% eliminate enzyme complex information from grRules (i.e., all rules will
% become "OR" expressions).
%
% This function will replace the input MODEL grRules with the corresponding
% original grRule from the reference model (refModel). In addition, the
% function will regenerate the "genes" and "rxnGeneMat" fields based on the
% restored grRules.
%
% USAGE:
%
%   restoredModel = restoreModelGrRules(model,refModel);
%
% INPUTS:
%
%   model      Model structure with grRules needing restoration.
%
%   refModel   Reference model from which the input MODEL was derived.
%
% OUTPUTS:
%
%   restoredModel  Input MODEL returned with grRules retrieved from
%                  refModel, as well as updated genes and geneRxnMat
%                  fields.
%


% restore model grRules with those in refModel
[~,ind] = ismember(model.rxns,refModel.rxns);
model.grRules = refModel.grRules(ind);

% regenerate "genes" and "rxnGeneMat" fields
[genes,rxnGeneMat] = getGenesFromGrRules(model.grRules);
model.genes = genes;
model.rxnGeneMat = rxnGeneMat;

% assign output
restoredModel = model;
