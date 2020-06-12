function charge_diff = crossCompChargeDiff(model,rxns,reqFormula)
%crossCompChargeDiff  Determine compartment charge changes caused by rxns.
%
%   This function calculates the net charge change in each compartment
%   caused by each given reaction based on the charges and compartments of 
%   metabolites involved in the reaction.
%
% USAGE:
%
%   charge_diff = crossCompChargeDiff(model,rxns,reqFormula);
%
%
% INPUTS:
%
%   model   Model structure.
%
%   rxns    (Optional, default is all model.rxns) A list of model reaction
%           IDs (or rxn indices) to be evaluated.
%
%   reqFormula   (Optional, default FALSE) If TRUE, reactions involving at
%                least one metabolite that does not have an associated
%                formula (in .metFormulas) will not be evaluated (i.e., the
%                result for that reaction will contain NaNs).
%                If FALSE, then all reactions provided will be evaluated.
%
% OUTPUTS:
%
%   charge_diff   Structure containing the following fields:
%
%                 .mat    An R x C matrix, where R = numer of input rxns
%                         and C = number of model compartments. For each
%                         reaction, the net change in charge of each
%                         compartment is given in the corresponding column.
%                        
%                      NOTE: If reqFormula is TRUE, then rows corresponding
%                            to rxns involving met(s) without formulas will
%                            contain all NaN's.
%
%                 .rxns   A list of reaction IDs corresponding to the rows
%                         of charge_diff.mat
%
%                 .comps  A list of compartment abbreviations corresponding
%                         to the columns of charge_diff.mat.
%
%                 .compNames   A list of compartment names corresponding
%                              to the columns of charge_diff.mat
%


% handle input args
if nargin < 3
    reqFormula = false;
end
if nargin < 2 || isempty(rxns)
    rxns = model.rxns;
end

if isnumeric(rxns)
    % if rxns are supplied as indices
    rxns = model.rxns(rxns);
elseif ~iscell(rxns)
    % if reaction is supplied as a string
    rxns = {rxns};
end

% get the reaction indices in the model
[~,rxn_ind] = ismember(rxns,model.rxns);
if any(rxn_ind == 0)
    error('One or more rxns in the provided list were not found in the model.');
end

% identify all reactions involving mets without formulas
if ( reqFormula )
    exclude_rxns = find(any(model.S(cellfun(@isempty,model.metFormulas),:) ~= 0,1)');
else
    exclude_rxns = [];
end

% generate metCompChargeMat
% This produces a matrix where each row corresponds to a metabolite and 
% each column corresponds to a compartment. For each metabolite, the
% metCharge is entered in the row corresponding to that metabolite's index,
% and the column corresponding to that metabolite's compartment.
metCompChargeMat = sparse(1:length(model.mets),model.metComps,double(model.metCharges));

% calculate net compartment charge changes for each reaction
rxnChargeDiffMat = model.S(:,rxn_ind)' * metCompChargeMat;
rxnChargeDiffMat(ismember(rxn_ind,exclude_rxns),:) = NaN;  % remove values for rxns containing mets without formulas (if specified)

% organize output
charge_diff.mat = full(rxnChargeDiffMat);
charge_diff.rxns = rxns;
charge_diff.comps = model.comps;
charge_diff.compNames = model.compNames;




