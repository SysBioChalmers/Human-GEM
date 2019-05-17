function rxnChanges = docRxnChanges(model1,model2,rxnNotes)
%docRxnChanges  Document changes to model reactions.
%
%   docRxnChanges identifies reaction-related differences between two
%   models. The function reports the original and new reaction equations, 
%   bounds, and grRules, and any user-supplied notes for each reaction.
%
% USAGE:
%
%   rxnChanges = docRxnChanges(model1,model2,rxnNotes,fileName);
%
% INPUT:
%
%   model1     Model structure before changing reaction(s).
%
%   model2     Model structure after changing reaction(s).
%
%   rxnNotes   (Optional) an Nx2 cell array, where the first column
%              contains the IDs of reactions that were changed, and the
%              second column contains notes (strings) associated with that
%              change. 
%              
%              Reactions that were changed but are not present in
%              the rxnNotes array, or if no rxnNotes array is provided,
%              will simply be given an empty note {''}.
%
%              Reactions that were not identified as changed by the
%              function, but were included in the rxnNotes array, will be
%              included in the documented changes, with the supplied note.
%
% OUTPUT:
%
%   rxnChanges  A cell structure containing the original and new reaction
%               properties, with the following fields:
%                 rxn         reaction ID
%                 eqnOrig     original reaction equation
%                 eqnNew      new reaction equation
%                 lbOrig      original lower bound
%                 lbNew       new lower bound
%                 ubOrig      original upper bound
%                 ubNew       new upper bound
%                 grRuleOrig  original grRules
%                 grRuleNew   new grRules
%                 notes       optional notes associated with the change
%
%
% Jonathan Robinson, 2019-05-17


% handle input arguments
if nargin < 3
    rxnNotes = [];
end


% get rxns from rxnNotes, and include them in the list of changed rxns
if ~isempty(rxnNotes)
    % check if any rxn IDs are repeated in the rxnNotes array, and if so,
    % merge those notes using a semicolon separator (;)
    [uniq_rxn,uniq_ind] = unique(rxnNotes(:,1));
    if length(uniq_rxn) < length(rxnNotes(:,1))
        for i = 1:length(uniq_rxn)
            ind = ismember(rxnNotes(:,1),uniq_rxn(i));
            if sum(ind) > 1
                rxnNotes{uniq_ind(i),2} = strjoin(rxnNotes(ind,2),'; ');
            end
        end
        rxnNotes = rxnNotes(uniq_ind,:);  % remove duplicated entries after merging notes
    end
    chg_rxn = rxnNotes(:,1);
else
    chg_rxn = [];
end

% get indices mapping original model to new model reactions
[is_present,rxn_ind] = ismember(model1.rxns,model2.rxns);
del_rxn = model1.rxns(~is_present);  % record any deleted reactions
nondel_rxn = model1.rxns(is_present);  % get list of non-deleted reactions

% generate reaction equations for each model
model1.eqns = constructEquations(model1);
model2.eqns = constructEquations(model2);

% check for added reactions
add_rxn = setdiff(model2.rxns,model1.rxns);

% check which reactions have a changed equation or grRule
ind = (1:length(model1.rxns));
chg_eqn = ~arrayfun(@(i) isequal(model1.eqns(i),model2.eqns(rxn_ind(i))),ind(is_present));
chg_rule = ~arrayfun(@(i) isequal(model1.grRules(i),model2.grRules(rxn_ind(i))),ind(is_present));

% check which reactions have a changed bound
chg_bnd = model1.lb(is_present) ~= model2.lb(rxn_ind(is_present)) | ...
          model1.ub(is_present) ~= model2.ub(rxn_ind(is_present));

% collect unique list of all changed reactions
model1.rxns = [model1.rxns; add_rxn];
chg_rxn = unique([chg_rxn; nondel_rxn(chg_eqn); nondel_rxn(chg_bnd); nondel_rxn(chg_rule); del_rxn; add_rxn]);
[~,chg_rxn_ind_orig] = ismember(chg_rxn,model1.rxns);
[~,chg_rxn_ind_new] = ismember(chg_rxn,model2.rxns);

% add a dummy entry to model2 to handle deleted reactions
chg_rxn_ind_new(chg_rxn_ind_new == 0) = length(model2.rxns)+1;
model2.eqns(end+1) = {'(DELETED)'};
model2.lb(end+1) = 0;
model2.ub(end+1) = 0;
model2.grRules(end+1) = {''};

% add dummy entries to model1 to handle added reactions
model1.eqns(end+1:end+numel(add_rxn)) = {'(ADDED)'};
model1.lb(end+1:end+numel(add_rxn)) = 0;
model1.ub(end+1:end+numel(add_rxn)) = 0;
model1.grRules(end+1:end+numel(add_rxn)) = {''};


% organize information into output structure
rxnChanges = {};
rxnChanges.rxns       = model1.rxns(chg_rxn_ind_orig);
rxnChanges.eqnOrig    = model1.eqns(chg_rxn_ind_orig);
rxnChanges.eqnNew     = model2.eqns(chg_rxn_ind_new);
rxnChanges.lbOrig     = model1.lb(chg_rxn_ind_orig);
rxnChanges.lbNew      = model2.lb(chg_rxn_ind_new);
rxnChanges.ubOrig     = model1.ub(chg_rxn_ind_orig);
rxnChanges.ubNew      = model2.ub(chg_rxn_ind_new);
rxnChanges.grRuleOrig = model1.grRules(chg_rxn_ind_orig);
rxnChanges.grRuleNew  = model2.grRules(chg_rxn_ind_new);


% add reaction notes
rxnChanges.notes = repmat({''},size(rxnChanges.rxns));
if ~isempty(rxnNotes)
    [~,rxn_ind] = ismember(rxnNotes(:,1),rxnChanges.rxns);
    if any(rxn_ind == 0)
        error('One or more reaction IDs in rxnNotes were not found in either model.');
    end
    rxnChanges.notes(rxn_ind) = rxnNotes(:,2);
end

