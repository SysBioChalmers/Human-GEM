function modelChanges = docModelChanges(model1,model2,changeNotes)
%docModelChanges  Document changes to model reactions and/or metabolites.
%
%   docModelChanges identifies differences in reaction- and/or metabolite-
%   related fields between two models.
%
% USAGE:
%
%   modelChanges = docModelChanges(model1,model2,changeNotes);
%
% INPUT:
%
%   model1       Model structure before changes.
%
%   model2       Model structure after changes.
%
%   changeNotes  (Optional) an Nx2 cell array, where the first column
%                contains the IDs of reactions and/or metabolites that were
%                changed, and the second column contains notes associated
%                with that change.
%              
%                Reactions or metabolites that were changed but are not
%                present in the changeNotes array, or if no changeNotes
%                array is provided, will simply have an empty note {''} 
%                associated with them.
%
%                Reactions or metabolites that were not identified as
%                changed, but were included in the changeNotes array, will
%                be included in the documented changes with the supplied
%                note.
%
% OUTPUT:
%
%   modelChanges  A cell structure containing the original and new reaction
%                 and/or metabolite properties. Reaction-related changes,
%                 if found, will be reported in the modelChanges.rxns
%                 field, whereas metabolite-related changes are reported in
%                 the modelChanges.mets field.
%
%                 The modelChanges.rxns field is further subdivided into
%                 the following fields:
%                    rxns         reaction ID
%                    eqnOrig      original reaction equation
%                    eqnNew       new reaction equation
%                    lbOrig       original lower bound
%                    lbNew        new lower bound
%                    ubOrig       original upper bound
%                    ubNew        new upper bound
%                    grRuleOrig   original grRules
%                    grRuleNew    new grRules
%                    notes        optional notes associated with the change
%
%                 The modelChanges.mets field is further subdivided into
%                 the following fields:
%                    mets         metabolite ID
%                    nameOrig     original metabolite name
%                    nameNew      new metabolite name
%                    formulaOrig  original metabolite formula
%                    formulaNew   new metabolite formula
%                    chargeOrig   original metabolite charge
%                    chargeNew    new metabolite charge
%                    notes        optional notes associated with the change
%


%% Preprocessing

% handle input arguments
if nargin < 3
    changeNotes = [];
end

if ~isempty(changeNotes)
    % check if any rxn or met IDs are repeated in the changeNotes array, 
    % and if so, merge those notes using a semicolon separator (;)
    [uniq_id,uniq_index] = unique(changeNotes(:,1));
    if length(uniq_id) < length(changeNotes(:,1))
        for i = 1:length(uniq_id)
            ind = ismember(changeNotes(:,1),uniq_id(i));
            if sum(ind) > 1
                changeNotes{uniq_index(i),2} = strjoin(changeNotes(ind,2),'; ');
            end
        end
        changeNotes = changeNotes(uniq_index,:);  % remove duplicated entries after merging notes
    end

    % get rxn and met IDs from changeNotes, and include them in the list of
    % changed fields
    chg_rxn = intersect(changeNotes(:,1),union(model1.rxns,model2.rxns));
    chg_met = intersect(changeNotes(:,1),union(model1.mets,model2.mets));
    if any(~ismember(changeNotes(:,1),[model1.rxns; model2.rxns; model1.mets; model2.mets]))
        error('One or more reaction or metabolite IDs in changeNotes were not found in either model.');
    end
else
    chg_rxn = [];
    chg_met = [];
end


%% Assess changes in reaction-related fields

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

% organize reaction change information into structure
rxnChanges = {};
if ~isempty(chg_rxn)
    rxnChanges.rxns       = model1.rxns(chg_rxn_ind_orig);
    rxnChanges.eqnOrig    = model1.eqns(chg_rxn_ind_orig);
    rxnChanges.eqnNew     = model2.eqns(chg_rxn_ind_new);
    rxnChanges.lbOrig     = model1.lb(chg_rxn_ind_orig);
    rxnChanges.lbNew      = model2.lb(chg_rxn_ind_new);
    rxnChanges.ubOrig     = model1.ub(chg_rxn_ind_orig);
    rxnChanges.ubNew      = model2.ub(chg_rxn_ind_new);
    rxnChanges.grRuleOrig = model1.grRules(chg_rxn_ind_orig);
    rxnChanges.grRuleNew  = model2.grRules(chg_rxn_ind_new);
    
    rxnChanges.notes = repmat({''},size(rxnChanges.rxns));
    if ~isempty(changeNotes)
        [is_rxn,rxn_ind] = ismember(changeNotes(:,1),rxnChanges.rxns);
        rxnChanges.notes(rxn_ind(is_rxn)) = changeNotes(is_rxn,2);
    end
end


%% Assess changes in metabolite-related fields

% get indices mapping original model to new model mets
[is_present,met_ind] = ismember(model1.mets,model2.mets);
del_met = model1.mets(~is_present);  % record any deleted mets
nondel_met = model1.mets(is_present);  % get list of non-deleted mets

% check for added mets
add_met = setdiff(model2.mets,model1.mets);

% check which metabolites have a changed name, metFormula, or metCharge
ind = (1:length(model1.mets));
chg_name = ~arrayfun(@(i) isequal(model1.metNames(i),model2.metNames(met_ind(i))),ind(is_present));
chg_formula = ~arrayfun(@(i) isequal(model1.metFormulas(i),model2.metFormulas(met_ind(i))),ind(is_present));
chg_charge = model1.metCharges(is_present) ~= model2.metCharges(met_ind(is_present));

% collect unique list of all changed metabolites
model1.mets = [model1.mets; add_met];
chg_met = unique([chg_met; nondel_met(chg_name); nondel_met(chg_formula); nondel_met(chg_charge); del_met; add_met]);
[~,chg_met_ind_orig] = ismember(chg_met,model1.mets);
[~,chg_met_ind_new] = ismember(chg_met,model2.mets);

% add a dummy entry to model2 to handle deleted metabolites
chg_met_ind_new(chg_met_ind_new == 0) = length(model2.mets)+1;
model2.metNames(end+1) = {'(DELETED)'};
model2.metFormulas(end+1) = {''};
model2.metCharges(end+1) = 0;

% add dummy entries to model1 to handle added metabolites
model1.metNames(end+1:end+numel(add_met)) = {'(ADDED)'};
model1.metFormulas(end+1:end+numel(add_met)) = {''};
model1.metCharges(end+1:end+numel(add_met)) = 0;

% organize metabolite change information into structure
metChanges = {};
if ~isempty(chg_met)
    metChanges.mets        = model1.mets(chg_met_ind_orig);
    metChanges.nameOrig    = model1.metNames(chg_met_ind_orig);
    metChanges.nameNew     = model2.metNames(chg_met_ind_new);
    metChanges.formulaOrig = model1.metFormulas(chg_met_ind_orig);
    metChanges.formulaNew  = model2.metFormulas(chg_met_ind_new);
    metChanges.chargeOrig  = model1.metCharges(chg_met_ind_orig);
    metChanges.chargeNew   = model2.metCharges(chg_met_ind_new);
    
    metChanges.notes = repmat({''},size(metChanges.mets));
    if ~isempty(changeNotes)
        [is_met,met_ind] = ismember(changeNotes(:,1),metChanges.mets);
        metChanges.notes(met_ind(is_met)) = changeNotes(is_met,2);
    end
end


%% Organize combined results structure

modelChanges.rxns = rxnChanges;
modelChanges.mets = metChanges;



