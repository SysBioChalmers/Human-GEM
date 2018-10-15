function rxnChanges = docRxnChanges(model1,model2,rxnNotes,fileName)
%docRxnChanges  Document changes to model reactions.
%
%   docRxnChanges looks for differences between two models, which should be
%   immediately before and after reaction(s) have been modified. It will
%   then report the original and new reaction equations and bounds, as well
%   as any user-supplied notes for each reaction.
%
% USAGE:
%
%   rxnChanges = docRxnChanges(model1,model2,rxnNotes,fileName);
%
% INPUT:
%
%   model1     Model structure prior to changing reaction(s).
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
%   fileName   (Optional) write the documented changes to a .tsv file with
%              the given name. If no name is provided, a file will not be 
%              written.
%
% OUTPUT:
%
%   rxnChanges  A cell structure containing the original and new reaction
%               properties, with the following fields:
%                 rxn       reaction ID
%                 eqnOrig   original reaction equation
%                 eqnNew    new reaction equation
%                 lbOrig    original lower bound
%                 lbNew     new lower bound
%                 ubOrig    original upper bound
%                 ubNew     new upper bound
%                 notes     notes associated with the change
%
% Jonathan Robinson, 2018-10-15


% handle input arguments
if nargin < 4
    fileName = [];
elseif ~contains(fileName,'.')
    % append with .tsv if no extension provided
    fileName = strcat(fileName,'.tsv');
end
if nargin < 3
    rxnNotes = [];
end


% get rxns from rxnNotes, and include them in the list of changed rxns
model1.rxnNotes = repmat({''},length(model1.rxns),1);  % intialize field
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
    [~,ind] = ismember(chg_rxn,model1.rxns);
    model1.rxnNotes(ind) = rxnNotes(:,2);
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

% check which reactions have a changed equation
ind = (1:length(model1.rxns));
chg_eqn = ~arrayfun(@(i) isequal(model1.eqns(i),model2.eqns(rxn_ind(i))),ind(is_present));

% check which reactions have a changed bound
chg_bnd = model1.lb(is_present) ~= model2.lb(rxn_ind(is_present)) | ...
          model1.ub(is_present) ~= model2.ub(rxn_ind(is_present));

% collect unique list of all changed reactions
chg_rxn = unique([chg_rxn; nondel_rxn(chg_eqn); nondel_rxn(chg_bnd); del_rxn]);
[~,chg_rxn_ind_orig] = ismember(chg_rxn,model1.rxns);
[~,chg_rxn_ind_new] = ismember(chg_rxn,model2.rxns);

% add a dummy entry to model2 to handle deleted reactions
chg_rxn_ind_new(chg_rxn_ind_new == 0) = length(model2.rxns)+1;
model2.eqns(end+1) = {'(DELETED)'};
model2.lb(end+1) = 0;
model2.ub(end+1) = 0;


% organize information into output structure
rxnChanges = {};
rxnChanges.rxns    = model1.rxns(chg_rxn_ind_orig);
rxnChanges.eqnOrig = model1.eqns(chg_rxn_ind_orig);
rxnChanges.eqnNew  = model2.eqns(chg_rxn_ind_new);
rxnChanges.lbOrig  = model1.lb(chg_rxn_ind_orig);
rxnChanges.lbNew   = model2.lb(chg_rxn_ind_new);
rxnChanges.ubOrig  = model1.ub(chg_rxn_ind_orig);
rxnChanges.ubNew   = model2.ub(chg_rxn_ind_new);
rxnChanges.notes   = model1.rxnNotes(chg_rxn_ind_orig);


% write to file, if specified
if ~isempty(fileName)
    % organize information for changed reactions
    docArray = [fieldnames(rxnChanges)';
                [rxnChanges.rxns, rxnChanges.eqnOrig, rxnChanges.eqnNew,...
                 num2cell([rxnChanges.lbOrig, rxnChanges.lbNew, ...
                 rxnChanges.ubOrig, rxnChanges.ubNew]), rxnChanges.notes]];
    % write to file
    writecell(docArray,fileName,true,'\t','%s\t%s\t%s\t%f\t%f\t%f\t%f\t%s\n');
end




