%
% FILE NAME:    curateMetFormula4PAPS.m
%
% PURPOSES: 1) Update the molecular formula of the PAPS (m02682) in the 
%              model based on the curation results in #58
%           2) Check if the updated molecular formula of PAPs improves the 
%              mass balancing status of the model
%           3) There are a total of 93 reactions are rebalanced in this 
%              curation through updating the forumulas of PAPS and balancing
%              the mismatched protons using function protonBalance4Rxns
%


%% 1. Load the model
if ~exist('ihuman','var')
    load('humanGEM.mat');  % version 0.8.1
end
ihuman_orig = ihuman;  % to keep track of changes made by the script
Eqns_before=constructEquations(ihuman);


%% 2. Compare the status of reaction balancing before and after changing
% the molecular formulas of PAPS 
[~,imbalancedMass_before,imbalancedCharge_before,~,~,~,~] = checkMassChargeBalance(ihuman);

% get the index PAPs and update their formulas
metsInd = find(startsWith(ihuman.mets, 'm02682'));
ihuman.metFormulas(metsInd)={'C10H11N5O13P2S'};

% check balancing status again
[~,imbalancedMass_after,imbalancedCharge_after,~,~,~,~] = checkMassChargeBalance(ihuman);

% the charge status remains the same
if isequal(imbalancedCharge_before, imbalancedCharge_after)
    disp('The formula change has no effect to the status of charge balancing.');
end

% find out the reactions that have PAPS involved (i.e. with changed mass status)
check = strcmp(imbalancedMass_before, imbalancedMass_after);
index_changedMass = find(check == 0);  % a total of 100 rxns with changes mass status
balanceStatus_before = imbalancedMass_before(index_changedMass);
balanceStatus_after  = imbalancedMass_after(index_changedMass);


%% 3. Group the reactions

% Through analyzing the status before and after updating PAPS formulas,
% these affected reactions can be classifed into 3 group:

% groupA: 28 reactions are now mass balanced
tmpInd=find(cellfun(@isempty,imbalancedMass_after(index_changedMass)));
ind_GroupA = index_changedMass(tmpInd);

% groupB: 63 reactions are not balanced solely due to mismatch of H+ ion
tmp=find(~cellfun(@isempty,regexp(imbalancedMass_after(index_changedMass),'^-\d+\sH$')));
ind_GroupB = index_changedMass(tmp);

% groupC: 9	reactions are not balanced due to other reasons, and thus
% documented for future refinement
ind_GroupC = setdiff(index_changedMass,[ind_GroupA;ind_GroupB]);
groupC(:,1)=ihuman.rxns(ind_GroupC);
groupC(:,2)=Eqns_before(ind_GroupC);
groupC(:,3)=imbalancedMass_after(ind_GroupC);
% Output group C reactions for manual curation


%% 4. fix the reactions in groupB by using proton balancing function
protonMetId = 'm02039';
model=ihuman;

% get the updated model by function protonBalance4Rxns.m
[new_model, ~, IndBalancedRxns] = protonBalance4Rxns(model, protonMetId);
length(IndBalancedRxns)
% a total 65 reactions are balanced, but it supposed to be 63!

% analyzing results for identifying exceptions
Eqns_after=constructEquations(new_model);
outliers = setdiff(IndBalancedRxns, index_changedMass);
ihuman.rxns(outliers)

% By checking these two outlier reactions (HMR_8824 and HMR_8825) before
% and after the conversion, it seems the rebalancing is a proper fixation
%
% Before:
% HMR_8824: 1-phosphatidyl-1D-myo-inositol-5-phosphate[c] + ATP[c] => ADP[c] + phosphatidylinositol-3,5-bisphosphate[c]
% HMR_8825: 1-phosphatidyl-1D-myo-inositol-5-phosphate[r] + ATP[r] => ADP[r] + phosphatidylinositol-3,5-bisphosphate[r]
% After:
% HMR_8824: 1-phosphatidyl-1D-myo-inositol-5-phosphate[c] + ATP[c] => ADP[c] + H+[c] + phosphatidylinositol-3,5-bisphosphate[c]
% HMR_8825: 1-phosphatidyl-1D-myo-inositol-5-phosphate[r] + ATP[r] => ADP[r] + H+[r] + phosphatidylinositol-3,5-bisphosphate[r]
%

% log information in rxnNotes array
rxnNotes = [ihuman.rxns(ind_GroupB), repmat({'proton balancing for the reactions with updated PAPS formulas'},length(ind_GroupB),1)];
rxnNotes = [rxnNotes; [ihuman.rxns(outliers), repmat({'proton balancing for the reactions with corrected formulas by #52'},length(outliers),1)]];


%% 5. Evaluate balancing results by comparing current status to the initial one
[~,imbalancedMass_balanced,imbalancedCharge_balanced,~,~,~,~] = checkMassChargeBalance(new_model);

% check charge status
ind_chargeChange = find(imbalancedCharge_before~=imbalancedCharge_balanced);
if isequal(ind_chargeChange, IndBalancedRxns) &&...        % charges are only changed by protonBalance4Rxns
    all(imbalancedCharge_balanced(ind_chargeChange) == 0)  % and they are now balanced
    fprintf('\nThe charge rebalanced reactions are just the ones fixed by protonBalance4Rxns.\n\n');
end

% check mass status
ind_massChange = find(strcmp(imbalancedMass_before, imbalancedMass_balanced) == 0);
if isequal(ind_massChange, sort([index_changedMass; outliers])) &&...  % mass are changed from changing PAPS formula and 2 outliers
    isequal(ind_massChange(getNonEmptyList(imbalancedMass_balanced(ind_massChange))), ind_GroupC) % only group C remain imbalanced in mass
    fprintf('\nThe mass rebalanced reactions are only from group A and B.\n\n');
end


%% 6. Document reaction changes, clear intermediate values and save results

ihuman = new_model;
rxnChanges = docRxnChanges(ihuman_orig,ihuman,rxnNotes);
writeRxnChanges(rxnChanges,'curateFormulas4PAPS_rxnChanges',1);
movefile('curateFormulas4PAPS_rxnChanges.tsv','../../ComplementaryData/modelCuration/');

% clear intermediate vars
clearvars -except ihuman groupC

% save model file
save('../../model/Human-GEM.mat','ihuman');
