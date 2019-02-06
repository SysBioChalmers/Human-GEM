%
% FILE NAME:    curateMetFormula4PAPS.m
%
% DATE STARTED: 2019-01-16
%     MODIFIED: 2019-02-06
%
% PROGRAMMERS:  Hao Wang, Pinar Kocabas 
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
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
load('humanGEM.mat'); 
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
    disp('The formula change has no effect to the status of charge balancing');
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
indGroupA=ind_changedMass(tmpInd);

% groupB: 63 reactions are not balanced solely due to mismatch of H+ ion
tmp=find(~cellfun(@isempty,regexp(imbalancedMass_after(ind_changedMass),'^-\d+\sH$')));
indGroupB=index_changedMass(tmp);

% groupC: 9	reactions are not balanced due to other reasons, and should be
% documented for future refinement
indGroupC=setdiff(index_changedMass,[indexA;indexB]);


%% 4. fix the reactions in groupB by using proton balancing function
protonMetId = 'm02039';
model=ihuman;

% get the updated model by function protonBalance4Rxns.m
[new_model, ~, IndBalancedRxns] = protonBalance4Rxns(model, protonMetId);
length(IndBalancedRxns)
% a total 65 reactions are balanced, but it supposed to be 63!

% analyzing results for identifying exceptions
Eqns_after=constructEquations(new_model);
outliers = setdiff(IndBalancedRxns, ind_changedMass);
ihuman.rxns(outliers)

% By checking these two outlier reactions before and after the conversion, 
% it seems the rebalancing is a proper fixation
%
% Before:
% 1-phosphatidyl-1D-myo-inositol-5-phosphate[c] + ATP[c] => ADP[c] + phosphatidylinositol-3,5-bisphosphate[c]
% 1-phosphatidyl-1D-myo-inositol-5-phosphate[r] + ATP[r] => ADP[r] + phosphatidylinositol-3,5-bisphosphate[r]
% After:
% 1-phosphatidyl-1D-myo-inositol-5-phosphate[c] + ATP[c] => ADP[c] + H+[c] + phosphatidylinositol-3,5-bisphosphate[c]
% 1-phosphatidyl-1D-myo-inositol-5-phosphate[r] + ATP[r] => ADP[r] + H+[r] + phosphatidylinositol-3,5-bisphosphate[r]
%


%% 5. Save updated model
ihuman = new_model;
save('humanGEM.mat','ihuman');

