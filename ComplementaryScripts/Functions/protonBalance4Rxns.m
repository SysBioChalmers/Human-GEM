function [new_model, new_S]=protonBalance4Rxns(model,protonMetId)
% 
%   protonBalance4Rxns aims to detect and rebalance the reactions that are
%   not balanced solely due to mismatched proton(s)
%
%   model           a model structure
%
%   protonMetId     compartment-free met id of proton
%
%   NOTE: The COBRA function checkMassChargeBalance is used for getting
%   imbalanced mass and charge values for each reaction. Rebalancing is 
%   restricted to reactions whose mets are all in the same compartment
%
%   Usage: [new_model, new_S]=protonBalance4Rxns(model,protonMetId)
%
%   Hao Wang, 2018-09-26
%

% handel input

% get imbalanced mass and charge values using checkMassChargeBalance
[~,imbalancedMass,imbalancedCharge,~,~,~,~] = checkMassChargeBalance(model);

% intermediate results
rxnInd = [];

% focus on the reactions with imbalanced mass
indImbalanceMass = getNonEmptyList(imbalancedMass);

fullS = full(model.S);
for i=1:length(indImbalanceMass)
    m = indImbalanceMass(i);
    protonCoeff = regexprep(imbalancedMass{m},' H$','');  % imbalanced mass
    
    % make sure that the imbalance is only contributed by proton
    % if the mass differences equals the charge differences
    if isequal(str2num(protonCoeff),imbalancedCharge(m))
        
        % check if all mets in a reaction are in the same compartment
        % here only RAVEN format model structure is considered
        checkComps=num2cell(model.metComps(find(model.S(:,m))));
        if isequal(checkComps{:})
            
            % compartment id is appended here, this only suits for HMR
            % met id so far and need to be adjusted for general usage
            proton=strcat(protonMetId,model.comps{checkComps{1}});
            
            % make sure proton is not present in this reaction            
            protonMetIndex = find(strcmp(model.mets,proton));
            if fullS(protonMetIndex,m) == 0
                fullS(protonMetIndex,m) = -1*imbalancedCharge(m);
            end
        end
    end
end

new_S = sparse(fullS);

% generating output
new_model = model;
new_model.S = new_S;

end


