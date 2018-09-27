function [newModel, newS]=protonBalance4Rxns(model,protonMetId)
% 
%   protonBalance4Rxns aims to detect and rebalance the reactions that are
%   not balanced solely due to mismatched proton(s)
%
%   model           a model structure
%   protonMetId     compartment-free met id of proton
% 
%   newModel        an updated model structure
%   newS            an updated S matrix with balanced proton for some
%                   reactions
%
%   NOTE: The COBRA function checkMassChargeBalance is used for getting
%   imbalanced mass and charge values for each reaction. Rebalancing is 
%   restricted to reactions whose mets are all in the same compartment.
%   Only the S matrix is modified in the update model structure.
%
%   Usage: [newModel, newS]=protonBalance4Rxns(model,protonMetId)
%
%   Hao Wang, 2018-09-27
%

% handel input

% get imbalanced mass and charge values using checkMassChargeBalance
[~,imbalancedMass,imbalancedCharge,~,~,~,~] = checkMassChargeBalance(model);

% intermediate results
indRxnImbalanceProton = [];

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
            % met id so far and need to be adjusted later for general usage
            proton=strcat(protonMetId,model.comps{checkComps{1}});
            
            protonMetIndex = find(strcmp(model.mets,proton));
            % when proton is NOT present in this reaction   
            if fullS(protonMetIndex,m) == 0
                fullS(protonMetIndex,m) = -1*imbalancedCharge(m);
            % when proton is present
            else
                fullS(protonMetIndex,m) = -1*imbalancedCharge(m) + fullS(protonMetIndex,m);
            end
            
            % record the rxn ids
            indRxnImbalanceProton = [indRxnImbalanceProton;m];
        end
    end
end
new_S = sparse(fullS);

% generating output
newModel = model;
newModel.S = newS;

end


