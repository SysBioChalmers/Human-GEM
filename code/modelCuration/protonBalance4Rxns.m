function [newModel, newS, indRxnImbalance]=protonBalance4Rxns(model,protonMetId)
% 
%   protonBalance4Rxns aims to detect and rebalance the reactions that are
%   not balanced solely due to mismatched proton(s)
%
%   model             a model structure
%   protonMetId       compartment-free met id of proton
% 
%   newModel          an updated model structure
%   newS              an updated S matrix with balanced proton for some
%                     reactions
%   indRxnImbalance   index of reactions with imbalanced proton(s)
%
%   NOTE: The COBRA function checkMassChargeBalance is used for getting
%   imbalanced mass and charge values for each reaction. Rebalancing is 
%   restricted to reactions whose mets are all in the same compartment.
%   Only the S matrix is modified in the update model structure.
%
%   Usage: [newModel, newS, indRxnImbalance]=protonBalance4Rxns(model,protonMetId)
%


if nargin<2
    EM='Missing input arguments';
    disp(EM);
end

% get imbalanced mass and charge values using checkMassChargeBalance
[~,imbalancedMass,imbalancedCharge,~,~,~,~] = checkMassChargeBalance(model);

% log the index of modified reactions
indRxnImbalance = [];

% focus on the reactions with imbalanced mass
indImbalanceMass = find(~cellfun(@isempty, imbalancedMass));


fullS = full(model.S);
for i=1:length(indImbalanceMass)
    m = indImbalanceMass(i);
    protonDiff = regexprep(imbalancedMass{m},' H$','');
    
    % make sure that the imbalance is only contributed by proton, this is
    % determined if the mass differences equals the charge differences
    if isequal(str2double(protonDiff),imbalancedCharge(m))
        
        % here only RAVEN format model structure is allowed
        if ~isfield(model,'metComps')
            error('model has to be in RAVEN format with "metComps" field.');
        else
            checkComps=num2cell(model.metComps(find(model.S(:,m))));    
        end
        
        % check if all mets in a reaction are in the same compartment
        if isequal(checkComps{:})
            
            % compartment id is appended here, this only suits for HMR
            % met id so far
            proton=strcat(protonMetId,model.comps{checkComps{1}});
            % this part needs to be adjusted later for general usage
            
            protonMetIndex = find(strcmp(model.mets,proton));
            % when proton is NOT present in this reaction   
            if fullS(protonMetIndex,m) == 0
                fullS(protonMetIndex,m) = -1*imbalancedCharge(m);
            % when proton is present
            else
                fullS(protonMetIndex,m) = -1*imbalancedCharge(m) + fullS(protonMetIndex,m);
            end
            
            % record the ids of modified rxns
            indRxnImbalance = [indRxnImbalance;m];
        end
    end
end

% generating output
newModel = model;
newS = sparse(fullS);
newModel.S = newS;

end


