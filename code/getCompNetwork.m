function [compNetwork, I]=getCompNetwork(model,comp,includePartial)
% getCompNetwork
%   Gets the metabolic network for a specified compartment
%
% Input:
%   model           a model structure
%   comp            string with the compartment id
%   includePartial  if true, include reactions with metabolites partially
%                   present in the specified compartment (opt, default false)
%
% Output:
%   compNetwork     a model structure for the specified compartment
%   I               index of reactions in the specified compartment
%
% Usage: [compNetwork, I]=getCompNetwork(model,comp,includePartial)
%


if ischar(comp)
    comp={comp};
else
    error('Incorrect compartment id!');
end
if nargin<3
    includePartial=false;
end

% get the reaction list for the specified compartment
I=getRxnsInComp(model,comp,includePartial);

rxnToRemove=setdiff(transpose(1:length(model.rxns)), I);
compNetwork=removeReactions(model,rxnToRemove,true,true,true);

end
