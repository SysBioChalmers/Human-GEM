function [subNetworks]=getSubNetworksInComp(model,comp,metsToRemove,includePartial)
% getSubNetworksInComp
%	Check network connectivity in specified compartment and output
%   the sub-networks in JSON fommat for map generation purpose
%
% Input:
%   model           a model structure
%   comp            string with the compartment id
%   metsToRemove    cell array of mets to be excluded from the network
%   includePartial  if true, include reactions with metabolites partially
%                   present in the specified compartment (opt, default false)
%
% Output:
%   subNetworks     a structure for the sub-networks
%      id           the id for identifying the sub-networks
%      number       number of reactions in each sub-network
%      rxns         cell array of reaction identifiers in each sub-network
%      subSystems   cell array of subSystems for reactions in each sub-network
%
%
% NOTE: The output filename "SubNetworks.json" and in defined JSON format 
%
% Usage: [subNetworks]=getSubNetworksInComp(model,comp,metsToRemove,includePartial)
%


% handle input arguments
if nargin < 2
    error('Missing input!');
end
if nargin < 3
    metsToRemove={''};
end
if nargin < 4
    includePartial=false;
end

% get the subnetworks in specified compartment
compNetwork = getCompNetwork(model,comp,includePartial);

% refine the network by excluding any metabolites to be removed and
% identifiy any ractions that have only currency metabolites
excludeMets = intersect(metsToRemove, compNetwork.mets);
if isempty(excludeMets)
    frintf('The provided metabolites are not found in the model\n');
    reducedNetwork = compNetwork;
    currencyRxns = '';
else
    reducedNetwork = removeMets(compNetwork,excludeMets,0,1);
    currencyRxns = setdiff(compNetwork.rxns, reducedNetwork.rxns);
    [~, currencyRxnInd] = ismember(currencyRxns, model.rxns);
end

% get the sub-graphs
subGraphs = getAllSubGraphs(reducedNetwork);

% generate output
graphNum = size(subGraphs,2);
subNetworks.id = transpose(1:graphNum);

% write in JSON format
fid = fopen('Subnetworks.json','w');
fprintf(fid, '{\n');

for i=1:graphNum
    index = find(subGraphs(:,i));
    [~, rxnList]=find(reducedNetwork.S(index,:));
    I = unique(rxnList);
    subNetworks.number(i,1) = length(I);
    subNetworks.rxns{i,1} = reducedNetwork.rxns(I);
    subNetworks.subSystems{i,1} = reducedNetwork.subSystems(I);
    
    fprintf(fid,['\t"' num2str(i) '":[\n']);
    writeSubNetworks(fid, reducedNetwork.rxns(I), reducedNetwork.subSystems(I));
    if i==graphNum && isempty(currencyRxns)
        fprintf(fid, '\t]\n');
    else
        fprintf(fid, '\t],\n');
    end
end

if ~isempty(currencyRxns)
    fprintf(fid,['\t"currencyRxns":[\n']);
    writeSubNetworks(fid, currencyRxns, model.subSystems(currencyRxnInd));
    fprintf(fid, '\t]\n');
end
fprintf(fid, '}\n');
fclose(fid);

end


function writeSubNetworks(file, rxnList, subSystemList)
if ~isequal(numel(rxnList), numel(subSystemList))
    error('Input reaction and subSystems lists do not match!');
else
    for i=1:numel(rxnList)
        if i==numel(rxnList)
            fprintf(file,['\t\t["' rxnList{i} '", "' subSystemList{i} '"]\n']);
        else
            fprintf(file,['\t\t["' rxnList{i} '", "' subSystemList{i} '"],\n']);
        end
    end
end

end