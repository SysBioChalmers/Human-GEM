function [subNetworks]=getSubNetworksInComp(model,comp,metsToRemove,includePartial)
% getSubNetworksInComp
%	Check network connectivity in specified compartment and output
%   the sub-graphs in JSON fommat for map generation purpose
%
%   Input:
%   model           a model structure
%   comp            string with the compartment id
%   metsToRemove    cell array of mets to be excluded from the network
%   includePartial  if true, include reactions with metabolites partially
%                   present in the specified compartment (opt, default false)
%
%	NOTE: The output filename is "SubNetworks.json" in defined JSON format 
%
%	Usage: [subNetworks]=getSubNetworksInComp(model,comp,metsToRemove,includePartial)
%
%	Hao Wang, 2018-11-12
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

% refine the network by excluding metabolites metsToRemove
excludeMets = intersect(metsToRemove, compNetwork.mets);
if isempty(excludeMets)
    frintf('The provided metabolites are not found in the model\n');
    reducedNetwork = compNetwork;
    currencyRxns = '';
else
    reducedNetwork = removeMets(compNetwork,excludeMets,0,1);

    % identifiy ractions that have only currency metabolites
    currencyRxns = setdiff(compNetwork.rxns, reducedNetwork.rxns);
    [~, currencyRxnInd] = ismember(currencyRxns, model.rxns);
end

% get the sub-graphs
subGraphs = getAllSubGraphs(reducedNetwork);
num = size(subGraphs,2);
subNetworks.id = transpose(1:num);

% generate output in JSON format
fid = fopen('Subnetworks.json','w');
fprintf(fid, '{\n');

for i=1:num
    fprintf(fid,['\t"' num2str(i) '":[\n']);
    index = find(subGraphs(:,i));
    [~, rxnList]=find(reducedNetwork.S(index,:));
    subNetworks.number(i,1) = length(unique(rxnList));
    subNetworks.rxnList{i,1} = unique(rxnList);
    for j=1:numel(subNetworks.rxnList{i,1})
        ind = subNetworks.rxnList{i,1}(j);
        if j==numel(subNetworks.rxnList{i,1})
            fprintf(fid,['\t\t["' reducedNetwork.rxns{ind} '", "' reducedNetwork.subSystems{ind} '"]\n']);
        else
            fprintf(fid,['\t\t["' reducedNetwork.rxns{ind} '", "' reducedNetwork.subSystems{ind} '"],\n']);
        end
    end
    if i==num && isempty(currencyRxns)
        fprintf(fid, '\t]\n');
    else
        fprintf(fid, '\t],\n');
    end
end

fprintf(fid,['\t"currencyRxns":[\n']);
for k=1:numel(currencyRxns)
    if k==numel(currencyRxns)
        fprintf(fid,['\t\t["' currencyRxns{k} '", "' model.subSystems{currencyRxnInd(k)} '"]\n']);
    else
        fprintf(fid,['\t\t["' currencyRxns{k} '", "' model.subSystems{currencyRxnInd(k)} '"],\n']);
    end
end
fprintf(fid, '\t]\n');

fprintf(fid, '}\n');

fclose(fid);


