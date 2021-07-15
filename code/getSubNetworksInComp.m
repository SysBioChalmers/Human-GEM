function subNetworks=getSubNetworksInComp(model,comp,metsToRemove,outFile,includePartial)
% getSubNetworksInComp
%	Check network connectivity in specified compartment and output
%   the sub-networks after excluding provided currency metabolites
%
% Input:
%   model           a model structure
%   comp            string with the compartment id
%   metsToRemove    cell array of met names to be excluded from the network
%   outFile         output results into a JSON format file (opt, default false)
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
% NOTE: The output file is "SubNetworks.json" that is in JSON format 
%
% Usage: subNetworks=getSubNetworksInComp(model,comp,metsToRemove,outFile,includePartial)
%


% handle input arguments
if nargin < 2
    error('Missing input!');
end
if nargin < 3
    metsToRemove={''};
else
    if ischar(metsToRemove)
        metsToRemove = (metsToRemove);
    end
end
if nargin < 4
    outFile=false;
end
if nargin < 5
    includePartial=false;
end


% get the sub-model of the specified compartment
compNetwork = getCompNetwork(model,comp,includePartial);


% refine the network by excluding any metabolites to be removed 
excludeMets = intersect(metsToRemove, compNetwork.metNames);
if ~isempty(excludeMets)
    reducedNetwork = removeMets(compNetwork,excludeMets,true,true);
else    
    fprintf('Note: the provided metabolites are not found in the model.\n');
    fprintf('      make sure they are metabolite names in cell array.\n');
    reducedNetwork = compNetwork;
end

% get the subnetworks
subGraphs = getAllSubGraphs(reducedNetwork);

% generate output
graphNum = size(subGraphs,2);
subNetworks.id = transpose(1:graphNum);

% write to a JSON file
if outFile
    fid = fopen('Subnetworks.json','w');
    fprintf(fid, '{\n');
end

% output rxns and subsystmes for each subnetwork
for i=1:graphNum
    index = find(subGraphs(:,i));
    [~, rxnList]=find(reducedNetwork.S(index,:));
    I = unique(rxnList);
    subNetworks.number(i,1) = length(I);
    subNetworks.rxns{i,1} = reducedNetwork.rxns(I);
    subNetworks.subSystems{i,1} = reducedNetwork.subSystems(I);

    if outFile
        fprintf(fid,['\t"' num2str(i) '":[\n']);
        writeSubNetworks(fid, reducedNetwork.rxns(I), reducedNetwork.subSystems(I));
        if i==graphNum
            fprintf(fid, '\t]\n');
        else
            fprintf(fid, '\t],\n');
        end
    end
    
end

% close file handle
if outFile
    fprintf(fid, '}\n');
    fclose(fid);
end

end


function writeSubNetworks(file, rxnList, subSystemList)
if ~isequal(numel(rxnList), numel(subSystemList))
    error('Input reaction and subSystems lists do not match!');
else
    for i=1:numel(rxnList)
        if i==numel(rxnList)
            fprintf(file,('\t\t["%s", "%s"]\n'),rxnList{i},subSystemList{i}{1});
        else
            fprintf(file,('\t\t["%s", "%s"],\n'),rxnList{i},subSystemList{i}{1});
        end
    end
end

end

