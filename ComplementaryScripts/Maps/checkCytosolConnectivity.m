%
% FILE NAME:    checkCytosolConnectivity.m
%
% DATE CREATED: 2018-11-06
%     MODIFIED: 2018-11-08
%
% PROGRAMMER:   Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
%
% PURPOSE: check network connectivity in Cytosol compartment, and output
%          the sub-graphs in JSON fommat
%

%% Generate the network of Cytosol compartment
load('humanGEM.mat');  %v0.5.2
cytosolNetwork = getCompNetwork(ihuman, 'c');

% load currency metabolites
fid = fopen('currencyMets.tsv','r');
input=textscan(fid,'%s');
fclose(fid);
currencyMets = input{1};

% refine the network by excluding currency metabolites
metsToRemove = intersect(currencyMets, cytosolNetwork.mets);
model = removeMets(cytosolNetwork,metsToRemove,0,1);

% identifiy ractions that have only currency metabolites
currencyRxns = setdiff(cytosolNetwork.rxns, model.rxns);

% get the sub-graphs
subGraphs = getAllSubGraphs(model);

% generate output in JSON format
fid = fopen('cytosolSubNetworks.txt','w');
fprintf(fid, '{\n');
for i=1:size(subGraphs,2)
    subNetworks.id(i,1) = i;
    fprintf(fid,['\t"' num2str(i) '":[\n']);
    index = find(subGraphs(:,i));
    [~, rxnList]=find(model.S(index,:));
    subNetworks.num(i,1) = length(unique(rxnList));
    subNetworks.rxnList{i,1} = unique(rxnList);
    for j=1:numel(subNetworks.rxnList{i,1})
        ind = subNetworks.rxnList{i,1}(j);
        if j==numel(subNetworks.rxnList{i,1})
            fprintf(fid,['\t\t["' model.rxns{ind} '", "' model.subSystems{ind} '"]\n']);
        else
            fprintf(fid,['\t\t["' model.rxns{ind} '", "' model.subSystems{ind} '"],\n']);
        end
    end
    if i==graphNum
        fprintf(fid, '\t]\n');
    else
        fprintf(fid, '\t],\n');
    end
end
fprintf(fid, '}\n');
fclose(fid);

