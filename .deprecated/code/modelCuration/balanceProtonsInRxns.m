%
% FILE NAME:    balanceProtonsInRxns.m
%
% PURPOSE: This script is to balance humanGEM (v0.3.1) targeting for the
% reactions that are imbalanced solely due to mismatch of proton(s)
%


%% Load model
load('humanGEM.mat');  % v0.3.1



%% Convert demand and sink type pseudoreactions to RAVEN fromat
model=addBoundaryMets(ihuman,false);
% A total of 49 boundary version metabolites were added to 253 reactions
% that are either demand or sink types, for reaction balancing and
% complying with RAVEN format



%% Fix imbalance reactions solely caused by mismatch of proton(s)

% compartmen-free met id of proton in humanGEM
protonMetId = 'm02039';

% get the updated model by function protonBalance4Rxns.m
[new_model, ~, modifiedRxns] = protonBalance4Rxns(model, protonMetId);
length(modifiedRxns)  % A total of 1419 reactons are rebalanced



%% Save updated model file
ihuman = new_model;
save('humanGEM.mat','ihuman');

