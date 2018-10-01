%
% FILE NAME:    balance_RGroup_CoA_rxns.m
% 
% DATE CREATED: 2018-09-27
%     MODIFIED: 2018-09-27
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script to balance reactions involving the "R group Coenzyme A"
%          metabolite.
%
%          A set of reactions imported from Recon3D involve metabolites
%          that are various forms of "R Group"s, which are used in pooling
%          reactions in e.g. lipid metabolism. However, some of these
%          reactions have improperly balanced their CoA component.
%          Therefore, additional equivalents of CoA need to be added to
%          some of these reactions in order to maintain mass balance.
%
%          For example, consider the following three reactions:
%
%           ARTFR207: eicosanoyl-CoA[c] => 1.25 R Group 2 Coenzyme A[c]
%              RTOT2: R Group 2 Coenzyme A[c] => R Total Coenzyme A[c]
%           ARTCOAL1: H2O[c] + R Total Coenzyme A[c] => CoA[c] + H+[c] + R Total[c]
%           
%          Combining these reactions consumes one equivalent of CoA, but 
%          generates 1.25 equivalents. Therefore, the first reaction needs
%          to be balanced by adding CoA consumption:
%
%            eicosanoyl-CoA[c] + 0.25 CoA => 1.25 R Group 2 Coenzyme A[c]
%
%          This script locates and corrects all cases of such imbalance in
%          the model.
%


% load model if not already present
if ~exist('ihuman','var')
    load('humanGEM.mat');
end

% identify all reactions requiring CoA re-balancing
met_ind = find(startsWith(ihuman.metNames,'R Group '));
rxn_ind = find(any(ihuman.S(met_ind,:) ~= 0 & abs(ihuman.S(met_ind,:)) ~= 1,1)');
met_ind(all(ihuman.S(met_ind,rxn_ind) == 0,2)) = [];  % discard mets that were not imbalanced

% get the stoich coeff for the "R Group # Coenzyme" met in each rxn
R_coeffs = ihuman.S(met_ind,rxn_ind);
if any(sum(R_coeffs ~= 0,1) ~= 1)
    % This is a check to make sure that there is only one of these types of
    % metabolites in each reaction being fixed. If this error occurs, those
    % reactions need to be corrected manually.
    error('Some reactions contain multiple R Group Coenzyme mets!');
end
R_coeffs = full(sum(R_coeffs,1)');  % sum the coeffs to get the only non-zero coeff for each reaction

% calculate the appropriate amount of CoA needed to balance each reaction
CoA_coeffs = 1 - R_coeffs;

% use the same compartment as the R Group
R_comps = unique(ihuman.comps(ihuman.metComps(met_ind)));
if length(R_comps) > 1
    % the current way in which the script is written does not allow for
    % different compartments among the imbalanced R-groups
    error('Some mets were in different compartments - this script will not work for these cases.');
end

% get CoA index
CoA_ind = getIndexes(ihuman,strcat('CoA[',R_comps{1},']'),'metscomps');

% update stoich matrix
ihuman.S(CoA_ind,rxn_ind) = CoA_coeffs;
fprintf('\nA total of %u reactions were re-balanced for CoA.\n\n',length(rxn_ind));


% clear intermediate variables
clear('CoA_coeffs','CoA_ind','met_ind','R_coeffs','R_comps','rxn_ind');




