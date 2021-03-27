function [grRules,new_genes] = integrateGeneRules(HMR,iHsa,Recon3D,rxnHMR2Recon3D)
%integrateGeneRules  Combine grRules from HMR, iHsa, and Recon3D models.
%
%
% USAGE:
%
%   grRules = integrateGeneRules(HMR,iHsa,Recon3D,rxnHMR2Recon3D)
%
% INPUT:
%
%   HMR         Human Metabolic Reaction (HMR) model structure, or some
%               other model structure with which the grRules will be
%               integrated.
%
%   iHsa        iHsa (EM Blais, JA Papin, et al. 2017) model structure.
%
%   Recon3D     Recon3D model structure.
%
%   rxnHMR2Recon3D    Cell array to convert between HMR rxn IDs (first 
%                     column) and Recon3D rxn IDs (second column).
%
% OUTPUT:
%
%   grRules     An updated grRules vector, where grRules from iHsa and
%               Recon3D have been incorporated into the grRules of HMR.
%
%   new_genes   A cell array listing new genes (if any) that were added to
%               each grRule, and were already present in some other
%               grRule(s) in the original model.
%


% Note: Loaded rxnHMR2Recon3D array using following commands:
% tmp = readtable('/Users/jonrob/Documents/PostDoc/HMR3/CurationFiles/reactions/countGeneNum4AssociatedRxns_Hao_20180726.xlsx');
% rxnHMR2Recon3D = [tmp.HMR,tmp.Recon3D];


% handle input arguments
if isempty(HMR)
    load('ComplementaryData/HMR2/HMRdatabase2_02.mat');  % loads as variable "ihuman"
    HMR = ihuman;
end
if isempty(iHsa)
    load('ComplementaryData/iHsa/iHsa.mat');  % loads as variable "iHsa"
end
if isempty(Recon3D)
    load('ComplementaryData/Recon3D/Recon3D_301.mat');  % loads as variable "Recon3D"
end

% initialize output
grRules = HMR.grRules;
new_genes = repmat({''},size(grRules));

% clean HMR.grRules if not yet done
fprintf('Cleaning HMR grRules... ');
HMR.grRules = cleanModelGeneRules(HMR.grRules);
fprintf('Done.\n');

% check if iHsa contains HMR rxn associations
if ~isfield(iHsa,'rxnHMRID')
    iHsa = addHMRrxnIDsToiHsa(iHsa);
end

% add iHsa rxn IDs to HMR model
HMR.rxniHsaID = repmat({''},size(HMR.rxns));
[hasmatch,ind] = ismember(HMR.rxns,iHsa.rxnHMRID);
HMR.rxniHsaID(hasmatch) = iHsa.rxns(ind(hasmatch));

% convert iHsa grRules to Ensembl IDs
ihsa_grRule = repmat({''},size(HMR.rxns));
ihsa_grRule(hasmatch) = iHsa.grRules(ind(hasmatch));
ihsa_grRule_ensg = translateGrRules(ihsa_grRule,'ENSG');


% extract grRules from Recon3D, translate, and align with HMR rules
[hasmatch,ind] = ismember(HMR.rxns,rxnHMR2Recon3D(:,1));
r3_rxn = repmat({''},size(HMR.rxns));
r3_rxn(hasmatch) = rxnHMR2Recon3D(ind(hasmatch),2);

[hasmatch,ind] = ismember(r3_rxn,Recon3D.rxns);
hasmatch(cellfun(@isempty,r3_rxn)) = false;
r3_grRule = repmat({''},size(HMR.rxns));
r3_grRule(hasmatch) = Recon3D.grRules(ind(hasmatch));
r3_grRule_ensg = translateGrRules(r3_grRule,'ENSG');


% get list of genes in each grRule
hmr_genes = genesInRxn(HMR.grRules);
ihsa_genes = genesInRxn(ihsa_grRule_ensg);
r3_genes = genesInRxn(r3_grRule_ensg);

% get list of all genes among all HMR rules
all_hmr_genes = unique(horzcat(hmr_genes{:}))';

% compare and integrate grRules
for i = 1:length(grRules)
    if isempty(hmr_genes{i})
        if isempty(ihsa_genes{i})
            if ~isempty(r3_genes{i})
                grRules{i} = r3_grRule_ensg{i};
                new_genes{i} = strjoin(intersect(r3_genes{i},all_hmr_genes),'; ');
            end
        else
            grRules{i} = ihsa_grRule_ensg{i};
            new_genes{i} = strjoin(ihsa_genes{i},'; ');
        end
    else
        if all(ismember(hmr_genes{i},ihsa_genes{i}))
            grRules{i} = ihsa_grRule_ensg{i};
            new_genes{i} = strjoin(intersect(ihsa_genes{i}(~ismember(ihsa_genes{i},hmr_genes{i})),all_hmr_genes),'; ');
        elseif all(ismember(hmr_genes{i},r3_genes{i}))
            grRules{i} = r3_grRule_ensg{i};
            new_genes{i} = strjoin(intersect(r3_genes{i}(~ismember(r3_genes{i},hmr_genes{i})),all_hmr_genes),'; ');
        end
    end
end



end  % function end



% function to get the list of unique genes in each grRule
function rxn_genes = genesInRxn(grRules)
    
% convert AND and OR to & and |, respectively
grRules = regexprep(grRules,' or ','|');
grRules = regexprep(grRules,' and ','&');

% identify genes associated with each reaction
rxn_genes = cellfun(@(r) unique(regexp(r,'[^&|\(\) ]+','match')),grRules,'UniformOutput',false);

end





