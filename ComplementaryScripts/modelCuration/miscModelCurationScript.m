%
% FILE NAME:    miscModelCurationScript.m
% 
% DATE CREATED: 2018-09-19
%     MODIFIED: 2018-09-19
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script for performing many miscellaneous, often small/minor
%          curations, updates, and corrections to HumanGEM.
%



%% load humanGEM
load('ModelFiles/mat/humanGEM.mat');


%% Correct the size of some non-standard model fields
% Some of the non-standard model fields (e.g., "rxnKEGGID") are not the
% correct size, and thus may cause confusion and loss of alignment with the
% list of reactions (model.rxns) to which they correspond. Therefore, they
% will be padded with blank entries, just to make them the proper
% dimensions. These will be updated properly elsewhere.
f = {'rxnKEGGID';'rxnEHMNID';'rxnBiGGID';'rxnHepatoNET1ID';'rxnREACTOMEID'};
for i = 1:length(f)
    ihuman.(f{i})(end+1:length(ihuman.rxns)) = {''};
end



%% Treatment of ubiquinone and FAD+
% The model essentially treats FAD/FADH2 the same as ubiquinone/ubiquinol,
% as exhibited by the presence of the following reaction:
%  HMR_6911: FADH2[m] + ubiquinone[m] <=> FAD[m] + ubiquinol[m]
%
% Some reactions in the model are duplicated such that one version of the
% reaction uses FAD, whereas the other uses ubiquinone. These are
% effectively identical reactions, and one should be removed (unless
% evidence suggests otherwise). 

% (Not yet implemented)


%% Update of EC number assignment

% The following reaction needs to have its EC number updated:
%  HMR_4365: 2-phospho-D-glycerate[c] <=> 3-phospho-D-glycerate[c]
%   EC orig: EC:5.4.2.1;EC:5.4.2.4
%   EC new:  EC:5.4.2.11
ihuman.eccodes(ismember(ihuman.rxns,'HMR_4365')) = {'EC:5.4.2.11'};



%% Update of grRules

% The reaction for Complex I (HMR_6921) should include ENSG00000283447
% (NDUFS1) as part of its grRule.
[~,rxn_ind] = ismember('HMR_6921',ihuman.rxns);
if ~contains(ihuman.grRules(rxn_ind),'ENSG00000283447')
    if contains(ihuman.grRules{rxn_ind},' and ')
        error('grRule contains AND expressions. The gene should be added manually to the rule.');
    end
    
    % add gene to end of the rxn grRule
    ihuman.grRules(rxn_ind) = strcat(ihuman.grRules(rxn_ind),' or ENSG00000283447');
    
    % find gene index; if it doesn't exist, add it to the model
    [~,gene_ind] = ismember('ENSG00000283447',ihuman.genes);
    if gene_ind == 0
        ihuman.genes(end + 1) = {'ENSG00000283447'};
        gene_ind = length(ihuman.genes);
    end
    
    % add gene to rxnGeneMat
    ihuman.rxnGeneMat(rxn_ind,gene_ind) = 1;
end


% The reaction for Complex III (HMR_6918) should include ENSG00000284493
% (UQCRC2) as part of its grRule.
[~,rxn_ind] = ismember('HMR_6918',ihuman.rxns);
if ~contains(ihuman.grRules(rxn_ind),'ENSG00000284493')
    if contains(ihuman.grRules{rxn_ind},' and ')
        error('grRule contains AND expressions. The gene should be added manually to the rule.');
    end
    
    % add gene to end of the rxn grRule
    ihuman.grRules(rxn_ind) = strcat(ihuman.grRules(rxn_ind),' or ENSG00000284493');
    
    % find gene index; if it doesn't exist, add it to the model
    [~,gene_ind] = ismember('ENSG00000284493',ihuman.genes);
    if gene_ind == 0
        ihuman.genes(end + 1) = {'ENSG00000284493'};
        gene_ind = length(ihuman.genes);
    end
    
    % add gene to rxnGeneMat
    ihuman.rxnGeneMat(rxn_ind,gene_ind) = 1;
end


% The lactate dehydrogenase rxns (HMR_4281, HMR_4388, and HMR_4280) have
% differing grRules, but they should be the same.
%  HMR_4281: H+[p] + NADH[p] + pyruvate[p] <=> L-lactate[p] + NAD+[p]
%  HMR_4388: H+[c] + NADH[c] + pyruvate[c] <=> L-lactate[c] + NAD+[c]
%  HMR_4280: H+[m] + NADH[m] + pyruvate[m] <=> L-lactate[m] + NAD+[m]
[~,rxn_ind] = ismember({'HMR_4281';'HMR_4388';'HMR_4280'},ihuman.rxns);
if length(unique(ihuman.grRules(rxn_ind))) > 1  % check that the grRules differ
    
    if any(contains(ihuman.grRules(rxn_ind),' and '))
        error('grRule(s) contains AND expressions, and should be updated manually.');
    end
    
    % combine rxnGeneMat rows and grRules for the rxns
    ihuman.rxnGeneMat(rxn_ind,:) = repmat(max(ihuman.rxnGeneMat(rxn_ind,:),[],1),length(rxn_ind),1);
    gene_inds = find(ihuman.rxnGeneMat(rxn_ind(1),:) == 1);
    ihuman.grRules(rxn_ind) = join(ihuman.genes(gene_inds),' or ');

end


% The reaction HMR_4379 is associated with ENSG00000160226, which encodes
% a protein, but the function is somewhat unclear, or related to DNA
% damage/repair, and not the reaction with which it is currently
% associated:
%  HMR_4379: ATP[c] + fructose-6-phosphate[c] => ADP[c] + fructose-1,6-bisphosphate[c]
% Therefore, the gene should be removed from the associated grRule.
[~,rxn_ind] = ismember('HMR_4379',ihuman.rxns);
if contains(ihuman.grRules(rxn_ind),'ENSG00000160226')
    if contains(ihuman.grRules(rxn_ind),' and ')
        error('grRule(s) contains AND expressions, and should be updated manually.');
    end
    [~,gene_ind] = ismember('ENSG00000160226',ihuman.genes);
    ihuman.rxnGeneMat(rxn_ind,gene_ind) = 0;
    ihuman.grRules(rxn_ind) = join(ihuman.genes(ihuman.rxnGeneMat(rxn_ind,:) == 1),' or ');
end


% The reaction HMR_4147 is associated with ENSG00000107104, which encodes
% a protein, but the function does not appear to be relevant to that rxn:
%  HMR_4147: CoA[m] + GTP[m] + succinate[m] <=> GDP[m] + Pi[m] + succinyl-CoA[m]
% Therefore, the gene should be removed from the associated grRule.
[~,rxn_ind] = ismember('HMR_4147',ihuman.rxns);
if contains(ihuman.grRules(rxn_ind),'ENSG00000107104')
    if contains(ihuman.grRules(rxn_ind),' and ')
        error('grRule(s) contains AND expressions, and should be updated manually.');
    end
    [~,gene_ind] = ismember('ENSG00000107104',ihuman.genes);
    ihuman.rxnGeneMat(rxn_ind,gene_ind) = 0;
    ihuman.grRules(rxn_ind) = join(ihuman.genes(ihuman.rxnGeneMat(rxn_ind,:) == 1),' or ');
end


% ATP synthase (HMR_6916) is incorrectly associated with many genes, such 
% as V-ATPases. The associated grRule therefore needs to be updated such
% that it contains only genes associated with ATP synthase.
[~,rxn_ind] = ismember('HMR_6916',ihuman.rxns);
ATPgenes = {'ENSG00000099624';'ENSG00000110955';'ENSG00000116459';'ENSG00000124172';...
            'ENSG00000125375';'ENSG00000135390';'ENSG00000152234';'ENSG00000154518';...
            'ENSG00000154723';'ENSG00000159199';'ENSG00000165629';'ENSG00000167283';...
            'ENSG00000167863';'ENSG00000169020';'ENSG00000198899';'ENSG00000241468';...
            'ENSG00000241837';'ENSG00000123472';'ENSG00000171953';'ENSG00000249222';...
            'ENSG00000228253';'ENSG00000156411';'ENSG00000173915'};
if ~all(cellfun(@(g) contains(ihuman.grRules(rxn_ind),g), ATPgenes)) || (sum(ihuman.rxnGeneMat(rxn_ind,:)) ~= 23)
    
    % ensure the rule doesn't contain ANDs
    if contains(ihuman.grRules(rxn_ind),' and ')
        error('grRule(s) contains AND expressions, and should be updated manually.');
    end
    
    % get the indices of each gene
    [~,gene_ind] = ismember(ATPgenes,ihuman.genes);
    if any(gene_ind == 0)
        % add any new genes to the model
        ihuman.genes = [ihuman.genes; ATPgenes(gene_ind == 0)];
        [~,gene_ind] = ismember(ATPgenes,ihuman.genes);
        ihuman.rxnGeneMat(1,length(ihuman.genes)) = 0;  % add new columns to rxnGeneMat
    end
    
    % update the rxnGeneMat row corresponding to the rxn
    ihuman.rxnGeneMat(rxn_ind,:) = 0;
    ihuman.rxnGeneMat(rxn_ind,gene_ind) = 1;
    
    % update the grRule (all ORs/isozymes - should be curated to incorporate enzyme complexes later)
    ihuman.grRules(rxn_ind) = join(ihuman.genes(ihuman.rxnGeneMat(rxn_ind,:) == 1),' or ');
    
end



%% Update direction of citrate-malate antiporter
% The citrate-malate antiport reaction is backwards:
%  HMR_4964: citrate[c] + H+[c] + malate[m] => citrate[m] + H+[m] + malate[c]
% It should be transporting citrate from M to C, and malate from C to M.
% Therefore, the reaction will be reversed.
[~,rxn_ind] = ismember('HMR_4964',ihuman.rxns);
ihuman.S(:,rxn_ind) = -ihuman.S(:,rxn_ind);



%% Removal of HMG-CoA transporter
% The model allows transport of HMR-CoA (3?hydroxy?3?methylglutaryl) across
% the mitochondrial membrane, which should not be possible:
%   HMR_1572: HMG-CoA[c] <=> HMG-CoA[m]
% In addition, there is another reaction that transports the same compound
% across the peroxisomal membrane:
%   HMGCOAtx: HMG-CoA[c] <=> HMG-CoA[p]
% Both of these reactions will therefore be deleted from the model.

if any(ismember({'HMR_1572';'HMGCOAtx'},ihuman.rxns))
    
    nrxns_orig = length(ihuman.rxns);  % record original number of rxns
    del_ind = ismember(ihuman.rxns,{'HMR_1572';'HMGCOAtx'});  % get indices of rxns to delete
    ihuman = removeReactions(ihuman,{'HMR_1572';'HMGCOAtx'});  % delete rxns
    
    % need to update some non-standard fields of the structure that are missed
    if length(ihuman.mets) == nrxns_orig
        error('This code will cause errors if the mets and rxns field are equal in size.');
    end
    model_fields = fieldnames(ihuman);
    for i = 1:length(model_fields)
        % if one of the field dimensions matches the length of the original
        % rxn list, then delete the indices of that field that correspond
        % to the delete rxns.
        dim = find(ismember(size(ihuman.(model_fields{i})),nrxns_orig));
        if dim == 1
            ihuman.(model_fields{i})(del_ind,:) = [];
        elseif dim == 2
            ihuman.(model_fields{i})(:,del_ind) = [];
        end
    end
    
end



%% Update protein-related fields
% Since some of the grRules have changed, the associated protein fields
% need to be updated.

[ihuman.prRules,ihuman.proteins,ihuman.rxnProtMat] = translateGeneRules(ihuman.grRules,'UniProt');



%% Save updated model file
save('humanGEM.mat','ihuman');

% remove intermediate variables
clear('ATPgenes','del_ind','dim','f','gene_ind','gene_inds','i','model_fields','nrxns_orig','rxn_ind');





