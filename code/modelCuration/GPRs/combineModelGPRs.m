function humanGEM_new = combineModelGPRs(humanGEM,iHsa,Recon3D)
%combineModelGPRs  Combine the GPRs of HMR, iHsa, and Recon3D.
% 
%   This is the master function for combining and integrating genes and
%   grRules from HMR2, iHsa, and Recon3D into humanGEM. The function makes
%   use of the "integrateGeneRules.m" and "translateGrRules.m" functions,
%   among others.
%
%   The function also adds additional protein-specific fields to the model:
%
%   .proteins   Analogous to the .genes field, but corresponds to proteins,
%               and contains UniProt IDs.
%               NOTE: This field is NOT ALIGNED with the .genes field. The
%                     size/ordering of the .protiens field differs from
%                     .genes, and its entries/indices should not be assumed
%                     to correspond to those in .genes.
%                     
%
%   .prRules    Analogous to the .grRules field, but corresponds to
%               proteins, and contains UniProt IDs.
%
%   .rxnProtMat    Analogous to the .rxnGeneMat field, but corresponds to
%                  proteins.
%                  NOTE: Like the .proteins field, the .rxnProtMat field
%                        will differ in size/ordering (of the columns) with
%                        .rxnGeneMat, because its columns correspond to
%                        entries in .proteins, NOT entries in .genes.
%
%   NOTE: if these protein-related fields already exist in humanGEM, they
%         will be overwritten.
%
%
% USAGE:
%
%   humanGEM_new = combineModelGPRs(humanGEM,iHsa,Recon3D);
%
% INPUTS:
%   
%   humanGEM    The humanGEM model. If not supplied, the function will try
%               to load the humanGEM.mat file.
%
%   iHsa        The iHsa model from Blais, Papin, et. al. If not supplied,
%               the function will try to load the iHsa.mat file.
%
%   Recon3D     The Recon3D model. If not supplied, the function will try
%               to load the Recon3D.mat file.
%
% OUTPUTS:
%
%   humanGEM_new    The humanGEM model, with updated grRules, genes, and
%                   rxnGeneMat fields. Also, with new or updated proteins,
%                   prRules, and rxnProtMat fields.
%


%% load/prepare each model

fprintf('Loading and preparing GEMs... ');

% Human GEM
if nargin < 1 || isempty(humanGEM)
    tmp = load('model/Human-GEM.mat');  % loads as variable "ihuman"
    ihuman = tmp.ihuman;
else
    ihuman = humanGEM;  % just rename it to make things easier
end

% % HMR2
% tmp = load('ComplementaryData/HMR2/HMRdatabase2_02.mat');  % loads as variable "ihuman"
% HMR = tmp.ihuman;

% iHsa
if nargin < 2 || isempty(iHsa)
    tmp = load('ComplementaryData/iHsa/iHsa.mat');  % loads as variable "iHsa"
    iHsa = tmp.iHsa;
end

% An older version of the Recon3D model is loaded, because its grRules are
% slightly better than the current version (differs only in parentheses
% and spacing) because the use of parentheses and spacing is more accurate.
if nargin < 3 || isempty(Recon3D)
    tmp = load('ComplementaryScripts/modelIntegration/Recon3DRaven.mat');  % loads as variable "Recon3D"
    Recon3D = tmp.Recon3DRaven;
    Recon3D.grRules = Recon3D.originalGrRules;
    Recon3D.genes = Recon3D.originalGenes;
    Recon3D.rxnGeneMat = Recon3D.originalRxnGeneMat;
end

fprintf('Done.\n');


%% assemble/pre-process HumanGEM-Recon3D associations

fprintf('Assembling rxn association information... ');

% retrieve the associations from HumanGEM
rxnRecon3DID = ihuman.rxnRecon3DID;

% if necessary, fill in some of the empty associations, which are
% reactions directly imported from Recon3D
empty_ind = find(cellfun(@isempty,rxnRecon3DID));
[has_match,ind] = ismember(ihuman.rxns(empty_ind),Recon3D.rxns);
rxnRecon3DID(empty_ind(has_match)) = Recon3D.rxns(ind(has_match));

% asssemble 1-to-1 association array between HumanGEM and Recon3D rxns
rxnHuman2Recon3D = {};
for i = 1:length(ihuman.rxns)
    ids = strsplit(rxnRecon3DID{i},';')';
    rxnHuman2Recon3D = [rxnHuman2Recon3D; [repmat(ihuman.rxns(i),numel(ids),1),ids]];
end

fprintf('Done.\n');


%% Integrate model grRules, and update "genes" and "rxnGeneMat" fields

% merge grRules
grRules = integrateGeneRules(ihuman,iHsa,Recon3D,rxnHuman2Recon3D);


%*********** Manual changes to some of the grRules ***********

% Some of the grRules were identified as having "ambiguous" AND/OR
% combinations due to a lack of parentheses. These rules have been manually
% inspected and modified below, so that no grRules remain in this ambiguous
% format.


% Pyruvate dehydrogenase complex. 
% The Recon3D rule is ambiguous, and iHsa is lacking one of the subunits.
% It was unclear from literature, but it appears as if the enzyme complex
% requires presence of all six subunits to function properly.
% Original Recon3D rule: (1738.1 and 8050.1) and (5161.1 and 5162.1) and (1737.1) or (1738.1 and 8050.1) and (5160.1 and 5162.1) and (1737.1)
ind = ismember(ihuman.rxns,{'HMR_4137'});
grRules(ind) = {'ENSG00000091140 and ENSG00000110435 and ENSG00000131828 and ENSG00000150768 and ENSG00000163114 and ENSG00000168291'};


% Fatty Acyl Coenzyme A Synthase, and B-Ketoacyl Synthetase. 
% These reactions were imported directly from Recon3D, along with the
% associated GPRs. The original grRule is ambiguous, and is functionally
% equivalent to a rule with only gene 2194 (FASN). The provided references
% (PMIDs) do not help clarify the rule, and it is suspected that the rule
% is missing some pieces (e.g., other ELOVL subunits). Due to the lack of
% information/certainty, this rule will be deleted entirely from all
% reactions with which it is associated.
% Original Recon3D rule: (2194.1) or (2194.1 and 79071.1) and (60481.1) and (54898.1)
ind = ismember(ihuman.rxns,{'FAS100COA','FAS120COA','FAS140COA','FAS160COA','FAS180COA','KAS8'});
grRules(ind) = {''};


% Beta-Galactosidase, Lysosomal. 
% These reactions (and associated GPR) were imported directly from Recon3D.
% The grRule is ambiguous, and there is no references or evidence
% supporting the rule, so it is unclear how it should be arranged.
% Therefore, the rule will be modified as all OR's, so that the
% gene-reaction associations are maintained, but the logic (ANDs) is
% removed.
% Original Recon3D rule: (5476.1) and (2720.1) and (2588.1) and (4758.1) and (5660.1) or (5476.1) and (2720.1) and (2588.1) and (4758.1) and (2760.1) or (5476.1) and (2720.1) and (5660.1) or (5476.1) and (2720.1) and (2760.1)
ind = ismember(ihuman.rxns,{'BGAL1l','BGAL2l'});
grRules(ind) = {'ENSG00000184494 or ENSG00000204386 or ENSG00000223957 or ENSG00000227129 or ENSG00000227315 or ENSG00000228691 or ENSG00000234343 or ENSG00000234846 or ENSG00000064601 or ENSG00000141012 or ENSG00000170266 or ENSG00000196743 or ENSG00000197746'};

%*************************************************************


% generate updated gene list and rxnGeneMat field based on merged grRules
[genes,rxnGeneMat] = getGenesFromGrRules(grRules);

% keep the existing grRule for traceability. If the field already exists,
% then it should not be overwritten.
if ~isfield(ihuman,'priorCombiningGrRules')
    ihuman.priorCombiningGrRules = ihuman.grRules;
end

% replace existing model fields with new fields
ihuman.genes = genes;
ihuman.grRules = grRules;
ihuman.rxnGeneMat = rxnGeneMat;

% unfortunately, the "geneFrom" field is no longer accurate, and therefore
% should be removed
if isfield(ihuman,'geneFrom')
    ihuman = rmfield(ihuman,'geneFrom');
end

% instead of updating the "geneComps" field, it will be removed from the
% model because it is uninformative/irrelevant.
if isfield(ihuman,'geneComps')
    ihuman = rmfield(ihuman,'geneComps');
end


%% Generate and add corresponding protein fields to model

% generate protein list, rules, and rxn association matrix
[prRules,proteins,rxnProtMat] = translateGrRules(ihuman.grRules,'UniProt');

% add new protein-related fields to model structure
ihuman.proteins = proteins;
ihuman.prRules = prRules;
ihuman.rxnProtMat = rxnProtMat;


%% Return model

humanGEM_new = ihuman;





