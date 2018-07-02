% Master script for balancing Recon4 reactions.
%   1) Minor curation/bug fixes in model metabolite information
%   2) Update metabolite form to reflect the major species at pH 7.3
%   3) Identify mass-imbalanced reactions, and correct


%% Load latest version of model
x = load('ComplementaryScripts/modelIntegration/HMR3_02.mat');
f = fields(x);
model = x.(f{1});  % ugly method to load model without knowing variable name
clear x f


%% Minor curation/bug fixes in model metabolite information

% reorder metabolite formulas so elements are in alphabetical order
% NOTE: this does not re-order formulas containing special characters 
%       (i.e., those other than A-Z,0-9).
model.metFormulas = alphabetizeMetFormulas(model.metFormulas);

% the metabolite formula of lepidimoide (m02357c and m02357s) contains a
% period "." at the end, which needs to be removed
ind = ismember(model.mets,{'m02357c','m02357s'});
model.metFormulas(ind) = regexprep(model.metFormulas(ind),'\.','');

% update protonation of a few metabolites
ind = contains(model.mets,'m02751');
model.metFormulas(ind) = {'HO4P'};  % originally listed as "H3PO4"
ind = contains(model.mets,'m02949');
model.metFormulas(ind) = {'O3S'};  % originally listed as "H2SO3"

% rename 'temp006x' to 'm01451x' (both correspond to 'cholesterol-ester pool')
model.mets(ismember(model.mets,'temp006x')) = {'m01451x'};

% rename some BiGG IDs to be consistent with the current BiGG database
model.metBiGGID = regexprep(model.metBiGGID,',|-|/','_');  % replace commas, dashes, and slashes with underscores
model.metBiGGID(ismember(model.metBiGGID,'Xyl_L_Ser_(protein)')) = {'xser'};



%% Update metabolites to correct protonation state at physiological pH

% First focus on mets originating from HMR

% load HMR model with Recon3D and MNX met association information
load('ComplementaryScripts/MetAssociation/ihumanMets2MNX_v2.mat');  % loads as variable "ihuman"

% load Recon3D model
load('ComplementaryScripts/Recon3D_301.mat');  % loads as variable "Recon3D"

% extract unique (compartment-free) met information from HMR
m = {};
mets_nocomp = regexprep(ihuman.mets,'\w$','');
[m.mets,ind] = unique(mets_nocomp);
m.metNames = ihuman.metNames(ind);
m.metFormulas = ihuman.metFormulas(ind);
m.metRecon3DID = ihuman.metRecon3DID(ind);

% retrieve metFormulas from Recon3D model for uniquely mapped mets
m.metFormulasRecon3D = repmat({''},size(m.mets));
m.metRecon3DID(cellfun(@numel,m.metRecon3DID) > 1) = {''};  % remove non-unique mapping entries
m.metRecon3DID = flattenCell(m.metRecon3DID,true);  % convert to column cell array
Recon3D.mets = regexprep,.......
[hasmatch,ind] = ismember(m.metRecon3DID,Recon3D.mets);
m.metFormulasRecon3D(hasmatch) = Recon3D.metFormulas(ind(hasmatch));






% Begin by comparing all model met formulas to the MNX met formulas.

% if metMNXID field is missing, try to retrieve the information
if ~isfield(model,'metMNXID')
    
    % initialize
    model.metMNXID = repmat({''},size(model.mets));
    
    % retrieve met MNXIDs from HMR
    load('ComplementaryScripts/MetAssociation/ihumanMets2MNX_v2.mat');
    [hasmatch,ind] = ismember(model.mets,ihuman.mets);
    model.metMNXID(hasmatch) = ihuman.metMNXID(ind(hasmatch));
    
    % retrieve met MNXIDs from Recon3D
    load('ComplementaryScripts/MetAssociation/Recon3Mets2MNX.mat');
    % the metMNXIDs in this model are in a different format
    Recon3D.metMNXID = cellfun(@(m) strsplit(m,'; '),Recon3D.metMNXID,'UniformOutput',false);
    [hasmatch,ind] = ismember(model.mets,Recon3D.mets);
    model.metMNXID(hasmatch) = Recon3D.metMNXID(ind(hasmatch));
    
    % reformat metMNXID field to ensure uniform format
    model.metMNXID = nestCell(flattenCell(model.metMNXID,true),true);
    
end

% obtain compartment-free met-MNXID pairs, to avoid working with duplicates
m = {};
mets_nocomp = regexprep(model.mets,'_?\w$','');
[m.mets,ind] = unique(mets_nocomp);
m.metFormulas = model.metFormulas(ind);
















