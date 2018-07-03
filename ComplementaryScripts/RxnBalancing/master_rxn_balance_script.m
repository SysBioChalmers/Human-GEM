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

% add a new field to the model: "metsNoComp"
model.metsNoComp = regexprep(model.mets,'_?\w$','');

% clear unnecessary intermediate variables
clear ind


%% Update metabolites to correct protonation state at physiological pH


%............ Obtain information for mets originating from HMR ............

% load HMR model with Recon3D and MNX met association information
load('ComplementaryScripts/MetAssociation/ihumanMets2MNX_v2.mat');  % loads as variable "ihuman"

% load Recon3D model
load('ComplementaryScripts/Recon3D_301.mat');  % loads as variable "Recon3D"

% extract unique (compartment-free) met information from HMR
m = {};
ihuman.metsNoComp = regexprep(ihuman.mets,'\w$','');
[m.mets,ind] = unique(ihuman.metsNoComp);
m.metNames = ihuman.metNames(ind);
m.metFormulas = alphabetizeMetFormulas(ihuman.metFormulas(ind));
m.metRecon3DID = ihuman.metRecon3DID(ind);
m.metMNXID = ihuman.metMNXID(ind);
if isfield(ihuman,'metCharges')
    m.metCharges = ihuman.metCharges(ind);
else
    m.metCharges = repmat({''},size(m.mets));
end

% retrieve metFormulas from Recon3D model for uniquely mapped mets
m.metFormulasR3D = repmat({''},size(m.mets));
m.metRecon3DID(cellfun(@numel,m.metRecon3DID) > 1) = {''};  % remove non-uniquely mapped entries
m.metRecon3DID = flattenCell(m.metRecon3DID,true);  % convert to column cell array
Recon3D.metsNoComp = regexprep(Recon3D.mets,'\[\w\]$','');  % remove compartment abbrevs from Recon3D met IDs
[hasMatch,ind] = ismember(m.metRecon3DID,Recon3D.metsNoComp);
m.metFormulasR3D(hasMatch) = Recon3D.metFormulas(ind(hasMatch));
m.metFormulasR3D = alphabetizeMetFormulas(m.metFormulasR3D);  % alphabetize met formulas

% also retrieve Recon3D met charges
m.metChargesR3D = repmat({[]},size(m.mets));
m.metChargesR3D(hasMatch) = num2cell(Recon3D.metCharges(ind(hasMatch)));  % make it a cell, so we can have empty entries for unknown charges

% identify met formulas differing only by the number of protons
m.metFormulas_noH = regexprep(m.metFormulas,'H\d*','');
m.metFormulasR3D_noH = regexprep(m.metFormulasR3D,'H\d*','');
matchInd = strcmp(m.metFormulas_noH,m.metFormulasR3D_noH);

% take the met formula and charge from Recon3D for matches
m.metFormulas_new = m.metFormulas;
m.metFormulas_new(matchInd) = m.metFormulasR3D(matchInd);
m.metCharges_new = m.metCharges;
m.metCharges_new(matchInd) = m.metChargesR3D(matchInd);

% now update merged model (Recon4) with new met properties
[hasMatch,ind] = ismember(model.metsNoComp,m.mets);
model.metFormulas(hasMatch) = m.metFormulas_new(ind(hasMatch));
model.metCharges(hasMatch) = m.metCharges_new(ind(hasMatch));
if ~isfield(model,'metMNXID')
    model.metMNXID = repmat({''},size(model.mets));
end
model.metMNXID(hasMatch) = m.metMNXID(ind(hasMatch));

% clear intermediate variables
clear hasMatch ind matchInd Recon3D m

%..........................................................................



%.......... Obtain information for mets originating from Recon3D ..........

% load Recon3D model with metMNXID information
load('ComplementaryScripts/MetAssociation/Recon3Mets2MNX.mat');  % loads as variable "Recon3D"

% the metFormulas have problems in the above model, so restore the original
x = load('ComplementaryScripts/Recon3D_301.mat');
Recon3D.metFormulas = x.Recon3D.metFormulas;
clear x

% generate compartment-free mets
Recon3D.metsNoComp = regexprep(Recon3D.mets,'_\w$','');

% update merged model (Recon4) with met properties retrieved from Recon3D
% Note: even though each compartment-free met can match with multiple
% Recon3D compartment-free mets, it doesn't matter, because the information
% associated with the mets (formula, charge, MNXID) is the same for all
% compartment-versions of each metabolite.
[hasMatch,ind] = ismember(model.metsNoComp,Recon3D.metsNoComp);
model.metFormulas(hasMatch) = Recon3D.metFormulas(ind(hasMatch));
model.metCharges(hasMatch) = num2cell(Recon3D.metCharges(ind(hasMatch)));

% Use the metMNXIDs obtained from BiGG DB mapping, unless the MNXID is
% missing, in which case the MNXID(s) obtained via metName and external IDs
% will be used.
R3MNXIDs = Recon3D.metBiGGDB2MNX;
empty_ind = cellfun(@isempty,R3MNXIDs);
if any(contains(Recon3D.metMNXID,';'))
    % if multiple IDs are separated by semicolon, split into nested cells
    Recon3D.metMNXID = cellfun(@(m) strsplit(m,{';',' '}),Recon3D.metMNXID,'UniformOutput',false);
end
R3MNXIDs(empty_ind) = Recon3D.metMNXID(empty_ind);
model.metMNXID(hasMatch) = R3MNXIDs(ind(hasMatch));  

% alphabetize met formulas
model.metFormulas = alphabetizeMetFormulas(model.metFormulas);

% clear intermediate variables
clear hasMatch ind empty_ind R3MNXIDs

%..........................................................................



%.............. Obtain met information via MNX associations ...............

% regenerate compartment-free met information structure, for calc speed
m = {};
mets_nocomp = regexprep(model.mets,'_?\w$','');
[m.mets,ind] = unique(mets_nocomp);
m.metNames = model.metNames(ind);
m.metFormulas = model.metFormulas(ind);
m.metCharges = model.metCharges(ind);
m.metMNXID = model.metMNXID(ind);
m.metMNXID = nestCell(flattenCell(m.metMNXID,true),true);  % ugly way to ensure correct format

% load met info retrieved from MNX database (~2 min load time)
mnx_met = buildMNXmodel('met');

% subset the MNX data structure to greatly reduce computation time
allMNXID = unique(flattenCell(m.metMNXID,true));  % get all MNXIDs in model
allMNXID(ismember(allMNXID,{''})) = [];  % remove empty string entry
ind = ismember(mnx_met.mets,allMNXID);

mnx_sub = {};
mnx_sub.mets = mnx_met.mets(ind);
mnx_sub.metNames = mnx_met.metNames(ind);
mnx_sub.inchis = mnx_met.inchis(ind);
mnx_sub.metFormulas = alphabetizeMetFormulas(mnx_met.metFormulas(ind));
mnx_sub.metCharges = mnx_met.metCharges(ind);
mnx_sub.metSMILES = mnx_met.metSMILES(ind);


% organize MNX information for each model met
mnxFields = {'metNames';'inchis';'metFormulas';'metCharges';'metSMILES'};
for i = 1:length(mnxFields)
    m.(strcat('MNX',mnxFields{i})) = cell(size(m.mets));
end
for i = 1:length(m.mets)
    if isempty(m.metMNXID{i})
        continue
    end
    [~,ind] = ismember(m.metMNXID{i},mnx_sub.mets);
    for j = 1:length(mnxFields)
        m.(strcat('MNX',mnxFields{j})){i} = mnx_sub.(mnxFields{j})(ind)';
    end
end


% check if each met formula in model agrees with the MNX formula(s)
m.metFormulaCheck = repmat({''},size(m.mets));
for i = 1:length(m.mets)
    if isempty(m.metFormulas{i}) || isempty(m.MNXmetFormulas{i}) || all(cellfun(@isempty,m.MNXmetFormulas{i}))
        % skip cases where model or MNX formulas are missing
        continue
    end
    if any(ismember(m.MNXmetFormulas{i},m.metFormulas{i}))
        m.metFormulaCheck{i} = 'confirmed';
    else
        m.metFormulaCheck{i} = 'no match';
    end
end

%..........................................................................










