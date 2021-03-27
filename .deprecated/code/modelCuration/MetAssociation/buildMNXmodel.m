function model = buildMNXmodel(model_type)
%BUILDMNXMODEL  Construct a model structure from MNX database files.
%
%   MODEL = BUILDMNXMODEL(MODEL_TYPE)
%
%--------------------------------- INPUTS ---------------------------------
%
%   model_type      (OPTIONAL) Specify the type of model that will be
%                   built, which determines what information is included:
%
%                   'model' (DEFAULT) All metabolite, reaction, and 
%                           compartment information is included in model
%                           construction. However, note that some mets and
%                           rxns will be excluded due to issues with lack
%                           of uniqueness and metabolites that don't exist
%                           in any reactions. Therefore, the model will not
%                           contain every rxn and met in the MXN database.
%
%                   'met'   Include only metabolite-related fields. This
%                           will include ALL metabolites in the MNX
%                           database, and will therefore contain more met
%                           information than if the 'model' option is
%                           specified.
%
%                   'rxn'   Included only reaction-related fields. No
%                           stoichiometry matrix will be constructed, but
%                           the model will include ALL reactions from the
%                           MNX database.
%
%                  'both'   Include fields related to either reactions or
%                           metabolites, but do not associate them with
%                           each other (i.e., no stoichimetry matrix). This
%                           will contain all the fields that would be
%                           obtained by running both the 'met' and 'rxn'
%                           options. More information will be present in
%                           this format than if the 'model' option is
%                           specified, because reactions and/or metabolites
%                           will not be removed to construct a functioning
%                           stoichiometry matrix.
%


% Specify directory containing MNX database files.
% * Note that the original files downloaded from the database were changed
%   from .tsv to .txt, and all header lines (those beginning with #) were 
%   removed, except the last one, which was kept, but the leading '#' was 
%   removed.
MNXdir = '';
if isempty(MNXdir)
    error('The MNX directory needs to be specified within the buildMNXmodel function.')
end

% handle inputs
if nargin < 1
    model_type = 'model';
end


%% Initialize model
% use similar field order as those from SBML
model=[];
model.id='MNXdatabase';
model.description='';
model.rxns={};
model.mets={};
model.S=[];
model.lb=[];
model.ub=[];
model.rev=[];
model.c=[];
model.b=[];
model.comps={};
model.compNames={};
model.rxnNames={};
model.grRules={};
model.rxnGeneMat=[];
model.eccodes={};
model.genes={};
model.metNames={};
model.metComps=[];
model.inchis={};
model.metFormulas={};
% model.metMiriams={};
model.metCharges=[];


% remove unnecessary fields for met-only or rxn-only structures
if strcmpi(model_type,'model')
    fprintf('*** Building MODEL structure from MNX Database (reactions AND metabolites, with stoichiometry matrix) ***\n\n')
    model.description = 'MetaNetX Database Model';
elseif strcmpi(model_type,'met')
    fprintf('*** Building information structure for MNX METABOLITES ***\n\n')
    model.description = 'MetaNetX Metabolite Data';
    model = rmfield(model,{'rxns','S','lb','ub','rev','c','b','comps','compNames','rxnNames','grRules','rxnGeneMat','eccodes','genes','metComps'});
elseif strcmpi(model_type,'rxn')
    fprintf('*** Building information structure for MNX REACTIONS ***\n\n')
    model.description = 'MetaNetX Reaction Data';
    model = rmfield(model,{'mets','S','lb','ub','rev','c','b','comps','compNames','rxnNames','grRules','rxnGeneMat','genes','metNames','metComps','inchis','metFormulas','metCharges'});
elseif strcmpi(model_type,'both')
    fprintf('*** Building information structure for MNX REACTIONS AND METABOLITES ***\n\n')
    model.description = 'MetaNetX Reaction and Metabolite Data';
    model = rmfield(model,{'S','lb','ub','rev','c','b','comps','compNames','grRules','rxnGeneMat','genes','metComps'});
else
    error('Invalid MODEL_TYPE argument.');
end


%% Load compartment properties

if strcmpi(model_type,'model')
    fprintf('Loading model compartment data... ');
    mnx = readtable([MNXdir,'comp_prop.txt'],'Delimiter','\t');
    
    % these compartment abbreviations will probably cause problems for some
    % functions, as they are not a single letter, but a full string
    model.comps = mnx.MNX_ID;
    model.compNames = mnx.Description;
    fprintf('Done.\n');
end

%% Load reaction properties

if ismember(model_type,{'model','rxn','both'})
    
    fprintf('Loading reaction data... ');
    mnx = readtable([MNXdir,'reac_prop.txt'],'Delimiter','\t');
    fprintf('Done.\n');
    
    if strcmp(model_type,'model')
        
        % remove reactions with ambiguous stoich coeff; e.g., (n), (N), (n-1), (2n), etc.
        fprintf('Processing reaction data... ');
        del_rxns = contains(mnx.Equation,'(');
        mnx(del_rxns,:) = [];
        
        % also remove some reactions that are duplicates, or imbalanced
        % (these were identified through a separate analysis)
        del_rxns = ismember(mnx.MNX_ID,{'MNXR109831','MNXR126812','MNXR117558','MNXR107762','MNXR135187','MNXR127223'});
        mnx(del_rxns,:) = [];
        
        % split each rxn equations into reactants and products
        eqns = split(mnx.Equation,' = ');
        reactants = eqns(:,1);
        products = eqns(:,2);
        
        % Remove reactions that have the same met as both a reactant and a product,
        % and with the same stoich coeff. This happens in cases where an enzyme has
        % been included in the reaction equation, or there is a mistake in the
        % reaction, or MNX is assuming that the two forms of a compound are
        % identical, when maybe it should not. Many are the result of the database
        % selecting the form of the compound that is most abundant at pH 7.3;
        % for example, the rxn: NH3 + H -> NH4 is  NH4 -> NH4 in the database.
        smet_reac = cellfun(@(r) regexp(r,'[0-9.]+ \w+@\w+','match'),reactants,'UniformOutput',false);
        smet_prod = cellfun(@(r) regexp(r,'[0-9.]+ \w+@\w+','match'),products,'UniformOutput',false);
        del_rxns = cellfun(@(reac,prod) ~isempty(intersect(reac,prod)),smet_reac,smet_prod);
        mnx(del_rxns,:) = [];
        reactants(del_rxns) = [];
        products(del_rxns) = [];
        
    else
        
        % add a field that contains a list of mets that participate in each reaction
        rxnMets = cellfun(@(r) regexp(r,'(\w+)@\w+','tokens'),mnx.Equation,'UniformOutput',false);
        model.rxnMets = cellfun(@(r) unique([r{:}]),rxnMets,'UniformOutput',false);
        
    end
    
    % extract reaction MNX IDs and other information
    model.rxns = mnx.MNX_ID;
    model.rxnBalanced = mnx.Balance;
    model.eccodes = mnx.EC;
    model.rxnNames = repmat({''},size(model.rxns));  % no rxnNames present
    model.rxnEqnMNX = mnx.Equation;
    model.rxnEqnNames = mnx.Description;
    
    if strcmp(model_type,'model')
        
        % associate every reaction with each of its constituent metabolites
        reactant_mets = cellfun(@(r) regexp(r,'\w+@\w+','match'),reactants,'UniformOutput',false);
        reactant_rxn_names = cellfun(@(rname,r) repmat({rname},1,numel(r)),mnx.MNX_ID,reactant_mets,'UniformOutput',false);
        reactant_rxn_met_pairs = [[reactant_mets{:}]',[reactant_rxn_names{:}]'];
        
        product_mets = cellfun(@(r) regexp(r,'\w+@\w+','match'),products,'UniformOutput',false);
        product_rxn_names = cellfun(@(rname,r) repmat({rname},1,numel(r)),mnx.MNX_ID,product_mets,'UniformOutput',false);
        product_rxn_met_pairs = [[product_mets{:}]',[product_rxn_names{:}]'];
        
        rxn_met_pairs = [reactant_rxn_met_pairs; product_rxn_met_pairs];
        
        
        % get a unique list of all metabolite-compartment combinations
        model.mets = unique([[reactant_mets{:}]';[product_mets{:}]']);
        model.metMNXID = regexprep(model.mets,'@.+$','');  % get list of metabolites without compartment names
        uniq_mets_nocomp = unique(model.metMNXID);  % unique list without compartments
        
        % extract metabolite compartments
        metComps = regexprep(model.mets,'^.+@','');
        [~,model.metComps] = ismember(metComps,model.comps);
        
        % obtain stoich coeffs
        reactant_stoichs = cellfun(@(r) -str2double(regexp(r,'^[0-9.]+ | [0-9.]+ ','match')), reactants,'UniformOutput',false);
        product_stoichs = cellfun(@(r) str2double(regexp(r,'^[0-9.]+ | [0-9.]+ ','match')), products,'UniformOutput',false);
        stoich_coeffs = [[reactant_stoichs{:}]';[product_stoichs{:}]'];
        
        % construct stoichiometry matrix
        [~,met_ind] = ismember(rxn_met_pairs(:,1),model.mets);
        [~,rxn_ind] = ismember(rxn_met_pairs(:,2),model.rxns);
        model.S = sparse(met_ind,rxn_ind,stoich_coeffs);
        
    end
    
    % get other rxn IDs
    rxnSources = mnx.Source;
    rxnSources(~contains(rxnSources,':')) = {''};  % only keep entries with a colon, which indicates they have a source beyond MNX
    rxnSourceNames = regexprep(rxnSources,':.+$','');  % retrieve source name
    rxnSourceIDs = regexprep(rxnSources,'^.+:','');  % retrieve source ID
    
    % add model fields corresponding to available source names
    model.rxnMNXID = model.rxns;
    model.rxnBiGGID = repmat({''},size(model.rxns));
    model.rxnBiGGID(ismember(rxnSourceNames,'bigg')) = rxnSourceIDs(ismember(rxnSourceNames,'bigg'));
    model.rxnKEGGID = repmat({''},size(model.rxns));
    model.rxnKEGGID(ismember(rxnSourceNames,'kegg')) = rxnSourceIDs(ismember(rxnSourceNames,'kegg'));
    model.rxnMetaCycID = repmat({''},size(model.rxns));
    model.rxnMetaCycID(ismember(rxnSourceNames,'metacyc')) = rxnSourceIDs(ismember(rxnSourceNames,'metacyc'));
    model.rxnREACTOMEID = repmat({''},size(model.rxns));
    model.rxnREACTOMEID(ismember(rxnSourceNames,'reactome')) = rxnSourceIDs(ismember(rxnSourceNames,'reactome'));
    model.rxnRheaID = repmat({''},size(model.rxns));
    model.rxnRheaID(ismember(rxnSourceNames,'rhea')) = rxnSourceIDs(ismember(rxnSourceNames,'rhea'));
    model.rxnSABIORKID = repmat({''},size(model.rxns));
    model.rxnSABIORKID(ismember(rxnSourceNames,'sabiork')) = rxnSourceIDs(ismember(rxnSourceNames,'sabiork'));
    model.rxnSEEDID = repmat({''},size(model.rxns));
    model.rxnSEEDID(ismember(rxnSourceNames,'seed')) = rxnSourceIDs(ismember(rxnSourceNames,'seed'));
    
    fprintf('Done.\n');
    
end


%% Load metabolite properties

if ismember(model_type,{'model','met','both'})
    
    fprintf('Loading metabolite data... ');
    mnx = readtable([MNXdir,'chem_prop.txt'],'Delimiter','\t');  % ~45 sec load time
    mnx_xref = readtable([MNXdir,'chem_xref.txt'],'Delimiter','\t');
    fprintf('Done.\n');
    
    fprintf('Processing metabolite data... ');
    if strcmp(model_type,'model') 
        % retrieve row indices corresponding to model met list
        [~,met_ind] = ismember(model.metMNXID,mnx.MNX_ID);
        mnx = mnx(met_ind,:);
    else
        model.mets = mnx.MNX_ID;
        model.metMNXID = mnx.MNX_ID;
    end
    
    % extract information and add to model
    model.metNames = mnx.Description;
    model.metFormulas = mnx.Formula;
    model.metCharges = mnx.Charge;
    model.inchis = mnx.InChI;
    model.metSMILES = mnx.SMILES;
    
    % get other IDs
    metSources = mnx.Source;
    metSources(~contains(metSources,':')) = {''};  % only keep entries with a colon, which indicates they have a source beyond MNX
    metSourceNames = regexprep(metSources,':.+$','');  % retrieve source name
    metSourceIDs = regexprep(metSources,'^.+:','');  % retrieve source ID
    
    
    % remove rows of mnx_xref that don't contain external IDs
    mnx_xref(~contains(mnx_xref.XREF,':') | contains(mnx_xref.XREF,'deprecated:'),:) = [];
    
    % extract source and ID information from MNX xref data
    metSourceNamesX = regexprep(mnx_xref.XREF,':.+$','');  % retrieve source name
    metSourceIDsX = regexprep(mnx_xref.XREF,'^.+:','');  % retrieve source ID
    mnxID2extID = [mnx_xref.MNX_ID,metSourceNamesX,metSourceIDsX];
    
    % Split met name field by delimiters '|' and '; '. This takes a very
    % long time, and should only be done if necessary.
    splitnames = true;  % CHANGE TO FALSE FOR FASTER PROCESSING
    if ( splitnames )
        % change semicolon+space delimiters to vertical bar "|"
        regexprep(mnx_xref.Description,'; ','|');
        
        % break processing into 10 chunks
        chunk_ind = 1:round(length(mnx_xref.Description)/10):length(mnx_xref.Description);
        chunk_ind(end) = length(mnx_xref.Description)+1;
        
        % split met names
        mnxDescr = {};
        fprintf('\nSplitting MNX metabolite names: ');
        for i = 1:length(chunk_ind)-1
            fprintf('%u%% ',round(chunk_ind(i)/length(mnx_xref.Description)*100));
            mnxDescr = [mnxDescr; cellfun(@(x) strsplit(x,'|'),mnx_xref.Description(chunk_ind(i):chunk_ind(i+1)-1),'UniformOutput',false)];
        end
        fprintf('Done.\n');
        
        % flatten nested cells to obtain an Nx2 matrix of MNXID-name pairs
        fprintf('Organizing MNXID-name pairs... ');
        mnxIDs = cellfun(@(id,descr) repmat({id},1,numel(descr)),mnx_xref.MNX_ID,mnxDescr,'UniformOutput',false);
        model.mnxID2name = [[mnxIDs{:}]',[mnxDescr{:}]'];
        
        % remove repeated ID-name pairs
        [~,numericPair] = ismember(lower(model.mnxID2name),unique(lower(model.mnxID2name)));  % make numeric for faster processing
        [~,uniq_ind] = unique(numericPair,'rows');
        model.mnxID2name = model.mnxID2name(uniq_ind,:);
        fprintf('Done.\n');
        
    else
        % add unsplit name field as additional column
        mnxID2extID = [mnxID2extID,mnx_xref.Description];
    end
    
    
    % add model fields corresponding to available source names
    sNames = {'bigg','chebi','envipath','hmdb','kegg','lipidmaps','metacyc','reactome','sabiork','seed','slm'};
    fNames = {'metBiGGID','metChEBIID','metEnviPathID','metHMDBID','metKEGGID','metLIPIDMAPSID','metMetaCycID','metREACTOMEID','metSABIORKID','metSEEDID','metSLMID'};

    % This portion of code would be nice, but it takes forever to run. 
%     for i = 1:length(sNames)
%         fprintf('Retrieving %s met IDs... ',upper(sNames{i}));
%         id2id = mnxID2extID(ismember(metSourceNamesX,sNames{i}),:);  % subset ID data for speed
%         ind = ismember(model.mets,id2id(:,1));  % subset met list for speed
%         model.(fNames{i}) = repmat({''},size(model.mets));
%         model.(fNames{i})(ind) = cellfun(@(id) id2id(ismember(id2id(:,1),id),3),model.mets(ind),'UniformOutput',false);
%         fprintf('Done.\n');
%     end

    % instead, just save the unprocessed information
    [~,ind] = ismember(mnxID2extID(:,2),sNames);
    mnxID2extID(:,2) = fNames(ind);  % rename to match field names
    model.mnxID2extID = mnxID2extID;
    
    
    % fill in other fields
    model.metBiGGID = repmat({''},size(model.mets));
    model.metBiGGID(ismember(metSourceNames,'bigg')) = metSourceIDs(ismember(metSourceNames,'bigg'));
    model.metChEBIID = repmat({''},size(model.mets));
    model.metChEBIID(ismember(metSourceNames,'chebi')) = metSourceIDs(ismember(metSourceNames,'chebi'));
    model.metEnviPathID = repmat({''},size(model.mets));
    model.metEnviPathID(ismember(metSourceNames,'envipath')) = metSourceIDs(ismember(metSourceNames,'envipath'));
    model.metHMDBID = repmat({''},size(model.mets));
    model.metHMDBID(ismember(metSourceNames,'hmdb')) = metSourceIDs(ismember(metSourceNames,'hmdb'));
    model.metKEGGID = repmat({''},size(model.mets));
    model.metKEGGID(ismember(metSourceNames,'kegg')) = metSourceIDs(ismember(metSourceNames,'kegg'));
    model.metLIPIDMAPSID = repmat({''},size(model.mets));
    model.metLIPIDMAPSID(ismember(metSourceNames,'lipidmaps')) = metSourceIDs(ismember(metSourceNames,'lipidmaps'));
    model.metMetaCycID = repmat({''},size(model.mets));
    model.metMetaCycID(ismember(metSourceNames,'metacyc')) = metSourceIDs(ismember(metSourceNames,'metacyc'));
    model.metREACTOMEID = repmat({''},size(model.mets));
    model.metREACTOMEID(ismember(metSourceNames,'reactome')) = metSourceIDs(ismember(metSourceNames,'reactome'));
    model.metSABIORKID = repmat({''},size(model.mets));
    model.metSABIORKID(ismember(metSourceNames,'sabiork')) = metSourceIDs(ismember(metSourceNames,'sabiork'));
    model.metSEEDID = repmat({''},size(model.mets));
    model.metSEEDID(ismember(metSourceNames,'seed')) = metSourceIDs(ismember(metSourceNames,'seed'));
    model.metSLMID = repmat({''},size(model.mets));  % SwissLipids database
    model.metSLMID(ismember(metSourceNames,'slm')) = metSourceIDs(ismember(metSourceNames,'slm'));
    
    
    fprintf('Done.\n');
    
    % %% Merge reactions with identical stoichiometry
    %
    % % group reactions with identical stoich columns
    % [~,~,rxn_groups] = unique(model.S','rows');
    %
    % % determine which groups have more than one member
    % rep_inds = ismember(rxn_groups,find(histcounts(rxn_groups,'BinMethod','integers') > 1));
    % rep_groups = unique(rxn_groups(rep_inds));
    
end


%% Final model adjustments

if strcmp(model_type,'model')
    model.c = zeros(size(model.rxns));
    model.b = zeros(size(model.mets));
    model.lb = -1000*ones(size(model.rxns));
    model.ub = 1000*ones(size(model.rxns));
    model.rev = true(size(model.rxns));
    
    model.genes = {'dummy'};
    model.grRules = repmat({''},size(model.rxns));
    model.rxnGeneMat = zeros(size(model.rxns));
end

if ismember(model_type,{'model','met','both'})
    % remove NAs from metFormulas
    model.metFormulas(ismember(model.metFormulas,'NA')) = {''};
    if strcmp(model_type,'model')
        % change format of mets from "met@comp" to "met[comp]"
        model.mets = strcat(regexprep(model.mets,'@','['),']');
    end
end



