function [genes,grRules,rxnGeneMat] = translateGeneRules(model,noMatch)
%translateGeneRules  Translate model genes, grRules, rxnGeneMat to alt IDs
%
% translateGeneRules converts the genes, grRules, and rxnGeneMat from a 
% model to six different gene ID types, where the original model gene list
% must correspond to one of the six ID types. The recognized gene ID types
% are:
%       'ENSG'     Ensembl gene ID
%       'ENST'     Ensembl transcript ID
%       'ENSP'     Ensembl protein ID
%       'UniProt'  UniProt protein ID
%       'Name'     HUGO gene abbreviation
%       'Entrez'   Entrez (NCBI) gene ID
%
% NOTE: All trailing ".#" after Entrez or Ensembl ID (e.g., ENSG0000123.1)
%       will be REMOVED. The grRules may therefore change; e.g.:
%           '123.1 or 123.2 or 4365.1' --> '123 or 4365'
%
% NOTE: The function also "cleans" the grRules, meaning duplicate genes,
%       extra parentheses, etc. will be removed from each grRule. See
%       "cleanModelGeneRules.m" for more details.
%
%
% USAGE:
%
%   [genes,grRules,rxnGeneMat] = translateGeneRules(model,noMatch);
%
%
% INPUTS:
%
%   model       A genome-scale metabolic model structure.
%
%   noMatch     What to do in case a gene in the original grRules does not
%               have a corresponding match in the new gene IDs to which the
%               grRules are to be translated. Choose from the following:
%               
%               'original'  (Default) keep the original gene ID. Note, this
%                           will result in mixed gene ID types in grRules
%                           and the corresponding list of genes, which may
%                           cause problems in other applications.
%
%               'delete'    Delete (exclude) the gene from the new rule.
%
%
% OUTPUTS:
%
%   genes       A structure containing gene lists corresponding to each of
%               the different types of gene IDs.
%
%   grRules     A structure containing grRules corresponding to each of the
%               different types of gene IDs.
%
%   rxnGeneMat  A structure containing a rxnGeneMat correspodning to each 
%               of the different types of gene IDs
%
%
% Jonathan Robinson, 2018-07-20


% recognized ID types for genes, transcripts, and proteins
gene_types = {'ENSG';'ENST';'ENSP';'UniProt';'Name';'Entrez'};

% handle input arguments
if nargin < 2
    noMatch = 'orig';
end

% notify user how unmatched genes will be handled
if contains('original',lower(noMatch))
    fprintf('\nNOTE: in the event a new gene ID cannot be matched to an original gene ID,\n');
    fprintf('the original gene ID will be retained (in both the generated new gene list and new grRules).\n\n');
elseif contains('delete',lower(noMatch))
    fprintf('\nNOTE: in the event a new gene ID cannot be matched to an original gene ID,\n');
    fprintf('the gene will be excluded from the generated new gene list and new grRules.\n\n');
else
    error('Invalid input for noMatch.');
end

% initialize empty structure with gene information
tmp = [gene_types,cell(size(gene_types))]';
genes = struct(tmp{:});

% preprocess gene list and gene-reaction rules from model, if necessary
if any(ismember({'GPI','GAPDH'},model.genes))  % need a more robust method to check if it's gene symbols
    % check if the list is gene names (symbols), if so, do not modify
    genes_orig = model.genes;
    rules_orig = model.grRules;
    gene_type_orig = 'Name';
else
    % remove ".#" from gene IDs (e.g. ENSG00000198888.2 -> ENSG00000198888, or 2597.1 -> 2597)
    genes_orig = regexprep(model.genes,'\.\d+$','');
    rules_orig = regexprep(model.grRules,'\.\d+','');
    gene_type_orig = {};
end

% determine the original gene ID type, if not gene names
if isempty(gene_type_orig)
    if ismember(genes_orig{1}(1:4),{'ENSG','ENST','ENSP'})
        % ensembl ID type
        gene_type_orig = genes_orig{1}(1:4);
    elseif all(ismember(cellfun(@length,genes_orig),[6,10])) && all(cellfun(@(x) ~isempty(regexp(x,'^[A-Z]\d\w\w\w\d','once')),genes_orig))
        % Uniprot IDs are 6 or 10 characters, with format [A-Z][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9]...
        gene_type_orig = 'UniProt';
    elseif ~any(isnan(str2double(genes_orig)))
        % Entrez IDs are only digits
        gene_type_orig = 'Entrez';
    else
        error('Unable to recognize the type of gene IDs in the model.');
    end
end

% remove duplicate genes (if any)
genes_orig = unique(genes_orig,'stable');

% import gene-transcript-protein conversion key
tmp = readtable('ComplementaryScripts/GPRs/IDconversion/ensembl_ID_mapping.txt');
conv_key_head = tmp.Properties.VariableNames';  % read header

% change header names to match contents of GENE_TYPES
[~,ind] = ismember(conv_key_head,{'Gene_stable_ID','Transcript_stable_ID',...
    'Protein_stable_ID','UniProtKB_Swiss_Prot_ID','Gene_name','NCBI_gene_ID'});
conv_key_head = gene_types(ind);

% convert numeric Entrez ID data to strings
tmp.NCBI_gene_ID = arrayfun(@num2str,tmp.NCBI_gene_ID,'UniformOutput',false);
conv_key = table2array(tmp); clear tmp;

% replace 'NaN' in Entrez ID column of conversion key with empty cells ('')
conv_key(:,ismember(conv_key_head,'Entrez')) = regexprep(conv_key(:,ismember(conv_key_head,'Entrez')),'NaN','');

% convert rules to all other gene ID types
for i = 1:length(gene_types)
    
    fprintf('Translating rules to %s...\n',gene_types{i});
    
    % extract portion of conversion key required; resulting variable will
    % be two columns of IDs, where column 1 is the original gene ID type,
    % and column 2 is the new gene ID type.
    [~,ind] = ismember([{gene_type_orig}, gene_types(i)], conv_key_head);
    conv_key_sub = conv_key(:,ind);
    
    % remove rows with empty key entries or unneeded genes
    conv_key_sub(any(cellfun(@isempty,conv_key_sub),2),:) = [];
    conv_key_sub(~ismember(conv_key_sub(:,1),genes_orig),:) = [];
    
    % only keep unique rows of conversion key
    [~,tmp] = ismember(conv_key_sub,unique(conv_key_sub));
    [~,ind] = unique(tmp,'rows');
    conv_key_sub = conv_key_sub(ind,:);
    
    if strcmpi(noMatch,'orig')
        % Identify genes that do not exist in the conversion key. For these
        % genes, add additional rows to the conversion key such that if an
        % original gene ID has no corresponding new gene ID, the original
        % gene ID will be used instead.
        unmatched = genes_orig(~ismember(genes_orig,conv_key_sub(:,1)));
        conv_key_sub = [conv_key_sub;repmat(unmatched,1,2)];
    end

    % define inline function to convert between gene IDs
    convertGene = @(g) strjoin(conv_key_sub(ismember(conv_key_sub(:,1),g),2), ' | ');
    
    % convert "and" to "&" and "or" to "|" (easier to work with)
    rules_new = rules_orig;
    rules_new = regexprep(rules_new, ' or ', ' | ');
    rules_new = regexprep(rules_new, ' and ', ' & ');
    
    % check if the target format is the original format
    if strcmp(gene_types{i},gene_type_orig)
        % if current model rules are already in the target format, do not
        % convert gene IDs
        fprintf('\tRules already contain gene IDs of this type, skipping conversion.\n\n');
    else
        % convert gene IDs
        fprintf('\tConverting gene IDs... ');
        
        % This next line identifies gene IDs as collections of characters that
        % are not spaces, parentheses, or the symbols & or |. It then replaces
        % those gene IDs with the new gene ID type, which calls the convertGene
        % inline function to retrieve.
        rules_new = regexprep(rules_new, '[^&|\(\) ]+', '(${convertGene($0)})');
        fprintf('Done.\n');
    end
    
    % clean up rules (removes extra parentheses, repeated genes, etc.)
    fprintf('\tCleaning grRules... ');
    rules_new = cleanModelRules(rules_new);    
    fprintf('Done.\n');
    
    % generate new rxnGeneMat and gene list based on new grRules
    
    % identify genes associated with each reaction
    rxnGenes = cellfun(@(r) unique(regexp(r,'[^&|\(\) ]+','match')),rules_new,'UniformOutput',false);
    
    % construct new gene list
    nonEmpty = ~cellfun(@isempty,rxnGenes);
    genes.(gene_types{i}) = unique([rxnGenes{nonEmpty}]');
    model.(strcat('genes_',gene_types{i})) = genes.(gene_types{i});
    
    % construct new rxnGeneMat
    rxnGeneCell = cellfun(@(rg) ismember(genes.(gene_types{i}),rg),rxnGenes,'UniformOutput',false);
    rxnGeneMat.(gene_types{i}) = sparse(double(horzcat(rxnGeneCell{:})'));
    
    
    % restore "&" as "and" and "|" as "or"
    rules_new = regexprep(rules_new, ' \| ', ' or ');
    rules_new = regexprep(rules_new, ' & ', ' and ');
    
    % add new rules to model structure
    grRules.(gene_types{i}) = rules_new;
    
end




