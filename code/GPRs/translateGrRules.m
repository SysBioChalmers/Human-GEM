function [grRules_new,genes,rxnGeneMat] = translateGrRules(grRules,targetFormat,origFormat,noMatch)
%translateGrRules  Translate grRules to other other gene ID types.
%
% translateGrRules converts grRules to a chosen target gene ID type(s).
% The original and target gene type(s) must be one of the following:
%
%       'ENSG'     Ensembl gene ID
%       'ENST'     Ensembl transcript ID
%       'ENSP'     Ensembl protein ID
%       'UniProt'  UniProt protein ID
%       'Name'     HUGO gene abbreviation
%       'Entrez'   Entrez (NCBI) gene ID
%
% The function will try to determine the original gene ID format, but this 
% can be manually specified with the "origFormat" input. Furthermore, the 
% user can supply a custom conversion key in the "targetFormat" input that 
% uses gene IDs other than the above built-in ID types (see details below).
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
%   [grRules_new,genes,rxnGeneMat] = translateGrRules(grRules,targetFormat,origFormat,noMatch);
%
%
% INPUTS:
%
%   grRules     A cell vector of gene-reaction rules (e.g., model.grRules).
%
%   targetFormat    The desired output format of the grRules. This can be
%                   one or multiple gene ID types, and must be selected
%                   from the following: 'ENSG', 'ENSP', 'ENST', 'UniProt',
%                   'Name', 'Entrez'.
%                   DEFAULT: All 6 formats.
%
%               *** NOTE: As a special case, instead of entering a target
%                   format or list of formats, this input can be replaced
%                   with a custom conversion key. This conversion key must
%                   be an Nx2 cell array, where the first column contains
%                   gene IDs corresponding to the current grRules, and the
%                   second column contains new gene IDs to which the rules
%                   will be translated.
%
%   origFormat  (Optional) A string specifying the original gene ID type 
%               (e.g., 'UniProt'). If not provided, the function will try 
%               to automatically determine the type of IDs present.
%
%   noMatch     What to do in case a gene in the original grRules does not
%               have a corresponding match in the new gene IDs to which the
%               grRules are to be translated. Choose from the following:
%               
%               'delete'    (DEFAULT) Delete (exclude) the gene from the 
%                           new rule.
%
%               'original'  keep the original gene ID. Note, this will
%                           result in mixed gene ID types in grRules and 
%                           the corresponding list of genes, which may
%                           cause problems in other applications.
%
%
% OUTPUTS:
%
%   grRules_new  grRules with the gene IDs converted to the target format.
%                If targetFormat is more than one ID type, then grRules 
%                will be returned as a structure, with a different field
%                for each gene ID type.
%
%   genes        A list of all genes found in the converted grRules.
%                If targetFormat is more than one ID type, then genes will
%                be returned as a structure, with a different field for
%                each gene ID type.
%
%   rxnGeneMat   A matrix associating reactions to genes, build from the
%                converted grRules and genes list. If targetFormat is more
%                than one ID type, then rxnGeneMat will be returned as a 
%                structure, with a different field for each gene ID type.
%


% handle input arguments
if nargin < 4
    noMatch = 'delete';
elseif ~ismember(lower(noMatch),{'original','delete'})
    error('Valid inputs for noMatch are "original" or "delete".');
end

if nargin < 3
    gene_type_orig = {};
else
    gene_type_orig = origFormat;
end

custom_key = false;
if nargin < 2 || isempty(targetFormat)
    % recognized ID types for genes, transcripts, and proteins
    targetFormat = {'ENSG';'ENST';'ENSP';'UniProt';'Name';'Entrez'};
elseif ischar(targetFormat)
    % convert target format from string to cell
    targetFormat = {targetFormat};
elseif (size(targetFormat,2) == 2) && (size(targetFormat,1) > 10)  % check that it has some arbitrarly long length above 10
    % use custom conversion key
    custom_key = true;
    conv_key = targetFormat;
    targetFormat = {'NewFormat'};
    conv_key_head = {'OldFormat','NewFormat'};
    gene_type_orig = 'OldFormat';
else
    % ensure it is a column vector
    targetFormat = targetFormat(:);
end

% allow for non-exact matches to "Name" option
targetFormat(ismember(lower(targetFormat),{'name','names'})) = {'Name'};


% get original list of genes from the grRules
genes_orig = getGenesFromGrRules(grRules);
if ismember('and',genes_orig) || ismember('or',genes_orig)
    error('Problem reading grRules. Verify that all "and" and "or" elements are lowercase and surrounded by spaces.');
end

% preprocess gene list and gene-reaction rules, if necessary
if any(ismember({'GPI','GAPDH'},genes_orig))  % need a more robust method to check if it's gene symbols
    % check if the list is gene names (symbols), if so, do not modify
    rules_orig = grRules;
    if ~custom_key && isempty(gene_type_orig)
        gene_type_orig = 'Name';
    end
else
    % remove ".#" from gene IDs (e.g. ENSG00000198888.2 -> ENSG00000198888, or 2597.1 -> 2597)
    genes_orig = regexprep(genes_orig,'\.\d+$','');
    rules_orig = regexprep(grRules,'\.\d+','');
end

% determine the original gene ID type, if not gene names
if isempty(gene_type_orig)
    if startsWith(genes_orig{1},{'ENSG','ENST','ENSP'})
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

if ~custom_key
    % import gene-transcript-protein conversion key
    
    % get the path
    [ST, I] = dbstack('-completenames');
    path = fileparts(ST(I).file);
    tmpfile = fullfile(path,'../../data/Ensembl','ensembl_ID_mapping.tsv');
    
    fid = fopen(tmpfile);
    tmp = textscan(fid,'%s%s%s%s%s%s','Delimiter','\t','Headerlines',1);
    fclose(fid);
    
    tmp = horzcat(tmp{:});
    conv_key_head = tmp(1,:);  % header
    conv_key = tmp(2:end,:);  % gene IDs
    clear tmp
    
    % change header names to match contents of GENE_TYPES
    [~,ind] = ismember(conv_key_head,{'Gene_stable_ID','Transcript_stable_ID',...
        'Protein_stable_ID','UniProtKB_Swiss_Prot_ID','Gene_name','NCBI_gene_ID'});
    type_abbrevs = {'ENSG';'ENST';'ENSP';'UniProt';'Name';'Entrez'};
    conv_key_head = type_abbrevs(ind);
end

% begin by "cleaning" the original grRules
rules_orig = cleanGrRules(rules_orig);

% convert rules to all other gene ID types
for i = 1:length(targetFormat)
    
    % extract portion of conversion key required; resulting variable will
    % be two columns of IDs, where column 1 is the original gene ID type,
    % and column 2 is the new gene ID type.
    [~,ind] = ismember([{gene_type_orig}, targetFormat(i)], conv_key_head);
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
    if strcmp(targetFormat{i},gene_type_orig)
        % if current model rules are already in the target format, do not
        % convert gene IDs
    else
        % convert gene IDs
        % This next line identifies gene IDs as collections of characters that
        % are not spaces, parentheses, or the symbols & or |. It then replaces
        % those gene IDs with the new gene ID type, which calls the convertGene
        % inline function to retrieve.
        rules_new = regexprep(rules_new, '[^&|\(\) ]+', '(${convertGene($0)})');
    end
    
    % clean up rules (removes extra parentheses, repeated genes, etc.)
    rules_new = cleanGrRules(rules_new);    
    
    % generate new rxnGeneMat and gene list based on converted grRules
    [genes.(targetFormat{i}),rxnGeneMat.(targetFormat{i})] = getGenesFromGrRules(rules_new);
    
    % restore "&" as "and" and "|" as "or"
    rules_new = regexprep(rules_new, ' \| ', ' or ');
    rules_new = regexprep(rules_new, ' & ', ' and ');
    
    % add new rules to model structure
    grRules_new.(targetFormat{i}) = rules_new;
    
end


% if only one target format is specified, structures are not necessary
if length(targetFormat) == 1
    grRules_new = grRules_new.(targetFormat{1});
    if ~isempty(genes)
        genes = genes.(targetFormat{1});
        rxnGeneMat = rxnGeneMat.(targetFormat{1});
    end
end





