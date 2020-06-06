function [complex_rules,potential_assoc] = addComplexesToGeneRules(grRules,complex_genes,complex_mat)
%addComplexesToGeneRules  Integrate enzyme complex info into model grRules.
%
% addComplexesToGeneRules incorporates enzyme complex information into
% model grRules by adding "AND" relationships between genes that belong to
% the same complex.
% 
% The existing grRules will be used only to determine which genes are 
% associated with each reaction; i.e., any existing ANDs/ORs etc. will be 
% removed, and the grRule will be regenerated based solely on the enzyme
% complex information.
%
% If an enzyme complex contains some genes that are associated with a
% reaction in addition to some genes that are NOT associated with the
% reaction, it will only group the genes associated with the reaction, and
% will not add the other genes to the association. For example:
%
%      grRule orig: GENE1 or GENE2 or GENE3
%
%   enzyme complex: GENE2, GENE3, GENE4
%
%       grRule new: GENE1 or (GENE2 and GENE3)
%
% In the above example, GENE4 is not included in the new grRule, but it
% will be included in the "potential_assoc" output, which indicates genes
% that are potentially associated with each reaction based on the enzyme
% complex information.
%
%
% USAGE:
%
%   [complex_rules,potential_assoc] = addComplexesToGeneRules(grRules,complex_genes,complex_mat);
%
% INPUT:
%
%   grRules         Gene-reaction (or protein-rxn, transcript-rxn, etc.)
%                   rules from a genome-scale metabolic model structure.
%                   The gene IDs must not contain any of the following
%                   characters: spaces, &, |, #. Logical operators can be
%                   lower-case text ("and","or") or symbols ("&","|").
%
%   complex_genes   List of gene IDs corresponding to the rows of
%                   complex_mat, and appearing in grRules. complex_genes
%                   can also contain genes that are not present in grRules.
%
%   complex_mat     A binary GxC matrix, where G is the number of genes and
%                   C is the number of complexes. Rows correspond to
%                   complex_genes, and each column represents an enzyme
%                   complex, where a "1" indicates that a gene belongs to
%                   that complex, and "0" otherwise.
%
% OUTPUT:
%
%   complex_rules    Gene-reaction rules with enzyme complex information
%                    incorporated in the form of added AND relationships.
%
%   potential_assoc  A list of potential gene associations for each
%                    reaction. If an enzyme complex includes some genes
%                    that are present in the grRule for a reaction, but
%                    also some genes that are not present, the genes that
%                    are not present in the grRule will not be added to the
%                    grRule, but will instead be added as a potential_assoc
%                    for that reaction.
%


% save original grRules
grRules_orig = grRules;

% check if the grRules use written or symbolic boolean operators
if any(contains(grRules,{'&','|'}))
    % fix some potential missing spaces between parentheses and operators
    grRules = regexprep(grRules,'\)&',') &');   % ")&"  ->  ") &"
    grRules = regexprep(grRules,'&\(','& (');   % "&("  ->  "& ("
    grRules = regexprep(grRules,'\)\|',') |');  % ")|"  ->  ") |"
    grRules = regexprep(grRules,'\|\(','| (');  % "|("  ->  "| ("
    Btype = 'symbol';  % record rule type, to revert back at the end
else
    % fix some potential missing spaces between parentheses and operators
    grRules = regexprep(grRules,'\)and',') and');  % ")and" ->  ") and"
    grRules = regexprep(grRules,'and\(','and (');  % "and(" ->  "and ("
    grRules = regexprep(grRules,'\)or',') or');    % ")or"  ->  ") or"
    grRules = regexprep(grRules,'or\(','or (');    % "or("  ->  "or ("
    
    % convert "and" to "&" and "or" to "|" (easier to work with symbols)
    grRules = regexprep(grRules, ' or ', ' | ');
    grRules = regexprep(grRules, ' and ', ' & ');
    Btype = 'text';  % record rule type, to revert back at the end
end


% identify genes appearing in each rule
rule_genes = cellfun(@(r) unique(regexp(r,'[^&|\(\) ]+','match')),grRules,'UniformOutput',false);

% get list of all genes appearing in rules
nonEmpty = find(~cellfun(@isempty,rule_genes));
all_genes = unique([rule_genes{nonEmpty}]');

% quick check for incorrect gene IDs
if ~any(ismember(all_genes,complex_genes))
    fprintf('None of the complex_genes were found in grRules.\n');
    fprintf('Verify that the gene IDs/names are of the same type.\n');
    fprintf('The original grRules will be returned, unmodified.\n');
    complex_rules = grRules_orig;
    potential_assoc = {};
    return
end

% iterate through each (non-empty) rule
potential_assoc = repmat({''},size(grRules));
for i = 1:length(nonEmpty)
    
    % get grRule index
    rule_ind = nonEmpty(i);
    
    % determine which grRule genes are in complex_genes, and their indices
    [is_cpx,gene_ind] = ismember(rule_genes{rule_ind},complex_genes);
    gene_ind(~is_cpx) = [];  % remove zero-indices
    
    if any(is_cpx)
        % find all complexes that include any genes in the current grRule
        cpx_ind = any(complex_mat(gene_ind,:),1);
        
        % identify all genes that are members of any of these complexes but are
        % NOT included in the grRule, and add them as suggested gene associations
        tmp_mat = complex_mat(:,cpx_ind);
        tmp_mat(gene_ind,:) = 0;
        sugg_ind = any(tmp_mat,2);
        potential_assoc(rule_ind) = join(complex_genes(sugg_ind),'; ');
        
        % extract subset of complex_mat including only genes in the current
        % grRule and their associated complexes
        cpx_mat_sub = complex_mat(gene_ind,cpx_ind);
        cpx_genes_sub = complex_genes(gene_ind);
        
        % remove non-unique gene complexes from matrix subset
        cpx_mat_sub = unique(cpx_mat_sub','rows')';
        
        % collect grRule components
        rule_pieces = rule_genes{rule_ind}(~is_cpx);
        for j = 1:size(cpx_mat_sub,2)
            if sum(cpx_mat_sub(:,j)) > 1
                % join complex subunits with &s, and enclose in parentheses
                rule_pieces(end+1) = strcat('(',join(cpx_genes_sub(cpx_mat_sub(:,j) > 0),' & '),')');
            else
                % only one of the complex subunits was present in the grRule,
                % so only that subunit will be added by itself
                rule_pieces(end+1) = cpx_genes_sub(cpx_mat_sub(:,j) > 0);
            end
        end
        
    else
        
        % genes not in complexes
        rule_pieces = rule_genes{rule_ind}(~is_cpx);
        
    end
    
    % generate grRule by joining complexes and single genes with |s
    grRules(rule_ind) = join(rule_pieces,' | ');
    
end


% change boolean operators back to original type
if strcmpi(Btype,'text')
    grRules = regexprep(grRules, ' \| ', ' or ');
    grRules = regexprep(grRules, ' & ', ' and ');
end


% assign output
complex_rules = grRules;









