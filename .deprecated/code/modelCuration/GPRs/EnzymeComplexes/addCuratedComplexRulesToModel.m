function model_new = addCuratedComplexRulesToModel(model,curationFile,delRedundant)
%addCuratedComplexRulesToModel  Incorporate curated grRules into model.
%
%   The CORUM database provides information on proteins that are found to
%   interact (i.e., may form an enzyme complex). These complexes were
%   integrated automatically into HumanGEM using addComplexesToGeneRules.m,
%   but only for grRules with only "OR" relationships. This new set of 
%   complex-incorporated grRules was compared with the original grRules.
%   Those that were different were manually curated to verify which of the
%   newly added complexes (AND relationships) were consistent with
%   literature (e.g., UniProt, NCBI, etc.).
%
%   This function incorporates the new, manually-curated, CORUM-informed 
%   grRules into the model. The incoporation involves a cleaning step,
%   which may remove redundant genes.
%
% USAGE:
%
%   model_new = addCuratedComplexRulesToModel(model,curationFile,delRedundant);
%
% INPUTS:
%
%   model          Model structure.
%
%   curationFile   Name of .txt file containing the curation information.
%                  The file should contain the following three columns
%                  (with the same headers):
%
%                    'rxn' - Rxn identifiers corresponding to each grRule.
%
%                    'grRule_original' - List of original grRules.
%
%                    'grRule_curated' - List of curated grRules to be
%                                       incorporated into the model.
%
%   delRedundant   (Opt, default = TRUE) If TRUE, genes that appear in the 
%                  grRule both individually and as part of a complex will
%                  be revised by deleting the individual form of the gene.
%                  If FALSE, the redundant individual form of the gene will
%                  not be removed from the grRule.
%                  
%              For example:
%                     Curated rule: (G1 and G2) or G1 or G3
%               delRedundant=FALSE: (G1 and G2) or G1 or G3
%               delRedundant= TRUE: (G1 and G2) or G3
%                    
% OUTPUTS:
%
%   model_new    New model structure, where grRules have been updated with
%                the curated grRules containing CORUM complex information.
%


% handle input args
if nargin < 3
    delRedundant = true;
end


% load curated grRule changes (with added enzyme complexes)
fid = fopen(curationFile,'r');
curation_data = textscan(fid,'%s %s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

% check if any reactions are missing from the given model
missing_rxns = find(~ismember(curation_data{1},model.rxns));
if ~isempty(missing_rxns)
    fprintf('WARNING! The following %u rxns are present in the curation file, but not the model:\n',length(missing_rxns));
    fprintf('\t%s\n',curation_data.rxn{missing_rxns});
    fprintf('\n');
end
% remove missing reactions
curation_data{1}(missing_rxns) = [];  
curation_data{2}(missing_rxns) = [];
curation_data{3}(missing_rxns) = [];

% extract information from curation_data
rxns = curation_data{1};
grRule_orig = curation_data{2};
grRule_curated = curation_data{3};

% pre-clean curation_data grRules
grRule_orig = cleanModelGeneRules(grRule_orig);
grRule_curated = cleanModelGeneRules(grRule_curated);

% retrieve and clean model grRules
[~,rxn_ind] = ismember(rxns,model.rxns);
grRule_model = cleanModelGeneRules(model.grRules(rxn_ind));

% check which model rules differ from the "original" rule in curation_data
diff_rules = find(~strcmp(grRule_orig,grRule_model));
if ~isempty(diff_rules)
    fprintf('WARNING! The model.grRules of the following %u rxns do not match the "original" grRules in the curation file:\n',length(diff_rules));
    fprintf('(and will therefore NOT be changed in the model)\n');
    fprintf('\t%s\n',rxns{diff_rules});
    fprintf('\n');
end

% skip rules that differ between model and "original" in curation_data
rxns(diff_rules) = [];
rxn_ind(diff_rules) = [];
grRule_orig(diff_rules) = [];
grRule_curated(diff_rules) = [];
grRule_model(diff_rules) = [];
if isempty(rxns)
    fprintf('*** None of the grRules in the model were changed. ***');
    model_new = model;
    return
end


% remove redundant genes from curated grRules (if specified)
if (delRedundant)
    
    % replace "and" and "or" with symbols (& and |, respectively)
    grRule_curated_sym = grRule_curated;
    grRule_curated_sym = regexprep(grRule_curated_sym,' and ',' & ');
    grRule_curated_sym = regexprep(grRule_curated_sym,' or ',' | ');
    
    % iterate through each of the curated grRules
    for i = 1:length(grRule_curated_sym)
        
        r = grRule_curated_sym{i};
        if contains(r,{'(',')'})
            
            % Identify all genes participating in a complex. This is a very
            % specific case in which complex genes can be identified simply as
            % those enclosed by parentheses. This approach does NOT generalize,
            % and should not be applied to other situations.
            complex_expr = regexp(r,'\([^\)]+\)','match');
            complex_genes = unique(regexp(strjoin(complex_expr),'[^&|\(\) ]+','match'));
            
            % Identify all remaining genes that do not participate in any
            % enzyme complexes
            remaining_expr = regexprep(r,'\([^\)]+\)','');
            indiv_genes = unique(regexp(remaining_expr,'[^&|\(\) ]+','match'));
            
            % remove individual genes that participate in any complexes
            indiv_genes(ismember(indiv_genes,complex_genes)) = [];
            
            % reconstruct grRule
            grRule_curated_sym(i) = join([complex_expr, indiv_genes], ' | ');
            
        end
    end
    
    % restore "and" and "or" phrases
    grRule_curated_sym = regexprep(grRule_curated_sym,' & ',' and ');
    grRule_curated_sym = regexprep(grRule_curated_sym,' \| ',' or ');
    grRule_curated = grRule_curated_sym;
    
end


% replace model grRules with the curated grRules
model_new = model;
model_new.grRules(rxn_ind) = grRule_curated;


