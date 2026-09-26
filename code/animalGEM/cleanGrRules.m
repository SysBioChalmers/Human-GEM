function cleaned_rules = cleanGrRules(grRules)
%cleanGrRules  Clean and simplify model grRules.
%
% cleanGrRules removes unnecessary parentheses, trailing ANDs/ORs, and
% duplicate genes in model grRules. The function requires that gene
% names/IDs do not contain spaces, parentheses, "&", "|", or "#".
% Furthermore, grRules must separate gene names/IDs from ANDs and ORs by at
% least one space. For example,
%
%       OK: "GENE1 and GENE2", or "GENE1 & GENE2"
%   NOT OK: "GENE1andGENE2", or "GENE1&GENE2"
% 
% Note that grRules can use written ("and" "or") or symbolic ("&" "|")
% boolean operators.
%
%
% USAGE:
%
%   cleaned_rules = cleanGrRules(grRules);
%
%
% INPUT:
%
%   grRules         The grRules field from a genome-scale metabolic model.
%
%
% OUTPUT:
%
%   cleaned_rules   Updated/cleaned grRules.
%



% check if the grRules use written or symbolic boolean operators
if any(contains(grRules,{'&','|'}))
    % fix some potential missing spaces between parentheses and &/|
    grRules = regexprep(grRules,'\)&',') &');   % ")&"  ->  ") &"
    grRules = regexprep(grRules,'&\(','& (');   % "&("  ->  "& ("
    grRules = regexprep(grRules,'\)\|',') |');  % ")|"  ->  ") |"
    grRules = regexprep(grRules,'\|\(','| (');  % "|("  ->  "| ("
    Btype = 'symbol';  % record rule type, to revert back at the end
else
    % fix some potential missing spaces between parentheses and AND/OR
    grRules = regexprep(grRules,'\)and',') and');  % ")and" ->  ") and"
    grRules = regexprep(grRules,'and\(','and (');  % "and(" ->  "and ("
    grRules = regexprep(grRules,'\)or',') or');    % ")or"  ->  ") or"
    grRules = regexprep(grRules,'or\(','or (');    % "or("  ->  "or ("
    
    % convert "and" to "&" and "or" to "|" (easier to work with symbols)
    grRules = regexprep(grRules, ' or ', ' | ');
    grRules = regexprep(grRules, ' and ', ' & ');
    Btype = 'text';  % record rule type, to revert back at the end
end


% iterate through the cleaning process until the grRules no longer change
grRules_orig = grRules;
ambiguous_ind = [];  % initialize variable
for i = 1:50  % perform a max of 50 iterations (it shouldn't require nearly that many)
    
    % remove extra parentheses
    has_AND_OR = contains(grRules,'&') & contains(grRules,'|');
    grRules(~has_AND_OR) = regexprep(grRules(~has_AND_OR),'\(|\)','');  % remove all parentheses from rules that don't have both ANDs and ORs
    grRules = regexprep(grRules,'\((\S*)\)','$1');  % remove extra parentheses: orphan '()' or those surrounding only one gene '(GENE1)'
    grRules = regexprep(grRules,'\((\([^\(\)]*\))\)','$1');  % remove double parentheses around expressions: '((GENE1 & GENE2))' or '((GENE1 | GENE2))'
    
    % remove repeated ORs and ANDs
    grRules = regexprep(grRules,'\|(\s+\|)+','|');  % remove repeated ORs
    grRules = regexprep(grRules,'&(\s+&)+','&');  % remove repeated ANDs
    
    % remove trailing ORs and ANDs
    grRules = regexprep(grRules,'\s*\|\s*(\))|(\()\s*\|\s*','$1');  % remove orphan ORs
    grRules = regexprep(grRules,'^\s*\|\s*|\s*\|\s*$','');  % remove ORs at the beginning or end of rule
    grRules = regexprep(grRules,'\s*&\s*(\))|(\()\s*&\s*','$1');  % remove orphan ANDs
    grRules = regexprep(grRules,'^\s*&\s*|\s*&\s*$','');  % remove ANDs at the beginning or end of rule
    
    % Remove unnecessary parentheses grouping OR or AND expressions.
    % e.g.: '((GENE1 or GENE2) or GENE3)'   -> '(GENE1 or GENE2 or GENE3)
    grRules = regexprep(grRules,'(?<=^\s*\(*)\(([^&\(\)]+)\)(?=\)*\s*\|)','$1');  % ORs, beginning of rule
    grRules = regexprep(grRules,'(?<=\|\s*\(*)\(([^&\(\)]+)\)(?=\)*\s*$)','$1');  % ORs, end of rule
    grRules = regexprep(grRules,'(?<=\|\s*\(*)\(([^&\(\)]+)\)(?=\)*\s*\|)','$1'); % ORs, middle of rule
    % e.g.: '((GENE1 and GENE2) and GENE3)' -> '(GENE1 and GENE2 and GENE3)
    grRules = regexprep(grRules,'(?<=^\s*\(*)\(([^\|\(\)]+)\)(?=\)*\s*&)','$1');  % ANDs, beginning of rule
    grRules = regexprep(grRules,'(?<=&\s*\(*)\(([^\|\(\)]+)\)(?=\)*\s*$)','$1');  % ANDs, end of rule
    grRules = regexprep(grRules,'(?<=&\s*\(*)\(([^\|\(\)]+)\)(?=\)*\s*&)','$1');  % ANDs, middle of rule
    
    
    % Look for duplicated genes that are within a group of all ORs or all
    % ANDs, and keep only the unique set of genes. Two examples:
    % '(GENE1 or GENE2 or GENE1) and (GENE4 or GENE4)' -> '(GENE1 or GENE2) and GENE4'
    % 'GENE1 and GENE2 and GENE2 and GENE2' -> 'GENE1 and GENE2'
    
    % define functions for replacing rule piece with unique set of genes
    keepUniques_AND = @(str,sep) strjoin(unique(regexp(str,'[^&|\(\) ]+','match')),' & ');
    keepUniques_OR = @(str,sep) strjoin(unique(regexp(str,'[^&|\(\) ]+','match')),' | ');
    
    % remove duplicates from groups of genes separated by only ANDs or ORs,
    % with no parentheses in between
    grRules = regexprep(grRules,'[^&|\(\) ]+( & [^&|\(\) ]+)+','${keepUniques_AND($0)}');
    grRules = regexprep(grRules,'[^&|\(\) ]+( \| [^&|\(\) ]+)+','${keepUniques_OR($0)}');

    
    %......................................................................
    % Perform more extensive grRule cleaning. This process identifies and
    % removes duplicated groups of genes, e.g.:
    %   '(G1 and G2) or (G2 and G1) or G3' --> '(G1 and G2) or G3 or G4'
    % as well as duplicate genes that are separated by another group, e.g.:
    %   'G1 and G2 and (G3 or G4) and G1' --> '(G3 or G4) and G1 and G2'
    % This process should even work on rules with many layers of nesting.
    
    % identify all rules containing both ANDs and ORs
    AndOrInd = find(contains(grRules,'&') & contains(grRules,'|'));
    
    % iterate through all rules containig both ANDs and ORs
    for ii = 1:length(AndOrInd)
        
        % extract grRule
        r = grRules{AndOrInd(ii)};
        
        % Specify phrases to search for in the grRule. These phrases will find
        % genes grouped by all ANDs (first phrase) or all ORs (second phrase).
        search_phrases = {'\([^&|\(\) ]+( & [^&|\(\) ]+)+\)','\([^&|\(\) ]+( \| [^&|\(\) ]+)+\)'};
        
        % initialize some variables
        chunks = {};  % this will contain groups of genes (separated by ANDs or ORs) that will be extracted from the rule
        c = 1;  % this is a counter keeping track of the group (chunk) number
        r_orig = r;  % keep track of the original rule, to see when it stops changing
        for k = 1:100  % iterate some arbitrarily high number of times
            for j = 1:length(search_phrases)  % search both phrases
                new_chunk = regexp(r,search_phrases{j},'match')';  % extract chunks
                if ~isempty(new_chunk)
                    chunks = [chunks; new_chunk];  % append to list of chunks
                    group_nums = arrayfun(@num2str,(c:length(chunks))','UniformOutput',false);  % get group numbers to be assigned to the new chunks, and convert to strings
                    r = regexprep(r,search_phrases{j},strcat('#',group_nums,'#'),'once');  % replace the chunks in the expression with their group numbers (enclosed by "#"s)
                    c = c + length(new_chunk);  % increase counter
                end
            end
            if isequal(r,r_orig)
                % if the grRule is no longer changing, then we can stop iterating
                break;
            else
                r_orig = r;
            end
        end
        
        % append the final state of the grRule as the last chunk
        chunks{end+1} = r;
        
        % Check whether an AND and OR appear in the final state of the 
        % grRule in an ambiguous manner (i.e., not separated by a parenthesis)
        if ~isempty(regexp(r,'\|[^\(\)]+&|&[^\(\)]+\|','once'))
            ambiguous_ind = [ambiguous_ind; AndOrInd(ii)];
        end
        
        % remove duplicates from chunks, iterating until no more changes
        chunks_orig = chunks;
        for k = 1:100
            
            % remove duplicate genes from within each chunk
            chunks = regexprep(chunks,'[^&|\(\) ]+( & [^&|\(\) ]+)+','${keepUniques_AND($0)}');
            chunks = regexprep(chunks,'[^&|\(\) ]+( \| [^&|\(\) ]+)+','${keepUniques_OR($0)}');
            
            % check if any of the chunks themselves are duplicated
            [~,~,ind] = unique(chunks,'stable');
            if max(ind) < length(chunks)  % this happens when at least one chunk is non-unique
                rep_ind = find(histcounts(ind,'BinMethod','integers') > 1);  % find the indices of the repeated chunks
                for j = 1:length(rep_ind)
                    nums = strcat('#',arrayfun(@num2str,find(ind == rep_ind(j)),'UniformOutput',false),'#');
                    chunks = regexprep(chunks,nums,nums{1});
                end
                
                % remove duplicate genes from within each chunk again
                chunks = regexprep(chunks,'[^&|\(\) ]+( & [^&|\(\) ]+)+','${keepUniques_AND($0)}');
                chunks = regexprep(chunks,'[^&|\(\) ]+( \| [^&|\(\) ]+)+','${keepUniques_OR($0)}');
            end
            
            if isequal(chunks,chunks_orig)
                % if no more changes occur, then exit loop
                break;
            else
                chunks_orig = chunks;
            end
            
        end
        
        % now re-insert the chunks into their appropriate location in the rule
        r = chunks{end};  % obtain updated final state of grRule
        for j = c-1:-1:1
            r = regexprep(r,strcat('#',num2str(j),'#'),chunks{j});
        end
        
        % update grRule
        grRules{AndOrInd(ii)} = r;
        
    end
    %......................................................................
    
    
    % determine if the grRules changed (in length); if not, exit the loop
    if all(cellfun(@length,grRules) == cellfun(@length,grRules_orig))
        % Note: the length of the rules, instead of the rules themselves,
        % are compared, because sometimes the function will get stuck in a
        % loop of simply rearranging rule pieces without actually
        % simplifying them further.
        break;
    else
        grRules_orig = grRules;
        if i == 50
            warning('After 50 iterations, grRules still changed. There is likely a problem with grRules or this function.');
        end
    end
    
end

% notify user of any ambiguous grRules (AND and OR expressions that are not
% separated by parentheses)
if ~isempty(ambiguous_ind)
    fprintf('\n***The following grRules contain ambiguous AND/OR combination(s) due to insufficient parentheses:\n');
    fprintf('\t%u\n',unique(ambiguous_ind));
    fprintf('\n');
end

% change boolean operators back to original type
if strcmpi(Btype,'text')
    grRules = regexprep(grRules, ' \| ', ' or ');
    grRules = regexprep(grRules, ' & ', ' and ');
end

% assign output
cleaned_rules = grRules;




