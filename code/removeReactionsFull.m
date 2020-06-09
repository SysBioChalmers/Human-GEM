function reducedModel = removeReactionsFull(model,rxnsToRemove,removeUnusedMets,removeUnusedGenes,removeUnusedComps)
% removeReactionsFull
%   Deletes a set of reactions from a model. The function also tries to
%   predict which model fields are associated with reaction properties, and
%   will update all such fields accordingly.
%
%   model             a model structure
%   rxnsToRemove      either a cell array of reaction IDs, a logical vector
%                     with the same number of elements as reactions in the model,
%                     or a vector of indexes to remove
%   removeUnusedMets  remove metabolites that are no longer in use (opt,
%                     default false)
%   removeUnusedGenes remove genes that are no longer in use (opt, default
%                     false)
%   removeUnusedComps remove compartments that are no longer in use (opt,
%                     default false)
%
%   reducedModel      an updated model structure
%
%   Usage: reducedModel=removeReactions(model,rxnsToRemove,removeUnusedMets,...
%           removeUnusedGenes,removeUnusedComps)
%


if nargin<3
    removeUnusedMets=false;
end
if nargin<4
    removeUnusedGenes=false;
end
if nargin<5
    removeUnusedComps=false;
end

if ischar(rxnsToRemove)
    rxnsToRemove={rxnsToRemove};
end

% initialize output
reducedModel=model;

if ~isempty(rxnsToRemove) || removeUnusedMets || removeUnusedGenes
    indexesToDelete=getIndexes(reducedModel,rxnsToRemove,'rxns');
    
    % remove reactions
    if ~isempty(indexesToDelete)
        
        % delete entries from fields predicted to be associated with rxns
        reducedModel = delModelFields(reducedModel,'rxns',indexesToDelete);
        
    end
    
    % remove unused metabolites
    if removeUnusedMets
        if isfield(reducedModel,'S')
            
            % find unused mets
            delMetInds = all(reducedModel.S == 0,2);
            
            % delete unused mets and their corresponding entries in related fields
            reducedModel = delModelFields(reducedModel,'mets',delMetInds);
        else
            error('Could not remove unused metabolites without an S matrix.');
        end
    end
    
    % remove unused genes
    if removeUnusedGenes
        if isfield(reducedModel,'rxnGeneMat')
            
            % find all genes that are not used
            delGeneInds = all(reducedModel.rxnGeneMat == 0,1);
            
            % delete unused genes and their corresponding entries in related fields
            reducedModel = delModelFields(reducedModel,'genes',delGeneInds);
        else
            error('Could not remove unused genes without a "rxnGeneMat" field.');
        end
    end
    
    % remove unused comps
    if removeUnusedComps
        
        usedComps = [];  % initialize variable
        if isfield(reducedModel,'geneComps')
            usedComps = unique([usedComps; reducedModel.geneComps]);
            geneComps_char = reducedModel.comps(reducedModel.geneComps);
        end
        if isfield(reducedModel,'rxnComps')
            usedComps = unique([usedComps; reducedModel.rxnComps]);
            rxnComps_char = reducedModel.comps(reducedModel.rxnComps);
        end
        if isfield(reducedModel,'metComps')
            usedComps = unique([usedComps; reducedModel.metComps]);
            metComps_char = reducedModel.comps(reducedModel.metComps);
        end
        
        if isempty(usedComps)
            error('Could not remove unused compartments without a "geneComps", "rxnComps", or "metComps" field.');
        end
        
        delCompInds = ~ismember(1:length(reducedModel.comps),usedComps);
        if any(delCompInds)
            
            % delete unused compartments and their corresponding entries in related fields
            reducedModel = delModelFields(reducedModel,'comps',delCompInds);
            
            % update the "geneComps", "rxnComps", and/or "metComps" fields,
            % because they contain indices that will be erroneous if any
            % compartments were removed.
            if isfield(reducedModel,'geneComps')
                [~,reducedModel.geneComps] = ismember(geneComps_char,reducedModel.comps);
            end
            if isfield(reducedModel,'rxnComps')
                [~,reducedModel.rxnComps] = ismember(rxnComps_char,reducedModel.comps);
            end
            if isfield(reducedModel,'metComps')
                [~,reducedModel.metComps] = ismember(metComps_char,reducedModel.comps);
            end
            
        end
        
    end
    
end


end  % end removeReactionsFull




function delModel = delModelFields(model,delField,delInds)
%delModelFields  Delete entries of field and associated fields from model.
%
% The function finds all fields in the model related to the given delField,
% and deletes the entries corresponding to the provided delInds.
%
% INPUTS:
%
%   model        Model structure.
%
%   delField     Either 'rxns', 'mets', 'genes', or 'comps'.
%
%   delInds      Indices of field entries to delete.
%
%
% OUTPUTS:
%
%   delModel     A model with the indices of the given field deleted, as
%                well as entries of related fields deleted.
%

if isfield(model,'proteins')
    if unique([length(model.rxns),length(model.mets),length(model.genes),length(model.comps),length(model.proteins)]) < 5
        error('This algorithm will cause problems if rxns, mets, genes, comps, and proteins are not unequal in length.');
    end
else
    if unique([length(model.rxns),length(model.mets),length(model.genes),length(model.comps)]) < 4
        error('This algorithm will cause problems if rxns, mets, genes, and comps are not unequal in length.');
    end
end

switch delField
    case 'rxns'
        n = length(model.rxns);
    case 'mets'
        n = length(model.mets);
    case 'genes'
        n = length(model.genes);
    case 'comps'
        n = length(model.comps);
    otherwise
        error('Invalid field_type.');
end

f = fieldnames(model);
for i = 1:length(f)
    [d1,d2] = size(model.(f{i}));
    if (n == d1)
        model.(f{i})(delInds,:) = [];
    elseif (n == d2)
        model.(f{i})(:,delInds) = [];
    end  
end

delModel = model;

end  % end delModelFields



