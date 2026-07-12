function [newModel, modelChanges] = addMetabolicNetwork(model,rxnsToAdd,metsToAdd)
%addMetabolicNetwork  Integrate a new metabolic network into a model
%
% NOTE: The function is to incorporate a subnetwork that is provided by two
%       arrays, rxnsToAdd and metsToAdd, which should be in defined 
%       structure as below，and with reactions and metabolites that are
%       NOT included in the model
%
% Usage:
%   [newModel, modelChanges] = addMetabolicNetwork(model,rxnsToAdd,metsToAdd)
%
%
% Inputs:
%   model        A model structure
%
%   rxnsToAdd    the reaction structure can have the following fields:
%            rxns               cell array with unique strings that
%                               identifies each reaction
%            equations          cell array with equation strings. Decimal
%                               coefficients are expressed as "1.2".
%                               Reversibility is indicated by "<=>" or "=>"
%                               The metabolites have to be written as
%                               "metNames[comps]". Only compartments in
%                               model.comps are allowed. New metabolites to
%                               the model should be added to metsToAdd array
%            rxnNames           cell array with the names of each reaction
%                               (opt, default '')
%            lb                 vector with the lower bounds (opt, default
%                               to model.annotations.defaultLB or -inf for
%                               reversible reactions and 0 for irreversible
%                               when "equations" is used. When "mets" and
%                               "stoichCoeffs" are used it defaults for all
%                               reactions to model.annotations.defaultLB or
%                               -inf)
%            ub                 vector with the upper bounds (opt, default
%                               to model.annotations.defaultUB or inf)
%            eccodes            cell array with the EC-numbers for each
%                               reactions. Delimit several EC-numbers with
%                               ";" (opt, default '')
%            subSystems         cell array with the subsystems for each
%                               reaction (opt, default '')
%            grRules            cell array with the gene-reaction
%                               relationship for each reaction. For example
%                               "(A and B) or (C)" means that the reaction
%                               could be catalyzed by a complex between
%                               A & B or by C on its own. All the genes
%                               have to be present in model.genes. Add
%                               genes with addGenesRaven before calling
%                               this function if needed (opt, default '')
%            rxnReferences      cell array with reaction references (opt,
%                               default '')
%            rxnConfidenceScores   vector with reaction confidence scores
%                               (opt, default NaN)
%
%   metsToAdd   the metabolite structure can have the following fields:
%            mets              cell array with unique strings that
%                              identifies each metabolite (opt, default IDs
%                              of new metabolites are numbered with the
%                              prefix defined below)
%            metNames          cell array with the names of each
%                              metabolite
%            compartments      cell array with the compartment of each
%                              metabolite. Should match model.comps.
%                              If this is a string rather than a cell
%                              array it is assumed that all mets are in
%                              that compartment
%            inchis            cell array with InChI strings for each
%                              metabolite (opt, default '')
%            metFormulas       cell array with the formulas for each of
%                              the metabolites (opt, default '')
%            metCharges        metabolite charge (opt, default NaN)
%
%
% Outputs:
%   newModel      a new model with the added reactons and metabolites
%
%   modelChanges  a cell structure containing the original and new reaction
%                 and/or metabolite properties. Reaction-related changes
%                 are reported in the modelChanges.rxns field, whereas
%                 metabolite-related changes are reported in the
%                 modelChanges.mets field
%

if nargin<2
    error('Missing input!');
end

if nargin<3
    metsToAdd = [];
end



%% Add new metabolites

if isempty(metsToAdd)
    newModel = model;
else
    % the following metabolite fields are required to be included
    requiredMetFields = {'mets','metNames','metFormulas','metCharges','compartments'};
    if ~all(ismember(requiredMetFields, fieldnames(metsToAdd)))
        error('There is missing fields in the metabolite structure to be added.');
    end

    % verify that none of the metabolites exist in the current model
    if any(ismember(metsToAdd.mets, model.mets))
        error('One or more metabolite IDs to be added already exist in the model.');
    end

    % add metabolites to the model
    newModel = addMets(model,metsToAdd);
end



%% Check new reactions

% the following reaction fields are required to be included
requiredRxnFields = {'rxns','equations','subSystems','grRules'};
if ~all(ismember(requiredRxnFields, fieldnames(rxnsToAdd)))
    error('There is missing fields in the metabolite structure to be added.');
end

% verify that none of the reactions exist in the current model
if any(ismember(rxnsToAdd.rxns, model.rxns))
    error('One or more reactions to be added already exist in the model.');
end

% make sure all equations have metabolites in a format of "metNames[comps]"
if ~all(contains(rxnsToAdd.equations, '[') & contains(rxnsToAdd.equations, ']'))
    error('The metabolites in equations should be in a format of "metNames[comps]".');
end



%% Add new genes 

% extract genes from grRules
newGenes = getGenesFromGrRules(rxnsToAdd.grRules);

% remove genes that are already existed in model
if any(ismember(newGenes, newModel.genes))
    fprintf('One or more genes to be added already exist in the model，therefore they will be merged.\n');
    newGenes = setdiff(newGenes, newModel.genes);
end

% append new genes and their names to model
newModel.genes = [newModel.genes; newGenes];
emptyGeneNames = newGenes;
emptyGeneNames(:) = {''};
newModel.geneShortNames = [newModel.geneShortNames; emptyGeneNames];

% add new columns to rxnGeneMat will be updated after the new reactions are added below.
newModel.rxnGeneMat(:, end+1:end+numel(newGenes)) = 0;



%% Integrate new subnetwork to the model
newModel = addRxns(newModel, rxnsToAdd, 3);


% update genes and rxnGeneMat with new grRules
[newModel.genes, newModel.rxnGeneMat] = getGenesFromGrRules(newModel.grRules);



%% Document model changes
modelChanges = docModelChanges(model, newModel);


