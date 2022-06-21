function rxns=getAATripletReactions(model, onlyRxnsWithoutGPRs)
% getAATripletReactions
%   Finds all reactions that work with amino acid triplets (or doublets). These are 
%   not very useful in the model, and removing them simplifies the problem.
%
%   model               The model to be examined.
%   onlyRxnsWithoutGPRs If true, the function only gets the rxns without GPRs
%
%   rxns                The resulting reactions
%
%   Usage: rxns=getAATripletReactions(model, onlyRxnsWithoutGPRs)
%

aminoAcidsCap = {'Alanine';'Arginine';'Asparagine';'Aspartate';'Cysteine';'Glutamine';'Glutamate';'Glycine';'Histidine';'Isoleucine';'Leucine';'Lysine';'Methionine';'Metheonine';'Phenylalanine';'Proline';'Serine';'Threonine';'Tryptophan';'Tyrosine';'Valine'};
aminoAcidoyls = {'Alanyl';'Alaninyl';'Alanine';'Arginyl';'Asparaginyl';'Aspartyl';'Cystyl';'Cystinyl';'Cysteinyl';'Glutaminyl';'Glutamyl';'Glutamatsyl';'Glycyl';'Histidyl';'Histidinyl';'Isoleucyl';'Isolecyl';'Leucyl';'Lysyl';'Lysine';'Methionyl';'Methioninyl';'Phenylalanyl';'Phenylalanine';'Phenylalaninyl';'Prolyl';'Seryl';'Threonyl';'Tryptophanyl';'Tyrosyl';'Tyrosinyl';'Valyl'};
singleTripletNames = {'Argtyrval'};

res = cell(length(model.metNames), 1);
for i = 1:length(res)
    res{i} = strsplit(model.metNames{i}, '-');
end

matchingMets = false(length(res),1);
skippedMets = false(length(res),1);
for i = 1:length(matchingMets)
    if any(strcmp(res{i}, singleTripletNames))
        matchingMets(i) = true;
    end
    if length(res{i}) > 1 
        stringList = res{i};
        success = true;
        for j = 1:(length(stringList)-1)
            if ~any(strcmp(aminoAcidoyls, stringList{j}))
               success = false;
               break;
            end
        end
        success = success & any(strcmp(aminoAcidsCap, stringList{length(stringList)}));
        if (success)
            matchingMets(i) = true;
        else
            skippedMets(i) = true;
        end
    end
end

%sum(matchingMets)%490
%sum(skippedMets)%3810

%model.metNames(matchingMets) %looks good

%unique(model.metNames(skippedMets)) %looks good, didn't miss any as far as I could see

%now investigate which reactions these are
rxnSel = (sum(model.S(matchingMets,:) ~= 0,1) > 0).';
%sum(rxnSel)%735
%constructEquations(model, model.rxns(rxnSel))
%model.grRules(rxnSel)
emptyGrRules = cellfun(@isempty, model.grRules(rxnSel));
toIgnoreAminoAcidTriplets = rxnSel;
if onlyRxnsWithoutGPRs
    toIgnoreAminoAcidTriplets(rxnSel) = emptyGrRules;
end
%constructEquations(model, model.rxns(rxnSel & ~toIgnoreAminoAcidTriplets))

%sum(toIgnoreAminoAcidTriplets)%489, to a large extent overlapping with import and exch reactions

rxns = model.rxns(toIgnoreAminoAcidTriplets);
end
