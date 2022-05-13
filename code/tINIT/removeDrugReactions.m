function outModel=removeDrugReactions(inModel)
% outModel=removeDrugReactions(inModel)
%
% Removes reactions from Human-GEM related to drugs
% 
% inModel       The model, typically Human-GEM or Mouse-GEM
% outModel      The resulting model.

drugNames = {'pravastatin';'Gliclazide';'atorvastatin';'fluvastatin';'fluvastain';'fluvstatin';'simvastatin';'cyclosporine';...
             'acetaminophen';'cerivastatin';'Tacrolimus';'ibuprofen';'lovastatin';'Losartan';'nifedipine';...
             'pitavastatin';'rosuvastatin';'Torasemide';'Midazolam'};
drugMets = contains(inModel.metNames,drugNames,'IgnoreCase',true);
%sum(drugMets)%456
%inModel.metNames(drugMets)
drugRxns = (sum(inModel.S(drugMets,:) ~= 0,1) > 0).';
%inModel.rxns(drugRxns)
%constructEquations(inModel,inModel.rxns(drugRxns))
%sum(drugRxns)%705
%inModel.grRules(drugRxns) %many of these have GPRs (although most of them don't)

outModel = removeReactions(inModel,inModel.rxns(drugRxns));
end
