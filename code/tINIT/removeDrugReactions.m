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
drugRxns = (sum(inModel.S(drugMets,:) ~= 0,1) > 0).';
outModel = removeReactions(inModel,inModel.rxns(drugRxns));
end
