function steps = getHumanGEMINITSteps(series)
% getHumanGEMINITSteps
%   Wrapper for getINITSteps. Currently only passes on the call to RAVEN.
%
%   series        See getINITSteps in RAVEN (opt, default 'default')
%
%   steps         Cell array of steps, used as input to ftINIT
%
%   Usage: steps = getHumanGEMINITSteps(series)

if nargin < 1
    steps = getINITSteps();
else
    steps = getINITSteps([], series);
end

%Recent tests show that we don't really need metsToIgnore, but we leave them in the code
%if anyone wants to experiment with them:
%We ignore these mets in the two first steps, where the import and transport reactions without GPR is there anyway.
%We then force them on in the third step.
%metsToIgnore.simpleMets.mets = {'H2O';'Pi';'PPi';'H+';'O2';'CO2';'Na+'};
%metsToIgnore.simpleMets.compsToKeep = {'i'};%the H+ in the 'i' compartment may matter

%steps = getINITSteps([], series);
%steps{1}.MetsToIgnore = metsToIgnore;
%steps{2}.MetsToIgnore = metsToIgnore;
%To only add these to step 1 and 2 can lead to failure when running the third step,
%since one of these metabolites can hinder some fluxes.

end
