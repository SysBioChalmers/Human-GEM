%
% FILE NAME:    addtionalManualCuration.m
% 
% PURPOSE: This script is for manual curation for the reaction
%          association results that accumulate contineously.
%


% 2018-05-28
load('ihumanRxns2BiGG.mat');    %HMR rxn association to BiGG
ihuman.HMR2BiGG{find(strcmp('HMR_1637',ihuman.rxns))}='RE3247M';
ihuman.HMR2BiGG(53:64)={'RE0958C';'RE0958E';'RE0951C';'RE0951E';'RE0944C';'RE0944E';'RE0935C';'RE0935E';'RE0926C';'RE0926E';'RE0915C';'RE0915E'};
save('ihumanRxns2BiGG.mat','ihuman');  % 2018-05-28

load('ihumanRxns2MNX.mat');     %HMR rxn association to MNX
ihuman.HMR2BiGG{find(strcmp('HMR_1637',ihuman.rxns))}='RE3247M';
% Update following BiGG and MNX associations based on subgroup ppt 2018-4-4
ihuman.HMR2BiGG(53:64)={'RE0958C';'RE0958E';'RE0951C';'RE0951E';'RE0944C';'RE0944E';'RE0935C';'RE0935E';'RE0926C';'RE0926E';'RE0915C';'RE0915E'};
ihuman.rxnMNXID{53}{1}='MNXR103496';
ihuman.rxnMNXID{54}{1}='MNXR103497';
ihuman.rxnMNXID{57}{1}='MNXR103495';
ihuman.rxnMNXID{58}{1}='MNXR101623';
ihuman.rxnMNXID{59}{1}='MNXR103491';
ihuman.rxnMNXID{60}{1}='MNXR101622';
ihuman.rxnMNXID{61}{1}='MNXR103488';
ihuman.rxnMNXID{62}{1}='MNXR101621';
ihuman.rxnMNXID{63}{1}='MNXR103479';
ihuman.rxnMNXID{64}{1}='MNXR101620';
save('ihumanRxns2MNX.mat','ihuman');  % 2018-05-28

load('mergedModel.mat');       %merged model struture
mergedModel.confirmedMNXID{find(strcmp('HMR_1637',mergedModel.rxns))}....
=mergedModel.rxnAssocMNXID{find(strcmp('HMR_1637',mergedModel.rxns))};
mergedModel.confirmedFilteredMNXID{find(strcmp('HMR_1637',mergedModel.rxns))}....
=mergedModel.rxnAssocMNXID{find(strcmp('HMR_1637',mergedModel.rxns))};
% Update following MNX associations based on subgroup ppt 2018-4-4
mergedModel.confirmedMNXID{51}{1}='MNXR103496';
mergedModel.confirmedFilteredMNXID{51}{1}='MNXR103496';
mergedModel.confirmedMNXID{56}{1}='MNXR101623';
mergedModel.confirmedFilteredMNXID{56}{1}='MNXR101623';
mergedModel.confirmedMNXID{62}{1}='MNXR101620';
mergedModel.confirmedFilteredMNXID{62}{1}='MNXR101620';
save('mergedModel.mat','mergedModel');  % 2018-05-28
