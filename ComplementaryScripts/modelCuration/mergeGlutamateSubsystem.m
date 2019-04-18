%
% FILE NAME:    mergeGlutamateSubsystem.m
% 
% DATE CREATED: 2019-04-05
%     MODIFIED: 2019-04-11
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Currently, there is only one reaction in the "Glutamate
%          metabolism" subsystem. The subsystem for this reaction will be
%          changed to "Alanine, aspartate and glutamate metabolism", which
%          contains 39 reactions.
%


% load model
load('humanGEM.mat');

% convert subsystems to cell array
subSystems = cellfun(@(x) x, ihuman.subSystems);

% find reaction in Glutamate metabolism subsystem
rxnInd = ismember(subSystems,'Glutamate metabolism');

% change subsystem to Alanine, aspartate and glutamate metabolism
ihuman.subSystems{rxnInd} = {'Alanine, aspartate and glutamate metabolism'};

% export model
exportHumanGEM(ihuman,'humanGEM','../../',{'mat','yml'},false,false);


