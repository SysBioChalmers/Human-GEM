%
% FILE NAME:    miscModelCurationScript_new.m
% 
% DATE CREATED: 2018-09-21
%     MODIFIED: 2018-09-21
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script for performing a number of different curations, updates,
%          and corrections to HumanGEM.
%




%% Treatment of ubiquinone and FAD+
% The model essentially treats FAD/FADH2 the same as ubiquinone/ubiquinol,
% as exhibited by the presence of the following reaction:
%  HMR_6911: FADH2[m] + ubiquinone[m] <=> FAD[m] + ubiquinol[m]
%
% Some reactions in the model are duplicated such that one version of the
% reaction uses FAD, whereas the other uses ubiquinone. These are
% effectively identical reactions, and one should be removed (unless
% evidence suggests otherwise). 

% (Not yet implemented)

