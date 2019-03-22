% FILE NAME:    issue71.m
%
% DATE STARTED: 2019-03-04
%     MODIFIED: 2019-03-22
%
% PROGRAMMERS:  Hao Wang, Pinar Kocabas
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
%
% PURPOSES:     Add Metadata to the annotation field of Human1

%% 1. Load the model
if ~exist('ihuman','var')
    load('humanGEM.mat');  % version 1.0.0
end

%% 2. Add Metadata to the annotation field of Human1

% Generate the Metadata as a structure
annotation.defaultLB=-1000;
annotation.defaultUB=1000;
annotation.taxonomy='9606';
annotation.note='Human genome-scale metabolic models are important tools for the study of human health and diseases, by providing a scaffold upon which different types of data can be analyzed. This is the latest version of human-GEM, which is a genome-scale model of the generic human cell. The objective of human-GEM is to serve as a community model for enabling integrative and mechanistic studies of human metabolism.';
annotation.sourceUrl='https://github.com/SysBioChalmers/human-GEM';
annotation.authors='Jonathan Robinson, Hao Wang, Pierre-Etienne Cholley';
annotation.email='nielsenj@chalmers.se';
annotation.organization='Chalmers University of Technology';

ihuman.annotation = annotation;
ihuman.description = 'Genome-scale model of the generic human cell';

%% 3. Save the model file
save('../../ModelFiles/mat/humanGEM.mat','ihuman');