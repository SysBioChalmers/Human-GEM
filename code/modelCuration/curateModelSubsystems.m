%
% FILE NAME:    curateModelSubsystems.m
% 
% PURPOSE: Script to curate the subsystems of humanGEM. The model
%          subsystems were manually evaluated to identify those that should
%          be combined/renamed based on similar naming or function. The
%          results of this evaluation were converted into a .tsv file
%          (subsystem_name_curated.tsv), which contains three columns:
%
%               1. Original subsystem names
%               2. New subsystem names
%               3. Notes/reasoning for the suggested change (if any)
%
%          Some examples include cases where two subsystems have nearly 
%          identical names except for a difference in capitalization or 
%          exact wording, e.g.:
%
%               'Fructose and Mannose metabolism'
%               'Fructose and mannose metabolism'
%
%               'Bile acid synthesis'
%               'Bile acid biosynthesis'
%
%          A larger change involved the transport subsystems, such as:
%
%               'Transport, Golgi apparatus'
%               'Transport, Golgi to extracellular'
%               'Transport, Golgi to lysosome'
%        
%         These transport subsystems will instead be merged into a single
%         subsystem named 'Transport reactions'.
%


% load model
load('humanGEM.mat');

% verify that there are no empty subsystems
if any(cellfun(@isempty,ihuman.subSystems))
    error('There should be no empty subSystem entries!');
end

% load subsystem conversion data
fid = fopen('../../ComplementaryData/modelCuration/subsystem_name_curated.tsv');
subData = textscan(fid,'%s%s%s','Delimiter','\t','HeaderLines',2);
fclose(fid);

% extract columns and remove header
subsys_orig = subData{1}(2:end);
subsys_new = subData{2}(2:end);

% update subsystem names
[hasMatch,ind] = ismember(ihuman.subSystems, subsys_orig);
subSystems = ihuman.subSystems;
subSystems(hasMatch) = subsys_new(ind(hasMatch));

% manually update the subsystem for a set of reactions
% These reactions were originally assigned to:
%   'Beta oxidation of unsaturated fatty acids (n-7) (mitochondrial)'
% but they take place in the peroxisome. They will therefore be renamed:
%   'Beta oxidation of unsaturated fatty acids (n-7) (peroxisomal)'
rxns = {'HMR_3356';'HMR_3357';'HMR_3358';'HMR_3359';'HMR_3360';'HMR_3361';'HMR_3362';'HMR_3363'};
[~,rxnInd] = ismember(rxns,ihuman.rxns);
[~,ind] = ismember('',subsys_orig);
subSystems(rxnInd) = subsys_new(ind);

% update model structure
ihuman.subSystems = subSystems;

% save model
save('../../model/Human-GEM.mat','ihuman');



