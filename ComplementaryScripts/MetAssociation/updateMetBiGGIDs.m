%
% FILE NAME:    updateMetBiGGIDs.m
% 
% DATE CREATED: 2019-08-08
%      UPDATED: 2019-08-08
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: 


% load HumanGEM and met association file
load('humanGEM.mat');
metAssoc = jsondecode(fileread('humanGEMMetAssoc.JSON'));

% retrieve metabolite data from BiGG database and convert to structure
webdata = webread('http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt');
webdata = textscan(webdata,repmat('%s',1,6));  % file contains 6 columns of text
bigg = {};
for i = 1:numel(webdata)
    bigg.(webdata{i}{1}) = webdata{i}(2:end);
end

% Some BiGG IDs contain the compartment abbrev, and therefore match the
% "bigg_id", but not the "universal_bigg_id". These IDs need to be updated
% to the universal IDs.
[hasMatch, matchInd] = ismember(metAssoc.metBiGGID, bigg.bigg_id);
metAssoc.metBiGGID(hasMatch) = bigg.universal_bigg_id(matchInd(hasMatch));

% determine which BiGG IDs are not found in "universal BiGG IDs"
badInd = find(~ismember(metAssoc.metBiGGID, bigg.universal_bigg_id) & ~cellfun(@isempty, metAssoc.metBiGGID));

% most of these can be fixed by adding an additional underscore to the ID
modifiedID = regexprep(metAssoc.metBiGGID(badInd), '_', '__');
newMatch = ismember(modifiedID, bigg.universal_bigg_id);
metAssoc.metBiGGID(badInd(newMatch)) = modifiedID(newMatch);

% the remaining unmatched IDs are manually updated as below:
IDconv = {'(alpha_D_Mannosyl)2_beta_D_mannosyl_diacetylchitobiosyldiphosphodolichol', 'm2mpdol'
          '(alpha_D_Mannosyl)3_beta_D_mannosyl_diacetylchitobiosyldiphosphodolichol', 'm3mpdol'
          '(alpha_D_Mannosyl)5_beta_D_mannosyl_diacetylchitobiosyldiphosphodolichol', 'm5mpdol'
          '(alpha_D_Mannosyl)6_beta_D_mannosyl_diacetylchitobiosyldiphosphodolichol', 'm6mpdol'
          '(alpha_D_Mannosyl)7_beta_D_mannosyl_diacetylchitobiosyldiphosphodolichol', 'm7mpdol'
          '(alpha_D_Mannosyl)8_beta_D_mannosyl_diacetylchitobiosyldiphosphodolichol', 'm8mpdol'
          'formcoa', 'forcoa'
          'guln', 'guln__L'
          'Lcystin', 'cysi__L'
          'phyQ', 'phllqne'
          'tagat_D', 'tag__D'
          'yvite', 'gtocophe'
          'eumelanin', ''};  % NOTE: "eumelanin" does not exist in BiGG!
[hasMatch, matchInd] = ismember(metAssoc.metBiGGID, IDconv(:,1));
metAssoc.metBiGGID(hasMatch) = IDconv(matchInd(hasMatch),2);

% verify that all BiGG IDs are now found in the BiGG database
if ~all(ismember(metAssoc.metBiGGID, bigg.universal_bigg_id) | cellfun(@isempty, metAssoc.metBiGGID))
    fprintf('FAIL: Some metBiGGIDs in metAssoc were NOT found in the BiGG Database!\n');
else
    fprintf('SUCCESS: All non-empty metBiGGIDs in metAssoc are found in the BiGG Database!\n\n');
    
    % export updated metAssoc structure to JSON
    jsonStr = jsonencode(metAssoc);
    fid = fopen('../../ComplementaryData/annotation/humanGEMMetAssoc.JSON', 'w');
    fwrite(fid, prettyJson(jsonStr));
    fclose(fid);
    
    fprintf('New metAssoc structure written to humanGEMMetAssoc.JSON.\n');
end














