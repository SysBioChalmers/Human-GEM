%
% FILE NAME:    updateBiGGIDs_issue124.m
% 
% PURPOSE: [Addresses Issue #124]
%
%          Many of the BiGG IDs in the humanGEMMetAssoc.JSON and 
%          humanGEMRxnAssoc.JSON annotation files are invalid.
%          For example, metabolite 'ala_L' was retrieved from Recon3D, but 
%          the valid BiGG ID is 'ala__L'.
%
%          This script updates the metBiGGIDs in the humanGEMMetAssoc.JSON
%          annotation file, and the rxnBiGGIDs in the humanGEMRxnAssoc.JSON
%          file, and confirms that they are all valid by comparing to the
%          corresponding data files retrieved from the BiGG database:
%
%          http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt
%          http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt
%
%          The script writes the updated IDs to the humanGEMMetAssoc.JSON
%          and humanGEMRxnAssoc.JSON files.
%
%
%          *** Additionally, this script updates ChEBI IDs in the
%          humanGEMMetAssoc.JSON file to include "CHEBI" in the ID, which
%          is consistent with the identifiers.org nomenclature.
%          For example, a ChEBI ID of "1234" should instead be written as
%          "CHEBI:1234".
%


%% Update metabolite BiGG IDs

% load met association file
metAssoc = jsondecode(fileread('../../ComplementaryData/annotation/humanGEMMetAssoc.JSON'));

% retrieve metabolite data from BiGG database and convert to structure
webdata = webread('http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt');
webdata = textscan(webdata, repmat('%s',1,6), 'Delimiter', '\t');  % file contains 6 columns of text
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
    fprintf('SUCCESS: All non-empty metBiGGIDs in metAssoc are found in the BiGG Database!\n');
    
    % correct met ChEBI ID format
    metAssoc.metChEBIID = regexprep(metAssoc.metChEBIID, 'chebi', 'CHEBI', 'ignorecase');  % make it all uppercase
    metAssoc.metChEBIID = regexprep(metAssoc.metChEBIID,'(^|\s)(\d+)','$1CHEBI:$2');  % add "CHEBI:" before numbers that don't already have it
    
    % export updated metAssoc structure to JSON
    jsonStr = jsonencode(metAssoc);
    fid = fopen('../../ComplementaryData/annotation/humanGEMMetAssoc.JSON', 'w');
    fwrite(fid, prettyJson(jsonStr));
    fclose(fid);
    
    fprintf('New metAssoc structure was written to humanGEMMetAssoc.JSON.\n');
end



%% Update reaction BiGG IDs

% load rxn association file
rxnAssoc = jsondecode(fileread('../../ComplementaryData/annotation/humanGEMRxnAssoc.JSON'));

% retrieve reaction data from BiGG database and convert to structure
webdata = webread('http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt');
webdata = textscan(webdata, repmat('%s',1,6), 'Delimiter', '\t');  % file contains 6 columns of text
bigg = {};
for i = 1:numel(webdata)
    bigg.(webdata{i}{1}) = webdata{i}(2:end);
end

% determine which BiGG IDs are not found in the database
badInd = find(~ismember(rxnAssoc.rxnBiGGID, bigg.bigg_id) & ~cellfun(@isempty, rxnAssoc.rxnBiGGID));

% Now check which of these are "old" IDs. First we need to reformat the ID
% mapping by splitting "old" ID lists by the "; " delimiter.
oldIDs_reformat = cellfun(@(id) strsplit(id, '; ')', bigg.old_bigg_ids, 'UniformOutput', false);
newIDs_reformat = arrayfun(@(i) repmat(bigg.bigg_id(i),numel(oldIDs_reformat{i}),1), (1:numel(bigg.bigg_id))', 'UniformOutput', false);
biggIDs_reformat = [vertcat(oldIDs_reformat{:}), vertcat(newIDs_reformat{:})];

% update "old" BiGG IDs to current BiGG IDs
[hasMatch, matchInd] = ismember(rxnAssoc.rxnBiGGID(badInd), biggIDs_reformat(:,1));
rxnAssoc.rxnBiGGID(badInd(hasMatch)) = biggIDs_reformat(matchInd(hasMatch),2);


% determine which BiGG IDs are still not found in the database
badInd = find(~ismember(rxnAssoc.rxnBiGGID, bigg.bigg_id) & ~cellfun(@isempty, rxnAssoc.rxnBiGGID));

% some reaction IDs can be fixed by changing a dash '-' to underscore '_'
modifiedID = regexprep(rxnAssoc.rxnBiGGID(badInd), '-', '_');
newMatch = ismember(modifiedID, bigg.bigg_id);
rxnAssoc.rxnBiGGID(badInd(newMatch)) = modifiedID(newMatch);

% some reaction IDs can be fixed by replacing compartment parentheses with
% underscores
modifiedID = regexprep(rxnAssoc.rxnBiGGID(badInd), '(\(|\[)(e|s|bl)(\)|\])$', '_e');
newMatch = ismember(modifiedID, bigg.bigg_id);
rxnAssoc.rxnBiGGID(badInd(newMatch)) = modifiedID(newMatch);

% some IDs need both modifications
modifiedID = regexprep(rxnAssoc.rxnBiGGID(badInd), '-', '_');
modifiedID = regexprep(modifiedID, '(\(|\[)(e|s|bl)(\)|\])$', '_e');
newMatch = ismember(modifiedID, bigg.bigg_id);
rxnAssoc.rxnBiGGID(badInd(newMatch)) = modifiedID(newMatch);

modifiedID = regexprep(rxnAssoc.rxnBiGGID(badInd), '-', '__');  % double underscore
modifiedID = regexprep(modifiedID, '(\(|\[)(e|s|bl)(\)|\])$', '_e');
newMatch = ismember(modifiedID, bigg.bigg_id);
rxnAssoc.rxnBiGGID(badInd(newMatch)) = modifiedID(newMatch);
 

% determine which BiGG IDs are still not found in the database
badInd = find(~ismember(rxnAssoc.rxnBiGGID, bigg.bigg_id) & ~cellfun(@isempty, rxnAssoc.rxnBiGGID));

% a lot of reactions can be matched using their rxnRecon3DID
for i = 1:numel(badInd)
    id = strsplit(rxnAssoc.rxnRecon3DID{badInd(i)},'; ');
    for j = 1:numel(id)
        if ismember(id{j}, bigg.bigg_id)
            rxnAssoc.rxnBiGGID{badInd(i)} = id{j};
            break
        else
            % check for matches to "old" BiGG IDs
            [~,ind] = ismember(id{j}, biggIDs_reformat(:,1));
            if ind > 0
                rxnAssoc.rxnBiGGID{badInd(i)} = biggIDs_reformat{ind,2};
                break
            end
        end
    end
end

% the remaining unmatched IDs are manually updated as below:
IDconv = {'CBPSAm', 'CPS'
          'MMCOAHm', 'r0571'
          'STRRer', 'CHLSTR'
          'STRR2er', 'ZYMSTR'
          'DM_Asn-X-Ser/Thr(ly)', ''  % rxn not in BiGG DB
          'sink_pre_prot(er)', ''};   % rxn not in BiGG DB
[hasMatch, matchInd] = ismember(rxnAssoc.rxnBiGGID, IDconv(:,1));
rxnAssoc.rxnBiGGID(hasMatch) = IDconv(matchInd(hasMatch),2);

% verify that all rxn BiGG IDs are now found in the BiGG database
if ~all(ismember(rxnAssoc.rxnBiGGID, bigg.bigg_id) | cellfun(@isempty, rxnAssoc.rxnBiGGID))
    fprintf('FAIL: Some rxnBiGGIDs in rxnAssoc were NOT found in the BiGG Database!\n');
else
    fprintf('SUCCESS: All non-empty rxnBiGGIDs in rxnAssoc are found in the BiGG Database!\n');
    
    % export updated rxnAssoc structure to JSON
    jsonStr = jsonencode(rxnAssoc);
    fid = fopen('../../ComplementaryData/annotation/humanGEMRxnAssoc.JSON', 'w');
    fwrite(fid, prettyJson(jsonStr));
    fclose(fid);
    
    fprintf('New rxnAssoc structure was written to humanGEMRxnAssoc.JSON.\n');
end









