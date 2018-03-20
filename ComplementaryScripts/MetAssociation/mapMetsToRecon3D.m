function model_mapped = mapMetsToRecon3D(model,Recon3D)
%mapMetsToRecon3D  Map metabolites to Recon3D via metabolite names.
%
% USAGE:
%
%   model_mapped = MapMetsToRecon3D(model,Recon3D);
%
%
% INPUTS:
%
%   model     Model structure containing metabolites that are to be mapped
%             to Recon3D metabolite identifiers.
%
%   Recon3D   Recon3D model structure.
%
%
% OUTPUTS:
%
%   model_mapped  Model structure with added field "metRecon3DID", which
%                 contains the Recon3D metIDs that were matched to the
%                 model metabolites.
%
%
% Jonathan Robinson 2018-03-20

% initialize output
model_mapped = model;

% strip compartments from Recon3D mets to obtain compartment-free
% met2metName pairs
Recon3D.mets = regexprep(Recon3D.mets,'\[.\]$','');
[~,ind] = unique(Recon3D.mets);
reconMet2Name = [Recon3D.mets(ind),Recon3D.metNames(ind)];

% split Recon3D metNames containing semicolon, and add as additional pairs
reconMet2Name(:,2) = cellfun(@(s) strsplit(s,'; '),reconMet2Name(:,2),'UniformOutput',false);
reconMet2Name = arrayfun(@(i) [repmat(reconMet2Name(i,1),numel(reconMet2Name{i,2}),1),vertcat(reconMet2Name{i,2})'],[1:length(reconMet2Name)]','UniformOutput',false);
reconMet2Name = vertcat(reconMet2Name{:});

% map metabolites
metRecon3DID = matchIDs(model.metNames,reconMet2Name);

% find unmapped metabolites
unmapped = cellfun(@isempty,metRecon3DID);


% use the metNamesAlt field to match unmapped metabolites
model.metNames = [model.metNames,model.metNamesAlt];
metRecon3DID(unmapped) = matchIDs(model.metNames(unmapped,:),reconMet2Name);

% find unmapped metabolites
unmapped = cellfun(@isempty,metRecon3DID);


% convert 'Coenzyme A' into 'CoA' to try and match remaining metabolites
reconMet2Name(:,2) = regexprep(reconMet2Name(:,2),'Coenzyme A','CoA');
metRecon3DID(unmapped) = matchIDs(model.metNames(unmapped,:),reconMet2Name);

% find unmapped metabolites
unmapped = cellfun(@isempty,metRecon3DID);


%-------------- specific targeted changes based on analysis ---------------

reconMet2Name(:,2) = regexprep(reconMet2Name(:,2),'\(.*-Density Lipoprotein\)','');
reconMet2Name(:,2) = regexprep(reconMet2Name(:,2),'13-Cis-Retinoyl Glucuronide','13-Cis-Retinoyl-beta-D-Glucuronide');
reconMet2Name(:,2) = regexprep(reconMet2Name(:,2),"Cytidine-5'-Diphosphate",'CDP');
reconMet2Name(:,2) = regexprep(reconMet2Name(:,2),'Human Liver Homolog','');
reconMet2Name(:,2) = regexprep(reconMet2Name(:,2),'Glutathionyl-Leuc4','glutathionyl-leukotriene C4');
reconMet2Name(:,2) = regexprep(reconMet2Name(:,2),'Omega-Cooh-Tetranor-Leukotriene E3','Omega-Cooh-Tetranor-LTE3');

%--------------------------------------------------------------------------

metRecon3DID(unmapped) = matchIDs(model.metNames(unmapped,:),reconMet2Name);


% flatten cell column into 2D cell array
metRecon3DID = flattenCell(metRecon3DID,true);

% assign output
model_mapped.metRecon3DID = metRecon3DID;


end


function metRecon3DID = matchIDs(metNames,reconMet2Name)

% make metNames in both models lowercase, and remove special characters
% (non-word or digit characters, e.g., dashes, parentheses, spaces, etc.)
reconMet2Name(:,2) = lower(regexprep(reconMet2Name(:,2),'[^a-zA-Z0-9]',''));
metNames = lower(regexprep(metNames,'[^a-zA-Z0-9]',''));

% compress model metNames into single column of nested entries
if size(metNames,2) > 1
    metNames = nestCell(metNames,true);
end

% find empty name indices to ignore
no_names = cellfun(@isempty,metNames);

% map model metNames to Recon3D metNames
metRecon3DID = repmat({''},size(metNames,1),1);
metRecon3DID(~no_names) = cellfun(@(x) reconMet2Name(ismember(reconMet2Name(:,2),x),1),metNames(~no_names),'UniformOutput',false);

end

