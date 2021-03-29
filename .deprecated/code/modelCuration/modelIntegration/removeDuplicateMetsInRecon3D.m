%
% FILE NAME:    removeDuplicateMetsInRecon3D.m
% 
% PURPOSE: Remove duplicate mets in Recon3D model after incorporating
%          metaboite association with HMR2
%


load('Recon3DRaven.mat');
model=Recon3DRaven;
metFreq=countFrequency(model.mets);        % count occurrence
multi_ind=find([metFreq.frequency{:}]>1);  % index to duplicate mets

idxDelete=[];  % store redundant met index

% Go through occurrence count results
for i=1:length(multi_ind)
    m=multi_ind(i);  % Get mets with multiple occurrences
    
    % Find the index for each these mets in the S matrix
    index=find(strcmp(model.mets,metFreq.uniqueList{m}));
    
    % Merging rows of duplicate mets into one
    newrow=model.S(index(1),:);   % the merged row
    for j=2:length(index)
        oldind=find(newrow);                % rxn index
        newind=find(model.S(index(j),:));   % rxn index

        % Make sure the duplicate mets do not happen in the same rxns
        if isempty(intersect(oldind,newind))
            newrow(newind)=model.S(index(j),newind);
            idxDelete=[idxDelete; index(j)];  % Save met index to delete 
        else
            fprintf('Conflicts detected when merging mets\n');
            return;   % Exit the function
            % all duplicate mets do not happen in the same rxns
        end
    end
    
    % Use the merged row
    model.S(index(1),:)=newrow;
end

% For these duplicate mets,the relevant fields (metNames, metComps, metMiriams)
% were confirmed with identical values; Others (inchis, metFormulas, metCharges) 
% are not and should be manually checked later

% Clear other redundant mets
if ~isempty(idxDelete)
    model.S(idxDelete,:) =[];
    model.mets(idxDelete) = [];
    model.metNames(idxDelete) = [];
    model.metComps(idxDelete) = [];
    model.b(idxDelete) = [];
    if isfield(model,'metFormulas')
        model.metFormulas(idxDelete) = [];
    end
    if isfield(model,'unconstrained')
        model.unconstrained(idxDelete) = [];
    end
    if isfield(model,'metMiriams')
        model.metMiriams(idxDelete) = [];
    end
    if isfield(model,'metCharges')
        model.metCharges(idxDelete) = [];
    end
    if isfield(model,'inchis')
        model.inchis(idxDelete) = [];
    end
    if isfield(model,'metFrom')
        model.metFrom(idxDelete) = [];
    end
end

% Backup S matrix and mets
model.oldS=Recon3DRaven.S;
model.oldMets=Recon3DRaven.mets;
Recon3DRaven=model;

% To assist model merge, modify the metNames for metabolites
% 'phacgly_s' and 'phacgly_c' by changing from 'Phenylacetylglycine'
% to 'Phenylacetylglycine_phacgly'
[~, index]=ismember({'phacgly_c','phacgly_s'},Recon3DRaven.mets);
Recon3DRaven.metNames(index)={'Phenylacetylglycine_phacgly'};

save('Recon3DRaven.mat','Recon3DRaven');
