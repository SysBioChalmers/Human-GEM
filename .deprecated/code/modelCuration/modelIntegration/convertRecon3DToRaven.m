%
% FILE NAME:    convertRecon3DToRaven.m
% 
% PURPOSE: Convert Recon3D model to RAVEN format
%          replace met with HMR info
%


% Load Recon3D
load('/Users/haowa/Box Sync/HMR3/Recon3D/Published/ModelFiles/Recon3D_301/Recon3D_301.mat');

% Load HMR2
load('HMRdatabase2_00.mat');

% Compartment information is missing, and is added here
Recon3D.comps=[ihuman.comps;'i'];
Recon3D.compNames=[ihuman.compNames;'Inner mitochondria'];

% Modify compartment system according to HMR2: e,x -> s,p
Recon3D.originalMets=Recon3D.mets;   % Make a backup
Recon3D.mets=regexprep(Recon3D.mets,'\[e\]','[s]');
Recon3D.mets=regexprep(Recon3D.mets,'\[x\]','[p]');

% The metPubChemID field is probablmatic and remove it
Recon3D=rmfield(Recon3D,'metPubChemID');

% Convert to RAVEN format
Recon3DRaven=ravenCobraWrapper(Recon3D);

% Add equations and orginal mets
Recon3DRaven.rxnEquations=constructEquations(Recon3DRaven);
Recon3DRaven.originalMets=Recon3D.originalMets;

% Converte met info from Recon3D to HMR
load('metAssoc.mat');
Recon3DRaven.metsBeforeConv=Recon3DRaven.mets;
metIDs=regexprep(Recon3DRaven.mets,'_\w$','');
for i=1:numel(metIDs)
		ind=find(strcmp(metIDs{i},metAssoc.metRecon3DID));
		compID=Recon3DRaven.comps{Recon3DRaven.metComps(i)};
		if length(ind)==1      % This met is unique in association
				Recon3DRaven.mets{i}=strcat(metAssoc.metHMRID{ind},compID);
				Recon3DRaven.metNames{i}=metAssoc.metNames{ind};
		elseif length(ind)>1   % This met has multiple associations
				metsWithComp=strcat(metAssoc.metHMRID(ind),compID);
				temp=intersect(ihuman.mets,metsWithComp);   % See the occurrence after adding comp id
				if numel(temp)==1  % If narrow down to unique association
						Recon3DRaven.mets{i}=temp;
						Recon3DRaven.metNames{i}=ihuman.metNames{find(strcmp(temp,ihuman.mets))};
				else               % If still have multiple association, then take the first match
						Recon3DRaven.mets{i}=strcat(metAssoc.metHMRID{ind(1)},compID);;
						Recon3DRaven.metNames{i}=metAssoc.metNames{ind(1)};
				end
		end
end
Recon3DRaven.mets=cellfun(@char,Recon3DRaven.mets,'un',0);  % uniform met id format

% Convert Recon3D grRules by replacing with Ensembl ids and
% converging multiple suffix versions into one in the field
[genes,grRules,rxnGeneMat] = translateGeneRules(Recon3DRaven);

% backup genes, grRules and rxnGeneMat
Recon3DRaven.originalGrRules=Recon3DRaven.grRules;
Recon3DRaven.originalGenes=Recon3DRaven.genes;
Recon3DRaven.originalRxnGeneMat=Recon3DRaven.rxnGeneMat;

% update genes, grRules and rxnGeneMat
Recon3DRaven.grRules=grRules.ENSG;
Recon3DRaven.genes=genes.ENSG;
Recon3DRaven.rxnGeneMat=rxnGeneMat.ENSG;

save('Recon3DRaven.mat','Recon3DRaven');  % 2018-08-03

