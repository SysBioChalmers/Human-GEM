%
% FILE NAME:    integrateRecon3DToHMR.m
% 
% PURPOSE: The master scrit for integrate Recon3D into HMR2 based on
%          reaction/metabolite association to MetaNetX/BiGG databases
%


% 1. Generate the comprehensive reaction associaiton between HMR2 and Recon3D

% a. Load the manually curated met association info
load('Recon3Rxns2HMR.mat');

% Directly save the one-to-one associations
tmp=reformatElements(Recon3D.rxnHMRID,'str2cell');
single_ind=find(cellfun(@numel,tmp)==1);
rxnAssoc.rxnHMRID=Recon3D.rxnHMRID(single_ind);
rxnAssoc.rxnRecon3DID=Recon3D.rxns(single_ind);

% Associate between one Recon3D id to multiple HMR ids
multi_ind=find(cellfun(@numel,tmp)>1);
for i=1:length(multi_ind)
		m=multi_ind(i);
		num=numel(tmp{m});
		temp=cell(num,1);
		temp(:)={Recon3D.rxns{m}};
		rxnAssoc.rxnRecon3DID=[rxnAssoc.rxnRecon3DID;temp];
		rxnAssoc.rxnHMRID=[rxnAssoc.rxnHMRID;transpose(tmp{m})];
end

% b. Get reaction pairs from overlapRxnDetection function
load('Recon3DRaven.mat');
load('humanGEM.mat');
overlapRxns=overlapRxnDetection(ihuman, Recon3DRaven);

% c. Combine automatically detected and manually confirmed association together
rxnAssoc.rxnRecon3DID=[rxnAssoc.rxnRecon3DID;overlapRxns.rxnModelB];
rxnAssoc.rxnHMRID=[rxnAssoc.rxnHMRID;overlapRxns.rxnModelA];
%rxnAssoc=uniqueArray(rxnAssoc);
check=strcat(rxnAssoc.rxnRecon3DID,';',rxnAssoc.rxnHMRID);
check=unique(check);
output=split(check,';');
rxnAssoc.rxnRecon3DID=cellstr(output(:,1));
rxnAssoc.rxnHMRID=cellstr(output(:,2));
save('rxnAssoc.mat','rxnAssoc');


% 2. Get the reduced Recon3D model by removing the duplicate reactions
rxnToRemove=unique(rxnAssoc.rxnRecon3DID);
reducedRecon3D=removeReactions(Recon3DRaven,rxnToRemove,1,1,1);


% 3. Merge models
reducedRecon3D.id='reducedRecon3D';
model=mergeModels({ihuman reducedRecon3D});
model.id='humanGEM';
model.description='Integrated from HMR2 and Recon3D';


% 4. Add external reaction identifiers
model.rxnRecon3DID=cell(numel(model.rxns),1);
model.rxnRecon3DID(:)={''};
% Recon3D rxn ids
for i=1:numel(ihuman.rxns)
		%[a, b]=ismember(model.rxns{i},rxnAssoc.rxnHMRID);
		ind=find(strcmp(rxnAssoc.rxnHMRID,model.rxns{i}));
		if length(ind)>0
				model.rxnRecon3DID{i}=rxnAssoc.rxnRecon3DID(ind);
		end
end
model.rxnRecon3DID=reformatElements(model.rxnRecon3DID,'cell2str');

% 5. Remove unused fields and save model
model=rmfield(model,{'geneMiriams','annotation'});
ihuman=model;
ihuman.version='0.2.0';
save('../../model/Human-GEM.mat','ihuman');
