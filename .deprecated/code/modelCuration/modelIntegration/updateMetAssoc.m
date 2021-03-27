%
% FILE NAME:    updateMetAssoc.m
% 
% PURPOSE: Create data structure of one-to-one met assocation between
%          HMR2 and Recon3D 
% 
% Note: The generated metAssoc.mat file is for convenient model integration,
% and this data structure should be updated everytime once metAssocHMR2Recon3.mat
% is changed!


% Get the manually curated met association info
load('metAssocHMR2Recon3.mat');

% Directly save the unique associations
single_ind=find(cellfun(@numel,metAssocHMR2Recon3.metR3DID)==1);
metAssoc.metHMRID=metAssocHMR2Recon3.metHMRID(single_ind);
metAssoc.metRecon3DID=reformatElements(metAssocHMR2Recon3.metR3DID(single_ind),'cell2str');
metAssoc.metNames=metAssocHMR2Recon3.metNames(single_ind);


% Associate between one HMR id to multiple Recon3D ids
multi_ind=find(cellfun(@numel,metAssocHMR2Recon3.metR3DID)>1);
for i=1:length(multi_ind)
		m=multi_ind(i);
		num=numel(metAssocHMR2Recon3.metR3DID{m});
		temp=cell(num,1);
		temp(:)={metAssocHMR2Recon3.metHMRID{m}};
		names=cell(num,1);
		names(:)={metAssocHMR2Recon3.metNames{m}};
		metAssoc.metHMRID=[metAssoc.metHMRID;temp];
		metAssoc.metRecon3DID=[metAssoc.metRecon3DID;transpose(metAssocHMR2Recon3.metR3DID{m})];
		metAssoc.metNames=[metAssoc.metNames;names];
end

% Save to modelIntegration subfolder
save('metAssoc.mat','metAssoc');         % 2018-09-20
