%
%   FILE NAME:    addManuallyCuratedMetAssoc.m
% 
%   PURPOSE: Incorporate manual curation results of met association 
%            with Recon/BiGG/MetaNetX ids to ihumanMets2MNX_v2.mat
%

% Load met association to MNX/BiGG
load('ihumanMets2MNX_v2.mat');

% Get the list of empty elements of Recon3D/BiGG
emptyList=find(cellfun(@isempty, ihuman.metRecon3DID));
MNXID=reformatElements(ihuman.metMNXID,'cell2str');
for i=1:length(emptyList)
		m=emptyList(i);
		fprintf('%s\t%s\t%s\n',ihuman.mets{m},ihuman.metNames{m},MNXID{m});
end

% Locate the ones with multiple Recon3D associations
HMRmetsNoComp=regexprep(ihuman.mets,'\w$','');
multi_index=find(cellfun(@numel, ihuman.metRecon3DID)>1);
% Obtain the non-unique assoc
metString=reformatElements(ihuman.metRecon3DID,'cell2str');
[multiMetAssoc.Recon3DID, I, ~]=unique(metString(multi_index));
multiMetAssoc.HMRID=HMRmetsNoComp(multi_index(I));
multiMetAssoc.MNXID=reformatElements(ihuman.metMNXID(multi_index(I)),'cell2str');
for i=1:length(multiMetAssoc.HMRID)
		fprintf('%s\t%s\t%s\n',multiMetAssoc.HMRID{i},multiMetAssoc.MNXID{i},multiMetAssoc.Recon3DID{i});
end

% These ones with missing and multiple Recon3D/BiGG ids were organized into
% an excel file metaboliteCuration_20180605.xlsx for manual cuartion

% The following updates are based on the manual curation results
% 1. For the elements without Recon3D/BiGG assoc
ihuman.metRecon3DID{find(strcmp('m00077p',ihuman.mets))}={'CE2416'};
ihuman.metRecon3DID{find(strcmp('m00678m',ihuman.mets))}={'dec24dicoa'};
ihuman.metRecon3DID{find(strcmp('m00678p',ihuman.mets))}={'dec24dicoa'};
ihuman.metRecon3DID{find(strcmp('m03035m',ihuman.mets))}={'dece3coa'};
ihuman.metRecon3DID{find(strcmp('m03035p',ihuman.mets))}={'dece3coa'};
ihuman.metRecon3DID{find(strcmp('m00980m',ihuman.mets))}={'dece4coa'};
ihuman.metRecon3DID{find(strcmp('m00980p',ihuman.mets))}={'dece4coa'};
ihuman.metRecon3DID{find(strcmp('m01422c',ihuman.mets))}={'cbtnCCP'};
ihuman.metRecon3DID{find(strcmp('m01942g',ihuman.mets))}={'gd1a_hs'};

ihuman.metMNXID{find(strcmp('m01942g',ihuman.mets))}={'MNXM11644','MNXM8637'};
ihuman.metMNXID{find(strcmp('m02149c',ihuman.mets))}={'MNXM56888'};
ihuman.metMNXID{find(strcmp('m02147m',ihuman.mets))}={'MNXM527231'};
ihuman.metMNXID{find(strcmp('m02147x',ihuman.mets))}={'MNXM527231'};
ihuman.metMNXID{find(strcmp('m02147s',ihuman.mets))}={'MNXM527231'};
ihuman.metMNXID{find(strcmp('m02147c',ihuman.mets))}={'MNXM527231'};


% 2. For the elements with multiple Recon3D/BiGG association
ihuman.metRecon3DID{find(strcmp('m00490c',ihuman.mets))}={'mag_hs'};
ihuman.metRecon3DID{find(strcmp('m02733c',ihuman.mets))}={'pa_hs'};
ihuman.metRecon3DID{find(strcmp('m02733r',ihuman.mets))}={'pa_hs'};
ihuman.metRecon3DID{find(strcmp('m00240c',ihuman.mets))}={'dag_hs'};
ihuman.metRecon3DID{find(strcmp('m00240g',ihuman.mets))}={'dag_hs'};
ihuman.metRecon3DID{find(strcmp('m00240n',ihuman.mets))}={'dag_hs'};
ihuman.metRecon3DID{find(strcmp('m01818c',ihuman.mets))}={'11_cis_retfa'};
ihuman.metRecon3DID{find(strcmp('m01818s',ihuman.mets))}={'11_cis_retfa'};
ihuman.metRecon3DID{find(strcmp('m01818x',ihuman.mets))}={'11_cis_retfa'};
ihuman.metRecon3DID{find(strcmp('m02959s',ihuman.mets))}={'tag_hs'};

save('ihumanMets2MNX_v2.mat','ihuman');  % 2018-06-12
% Some MNX IDs also need refinement (to be continued)

% 3. Generate the array structure (metAssocHMR2Recon3D) of metabolite
% association for further curation and model integration, by trimming
% off duplicate mets in multiple compartments
HMRmets=regexprep(ihuman.mets,'\w$','');   % HMR met ids without comp id
[metAssocHMR2Recon3.metHMRID, I, ~]=unique(HMRmets,'stable');
metAssocHMR2Recon3.metBiGGID=ihuman.metBiGGID(I);
metAssocHMR2Recon3.metMNXID=ihuman.metMNXID(I);
metAssocHMR2Recon3.metRecon3DID=ihuman.metRecon3DID(I);
metAssocHMR2Recon3.metNames=ihuman.metNames(I);
save('metAssocHMR2Recon3.mat','metAssocHMR2Recon3');  % 2018-06-17


% 4. Detect mets associated from one Recon3D id to multiple HMR ids
% get the array of non-empty Recon3D met assocations
Recon3DID=reformatElements(metAssocHMR2Recon3.metRecon3DID,'cell2str');
Recon3DID_nonEmpty=Recon3DID(getNonEmptyList(Recon3DID));

% get the array of unique Recon3D met id and the occurrences
check=countFrequency(Recon3DID_nonEmpty);
list=find([check.frequency{:}] > 1);   % Recon3D mets with multiple occurrences
for i=1:numel(list)
		m=list(i);
		ind=find(strcmp(check.uniqueList{m},Recon3DID));
		fprintf('%s\t%s\n',check.uniqueList{m},strjoin(metAssocHMR2Recon3.metHMRID(ind),';'));
end
% These Recon3D mets associated to multiple HMR ids were subjected to
% manual cuartion 2018-06-18


% 5. The following updates of HMR-Recon3D assoc are based on the curation results
% in excel file metaboliteCuration_20180618_HW.xlsx
metAssocHMR2Recon3.metRecon3DID{find(strcmp('m00591',metAssocHMR2Recon3.metHMRID))}='';
metAssocHMR2Recon3.metRecon3DID{find(strcmp('m00352',metAssocHMR2Recon3.metHMRID))}='';
metAssocHMR2Recon3.metRecon3DID{find(strcmp('m00379',metAssocHMR2Recon3.metHMRID))}='';
metAssocHMR2Recon3.metRecon3DID{find(strcmp('m00555',metAssocHMR2Recon3.metHMRID))}='';
metAssocHMR2Recon3.metRecon3DID{find(strcmp('m01123',metAssocHMR2Recon3.metHMRID))}={'CE7097'};
metAssocHMR2Recon3.metRecon3DID{find(strcmp('m01911',metAssocHMR2Recon3.metHMRID))}='';
metAssocHMR2Recon3.metRecon3DID{find(strcmp('m01942',metAssocHMR2Recon3.metHMRID))}='';
metAssocHMR2Recon3.metRecon3DID{find(strcmp('m02410',metAssocHMR2Recon3.metHMRID))}='';
metAssocHMR2Recon3.metRecon3DID{find(strcmp('m02487',metAssocHMR2Recon3.metHMRID))}='';
metAssocHMR2Recon3.metRecon3DID{find(strcmp('m02578',metAssocHMR2Recon3.metHMRID))}='';


% 6. Associate Recon3D mets 'M00196' and 'protein' to HMR2 met 'm00196'
ind=find(strcmp(metAssocHMR2Recon3.metHMRID,'m00196'));
metAssocHMR2Recon3.metRecon3DID{ind}{2}='protein';
save('metAssocHMR2Recon3.mat','metAssocHMR2Recon3');  % 2018-08-03


% 7. Synchronize curtated metabolite association previously in % 5 and % 6
% from 'metAssocHMR2Recon3.mat' to 'ihumanMets2MNX_v2.mat'
load('ihumanMets2MNX_v2.mat');
ihuman.metRecon3DID{find(strcmp('m00591c',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m00591s',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m00352r',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m00379c',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m00555c',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m00555g',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m00555r',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m01123c',ihuman.mets))}{1}='CE7097';
ihuman.metRecon3DID{find(strcmp('m01911c',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m01942g',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m02410c',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m02410m',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m02410r',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m02487c',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m02487m',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m02578c',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m02578m',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m02578n',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m02578p',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m02578r',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m02578s',ihuman.mets))}='';
ihuman.metRecon3DID{find(strcmp('m02578x',ihuman.mets))}='';
% Associate Recon3D mets 'M00196' and 'protein' to HMR2 met 'm00196'
ihuman.metRecon3DID{find(strcmp('m00196c',ihuman.mets))}{2}='protein';
save('ihumanMets2MNX_v2.mat','ihuman');  % 2018-09-03
