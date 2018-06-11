%
%   FILE NAME:    addManuallyCuratedMetAssoc.m
% 
%   DATE CREATED: 2018-06-05
%       MODIFIED: 2018-06-11
%        
%   PROGRAMMER:   Hao Wang
%                 Department of Biology and Biological Engineering
%                 Chalmers University of Technology
% 
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
		m=list(i);
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
% an excel file metaboliteCuration_20180605.xlsx for manula cuartion

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
ihuman.metMNXID{find(strcmp('m02149c',ihuman.mets))}={'MNXM527231'};

