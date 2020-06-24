%
%   FILE NAME: getMetsFromBiGG.m
% 
%   PURPOSE: Generate data structure for BiGG metabolites
%


% Move to the target folder
cd('/Users/haowa/Box Sync/HMR3/BiGG');

% Load the text format BiGG reactions, after fixing errors with M01870
T=readtable('bigg_models_metabolites_20180424.txt','ReadVariableNames',1);
BiGGMets=table2struct(T,'ToScalar',true);

% Rename the fields according to RAVEN specification
BiGGMets.mets=BiGGMets.bigg_id;
BiGGMets.metNames=BiGGMets.name;
BiGGMets.universalID=BiGGMets.universal_bigg_id;

% Add oldids field as cell array based on info from old_bigg_ids
num=numel(BiGGMets.bigg_id);
BiGGMets.oldids=cell(num,1);
BiGGMets.oldids(:)={''};
% Add MNX field
metMNXID=regexp(BiGGMets.database_links,'MNXM\d+','match');
BiGGMets.metMNXID=cell(num,1);
BiGGMets.metMNXID(:)={''};
count=0;
for i=1:num
		if ~isempty(BiGGMets.old_bigg_ids{i})
				% Remove existing BiGG ids and universal ids from oldids
				BiGGMets.oldids{i}=setdiff(transpose(strsplit(BiGGMets.old_bigg_ids{i},'; ')),BiGGMets.mets{i});
				BiGGMets.oldids{i}=setdiff(BiGGMets.oldids{i},BiGGMets.universalID{i});
		end
		if ~isempty(metMNXID{i})
				if numel(metMNXID{i})==1
						BiGGMets.metMNXID{i}=metMNXID{i}{1};
				else
						BiGGMets.metMNXID{i}=metMNXID{i};
						count=count+1;
				end
		end
end
%count=0
numel(find(~cellfun(@isempty,BiGGMets.metMNXID)))  %ans = 9915

% Remove some fields
BiGGMets=rmfield(BiGGMets,{'bigg_id','name','universal_bigg_id','old_bigg_ids'});

save('BiGGMets.mat','BiGGMets');

