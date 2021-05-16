%
%   FILE NAME: getRxnsFromBiGG.m
% 
%   PURPOSE: Generate Matlab structure for BiGG reactions
%


% Move to the target folder
cd('/Users/haowa/Box Sync/HMR3/BiGG');

% Load the text format BiGG reactions
T=readtable('bigg_models_reactions_20180424.txt','ReadVariableNames',1);
BiGGRxns=table2struct(T,'ToScalar',true);

% Rename the fields according to RAVEN specification
BiGGRxns.rxns=BiGGRxns.bigg_id;
BiGGRxns.rxnNames=BiGGRxns.name;
BiGGRxns.rxnEquations=BiGGRxns.reaction_string;

% Add oldids field as cell array based on info from old_bigg_ids
num=numel(BiGGRxns.bigg_id);
BiGGRxns.oldids=cell(num,1);
BiGGRxns.oldids(:)={''};
% Add MNX field
rxnMNXID=regexp(BiGGRxns.database_links,'MNXR\d+','match');
BiGGRxns.rxnMNXID=cell(num,1);
BiGGRxns.rxnMNXID(:)={''};
count=0;
for i=1:num
		% Ignore old_bigg_ids identical to bigg_id
		if ~isequal(BiGGRxns.rxns{i},BiGGRxns.old_bigg_ids{i})
				% Remove existing BiGG ids from oldids
				BiGGRxns.oldids{i}=setdiff(transpose(strsplit(BiGGRxns.old_bigg_ids{i},'; ')),BiGGRxns.rxns{i});
		end
		if ~isempty(rxnMNXID{i})
				if numel(rxnMNXID{i})==1
						BiGGRxns.rxnMNXID{i}=rxnMNXID{i}{1};
				else
						BiGGRxns.rxnMNXID{i}=rxnMNXID{i};
						count=count+1;
				end
		end
end
%count=0
numel(find(~cellfun(@isempty,BiGGRxns.rxnMNXID)))  %ans = 15904

% Remove some fields
BiGGRxns=rmfield(BiGGRxns,{'bigg_id','name','reaction_string','old_bigg_ids'});

save('BiGGRxns.mat','BiGGRxns');

