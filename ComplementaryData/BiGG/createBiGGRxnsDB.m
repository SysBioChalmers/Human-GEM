%
%   FILE NAME: createBiGGRxnsDB.m
% 
%   DATE CREATED: 2018-01-17
%        
%	
%   PROGRAMMER:   Hao Wang
%                 Department of Biology and Biological Engineering
%                 Chalmers University of Technology
% 
% 
%   PURPOSE: Generate BiGG reaction database to Matlab structure
% 

% Move to the target folder
cd('/Users/haowa/Box Sync/HMR3/BiGG');

% Load the text format BiGG reactions
T=readtable('bigg_models_reactions_20180117.txt','ReadVariableNames',1);
BiGGRxns=table2struct(T,'ToScalar',true);

% Add oldids field as cell array based on string elements from old_bigg_ids
num=numel(BiGGRxns.bigg_id);
BiGGRxns.oldids=cell(num,1);
BiGGRxns.oldids(:,1)={''};
for i=1:num
	% Ignore old_bigg_ids identical to bigg_id
	if ~isequal(BiGGRxns.bigg_id{i},BiGGRxns.old_bigg_ids{i})
		BiGGRxns.oldids{i}=transpose(strsplit(BiGGRxns.old_bigg_ids{i},'; '));
	end
end

% Extract associated MNXref reaction ids
MNXrefid=regexp(BiGGRxns.database_links,'MNXR\d+','match');
% Add MNX reaction field to BiGGRxns structure
BiGGRxns.MNXrefid=cell(num,1);
BiGGRxns.MNXrefid(:,1)={''};
for i=1:num
	if ~isempty(MNXrefid{i})
		BiGGRxns.MNXrefid{i}=MNXrefid{i}{1};
	end
end

numel(find(~cellfun(@isempty,BiGGRxns.MNXrefid)))  %ans = 15094
save('BiGGRxns.mat','BiGGRxns');

