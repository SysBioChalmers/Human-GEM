function [model_new,results] = mapMetsViaRxns(model,mnx)
%mapMetsViaRxns  Add/remove met MNXID associations based on involved rxns.
%
% mapMetsViaRxns analyzes model metabolite and reaction MNXID associations
% to determine new metMNXIDs that should be added to the model based on the
% reactions in which the corresponding metabolites are involved.
%
% In addition, the function will then determine which IDs should be removed
% in situations where a met is associated with multiple MNXIDs. This is 
% accomplished by identifying all model reactions that involve the met, 
% obtaining their associated MNX reaction IDs, and retrieving those rxns 
% from the MNX database. Any MNXIDs associated with the met that are not
% included in that set of rxns from the MNX database will be removed.
% 
%
% USAGE:
%
%   [model_new,results] = mapMetsViaRxns(model,mnx);
%
% INPUT:
%
%   model    A model structure containing reaction and metabolite MNX ID
%            association fields ("rxnMNXID" and "metMNXID", respectively).
%
%   mnx      (Optional) An MNX database structure, containing reaction-
%            related information retrieved from the MNX database, generated
%            using the following command: mnx = buildMNXmodel('rxn');
%            By default, the function will automatically run the above
%            command to regenerate the MNX database structure (slower).
%
% OUTPUT:
%
%   model_new  A model which has added met MNXIDs based on their associated
%              rxns, and/or removed met MNXIDs that were not found in any 
%              of the MNX rxns associated with that metabolite.
%
%   results    A results structure containing more detailed information on
%              the analysis and addition/removal results.
%
%
% Jonathan Robinson 2018-05-25


% handle input arguments
if nargin < 2
    mnx = buildMNXmodel('rxn');
end

% if metMNXID field contains multiple columns, consolidate these into a
% single column with nested cell entries
if size(model.metMNXID,2) > 1
    model.metMNXID = nestCell(model.metMNXID,true);
end


% Identify all metabolites with two or more MNXID associations. Metabolites
% with one or zero MNXIDs will be ignored.
multi_ind = find(cellfun(@numel,model.metMNXID) > 1);

% Iterate through mets with multiple MNXIDs, and determine which MNXIDs
% should be removed.
for i = 1:length(multi_ind)
    rxn_ind = model.S(multi_ind(i),:) ~= 0;
    
    
    
    
end









