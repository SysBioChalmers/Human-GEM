function [] = compareMNXmetFormulas(model,mnx,ignoreComp)

% first subset MNX database structure for quicker analysis
allids = unique(horzcat(model.metMNXID{:})');
rem_ind = ~ismember(mnx.mets,allids);
mnx.mets(rem_ind) = [];
mnx.metFormulas(rem_ind) = [];
mnx.metCharges(rem_ind) = [];
mnx.metNames(rem_ind) = [];

if ( ignoreComp )
    model.mets = regexprep(model.mets,'.$','');
    [model.mets,uniq_ind] = unique(model.mets);
    model.metMNXID = model.metMNXID(uniq_ind);
end

num_ids = cellfun(@numel,model.metMNXID);
num_missing = zeros(size(model.mets));
for i = 1:length(model.mets)
    if isempty(model.metMNXID{i})
        continue
    end
    ind = ismember(mnx.mets,model.metMNXID{i});
    num_missing(i) = sum(cellfun(@isempty,mnx.metFormulas(ind)) | cellfun(@isempty,mnx.metCharges(ind)));
end

% num_same = zeros(size(model.mets));
% for i = 1:length(model.mets)
%     if (num_ids(i) - num_missing(i)) > 1
%         
%         
%     end
% end










