function result=metConsistencyCheckViaMNX(queryList)
% 
%   metConsistencyCheckViaMNX aims to check the consistency when multiple
%   MNX metids are associated to one HMR metabolite, and report to the
%   result cell array with four conditions:
%          'Empty'  -  no assocaiton to this metabolite
%          'Single' -  with single MNX metid association
%          'Pass'   -  the multiple MNX metids share the same formula and charge
%          'Fail'   -  the multiple MNX metids have different formula and charge
%
%   queryList    cell array of associated MNX metabolite ids (nested array)
%
%   Usage: result=metConsistencyCheckViaMNX(queryList)
%


if nargin<1
		EM='Missing input arguments';
		disp(EM);
end

% Initilize output cell array
result=cell(numel(queryList),1);
result(:)={''};

% If there is No or Single MNX metID associated
empty_ind=find(cellfun(@isempty,queryList));
single_ind=find(cellfun(@numel,queryList)==1);
result(empty_ind)={'Empty'};
result(single_ind)={'Single'};

% Deal with mets with multiple MNX association
load('MNXMets.mat');
multi_ind=find(cellfun(@numel,queryList) > 1);
for i = 1:length(multi_ind)
		m=multi_ind(i);
		[hit, index]=ismember(queryList{m},MNXMets.mets);
		if all(hit)
				charges=MNXMets.metCharges(index);
				if isequal(MNXMets.metFormulas{index}) && all(charges==charges(1))
						result{m}='Pass';
				else
						result{m}='Fail';
				end
		else
				dispEM('There is mistakes in metabolite association to MNX id!');
		end
end

