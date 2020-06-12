function iHsa_new = addHMRrxnIDsToiHsa(iHsa)
%addHMRrxnIDsToiHsa  add HMR reaction IDs to iHsa model as a new field.
%
% USAGE:
%
%   iHsa_new = addHMRrxnIDsToiHsa(iHsa);
%
% INPUT:
%
%   iHsa       iHsa model structure.
%
% OUTPUT:
%
%   iHsa_new   iHsa model structure with additional field "rxnHMRID",
%              which contains the HMR rxn IDs that correspond to each of
%              the iHsa rxns. Note that some of the iHsa rxns do not have a
%              corresponding HMR rxn, and therefore will have a blank entry
%              in the rxnHMRID field.
%


% import rxn associations from supporting information dataset
supp_data = readtable('ComplementaryData/iHsa/iHsa_supp_data_3.xlsx','Range','A2:T8338');  % specify range to exclude first line
ihsa_id = supp_data.rxn_id;
hmr_id = supp_data.hmr2_id;
[~,ind] = ismember(iHsa.rxns,ihsa_id);
iHsa.rxnHMRID = hmr_id(ind);

% assign output
iHsa_new = iHsa;




