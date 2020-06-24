function GPRdata = compare_HMR_iHsa_Recon3D_GPRs(HMR,iHsa,Recon3D,rxnHMR2Recon3D,writefile)
% compare_HMR_iHsa_Recon3D_GPRs  Align HMR, iHsa, and Recon3D grRule data.
%
%   compare_HMR_iHsa_Recon3D_GPRs extracts grRule information from HMR,
%   iHsa, and Recon3D models, and aligns them in a single cell array. The
%   function will also retrieve additional information on iHsa grRules from
%   the supporting information associated with the iHsa publication.
%
%
% USAGE:
%
%   GPRdata = compare_HMR_iHsa_Recon3D_GPRs(HMR,iHsa,Recon3D,rxnHMR2Recon3D,writefile);
%
%
% INPUT:
%
%   HMR         Human Metabolic Reaction (HMR) model structure.
%
%   iHsa        iHsa (EM Blais, JA Papin, et al. 2017) model structure.
%
%   Recon3D     Recon3D model structure.
%
%   rxnHMR2Recon3D    Cell array to convert between HMR rxn IDs (first 
%                     column) and Recon3D rxn IDs (second column).
%
%   writefile   (Optional) Name of file to which results will be written.
%               If left blank, no file will be written.
%
%
% OUTPUT:
%
%   GPRdata     A cell array containing the aligned reaction and grRule
%               information corresponding to each of the three models. The
%               column headers are as follows:
%
%               'HMR rxn'           HMR rxn IDs
%               'HMR grRule'        HMR grRules
%               'HMR Ngenes'        number of genes in each HMR grRule
%               'iHsa rxn'          iHsa rxn IDs
%               'iHsa grRule'       iHsa grRules
%               'iHsa Ngenes'       number of genes in each iHsa grRule
%               'iHsa rule action'  changes made to HMR grRule to obtain iHsa grRule
%               'iHsa rule note'    
%               'iHsa rule comment' 
%               'Recon3D rxn'       
%               'Recon3D grRule'
%               'Recon3D Ngenes'
%


% handle input arguments
if nargin < 1 || isempty(HMR)
    load('ModelFiles/mat/HMRdatabase2_02.mat');  % loads as variable "ihuman"
    HMR = ihuman;
end
if nargin < 2 || isempty(iHsa)
    load('ComplementaryData/iHsa/iHsa.mat');  % loads as variable "iHsa"
end
if nargin < 3 || isempty(Recon3D)
    load('ComplementaryData/Recon3D/Recon3D_301.mat');  % loads as variable "Recon3D"
end
if nargin < 5
    writefile = [];
elseif isequal(writefile,true)
    % if a filename isn't provided, use a default filename
    writefile = 'GPRcomparison_output.txt';
end

% clean HMR.grRules if not yet done
fprintf('Cleaning HMR grRules... ');
HMR.grRules = cleanModelGeneRules(HMR.grRules);
fprintf('Done.\n');

% initialize outputs
GPRdata = [HMR.rxns,HMR.grRules,genesPerRule(HMR.grRules)];
GPRdata_head = {'HMR rxn'            % 1
                'HMR grRule'         % 2
                'HMR Ngenes'         % 3
                'iHsa rxn'           % 4
                'iHsa grRule'        % 5
                'iHsa Ngenes'        % 6
                'iHsa rule action'   % 7
                'iHsa rule note'     % 8
                'iHsa rule comment'  % 9
                'Recon3D rxn'        % 10
                'Recon3D grRule'     % 11
                'Recon3D Ngenes'}';  % 12
            
            
% check if iHsa contains HMR rxn associations
if ~isfield(iHsa,'rxnHMRID')
    iHsa = addHMRrxnIDsToiHsa(iHsa);
end

% add iHsa rxn IDs to HMR and GPRdata
HMR.rxniHsaID = repmat({''},size(HMR.rxns));
[hasmatch,ind] = ismember(HMR.rxns,iHsa.rxnHMRID);
HMR.rxniHsaID(hasmatch) = iHsa.rxns(ind(hasmatch));
GPRdata(:,ismember(GPRdata_head,'iHsa rxn')) = HMR.rxniHsaID;

% convert iHsa grRules to Ensembl IDs, and add to GPRdata
ihsa_grRule = repmat({''},size(HMR.rxns));
ihsa_grRule(hasmatch) = iHsa.grRules(ind(hasmatch));
ihsa_grRule_ensg = translateGrRules(ihsa_grRule,'ENSG');
GPRdata(:,ismember(GPRdata_head,'iHsa grRule')) = ihsa_grRule_ensg;
GPRdata(:,ismember(GPRdata_head,'iHsa Ngenes')) = genesPerRule(ihsa_grRule_ensg);

% retrieve grRule modification information from iHsa Supp Data #1
supp_data = readtable('ComplementaryData/iHsa/iHsa_supp_data_1.xlsx','Sheet','GPR Associations');
rat_ind = ismember(supp_data.variable,'gpr_rno');
supp_data(rat_ind,:) = [];  % remove rows corresponding to Rat model

% map supp data to HMR reactions
[hasmatch,ind] = ismember(HMR.rxniHsaID,supp_data.rxn_id);
hasmatch(cellfun(@isempty,HMR.rxniHsaID)) = false;

% add information to GPRdata
ihsa_rule_action = repmat({''},size(HMR.rxns));
ihsa_rule_action(hasmatch) = supp_data.action(ind(hasmatch));
GPRdata(:,ismember(GPRdata_head,'iHsa rule action')) = ihsa_rule_action;

ihsa_rule_note = repmat({''},size(HMR.rxns));
ihsa_rule_note(hasmatch) = supp_data.note(ind(hasmatch));
GPRdata(:,ismember(GPRdata_head,'iHsa rule note')) = ihsa_rule_note;

ihsa_rule_comment = repmat({''},size(HMR.rxns));
ihsa_rule_comment(hasmatch) = supp_data.comment(ind(hasmatch));
GPRdata(:,ismember(GPRdata_head,'iHsa rule comment')) = ihsa_rule_comment;


% extract information from Recon3D
[hasmatch,ind] = ismember(HMR.rxns,rxnHMR2Recon3D(:,1));
r3_rxn = repmat({''},size(HMR.rxns));
r3_rxn(hasmatch) = rxnHMR2Recon3D(ind(hasmatch),2);
GPRdata(:,ismember(GPRdata_head,'Recon3D rxn')) = r3_rxn;

[hasmatch,ind] = ismember(r3_rxn,Recon3D.rxns);
hasmatch(cellfun(@isempty,r3_rxn)) = false;
r3_rule = repmat({''},size(HMR.rxns));
r3_rule(hasmatch) = Recon3D.grRules(ind(hasmatch));
r3_rule_ensg = translateGrRules(r3_rule,'ENSG');
GPRdata(:,ismember(GPRdata_head,'Recon3D grRule')) = r3_rule_ensg;
GPRdata(:,ismember(GPRdata_head,'Recon3D Ngenes')) = genesPerRule(r3_rule_ensg);


% append GPRheader to GPRdata
GPRdata = [GPRdata_head;GPRdata];

% write results to file, if specified
if ~isempty(writefile)
    
    % convert all numbers to strings
    GPRdata_str = cellfun(@num2str,GPRdata,'UniformOutput',false);
    
    % write to file
    writecell2file(GPRdata_str,writefile,true,'\t');
end

end  % function end



% function to count the number of unique genes in each grRule
function ngenes = genesPerRule(grRules)
    
% convert AND and OR to & and |, respectively
grRules = regexprep(grRules,' or ','|');
grRules = regexprep(grRules,' and ','&');

% identify genes associated with each reaction
ngenes = cellfun(@(r) numel(unique(regexp(r,'[^&|\(\) ]+','match'))),grRules,'UniformOutput',false);

end



