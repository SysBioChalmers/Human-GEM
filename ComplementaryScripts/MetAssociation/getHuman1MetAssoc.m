% FILE NAME:    getHuman1MetAssoc.m
%
% DATE CREATED: 2019-05-23
%
% PROGRAMMERS:  Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
%
% PURPOSE: Sort out and refine the previous .mat files generated during
%          metabolite association/curation, then extensively extract
%          out exteranl identifiers and save to JSON format (#75).
%

%% Convert metAssocHMR2Recon3 file format from matlab to JSON

% Load met association information
load('metAssocHMR2Recon3.mat');

% change variable name and extract field names
m = metAssocHMR2Recon3;
fields = fieldnames(m);

% Convert elements of two fields (metChEBIID, metR3DID) from cell to string
for i=1:numel(fields)
    if iscell(m.(fields{i}){1})
        m.(fields{i}) = reformatElements(m.(fields{i}),'cell2str');
    end
end
m.metChEBIID = regexprep(m.metChEBIID, 'CHEBi:', '');  % small refinements

% some fields need to be reformated for either
% adding space after delimiter (metMNXID, metCuratedMNXID) or
% removing additional spaces (metR3DFormulas,metR3DKEGGID,metCuratedFormulas)
reformatFields = {'metMNXID';'metCuratedMNXID';'metR3DFormulas';'metR3DKEGGID';'metCuratedFormulas'};
for i=1:numel(reformatFields)
    convert2Cell = reformatElements(m.(reformatFields{i}),'str2cell');
    m.(reformatFields{i})  = reformatElements(convert2Cell,'cell2str');
end

% The metLIPIDMAPSID field has three elements that need to manually fixed
% m.metLIPIDMAPSID{336} ='LMFA01050113 LMFA01050349 LMFA01050359 LMFA02000035';
% m.metLIPIDMAPSID{1216}='LMFA01070018 LMFA02000037';
% m.metLIPIDMAPSID{606} ='LMST01010086;LMST01010144';
m.metLIPIDMAPSID{336} ='LMFA01050113; LMFA01050349; LMFA01050359; LMFA02000035';
m.metLIPIDMAPSID{1216}='LMFA01070018; LMFA02000037';
m.metLIPIDMAPSID{606} ='LMST01010086; LMST01010144';

% get the index PAPs and update it formula, which was fixed in #81
metsInd = find(strcmp(m.metHMRID, 'm02682'));
m.metFormulas(metsInd)={'C10H11N5O13P2S'};

% Convert elements of metCuratedCharges field from number to string
m.metCuratedCharges=cellfun(@num2str, m.metCuratedCharges, 'UniformOutput', false);



