%
% FILE NAME:    curateModelECcodes.m
% 
% PURPOSE: This script corrects some errors in the eccodes field of
%          Human-GEM, trims the "EC:" preceding many codes, and removes
%          non-standard codes (TCDB and "Spontaneous").


% load model
ihuman = importHumanYaml('HumanGEM.yml');

% clean EC codes
ec = ihuman.eccodes;
ec = regexprep(ec, '\s', '');  % remove all spaces
ec = regexprep(ec, ';$', '');  % remove any trailing semicolons

% extract all TCDB codes
tcdb = repmat({''}, size(ec));
for i = 1:numel(ec)
    code = ec{i};
    if isempty(code)
        continue
    end
    code_pieces = strsplit(code, ';');
    is_tcdb = startsWith(code_pieces, 'TCDB');
    
    ec{i} = strjoin(code_pieces(~is_tcdb), ';');
    tcdb{i} = strjoin(code_pieces(is_tcdb), ';');
end

% remove all 'EC' or 'EC:' prefixes, and 'TCDB' or 'TCDB:' prefixes
ec = regexprep(ec, 'EC:*', '');
tcdb = regexprep(tcdb, 'TCDB:*', '');

% fix some problematic EC codes
ec = regexprep(ec, '^:', '');
ec(ismember(ec, '3.5.4-')) = {'3.5.4.-'};
ec(ismember(ec, '1.1.1.141,1.1.1.196')) = {'1.1.1.141;1.1.1.196'};
ec(ismember(ec, '1.3.3.64.2.1.17')) = {'1.3.3.6;4.2.1.17'};

% remove eccodes of "spontaneous" or "nonenzymatic"
spont = ismember(lower(ec), {'spontaneous','nonenzymatic'});
ec(spont) = {''};

% verify that all EC codes are in the correct format
badcodes = repmat({''}, size(ec));
badind = false(size(ec));
for i = 1:numel(ec)
    code = ec{i};
    if isempty(code)
        continue
    end
    code_pieces = strsplit(code, ';');
    isbad = cellfun(@isempty, regexp(code_pieces, '^(\d+|-)\.(\d+|-)\.(\d+|-)\.(\d+|-)$'));
    badcodes{i} = strjoin(code_pieces(isbad), ';');
    if any(isbad)
        badind(i) = true;
    end
end
if any(badind)
    error('Investigate invalid EC codes!');
end

% load rxnAssoc JSON
rxnAssoc = jsondecode(fileread('../../data/annotation/humanGEMRxnAssoc.JSON'));
if ~isequal(rxnAssoc.rxns, ihuman.rxns)
    error('rxnAssoc and Human-GEM reaction fields are not aligned!');
end

% make new "rxnTCDBID" field
rxnAssoc.rxnTCDBID = tcdb;

% make new "spontaneous" field
rxnAssoc.spontaneous = double(spont);

% export rxnAssoc to JSON
jsonStr = jsonencode(rxnAssoc);
fid = fopen('../../data/annotation/humanGEMRxnAssoc.JSON', 'w');
fwrite(fid, prettyJson(jsonStr));
fclose(fid);

% update eccodes field in Human-GEM
ihuman.eccodes = ec;

% export to yml
writeHumanYaml(ihuman, '../../model/HumanGEM.yml');


