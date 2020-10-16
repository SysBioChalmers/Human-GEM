function [orthologPairs, orthologStructure] = extractAllianceGenomeOrthologs(homologFilename, countBest)
% extractAllianceGenomeOrthologs
%   Read the JSON format ortholog pairs downloaded from Alliance Genome
%   database (alliancegenome.org) into a structure that is futher filtered
%   according to following criteria:
%
% 1. Filter out the pairs that are neither bestForward nor bestReverse
% 2. Keep all single-hit pairs
% 3. For the rest, keep those bestForward and bestReverse are both 'Yes'
% 4. If step 3 excludes all hits, only keep the one(s) has the highest 
%    `methodCount`
%
% Input:
%   homologFilename    the homolog file downloaded from AllianceGenome
%                      through API
%
%   countBest          filter out the ones that are neither best forward
%                      nor best reverse orthologs (optional, default TRUE)
%
% Output:
%   orthologPairs      ortholog pairs in NX2 cell array
%
%   orthologStructure  the output structure with following fields: 
%              fromGeneId:       query gene id
%              fromSymbol        query gene symbol
%              toGeneId:         gene id in the targe organism
%              toSymbol:         gene symbol of targe organism
%              best:             if the best forward match (Yes/No)
%              bestReverse:      if the best reverse match (Yes/No)
%              methodCount:      number of supportive methods
%              totalMethodCount: total number of ortholog finding methods
%
% Usage: [orthologPairs, orthologStructure] = extractAllianceGenomeOrthologs(homologFilename, countBest)
%
%

if nargin<2
    countBest = true;
end

if ~(exist(homologFilename,'file')==2)
    error('Input file %s cannot be found',string(homologFilename));
end

% load input file
inputPairs = jsondecode(fileread(homologFilename));

% define output structure
fieldList = {'fromGeneId',
             'fromSymbol'
             'toGeneId',
             'toSymbol',
             'best',
             'bestReverse',
             'methodCount',
             'totalMethodCount'};


% initialize output structure
for i = 1:length(fieldList)
    orthologStructure.(fieldList{i}) = {};
end


% loop through the content to extract corresponding values
for i=1:length(inputPairs.results)
    orthologStructure.fromGeneId = [orthologStructure.fromGeneId; inputPairs.results(i).gene.id];
    orthologStructure.fromSymbol = [orthologStructure.fromSymbol; inputPairs.results(i).gene.symbol];
    orthologStructure.toGeneId   = [orthologStructure.toGeneId; inputPairs.results(i).homologGene.id];
    orthologStructure.toSymbol   = [orthologStructure.toSymbol; inputPairs.results(i).homologGene.symbol];
    for j = 5:length(fieldList)
        orthologStructure.(fieldList{j}) = [orthologStructure.(fieldList{j}); inputPairs.results(i).(fieldList{j})];
    end
end


% check countBest
if countBest
    % get the index of ortholog-pairs that are neigher best forward nor reverse
    indToRemove = intersect(find(ismember(orthologStructure.best,'No')),...
                            find(ismember(orthologStructure.bestReverse,'No')));
    for i = 1:length(fieldList)
        orthologStructure.(fieldList{i})(indToRemove) = [];
    end
end


%% filter out ortholog with multiple pairs

% count the number of mapped orthologs
countHits = countFrequency(orthologStructure.fromGeneId);


% keep the pairs where there is only 1 ortholog
indSingleHit = find(countHits.frequency == 1);
indSingle = find(ismember(orthologStructure.fromGeneId, countHits.uniqueList(indSingleHit)));

output.from = orthologStructure.fromSymbol(indSingle);
output.to   = orthologStructure.toSymbol(indSingle);


% process the pairs with more ortholog hits
indMoreHit = find(countHits.frequency > 1);

for i=1:length(indMoreHit)
    indMore = find(ismember(orthologStructure.fromGeneId, countHits.uniqueList{indMoreHit(i)}));
    
    % try to keep only those bestForward and bestReverse are both 'Yes'
    noBestForwardReverse = 1;
    for j=1:length(indMore)
        if isequal(orthologStructure.best{indMore(j)}, 'Yes') &&...
            isequal(orthologStructure.bestReverse{indMore(j)}, 'Yes')
            output.from = [output.from; orthologStructure.fromSymbol{indMore(j)}];
            output.to   = [output.to;   orthologStructure.toSymbol{indMore(j)}];
            noBestForwardReverse = 0;        
        end
    end
    
    % if the pairs are neither bestForward nor bestReverse, then only
    % reserve the one with highest number of methodCount
    if noBestForwardReverse
        [~, topMethodCount]=maxk(cell2mat(orthologStructure.methodCount(indMore)), 1);
        output.from = [output.from; orthologStructure.fromSymbol{indMore(topMethodCount)}];
        output.to   = [output.to;   orthologStructure.toSymbol{indMore(topMethodCount)}];
    end
end

orthologPairs = [output.from, output.to];


