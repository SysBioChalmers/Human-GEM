% This function is to extract the compartment prediction from DeepLoc2 tool
% to expand the gene.tsv file.
% Only when there is not such info from cellAtlas or SwissPort, we use the
% compartment prediction from the DeepLoc2.
% DeepLoc2 is used for the localization prediction, the protein fasta file
% for all metabolic proteins in the model was used as the input file, the
% default parameter was used in the prediction, output file is stored as
% DeepLoc2_compartment.tsv in ../../data/modelCuration

% Read the predicted compartment info
fileName = '../../data/modelCuration/DeepLoc2_compartment.csv';
fid  = fopen(fileName);
DeeplocData = textscan(fid,[repmat('%s ',[1,12]) '%s'],'Delimiter',',');
DeeplocData = [DeeplocData{1:end}];
DeeplocData = DeeplocData(2:end,:);
DeeplocData(:,[1,14]) = split(DeeplocData(:,1),';');
fclose(fid);

% Read gene.tsv
fileName = '../../model/genes.tsv';
fid  = fopen(fileName);
geneData = textscan(fid,[repmat('%s ',[1,9]) '%s'],'Delimiter','\t');
geneData = [geneData{1:end}];
fclose(fid);
geneData = strrep(geneData,'"','');

% index the DeeplocDdata by geneData order
DeeplocDdata_sorted(1:length(geneData(:,1)),1:length(DeeplocData(1,:))) = {''};
[~,idx] = ismember(geneData(:,1),DeeplocData(:,1));
DeeplocDdata_sorted(idx~=0,:) = DeeplocData(idx(idx~=0),:);


% build a dict for key word mapping of compartment assignment
DeepLoc_keywords = {'Cytoplasm','Cytosol';
'Nucleus','Nucleus';
'Extracellular','Extracellular';
'Cell membrane','Cell membrane';
'Mitochondrion','Mitochondria';
'Plastid','';
'Endoplasmic reticulum','Endoplasmic reticulum';
'Lysosome/Vacuole','Lysosome';
'Golgi apparatus','Golgi apparatus';
'Peroxisome','Peroxisome'
'|',';'};
for i = 1:length(DeepLoc_keywords(:,1))
    DeeplocDdata_sorted(:,2) = strrep(DeeplocDdata_sorted(:,2),DeepLoc_keywords{i,1},DeepLoc_keywords{i,2});
end
% fill the compartment for protein without localization info from CellAtlas
%and SwissPort
idx = cellfun(@isempty,geneData(:,9)) & ~(cellfun(@isempty,DeeplocDdata_sorted(:,1)));
geneData(idx,9) = DeeplocDdata_sorted(idx,2);
geneData(idx,10) = {'DeepLoc2'};

% write out gene.tsv
out = geneData';
fID = fopen('../../model/genes.tsv','w');
fprintf(fID,[repmat('%s\t',[1,9]) '%s\n'],out{1:10});
fprintf(fID,[repmat('"%s"\t',[1,9]) '"%s"\n'],out{11:end});
fclose(fID);

