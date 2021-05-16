%
% FILE NAME:    createComplexStructure.m
% 
% PURPOSE: Create enzymeComplex structure for HMR curation
%	


% 1. Dump in external database
% Move to the path
cd('/Users/haowa/Box Sync/HMR3/Complex-Subunit/CORUM/');
T=readtable('coreComplexes.txt','ReadVariableNames',1);
CORUM=table2struct(T,'ToScalar',true);
CORUM.ComplexID=num2cell(CORUM.ComplexID);
CORUM.ComplexID=cellfun(@num2str,CORUM.ComplexID,'un', 0);
CORUM.PubMedID=num2cell(CORUM.PubMedID);
CORUM.PubMedID=cellfun(@num2str,CORUM.PubMedID,'un', 0);


% 2. Add the field of subunits with UniProt ids
% Corum
subunitsUniProtIDs={};
for i=1:numel(CORUM.ComplexID)
		subunitsUniProtIDs{i}=strsplit(CORUM.subunits_UniProtIDs_{i},';');
		subunitsUniProtIDs{i}=regexprep(subunitsUniProtIDs{i},'-\d$','');
		subunitsUniProtIDs{i}=regexprep(subunitsUniProtIDs{i},'-$','');
		subunitsUniProtIDs{i}=unique(subunitsUniProtIDs{i});
end
CORUM.subunitsUniProtIDs=transpose(subunitsUniProtIDs);


% 3. Generate Human complexes for HMR GTRs curation
human_id=find(strcmp('Human', CORUM.Organism));
multi_ind=find(cellfun(@numel,CORUM.subunitsUniProtIDs)>1);
index=intersect(human_id,multi_ind);
CORUM.HumanCplxID=CORUM.ComplexID(index);
CORUM.HumanSubunits=CORUM.subunitsUniProtIDs(index);

UniProtID=reformatElements(CORUM.HumanSubunits,'cell2str');
UniProtID=strjoin(UniProtID,';');
UniProtID=unique(strsplit(UniProtID,';'));
CORUM.HumanUniProtID=transpose(UniProtID);


% 4. Generate complex-subunit matrix
% Corum
cplxSubMat=zeros(numel(CORUM.HumanUniProtID),numel(CORUM.HumanCplxID));
for i=1:numel(CORUM.HumanCplxID)
		[~, index]=ismember(CORUM.HumanSubunits{i},CORUM.HumanUniProtID);
		cplxSubMat(index,i)=1;
end
CORUM.HumanCplxSubMat=sparse(cplxSubMat);


% 5. Associate with Ensembl ids
load('Ensembl2Uniprot.mat');
CORUM.EnsemblID=cell(numel(CORUM.HumanUniProtID),1);
CORUM.EnsemblID(:)={''};
[a, b]=ismember(CORUM.HumanUniProtID,Ensembl2Uniprot.SwissProtID);
CORUM.EnsemblID(find(a))=Ensembl2Uniprot.genes(b(find(a)));
[a, b]=ismember(CORUM.HumanUniProtID,Ensembl2Uniprot.TrEMBLID);
CORUM.EnsemblID(find(a))=Ensembl2Uniprot.genes(b(find(a)));


% 6. Remove duplicate complexes and save the structure
[uniqueMat, I, ~]=unique(transpose(CORUM.HumanCplxSubMat),'stable','rows');
CORUM.HumanCplxID=CORUM.HumanCplxID(I);
CORUM.HumanSubunits=CORUM.HumanSubunits(I);
CORUM.HumanCplxSubMat=transpose(uniqueMat);
save('CORUM.mat','CORUM');   % 2018-06-10
