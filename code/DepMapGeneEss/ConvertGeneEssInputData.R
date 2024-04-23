library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

d = read_csv(paste0(rstudioapi::getActiveDocumentContext(), "/data/CCLE_expression_full.csv"))
dim(d)#1377 52055

#transpose
ds = as_tibble(cbind(gene = names(d)[-1], t(as.matrix(d[,-1]))))
colnames(ds)[-1] = d[[1]]

#make sure the columns are numeric instead of chr
ds2 = ds
for (i in 2:ncol(ds)) {
  ds2[[i]] = as.numeric(ds[[i]])

}

#convert the gene names

#To get ensembl: 
#pattern = ".*\\(([A-Z0-9]*)\\)"
#newGenes = str_match(ds$gene, pattern)
#ds2$gene = newGenes[,2]
#fix the ERCC genes
#ds2$gene[is.na(newGenes[,2])] = ds$gene[is.na(newGenes[,2])];

#To get gene symbols:
#It is a bit tricky, not all genes follow the pattern. Some are like LINC00328-2P (ENSG00000225016),
#some just an ensembl id (we then take the assembl id), some ERCC 
newGenes = ds2$gene
x  = strsplit(ds2$gene[!ERCCGenesSel], " ")
for(i in 1:length(x)) {
  newGenes[i] = x[[i]][1] #handles all cases
}
length(newGenes)
length(unique(newGenes)) #not the same, we need to merge a few rows, done later

ds2$gene = newGenes

ds2


#now convert the data. It is currently as log2(TPM + 1)
dsTPM = ds2
dsTPM[,-1] = 2^ds2[,-1] - 1
colSums(dsTPM[,-1])#very minor roundoff differences, ok

#Sum up the rows that have the same gene name
duplGenes = unique(dsTPM$gene[duplicated(dsTPM$gene)])
length(duplGenes)#9
dsTPM[dsTPM$gene %in% duplGenes,1:10]
#CCDC39 is a good example to test
#is 1.59 + 0.150 in the first row, those should be summed up

rowsToRem = rep(FALSE, nrow(dsTPM))
for (i in 1:length(duplGenes)) {
  inds = which(dsTPM$gene == duplGenes[i])
  dsTPM[inds[1], -1] = colSums(dsTPM[inds, -1]) #the first row now gets the sum of all rows
  rowsToRem[inds[-1]] = TRUE #the other rows are marked for deletion
}
sum(rowsToRem) #9, looks good
#test: 
dsTPM[[2]][dsTPM$gene == "CCDC39"] #1.74 0.15, as expected, ok
dsTPM$gene[rowsToRem] # CCDC39 is in there, ok
dsTPM2 = dsTPM[!rowsToRem,]

write_tsv(dsTPM2, paste0(dataFolder, "DepMap_tpm_gene_symbols.txt"))

