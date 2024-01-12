library(ggplot2)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


figPath = "figures/"


############################################
# Fig 1B - Gene essentiality
############################################

# load data
# fns = c("data/geneEss_model1.txt",
#         "data/geneEss_newalg.txt",
#         "data/geneEss_newalg2.txt"
#         )
fns = c("data/geneEss_model1.txt")


#names = c("tINIT","ftINIT 1+0", "ftINIT 1+1")
names = c("Human2")

gea_res = NULL
for (i in 1:length(fns)) {
  x = read.delim(file = fns[i], sep='\t', stringsAsFactors=F)
  x$model = names[i]
  gea_res = rbind(gea_res, x)
}
gea_res$model = factor(gea_res$model, as.character(names)[1:length(fns)])  # to enforce the model order


color_palette <- c('#B5D39B','#6B97BC','#E7B56C')  # light green, light blue, light yellow

p1B = ggplot(gea_res, aes(x = model, y = MCC, fill = model)) +
  geom_violin(trim=F, show.legend=F, scale='count') +
  scale_fill_manual(values=color_palette) +
  theme_classic() + 
  ylab('MCC') +
  xlab('') +
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,
                                   color='black', size=14),
        axis.text.y = element_text(color='black', size=14),
        axis.line.x = element_blank()) +
  ylim(c(0.08,0.40)) # + #manipulate these numbers to include all data
  #ylim(c(0,0.5)) # +
p1B


ggsave(
  paste0(figPath, "FigGeneEss.png"),
  plot = p1B,
  width = 3.5, height = 3.2, dpi = 300)

ggsave(
  paste0(figPath, "FigGeneEss.eps"),
  plot = p1B,
  width = 3.5, height = 3.2, dpi = 300)
  
  
