#Author Min Yang

library("msigdbr")
library("dplyr")
library("GSVA")
#check the supporting species
msigdbr_show_species()
#mouse gene set
m_df = msigdbr(species = "Mus musculus")
#checking gene set
a <- m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
#selecting gene set
m_df = msigdbr(species = "Mus musculus", category = "H")
#constract the reference dataset
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
#get the normalized data
Mac.data <- as.matrix(Mac@assays$RNA@data)
meta <- Mac@meta.data[,c("tdTomato")]
meta <- gsub("TRUE", "positive", meta)
meta <- gsub("FALSE", "negative", meta)
#run GSVA
mac.cp <- gsva(Mac.data, msigdbr_list, kcdf = "Gaussian", method = "gsva", parallel.sz = 20)
#differential pathways using limma
library(limma)
group <- factor(meta,levels = c("positive", "negative"),ordered = F)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(Mac.data)
fit <- lmFit(mac.cp, design)
con.matrx <- makeContrasts(contrasts = c("positive-negative"), levels = design)
fit2 <- contrasts.fit(fit, con.matrx)
fit2 <- eBayes(fit2)
mac.pathway <- topTable(fit2, coef = "positive-negative", number = 20)
mac.pathway$pathway <- rownames(mac.pathway)
mac.pathway$pathway <- factor(mac.pathway$pathway, levels = mac.pathway[order(mac.pathway$t), ]$pathway)
ggplot() + 
  geom_bar(mac.pathway, mapping = aes(x = t, y = pathway), stat = "identity") +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "black", fill = "white"))
