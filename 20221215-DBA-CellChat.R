library(CellChat)
library(ggplot2)
library(ggalluvial)
library(NMF)
library(ComplexHeatmap)
library(tidyverse)
cellsubset <- c("D1_eSF", "D2_Pre-DSC", "D3_Top2a", "D4_Ptn", "D5_Gatm", "D6_S100a8", "D7_Prl8a2", 
                "SMC", "Epithelial", "Mac", "NK", "Neutrophil", "Proliferating Mac", 
                "DC-1", "T cell", "Proliferating EC", "Angiogenic EC", "Arterial EC",  "Venous EC1", "Venous EC2", "EPC",
                "iDSC0", "iDSC1", "iDSC2",
                "tropho-0", "tropho-1", "tropho-2", "tropho-3", "tropho-4", "tropho-5", "Viseral endoderm")
load("20220622-Comments/05.DBA/DBA.1-merge-embryo-SCT.RData")
DBA.1.ST.inte$cellcluster <- gsub("Proliterating EC", "Proliferating EC", DBA.1.ST.inte$cellcluster)
DBA.1.cellchat <- subset(DBA.1.ST.inte, subset = cellcluster %in% cellsubset)
load("20210827-figures/files/DBA.1-iDSC-three-subclusters-20230306.RData")
DimPlot(DBA.1.iDSC, label = T)
DBA.1.iDSC <- RenameIdents(DBA.1.iDSC, "DBA.1.iDSC0" = "iDSC0", "DBA.1.iDSC1" = "iDSC1", "DBA.1.iDSC2" = "iDSC2")
DimPlot(DBA.1.iDSC)
DBA.1.iDSC$cellcluster <- DBA.1.iDSC@active.ident
DBA.1.cellchat <- merge(DBA.1.cellchat, DBA.1.iDSC)
DBA.1.cellchat@active.ident <- factor(DBA.1.cellchat$cellcluster)
DBA.1.cellchat <- NormalizeData(DBA.1.cellchat)
#subsample 1000 cells
DBA.1.cellchat <- subset(DBA.1.cellchat, downsample = 1000)
DBA.1.cellchat <- NormalizeData(DBA.1.cellchat)
DBA.1.cellchat$cellname <- DBA.1.cellchat@active.ident
options(stringsAsFactors = FALSE)
data.input <- DBA.1.cellchat@assays$RNA@data
meta <- DBA.1.cellchat@meta.data
meta <- meta[, c(1:4,7)]
#创建cellchat对象
DBA.1.cellchat <- createCellChat(object = DBA.1.cellchat, meta = meta, group.by = "cellcluster")
#set default idents
DBA.1.cellchat <- setIdent(DBA.1.cellchat, ident.use = "cellcluster", levels = cellsubset)
groupSize <- as.numeric(table(DBA.1.cellchat@idents))
#set the ligand-receptor interaction database
DBA.1.cellchatDB <- CellChatDB.mouse
DBA.1.cellchatDB.use <- subsetDB(DBA.1.cellchatDB, search = "Secreted Signaling")
DBA.1.cellchat@DB <- DBA.1.cellchatDB.use
#processing the expression data for cell-cell communication
DBA.1.cellchat <- subsetData(DBA.1.cellchat)
future::plan("multiprocess", workers = 20)#do parallel
DBA.1.cellchat <- identifyOverExpressedGenes(DBA.1.cellchat)
DBA.1.cellchat <- identifyOverExpressedInteractions(DBA.1.cellchat)
DBA.1.cellchat <- projectData(DBA.1.cellchat, PPI.mouse)
#inference of cell-cell communication network
#1) compute the communication probability and infer the cellular communication network
DBA.1.cellchat <- computeCommunProb(DBA.1.cellchat)
DBA.1.cellchat <- filterCommunication(DBA.1.cellchat, min.cells = 10)
#2) infer the cell-cell communication at a signaling pathway level
DBA.1.cellchat <- computeCommunProbPathway(DBA.1.cellchat)
#3) calculate the aggregated cell-cell communication network
DBA.1.cellchat <- aggregateNet(DBA.1.cellchat)
DBA.1.cellchat <- netAnalysis_computeCentrality(DBA.1.cellchat, slot.name = "netP")
save(DBA.1.cellchat, file = "20210827-figures/files/DBA1-cellchat-20230306.RData")
#global
ht1 <- netAnalysis_signalingRole_heatmap(DBA.1.cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(DBA.1.cellchat, pattern = "incoming")
ht1 + ht2
load("../20210201-integrate/UE-harmony-cellchat-subname-FB-subclusters-results-20211007.RData")
ht1 <- netAnalysis_signalingRole_heatmap(UE.harmony.cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(UE.harmony.cellchat, pattern = "incoming")
ht1 + ht2

#compare DBA and Normal cellchat results
load("20210827-figures/files/DBA1-cellchat-20221221.RData")
netAnalysis_signalingRole_heatmap(DBA.1.cellchat, pattern = "outgoing", width = 5, height = 6)
DBA.1.cellchat <- updateCellChat(DBA.1.cellchat)
#UE.harmony.cellchat
load("../20210201-integrate/UE-harmony-cellchat-subname-FB-subclusters-results-20221221.RData")
UE.harmony.cellchat <- updateCellChat(UE.harmony.cellchat)
netAnalysis_signalingRole_heatmap(UE.harmony.cellchat, pattern = "outgoing", width = 5, height = 6)

#DBA.1.cellchat
DBA.1.cellchat@meta$cellcluster <- gsub("Proliterating EC", "Proliferating EC", DBA.1.cellchat@meta$cellcluster)
DBA.1.cellchat@meta$cellcluster <- factor(DBA.1.cellchat@meta$cellcluster, 
                                          levels = c(levels(DBA.1.cellchat@idents)[1:15], "Proliferating EC", levels(DBA.1.cellchat@idents)[17:32]))
DBA.1.cellchat <- setIdent(DBA.1.cellchat, ident.use = "cellcluster")

group.new1 <- c(levels(UE.harmony.cellchat@idents), "EPC", "iDSC0-2")
group.new2 <- c(levels(DBA.1.cellchat@idents), "B", "DC-2", "EC-5")
UE.harmony.cellchat <- liftCellChat(UE.harmony.cellchat, group.new = group.new1)
DBA.1.cellchat <- liftCellChat(DBA.1.cellchat, group.new = group.new2)

DBA.1.cellchat <- updateCellChat(DBA.1.cellchat)
UE.harmony.cellchat <- updateCellChat(UE.harmony.cellchat)
object.list <- list(Normal = UE.harmony.cellchat, DBA = DBA.1.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
#differential interaction strength 
netVisual_heatmap(cellchat, measure = "weight")
pdf("20220622-Comments/05.DBA/00.figures/45.Normal-DBA-merge-cellchat-rank.pdf", width = 10, height = 10)
rankNet(cellchat, mode = "single", stacked = F, do.stat = F)
dev.off()

rankNet(DBA.1.cellchat, mode = "single", stacked = F, do.stat = F)

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
pdf("20220622-Comments/05.DBA/00.figures/40.Normal-DBA-cellchat-outgoing-heatmap.pdf", width = 10, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
pdf("20220622-Comments/05.DBA/00.figures/40.Normal-DBA-cellchat-incoming-heatmap.pdf", width = 10, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

netVisual_bubble(cellchat, sources.use = 24, targets.use = 18:23,  comparison = c(1, 2), angle.x = 45)

signaling.use <- c("CALCR", "ANGPTL", "ANGPT", "CHEMERIN", "SEMA3", "VEGF")
source.cell <- c("D3_Top2a", "D5_Gatm", "D6_S100a8", "D7_Prl8a2", "iDSC0")
target.cell <- c("Proliferating EC", "Angiogenic EC", "Arterial EC", "Venous EC1", "Venous EC2")
gg1 <- netVisual_bubble(cellchat, sources.use = source.cell, targets.use = target.cell, comparison = c(1, 2), max.dataset = 2, signaling = signaling.use[6], title.name = "Increased signaling in DBA", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = source.cell[c(3,5)], targets.use = target.cell, comparison = c(1, 2), signaling = signaling.use, title.name = "Decreased signaling in DBA", angle.x = 45, remove.isolate = F)
#> Comparing communications on a merged object
gg3 <- netVisual_bubble(cellchat, sources.use = c(3,5:7), targets.use = c(26),  comparison = c(1, 2), max.dataset = 1, signaling = "CD137", title.name = "Decreased signaling in DBA", angle.x = 45, remove.isolate = T)

gg2.test <- netVisual_bubble(cellchat, sources.use = source.cell[c(3,5)], targets.use = target.cell, comparison = c(1, 2), signaling = signaling.use[1:4], title.name = "Decreased signaling in DBA", angle.x = 45, remove.isolate = F)

#gg2 data
gg2.data <- gg2$data[, c(5, 8, 9, 12, 14, 15)]
gg2.data <- gg2.data[gg2.data$pathway_name != "VEGF", ]
pdf("20220622-Comments/05.DBA/00.figures/32.DBA-iDSC0-D6-EC-maturation-decrease-plot.pdf", width = 10, height = 10)
ggplot() + 
  geom_point(gg2.data, mapping = aes(x = source.target, y = interaction_name_2, color = prob, size = prob)) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral"))) +
  scale_size_continuous(range = c(1,10)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

#gg1 data
data <- gg1$data[, c(5, 8, 9, 12, 14, 15)]
gg1.data <- gg1.data[gg1.data$pathway_name == "VEGF", ]
st.levels <- levels(gg1.data$source.target)
gg1.data$source.target <- factor(gg1.data$source.target, levels = c(st.levels[1:30],
                                                                    "D7_Prl8a2 -> Proliferating EC (Normal)", st.levels[31], 
                                                                    "D7_Prl8a2 -> Angiogenic EC (Normal)", st.levels[32], 
                                                                    "D7_Prl8a2 -> Arterial EC (Normal)", st.levels[33], 
                                                                    "D7_Prl8a2 -> Venous EC1 (Normal)", st.levels[34], 
                                                                    "D7_Prl8a2 -> Venous EC2 (Normal)", st.levels[35:45]))

gg1.data <- rbind(gg1.data, c(NA, "Pgf  - Vegfr1", "VEGF", "D7_Prl8a2 -> Proliferating EC (Normal)", "D7_Prl8a2 -> Proliferating EC", "Normal"), 
                  c(NA, "Pgf  - Vegfr1", "VEGF", "D7_Prl8a2 -> Angiogenic EC (Normal)", "D7_Prl8a2 -> Angiogenic EC", "Normal"),
                  c(NA, "Pgf  - Vegfr1", "VEGF", "D7_Prl8a2 -> Arterial EC (Normal)", "D7_Prl8a2 -> Arterial EC", "Normal"), 
                  c(NA, "Pgf  - Vegfr1", "VEGF", "D7_Prl8a2 -> Venous EC1 (Normal)", "D7_Prl8a2 -> Venous EC1", "Normal"), 
                  c(NA, "Pgf  - Vegfr1", "VEGF", "D7_Prl8a2 -> Venous EC2 (Normal)", "D7_Prl8a2 -> Venous EC2", "Normal"))
gg1.data$prob <- as.numeric(gg1.data$prob)
pdf("20220622-Comments/05.DBA/00.figures/32.DBA-iDSC0-D6-EC-angiogenic-increase-plot.pdf", width = 16, height = 5)
ggplot() + 
  geom_point(gg1.data, mapping = aes(x = source.target, y = interaction_name_2, color = prob, size = prob)) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral"))) +
  scale_size_continuous(range = c(0,10)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

#gg3 data
gg3.data <- gg3$data[, c(5, 8, 9, 12, 14, 15)]
st.levels <- levels(gg3.data$source.target)
gg3.data$source.target <- factor(gg3.data$source.target, levels = c(st.levels[1], "D3_Top2a -> iDSC2 (DBA)",
                                                                    st.levels[2], "D5_Gatm -> iDSC2 (DBA)", 
                                                                    st.levels[3], "D6_S100a8 -> iDSC2 (DBA)"))
gg3.data <- rbind(gg3.data, c(NA, "Tnfsf9  - Tnfrsf9", "CD137", "D3_Top2a -> iDSC2 (DBA)", "D3_Top2a -> iDSC2", "DBA"), 
                  c(NA, "Tnfsf9  - Tnfrsf9", "CD137", "D5_Gatm -> iDSC2 (DBA)", "D5_Gatm -> iDSC2", "DBA"),
                  c(NA, "Tnfsf9  - Tnfrsf9", "CD137", "D6_S100a8 -> iDSC2 (DBA)", "D6_S100a8 -> iDSC2", "DBA"))
gg3.data$prob <- as.numeric(gg3.data$prob)
pdf("20220622-Comments/05.DBA/00.figures/32.DBA-iDSC0-D6-EC-angiogenic-increase-plot.pdf", width = 20, height = 5)
ggplot() + 
  geom_point(gg3.data, mapping = aes(x = source.target, y = interaction_name_2, color = prob, size = prob)) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral"))) +
  scale_size_continuous(range = c(0,10)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

#Gene expression
load("20210827-figures/files/DBA.1-DSC-20221212.Rdata")
DBA.1.DSC$cellcluster <- gsub("_.*", "", DBA.1.DSC$cellcluster)
DBA.1.DSC$sample <- "DBA"
load("20210827-figures/files/UE-DSC-web.RData")
UE.decidual.susbet$cellcluster <- gsub("_.*", "", UE.decidual.susbet$sub_clusters)
UE.decidual.susbet$sample <- "Normal"
UE.decidual.susbet <- subset(UE.decidual.susbet, subset = cellcluster %in% c("D3", "D4", "D5", "D6", "D7"))
load("20210827-figures/files/DBA.1-iDSC-subclusters-20221212.RData")
DBA.1.iDSC$cellcluster <- gsub("DBA.1.", "", DBA.1.iDSC$cellcluster)
DBA.1.iDSC$cellcluster <- gsub("iDSC0-1", "iDSC0", DBA.1.iDSC$cellcluster)
DBA.1.iDSC$sample <- "DBA"
load("20210827-figures/figure5/FB-immune-sc-subclusters-SCT.RData")
FB.immune.SCT$cellcluster <- gsub("Normal", "iDSC", FB.immune.SCT$subcluster)
FB.immune.SCT$sample <- "Normal"
Normal.DBA.DSC.iDSC <- merge(UE.decidual.susbet, y = list(DBA.1.DSC, FB.immune.SCT, DBA.1.iDSC))
Normal.DBA.DSC.iDSC <- NormalizeData(Normal.DBA.DSC.iDSC)

load("20210827-figures/files/UE-EC-web.RData")
UE.EC$sample <- "Normal"
UE.EC$sub_clusters <- gsub("Venous EC2", "AAAAAA", UE.EC$sub_clusters)
UE.EC$sub_clusters <- gsub("Venous EC", "Venous EC1", UE.EC$sub_clusters)
UE.EC$sub_clusters <- gsub("AAAAAA", "Venous EC2", UE.EC$sub_clusters)
UE.EC$cellcluster <- UE.EC$sub_clusters
UE.EC <- subset(UE.EC, )
load("20220622-Comments/05.DBA/DBA-EC-clusters-20221215.RData")
DimPlot(DBA.1.EC, label = T)
DBA.1.EC$sample <- "DBA"
DBA.1.EC$cellcluster <- gsub("Proliterating EC", "Proliferating EC", DBA.1.EC$cellcluster)
DBA.1.EC$cellcluster <- gsub("EPC", "Inflammation EC", DBA.1.EC$cellcluster)
Normal.DBA.EC <- merge(UE.EC, DBA.1.EC)
Normal.DBA.EC$cellcluster <- factor(Normal.DBA.EC$cellcluster, levels = c("Proliferating EC", "Angiogenic EC", "Arterial EC", "Venous EC1", "Venous EC2", "Lymphatic EC", "Inflammation EC"))
Normal.DBA.EC <- NormalizeData(Normal.DBA.EC)

Angiogenesis.ligand <- c("Vegfa", 'Vegfb', "Vegfd", "Adm", "Calca", "Rarres2", "Angptl2", "Angptl4", "Angpt2", "Sema3a", "Sema3b", "Sema3c", "Sema3e")
Angiogenesis.receptor <- c("Flt1", 'Kdr', "Flt4", "Calcrl", "Cmklr1", "Tlr4", "Cdh5", "Cdh11", "Sdc1", "Sdc2", "Sdc4", "Itga5", "Itgav", "Itgb1", "Itgb3", "Tek", "Plxnd1", "Nrp1", "Nrp2", "Plxna1", "Plxna2")
VlnPlot(Normal.DBA.DSC.iDSC, Angiogenesis.ligand[7:13], split.by = "sample", group.by = "cellcluster", pt.size = 0, combine = T, split.plot = F, log = T, ncol = 3)
plots <- VlnPlot(Normal.DBA.DSC.iDSC, Angiogenesis.ligand, split.by = 'sample', group.by = "cellcluster", pt.size = 0, combine = F, log = T, split.plot = T)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = c('#E87D72', '#56BCC2'))
)
pdf("20220622-Comments/05.DBA/00.figures/34.Normal-DBA-DSC-iDSC-ligand-expression-violin.pdf", width = 15, height = 10)
CombinePlots(plots = plots, legend = 'right')
dev.off()

pdf("20220622-Comments/05.DBA/00.figures/34.Normal-DBA-EC-expression-violin.pdf", width = 15, height = 15)
VlnPlot(Normal.DBA.EC, Angiogenesis.receptor, split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 5, cols = c('#E87D72', '#56BCC2'))
dev.off()


#CD137 pathway
load("20210827-figures/files/DBA.1-immune-20221215.RData")
DimPlot(DBA.1.immune.m, label = T)
DBA.1.immune.m$sample <- "DBA"
load("20210827-figures/files/UE-immune-harmony-web.RData")
DimPlot(UE.immune.harmony, label = T)
UE.immune.harmony$sample <- "Normal"
UE.immune.harmony$cellcluster <- gsub("Mono", "Monocyte", UE.immune.harmony$sub_clusters)
UE.immune.harmony$cellcluster <- gsub("Cycling Mac", "Proliferating Mac", UE.immune.harmony$cellcluster)
UE.immune.harmony$cellcluster <- gsub("T", "T cell", UE.immune.harmony$cellcluster)
Normal.DBA.immune <- merge(UE.immune.harmony, DBA.1.immune.m)
Normal.DBA.immune@active.ident <- factor(Normal.DBA.immune$cellcluster)
Normal.DBA.immune <- subset(Normal.DBA.immune, subset = cellcluster %in% c("Monocyte", "Mac", "proliferating Mac", "Neutrophil", "NK"))
pdf("20220622-Comments/05.DBA/00.figures/43.CD137-DSC-iDSC-expression-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.DSC.iDSC, features = c("Tnfsf9", "Tnfrsf9"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, cols = c('#E87D72', '#56BCC2'))
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/43.CD137-Immune-expression-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.immune, features = c("Tnfsf9", "Tnfrsf9"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, cols = c('#E87D72', '#56BCC2'))
dev.off()

#Immune 相关的通路
#1. iDSC中的CXCL和CCL的ligands表达上调
DBA.Normal.cellchat.df <- subsetCommunication(cellchat, slot.name = "net")
DBA.cellchat.df.chemokine <- DBA.Normal.cellchat.df$DBA[DBA.Normal.cellchat.df$DBA$pathway_name %in% c("CXCL", "CCL"), ]
Normal.cellchat.df.chemokine <- DBA.Normal.cellchat.df$Normal[DBA.Normal.cellchat.df$Normal$pathway_name %in% c("CXCL", "CCL"), ]
chemokine.ligands <- unique(c(DBA.cellchat.df.chemokine$ligand, Normal.cellchat.df.chemokine$ligand))

load("20210827-figures/files/DBA.1-iDSC-subclusters-20221212.RData")
DBA.1.iDSC$cellcluster <- gsub("DBA.1.", "", DBA.1.iDSC$cellcluster)
DBA.1.iDSC$cellcluster <- gsub("iDSC0-1", "iDSC0", DBA.1.iDSC$cellcluster)
DBA.1.iDSC$sample <- "DBA"
load("20210827-figures/figure5/FB-immune-sc-subclusters-SCT.RData")
FB.immune.SCT$cellcluster <- gsub("Normal", "iDSC", FB.immune.SCT$subcluster)
FB.immune.SCT$sample <- "Normal"
DBA.Normal.iDSC <- merge(DBA.1.iDSC, FB.immune.SCT)
DBA.Normal.iDSC@active.ident <- factor(DBA.Normal.iDSC$cellcluster)
DBA.Normal.iDSC <- NormalizeData(DBA.Normal.iDSC)
DefaultAssay(DBA.Normal.iDSC) <- "RNA"
DBA.Normal.iDSC$cellcluster <- factor(DBA.Normal.iDSC$cellcluster, levels = c("iDSC0", "iDSC1", "iDSC2", "iDSC0-2"))
pdf("20220622-Comments/05.DBA/00.figures/39.DBA-Normal-iDSC-chemokines-expression-violin.pdf", width = 10, height = 10)
VlnPlot(DBA.Normal.iDSC, features = chemokine.ligands, split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 5, cols = c('#E87D72', '#56BCC2') )
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/39.DBA-Normal-iDSC-angiogenesis-granzyme-expression-violin.pdf", width = 10, height = 10)
VlnPlot(DBA.Normal.iDSC, features = c("Vegfd", "Vegfb", "Pgf", "Angptl4", "Angptl2", "Calca", "Adm", "Gzma", "Gzmb", "Gzme", "Spp1", "Prf1", "Tnfrdf9"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 5, cols = c('#E87D72', '#56BCC2') )
dev.off()


#CXCL, CCL通路的变化
signaling.use <- c("CXCL", "CCL")
source.cell <- c("iDSC0", "iDSC1", "iDSC2")
target.cell <- c("Mac", "Proliferating Mac", "DC-1", "NK", "Neutrophil")
gg4 <- netVisual_bubble(cellchat, sources.use = source.cell, targets.use = target.cell, comparison = c(1, 2), max.dataset = 2, signaling = signaling.use, angle.x = 45, remove.isolate = F)
#gg4 data
gg4.data <- gg4$data[, c(5, 8, 9, 12, 14, 15)]
gg4.data$interaction_name_2 <- factor(gg4.data$interaction_name_2, 
                                      levels = rev(c("Cxcl2  - Cxcr2", "Cxcl1  - Cxcr2", "Ccl12  - Ccr2", "Ccl7  - Ccr2", "Ccl2  - Ccr2", "Ccl7  - Ccr1", "Ccl3  - Ccr1", 
                                                     "Ccl8  - Ccr2", "Ccl6  - Ccr2", "Ccl8  - Ccr5", "Ccl4  - Ccr5", "Ccl3  - Ccr5", "Ccl8  - Ccr1", "Ccl6  - Ccr1", "Ccl9  - Ccr1")))
pdf("20220622-Comments/05.DBA/00.figures/48.DBA-iDSC-immune-CXCL-CCL-increase-plot.pdf", width = 15, height = 6)
ggplot() + 
  geom_point(gg4.data, mapping = aes(x = source.target, y = interaction_name_2, color = prob, size = prob)) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral"))) +
  scale_size_continuous(range = c(1,10)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

#COMPLEMENT通路的变化
signaling.use <- c("COMPLEMENT")
source.cell <- c("Epithelial")
target.cell <- c("iDSC0", "iDSC1", "iDSC2", "iDSC0-2")
gg5 <- netVisual_bubble(cellchat, sources.use = source.cell, targets.use = target.cell, comparison = c(1, 2), max.dataset = 2, signaling = signaling.use, angle.x = 45, remove.isolate = F)
#gg5 data
gg5.data <- gg5$data[, c(5, 8, 9, 12, 14, 15)]
pdf("20220622-Comments/05.DBA/00.figures/49.DBA-Epithelial-iDSC-COMPLEMENT-increase-plot.pdf", width = 4.8, height = 2.1)
ggplot() + 
  geom_point(gg5.data, mapping = aes(x = source.target, y = interaction_name_2, color = prob, size = prob)) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral"))) +
  scale_size_continuous(range = c(1,10)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()


library(stringr)
wnt.gene <- read.table("20220622-Comments/05.DBA/WNT_from_PB.csv", sep = ",", header = F, check.names = F)
wnt.gene <- sort(str_to_title(as.character(wnt.gene$V1)))
#DSC+iDSC
pdf("20220622-Comments/05.DBA/00.figures/50.DSC-iDSC-Normal-DBA-WNT-expression-violinplot.pdf", width = 40, height = 50)
VlnPlot(Normal.DBA.DSC.iDSC, features = wnt.gene, split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 10, cols = c('#E87D72', '#56BCC2') )
dev.off()
#EC
pdf("20220622-Comments/05.DBA/00.figures/50.EC-Normal-DBA-WNT-expression-violinplot.pdf", width = 40, height = 70)
VlnPlot(Normal.DBA.EC, features = wnt.gene, split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 10, cols = c('#E87D72', '#56BCC2') )
dev.off()

#compare gene expression in CBA and Normal
DBA.all <- DBA.1.cellchat
save(DBA.all, file = "20220622-Comments/05.DBA/DBA-all-subclusters.RData")
Normal.all <- UE.harmony.pc30.cellchat
save(Normal.all, file = "20220622-Comments/Normal-all-subclusters.RData")

load("20220622-Comments/Normal-all-subclusters.RData")
load("20220622-Comments/05.DBA/DBA-all-subclusters.RData")
Normal.DBA.all <- merge(Normal.all, DBA.all)
Normal.DBA.all <- NormalizeData(Normal.DBA.all)
#Vegf expression
VlnPlot(Normal.DBA.all, features = c("Vegfa", "Vegfb", "Vegfd", "Pgf"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
#Cytokine
Cytokine.path <- read.table("20220622-Comments/Cytokines-GO_term_summary_20230112_200034.txt", sep = "\t", header = T, check.names = F, row.names = NULL)
for (p in unique(Cytokine.path$Qualifier)) {
  genes <- unique(Cytokine.path[Cytokine.path$Qualifier == p, ]$`MGI Gene/Marker ID`)
  Normal.DBA.all.sub <- AddModuleScore(Normal.DBA.all, features = list(genes), name = p)
}
Normal.DBA.all.iDSC <- subset(Normal.DBA.all, cellcluster %in% c("iDSC1","iDSC2","iDSC0","iDSC0-2"))
pdf("20220622-Comments/05.DBA/00.figures/51.Normal-DBA-iDSC-celltype-Cytokin-score.pdf", width = 100, height = 50)
VlnPlot(Normal.DBA.all.iDSC, features = colnames(Normal.DBA.all.sub@meta.data)[11:69], split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 10, cols = c('#E87D72', '#56BCC2'))
dev.off()

#pro-inflammatory gene
pro.inflammatory <- c("Il2", "Il12a", "Il12b", "Il17a", "Il18", "Ifng", "Tnf", "Il6", "Il1a", "Il1b", "Il36a", "Il36b", "Il36g")
anti.inflammatory <- c("Il1rn", "Il4", "Il10", "Il11", "Il13", "Tgfb1")

pdf("20220622-Comments/05.DBA/00.figures/52.Normal-DBA-all-celltype-pro-inflammatory.pdf", width = 35, height = 20)
VlnPlot(Normal.DBA.all, features = pro.inflammatory, split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 3, cols = c('#E87D72', '#56BCC2'))
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/52.Normal-DBA-all-celltype-anti-inflammatory.pdf", width = 20, height = 15)
VlnPlot(Normal.DBA.all, features = anti.inflammatory, split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
dev.off()

#> Comparing communications on a merged object
source.cells <- c("D3_Top2a", "D5_Gatm", "D6_S100a8", "iDSC0", "iDSC0-2")
target.cells <- c("Proliferating EC", "Angiogenic EC", "Arterial EC", "Venous EC1", "Venous EC2")
pathways <- c("VEGF", "ANGPTL", "CALCR")
gg6 <- netVisual_bubble(cellchat, sources.use = source.cells, targets.use = target.cells, comparison = c(1, 2), signaling = pathways, title.name = "Angiogenesis pathways", angle.x = 45, remove.isolate = F, thresh = 1)

#gg6 data
gg6.data <- gg6$data[, c(5, 8, 9, 12, 14, 15)]
gg6.data <- na.omit(gg6.data)
gg6.data$pathway_cellname_dataset <- paste(gg6.data$pathway_name, gg6.data$source.target, sep = ";")
gg6.data.agg <- aggregate(gg6.data$prob, by = list(gg6.data$pathway_cellname_dataset), max)
gg6.data.agg <- separate(gg6.data.agg, col = "Group.1", into = c("pathway", "source.target"), sep = ";")
gg6.data <- gg6.data[!gg6.data$interaction_name_2 %in% c("Vegfa  - Vegfr1", "Vegfa  - Vegfr1r2", "Vegfa  - Vegfr2", "Pdgfa  - Pdgfra", 
                                                         "Angptl4  - Sdc1", "Angptl4  - Cdh11", "Angptl4  - Sdc2", "Angptl2  - Tlr4", 
                                                         "Angptl4  - (Itgav+Itgb3)", "Angptl4  - Sdc4"), ]
pdf("20220622-Comments/05.DBA/00.figures/53.DBA-compare-normal-angiogenesis-change-pathway-plot.pdf", width = 18, height = 5)
ggplot() + 
  geom_point(gg6.data, mapping = aes(x = source.target, y = interaction_name_2, color = prob, size = prob)) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral"))) +
  scale_size_continuous(range = c(1,10)) + 
  theme_classic() + 
theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

#enhance VEGF expression
target.cells <- c("D3_Top2a", "D5_Gatm", "D6_S100a8", "iDSC0", "iDSC0-2")
source.cells <- c("iDSC0", "iDSC1", "iDSC2", "iDSC0-2")
pathways <- c("WNT", "PDGF")
gg7 <- netVisual_bubble(cellchat, sources.use = source.cells, targets.use = target.cells, comparison = c(1, 2), signaling = pathways, title.name = "Angiogenesis pathways", angle.x = 45, remove.isolate = F, thresh = 1)

#gg7 data
gg7.data <- gg7$data[, c(5, 8, 9, 12, 14, 15)]
gg7.data <- rbind(gg7.data, c(NA, "Pdgfc  - Pdgfra", "PDGF", "iDSC0 -> iDSC0-2 (Normal)", "iDSC0 -> iDSC0-2", "Normal"), 
                  c(NA, "Pdgfc  - Pdgfra", "PDGF", "iDSC1 -> iDSC0-2 (Normal)", "iDSC1 -> iDSC0-2", "Normal"), 
                  c(NA, "Pdgfc  - Pdgfra", "PDGF", "iDSC2 -> iDSC0-2 (Normal)", "iDSC2 -> iDSC0-2", "Normal"))
gg7.data <- gg7.data[!gg7.data$interaction_name_2 %in% c("Pdgfa  - Pdgfra", "Wnt10a  - (Fzd1+Lrp6)", "Wnt10a  - (Fzd2+Lrp6)"), ]
gg7.data$prob <- as.numeric(gg7.data$prob)
gg7.data$source.target <- factor(gg7.data$source.target, 
                                 levels = c("iDSC0 -> D3_Top2a (Normal)", "iDSC0 -> D3_Top2a (DBA)", "iDSC1 -> D3_Top2a (Normal)", "iDSC1 -> D3_Top2a (DBA)", 
                                            "iDSC2 -> D3_Top2a (Normal)", "iDSC2 -> D3_Top2a (DBA)", "iDSC0-2 -> D3_Top2a (Normal)", "iDSC0-2 -> D3_Top2a (DBA)", 
                                            "iDSC0 -> D5_Gatm (Normal)", "iDSC0 -> D5_Gatm (DBA)", "iDSC1 -> D5_Gatm (Normal)", "iDSC1 -> D5_Gatm (DBA)", 
                                            "iDSC2 -> D5_Gatm (Normal)", "iDSC2 -> D5_Gatm (DBA)", "iDSC0-2 -> D5_Gatm (Normal)", "iDSC0-2 -> D5_Gatm (DBA)", 
                                            "iDSC0 -> D6_S100a8 (Normal)", "iDSC0 -> D6_S100a8 (DBA)", "iDSC1 -> D6_S100a8 (Normal)", "iDSC1 -> D6_S100a8 (DBA)", 
                                            "iDSC2 -> D6_S100a8 (Normal)", "iDSC2 -> D6_S100a8 (DBA)", "iDSC0-2 -> D6_S100a8 (Normal)", "iDSC0-2 -> D6_S100a8 (DBA)",
                                            "iDSC0 -> iDSC0 (Normal)", "iDSC0 -> iDSC0 (DBA)", "iDSC1 -> iDSC0 (Normal)", "iDSC1 -> iDSC0 (DBA)", 
                                            "iDSC2 -> iDSC0 (Normal)", "iDSC2 -> iDSC0 (DBA)", "iDSC0-2 -> iDSC0 (Normal)", "iDSC0-2 -> iDSC0 (DBA)",
                                            "iDSC0 -> iDSC1 (Normal)", "iDSC0 -> iDSC1 (DBA)", "iDSC1 -> iDSC1 (Normal)", "iDSC1 -> iDSC1 (DBA)", 
                                            "iDSC2 -> iDSC1 (Normal)", "iDSC2 -> iDSC1 (DBA)", "iDSC0-2 -> iDSC1 (Normal)", "iDSC0-2 -> iDSC1 (DBA)",
                                            "iDSC0 -> iDSC2 (Normal)", "iDSC0 -> iDSC2 (DBA)", "iDSC1 -> iDSC2 (Normal)", "iDSC1 -> iDSC2 (DBA)", 
                                            "iDSC2 -> iDSC2 (Normal)", "iDSC2 -> iDSC2 (DBA)", "iDSC0-2 -> iDSC2 (Normal)", "iDSC0-2 -> iDSC2 (DBA)",
                                            "iDSC0 -> iDSC0-2 (Normal)", "iDSC0 -> iDSC0-2 (DBA)", "iDSC1 -> iDSC0-2 (Normal)", "iDSC1 -> iDSC0-2 (DBA)", 
                                            "iDSC2 -> iDSC0-2 (Normal)", "iDSC2 -> iDSC0-2 (DBA)", "iDSC0-2 -> iDSC0-2 (Normal)", "iDSC0-2 -> iDSC0-2 (DBA)"))
pdf("20220622-Comments/05.DBA/00.figures/54.DBA-compare-normal-pro-VEGF-change-pathway-plot.pdf", width = 15, height = 3.5)
ggplot() + 
  geom_point(gg7.data, mapping = aes(x = source.target, y = interaction_name_2, color = prob, size = prob)) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral"))) +
  scale_size_continuous(range = c(1,10)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

#plot gene expression of TNF, IL1, EGF, KIT, PROS, GAS, IGF pathways
plotGeneExpression(cellchat, signaling = "TNF", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "IL1", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "EGF", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "KIT", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "PROS", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "GAS", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "PDGF", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "ANGPTL", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "CALCR", split.by = "datasets", colors.ggplot = T)
Normal.DBA.all.s <- subset(Normal.DBA.all, cellcluster %in% c("D3_Top2a", "D5_Gatm", "D6_S100a8", "Epithelial", "Mac", "NK", "Neutrophil", "Proliferating Mac", "DC-1", 
                                                              "Proliferating EC", "Angiogenic EC", "Arterial EC", "Venous EC1", "Venous EC2", "EPC", "iDSC0", "iDSC1", "iDSC2"))
Normal.DBA.all.s$cellcluster <- factor(Normal.DBA.all.s$cellcluster, levels = c("D3_Top2a", "D5_Gatm", "D6_S100a8", "Epithelial", "Mac", "NK", "Neutrophil", "Proliferating Mac", "DC-1", 
                                                                                "Proliferating EC", "Angiogenic EC", "Arterial EC", "Venous EC1", "Venous EC2", "iDSC0", "iDSC1", "iDSC2", "EPC"))
#TNF pathway
pdf("20220622-Comments/05.DBA/00.figures/55.Normal-DBA-TNF-pathway-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.all.s, features = c("Tnf", "Tnfrsf1a", "Tnfrsf1b"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
dev.off()

#IL1 pathway
pdf("20220622-Comments/05.DBA/00.figures/55.Normal-DBA-IL1-pathway-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.all.s, features = c("Il1a", "Il1r2"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
dev.off()

#EGF pathway
pdf("20220622-Comments/05.DBA/00.figures/55.Normal-DBA-EGF-pathway-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.all.s, features = c("Hbegf", "Egfr", "Erbb2"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
dev.off()

#KIT pathway
pdf("20220622-Comments/05.DBA/00.figures/55.Normal-DBA-KIT-pathway-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.all.s, features = c("Kitl", "Kit"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
dev.off()

#PROS/GAS pathway
pdf("20220622-Comments/05.DBA/00.figures/55.Normal-DBA-PROS-GAS-pathway-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.all.s, features = c("Pros1", "Gas6", "Axl"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
dev.off()

#PDGF pathway
pdf("20220622-Comments/05.DBA/00.figures/55.Normal-DBA-PDGF-pathway-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.all.s, features = c("Pdgfa", "Pdgfc", "Pdgfra"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
dev.off()

#WNT pathway
pdf("20220622-Comments/05.DBA/00.figures/55.Normal-DBA-WNT-pathway-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.all.s, features = c("Wnt10a", "Wnt4", "Wnt6", "Fzd1", "Fzd2", "Fzd4", "Fzd6", "Lrp6"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
dev.off()

#CD137
pdf("20220622-Comments/05.DBA/00.figures/55.Normal-DBA-CD137-pathway-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.all.s, features = c("Tnfsf9", "Tnfrsf9"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
dev.off()

#VEGF
pdf("20220622-Comments/05.DBA/00.figures/55.Normal-DBA-VEGF-pathway-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.all.s, features = c("Vegfd", "Vegfb", "Pgf", "Flt4", "Kdr", "Flt1"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
dev.off()

#CALCR
pdf("20220622-Comments/05.DBA/00.figures/55.Normal-DBA-CALCR-pathway-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.all.s, features = c("Adm", "Calca", "Calcrl"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
dev.off()

#ANGPTL
pdf("20220622-Comments/05.DBA/00.figures/55.Normal-DBA-ANGPTL-pathway-violin.pdf", width = 10, height = 10)
VlnPlot(Normal.DBA.all.s, features = c("Angptl4", "Angptl2", "Cdh5", "Itga5", "Itgb1"), split.by = "sample", group.by = "cellcluster", 
        pt.size = 0, combine = T, split.plot = T, log = T, ncol = 2, cols = c('#E87D72', '#56BCC2'))
dev.off()


load("/Data3/dingli/RData/allGenes_WHY_hTBLC.data")
DimPlot(embryo.combined2)

embryo.combined2@active.ident <- factor(embryo.combined2$id)
DefaultAssay(embryo.combined2) <- "RNA"
embryo.combined2.markers <- FindAllMarkers(embryo.combined2, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
embryo.combined2 <- ScaleData(embryo.combined2, features = rownames(embryo.combined2))
embryo.combined2.markers.top5 <- embryo.combined2.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
pdf("test.pdf", width = 20, height = 20)
DoHeatmap(embryo.combined2, features = embryo.combined2.markers.top5$gene) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()

#20230227
pathways <- c("CXCL", "CCL")
source.cell <- c("iDSC0", "iDSC1", "iDSC2", "iDSC0-2")
target.cell <- c("Mac", "Proliferating Mac", "DC-1", "NK", "Neutrophil")
gg8 <- netVisual_bubble(cellchat, sources.use = source.cell, targets.use = target.cell, comparison = c(1, 2), signaling = pathways, title.name = "CXCL/CCL pathways", angle.x = 45, remove.isolate = F, thresh = 1)
gg8.data <- gg8$data[, c(5, 8, 9, 12, 14, 15)]
gg8.data.dcast <- dcast(gg8.data, interaction_name_2 ~ source.target, value.var = "prob", fill = 0)
rownames(gg8.data.dcast) <- gg8.data.dcast$interaction_name_2
gg8.data.dcast$interaction_name_2 <- NULL
pdf("20220622-Comments/05.DBA/00.figures/59.DBA-normal-iDSC-immune-cell-CXCL-CCL-pheatmap.pdf", width = 12, height = 8)
pheatmap(gg8.data.dcast, scale = "none", cluster_rows = T, cluster_cols = F, 
         color = colorRampPalette(colors = rev(brewer.pal(n = 10, name = "Spectral")))(100))
dev.off()

EC.sprouting <- c("VEGF", "CALCR", "ANGPTL", "ANGPT")
source.cell <- c("D6_S100a8", "iDSC0", "iDSC0-2")
target.cell <- c("Proliferating EC", "Angiogenic EC", "Arterial EC", "Venous EC1", "Venous EC2")
gg9 <- netVisual_bubble(cellchat, sources.use = source.cell, targets.use = target.cell, comparison = c(1, 2), signaling = EC.sprouting, title.name = "VEGF", angle.x = 45, remove.isolate = F, thresh = 1)
gg9.data <- gg9$data[, c(5, 8, 9, 12, 14, 15)]
gg9.data.dcast <- dcast(gg9.data, interaction_name_2 ~ source.target, value.var = "prob", fill = 0)
rownames(gg9.data.dcast) <- gg9.data.dcast$interaction_name_2
gg9.data.dcast$interaction_name_2 <- NULL
pdf("20220622-Comments/05.DBA/00.figures/59.DBA-normal-iDSC-EC-angiogenesis-pheatmap.pdf", width = 12, height = 8)
pheatmap(gg9.data.dcast, scale = "none", cluster_rows = T, cluster_cols = F, 
         color = colorRampPalette(colors = rev(brewer.pal(n = 10, name = "Spectral")))(100))
dev.off()





