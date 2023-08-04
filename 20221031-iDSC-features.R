load("20210827-figures/files/UE-immune-harmony-web.RData")
load("20210827-figures/files/UE-DSC-web.RData")
load("20210827-figures/figure5/FB-immune-sc-subclusters-SCT.RData")
decidual.color <- c("#6A5ACD","#9ACD32","#FB8072","#8DD3C7","#08519C","#FDB462","#228B22","#FCCDE5","#D9D9D9")
immune.color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#2E8B57", "#CAB2D6", "#BDB76B", "#FF7F00", "#6A3D9A")

decidual.levels <- c("Lum-D", "Sfrp4-D", "Top2a-D", "Ptn-D", "Gatm-D", "S100a8-D", "Prl8a2-D", "cluster9-1", "Prg-")
immune.levels <- c("Mono", "Mac", "Pro-Mac", "DC-1", "DC-2", "Neutrophil", "NK", "B", "T", "FB-immune")

DimPlot(UE.decidual.susbet, label = T, group.by = "sub_clusters")
DimPlot(UE.immune.harmony, label = T, group.by = "sub_clusters")
UE.immune.harmony <- subset(UE.immune.harmony, sub_clusters != "iDSC")
FB.immune.SCT$sub_clusters <- paste("iDSC_", FB.immune.SCT$seurat_clusters, sep = "")
DimPlot(FB.immune.SCT, label = T, group.by = "sub_clusters")
Normal.DSC.immune.iDSC <- merge(UE.decidual.susbet, y = list(UE.immune.harmony, FB.immune.SCT))
Normal.DSC.immune.iDSC <- NormalizeData(Normal.DSC.immune.iDSC)
save(Normal.DSC.immune.iDSC, file = "20220622-Comments/07.iDSC_features/Normal-DSC-immue-iDSC-merge-20210104.RData")
#add spliced and unspliced reads for immune and iDSC cells in velocity.R file.
load("20220622-Comments/07.iDSC_features/Normal-DSC-immue-iDSC-merge-splicing-20210104.RData")
Normal.DSC.immune.iDSC$sub_clusters <- factor(Normal.DSC.immune.iDSC$sub_clusters, 
                                              levels = c("D1_eSF", "D2_Pre-DSC", "D3_Proliferating DSC", "D4_Biomineral-regulatory DSC", "D5_Nourishing DSC", "D6_Angiogenesis DSC", "D7_Postmature DSC", 
                                                         "Mono", "Mac", "Cycling Mac", "DC-1", "DC-2", "Neutrophil", "NK", "T", "B", "iDSC_0", "iDSC_1", "iDSC_2"))
phagocytocis.gene <- read.table("20220622-Comments/07.iDSC_features/GO_0006909_phagocytosis.csv", sep = ",", header = T, check.names = F)
Normal.DSC.immune.iDSC <- AddModuleScore(Normal.DSC.immune.iDSC, features = list(phagocytocis.gene$Symbol), name = "phago_score")
Normal.DSC.immune.iDSC.color <- c(decidual.color[1:7], immune.color[c(1:10,10,10)])
Normal.DSC.immune.iDSC$subcluster <- gsub("iDSC_.*", "iDSC", Normal.DSC.immune.iDSC$sub_clusters)
Normal.DSC.immune.iDSC$subcluster <- factor(Normal.DSC.immune.iDSC$subcluster, 
                                            levels = c("D1_eSF", "D2_Pre-DSC", "D3_Proliferating DSC", "D4_Biomineral-regulatory DSC", "D5_Nourishing DSC", "D6_Angiogenesis DSC", "D7_Postmature DSC", 
                                                       "Mono", "Mac", "Cycling Mac", "DC-1", "DC-2", "Neutrophil", "NK", "T", "B", "iDSC"))
pdf("20220622-Comments/07.iDSC_features/DSC-iDSC-phagocytosis-transcripts-violin.pdf", width = 10, height = 4)
VlnPlot(Normal.DSC.immune.iDSC, features = c("phago_score1", "nFeature_RNA", "nCount_RNA"), pt.size = 0, group.by = "subcluster", cols = Normal.DSC.immune.iDSC.color)
dev.off()

UE.harmony.pc30.markers <- read.table("20210827-figures/files/UE-harmony-major-DEGs.csv", sep = ",", header = T, check.names = F)
DSC.genes <- UE.harmony.pc30.markers[UE.harmony.pc30.markers$cluster == "Decidual", ]$gene
immune.genes <- UE.harmony.pc30.markers[UE.harmony.pc30.markers$cluster == "Immune", ]$gene

Normal.DSC.immune.iDSC.unspliced.spliced <- as.matrix(Normal.DSC.immune.iDSC@assays$unspliced.spliced@data)
Normal.DSC.immune.iDSC.splicing <- data.frame(cellname = colnames(Normal.DSC.immune.iDSC), 
                                              celltype = Normal.DSC.immune.iDSC$sub_clusters)
Normal.DSC.immune.iDSC.unspliced.spliced.df <- data.frame(t(Normal.DSC.immune.iDSC.unspliced.spliced), check.names = F)
Normal.DSC.immune.iDSC.unspliced.spliced.df$cellname <- rownames(Normal.DSC.immune.iDSC.unspliced.spliced.df)
Normal.DSC.immune.iDSC.splicing <- merge(Normal.DSC.immune.iDSC.splicing, Normal.DSC.immune.iDSC.unspliced.spliced.df, by = "cellname", all = F)
Normal.DSC.immune.iDSC.splicing$celltype <- gsub("iDSC_.*", "iDSC", Normal.DSC.immune.iDSC.splicing$celltype)
Normal.DSC.immune.iDSC.splicing$celltype <- factor(Normal.DSC.immune.iDSC.splicing$celltype, 
                                                   levels = c("D1_eSF", "D2_Pre-DSC", "D3_Proliferating DSC", "D4_Biomineral-regulatory DSC", "D5_Nourishing DSC", "D6_Angiogenesis DSC", "D7_Postmature DSC", 
                                                              "Mono", "Mac", "Cycling Mac", "DC-1", "DC-2", "Neutrophil", "NK", "T", "B", "iDSC"))
Normal.DSC.immune.iDSC.splicing.celltype <- aggregate(Normal.DSC.immune.iDSC.splicing[, -c(1:2)], by = list(Normal.DSC.immune.iDSC.splicing$celltype), mean)
rownames(Normal.DSC.immune.iDSC.splicing.celltype) <- Normal.DSC.immune.iDSC.splicing.celltype$Group.1
Normal.DSC.immune.iDSC.splicing.celltype$Group.1 <- NULL
Normal.DSC.immune.iDSC.splicing.celltype.t <- data.frame(t(Normal.DSC.immune.iDSC.splicing.celltype), check.names = F)
Normal.DSC.immune.iDSC.splicing.celltype.DSC.markers <- Normal.DSC.immune.iDSC.splicing.celltype.t[rownames(Normal.DSC.immune.iDSC.splicing.celltype.t) %in% DSC.genes, ]
Normal.DSC.immune.iDSC.splicing.celltype.DSC.markers$gene_type <- "DSC"
Normal.DSC.immune.iDSC.splicing.celltype.immune.markers <- Normal.DSC.immune.iDSC.splicing.celltype.t[rownames(Normal.DSC.immune.iDSC.splicing.celltype.t) %in% immune.genes, ]
Normal.DSC.immune.iDSC.splicing.celltype.immune.markers$gene_type <- "Immune"

Normal.DSC.immune.iDSC.splicing.celltype.markers <- rbind(Normal.DSC.immune.iDSC.splicing.celltype.DSC.markers, Normal.DSC.immune.iDSC.splicing.celltype.immune.markers)
Normal.DSC.immune.iDSC.splicing.celltype.markers.m <- melt(Normal.DSC.immune.iDSC.splicing.celltype.markers, id.vars = "gene_type")

pdf("20220622-Comments/07.iDSC_features/DSC-iDSC-unspliced-spliced-boxplot.pdf", width = 10, height = 4)
ggplot() + 
  geom_boxplot(Normal.DSC.immune.iDSC.splicing.celltype.markers.m, mapping = aes(x = variable, y = value, fill = gene_type), outlier.shape = NA) + 
  theme_classic()
dev.off()



