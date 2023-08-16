#Author: Min Yang

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(NMF)
library(circlize)

options(stringsAsFactors = FALSE)

#Normal samples
#1. Merge all the DSCs/Immune cells/ECs/iDSCs subclusters and Epithelial cells in one object, named as: Normal.cellchat
Normal.cellchat$active.ident <- Normal.cellchat$subname
Normal.cellchat@active.ident <- factor(Normal.cellchat$subname)
normal.data.input <- Normal.cellchat@assays$RNA@data
normal.meta <- Normal.cellchat@meta.data

Normal.cellchat <- createCellChat(object = Normal.cellchat, meta = normal.meta, group.by = "active.ident")
#set default idents
Normal.cellchat <- setIdent(Normal.cellchat, ident.use = "active.ident")
groupSize <- as.numeric(table(Normal.cellchat@idents))
#set the ligand-receptor interaction database
Normal.cellchatDB <- CellChatDB.mouse
Normal.cellchatDB.use <- subsetDB(Normal.cellchatDB, search = "Secreted Signaling")
Normal.cellchat@DB <- Normal.cellchatDB.use
#processing the expression data for cell-cell communication
Normal.cellchat <- subsetData(Normal.cellchat)
future::plan("multiprocess", workers = 10)#do parallel
Normal.cellchat <- identifyOverExpressedGenes(Normal.cellchat)
Normal.cellchat <- identifyOverExpressedInteractions(Normal.cellchat)
Normal.cellchat <- projectData(Normal.cellchat, PPI.mouse)
#inference of cell-cell communication network
#1) compute the communication probability and infer the cellular communication network
Normal.cellchat <- computeCommunProb(Normal.cellchat)
#2) infer the cell-cell communication at a signaling pathway level
Normal.cellchat <- computeCommunProbPathway(Normal.cellchat)
#3) calculate the aggregated cell-cell communication network
Normal.cellchat <- aggregateNet(Normal.cellchat)
save(Normal.cellchat, file = "Normal.cellchat.RData")

#2. The same as CBA sample and saved as: CBA.cellchat

#3. compare Normal and CBA cellchat results
CBA.cellchat <- updateCellChat(CBA.cellchat)
Normal.cellchat <- updateCellChat(Normal.cellchat)
#extract the cell identity, make sure these two datasets are the same
group.new1 <- c(levels(CBA.cellchat@idents), setdiff(levels(Normal.cellchat@idents), levels(CBA.cellchat@idents)))
group.new2 <- c(levels(Normal.cellchat@idents), setdiff(levels(CBA.cellchat@idents), levels(Normal.cellchat@idents)))
Normal.cellchat <- liftCellChat(Normal.cellchat, group.new = group.new2)
CBA.cellchat <- liftCellChat(CBA.cellchat, group.new = group.new1)
#Check the level of cell identity, and they need to have the same level
object.list <- list(Normal = Normal.cellchat, DBA = CBA.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

#Get some values
df.net <- subsetCommunication(cellchat, slot.name = 'netP')
#Example: VEGF strength comparing
gg1 <- netVisual_bubble(cellchat, sources.use = source.cell, targets.use = target.cell, comparison = c(1, 2), max.dataset = 2, signaling = "VEGF", title.name = "Increased signaling in CBA", angle.x = 45, remove.isolate = T)
gg1.data <- gg1.data[gg1.data$pathway_name == "VEGF", ]
st.levels <- levels(gg1.data$source.target)
gg1.data$prob <- as.numeric(gg1.data$prob)
gg1.data.dcast <- dcast(gg1.data, interaction_name_2 ~ source.target, value.var = "prob", fill = 0)
rownames(gg1.data.dcast) <- gg1.data.dcast$interaction_name_2
gg1.data.dcast$interaction_name_2 <- NULL
pheatmap(gg1.data.dcast, scale = "none", cluster_rows = T, cluster_cols = F, 
         color = colorRampPalette(colors = rev(brewer.pal(n = 10, name = "Spectral")))(100))

#Chord plot
chordDiagram(gg1.data.dcast, 
             order = c(cell.order), 
             grid.col = grid.col, directional = 1, 
             direction.type = c("diffHeight", "arrows"), 
             link.arr.type = "big.arrow", 
             annotationTrack = c("name", "grid"), 
             scale = "TRUE", 
             link.visible = gg1.data.dcast[[3]] >= 0.01) 


