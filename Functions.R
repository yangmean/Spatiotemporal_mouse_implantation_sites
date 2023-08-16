#Author: Min Yang

Sptial_plot <- function(object, pt.size = 0.01, pt.shape = 16, features = NULL, min.cutoff = NA, max.cuttoff = NA,
                            slot = "data", group.by = Idents(object), ncol = 3, shape = 16,
                            color = rev(brewer.pal(n = 11, name = "Spectral"))){
  DefaultAssay(object) <- "SCT"
  image <- GetTissueCoordinates(object)
  image <- image[, c("x", "y", "cells")]
  imgId <- data.frame(image, group.by = group.by)
  data <- imgId
  plot <- ggplot(data, aes(x = x, y = y, color = group.by)) + 
    geom_point(size = pt.size, shape = pt.shape) + 
    theme(panel.grid = element_blank(), 
          panel.background = element_rect(color = "black", fill = "white"))
  if (is.null(features)) {
    return(plot)
  }
  else {
    features.data <- FetchData(object = object, vars = features, slot = slot)
    if (length(features) == 1){
      colnames(features.data) <- features
    }
    else { 
      features <- colnames(features.data)
    }
    plots <- list()
    for (i in 1:length(features)) {
      data <- cbind(imgId, features.data[, i])
      colnames(data)[5] <- "value" #important
      plot <- ggplot(data, aes(x = x, y = y, color = value)) + 
        geom_point(size = pt.size, shape = pt.shape) + 
        scale_color_gradientn(colours = color) + 
        theme(panel.grid = element_blank(), 
              panel.background = element_rect(color = "black", fill = "white"), 
              axis.text = element_blank()) + 
        ggtitle(features[i])
      plots[[i]] <- plot
    }
    return(wrap_plots(plots, ncol = ncol))
  }
}

enrich.GO.function <- function(geneset){
  enrich.GO <- enrichGO(gene = geneset,
                        OrgDb = "org.Mm.eg.db",
                        keyType = "SYMBOL",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)
  enrich.GO <- simplify(enrich.GO, cutoff = 0.7, by = "p.adjust", select_fun = min)
  enrich.GO <- enrich.GO@result[enrich.GO@result$pvalue<0.01, ]
  return(enrich.GO)
}