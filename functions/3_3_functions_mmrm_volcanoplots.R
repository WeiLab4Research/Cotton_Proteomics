
########## Volcano Plots 所需函数 ##########

volcano_plot <- function (results, inter = F, titles,
                          x.limits = c(-120,120),
                          x.breaks = seq(-120, 120, 20),
                          y.limits = c(0,7),
                          y.breaks = seq(0, 7, 1),
                          FDR = F,
                          gene_name = F) {
  
  if (inter) {
    volres <- ggplot(data = results) +
      geom_point(aes(x=beta, y= -log10(Pval), color = change, alpha = 0.4)) +
      # geom_hline(yintercept = -log10(0.05), linetype = 2) +
      labs(x="beta", y="-log10(P)") +
      scale_color_manual(values=c("Nosig" = "#d2dae2",
                                  "Up" = "#BC3C29FF",
                                  "Down" = "#0072B5FF")) +
      scale_x_continuous(limits = c(-3,3),
                         breaks = seq(-3, 3, 1),
                         expand = c(0,0)) +
      scale_y_continuous(limits = c(0,7),
                         breaks = seq(0, 7, 1),
                         expand = c(0,0)) +
      geom_text_repel(data = results[which(results$FDR < 0.05),],
                      aes(x= beta, y= -log10(Pval), label = gene.name),
                      min.segment.length = 0.05,
                      segment.alpha = 0.6,
                      max.overlaps = 30) +
      # nudge_x = 0.2, 
      # nudge_y = 0.2) +
      theme_classic() +
      theme(legend.position = "none",
            axis.title = element_text(face = "bold",
                                      size = 15),
            axis.text = element_text(size = 12)) +
      ggtitle(titles)
  } else {
    volres <- ggplot(data = results) +
      geom_point(aes(x=beta, y= -log10(Pval), color = change, alpha = 0.4)) +
      # geom_hline(yintercept = -log10(0.05), linetype = 2) +
      labs(x="beta", y="-log10(P)") +
      scale_color_manual(values=c("Nosig" = "#d2dae2",
                                  "Up" = "#BC3C29FF",
                                  "Down" = "#0072B5FF")) +
      scale_x_continuous(limits = x.limits,
                         breaks = x.breaks,
                         expand = c(0,0)) +
      scale_y_continuous(limits = y.limits,
                         breaks = y.breaks,
                         expand = c(0,0)) +
      # nudge_x = 0.2, 
      # nudge_y = 0.2) +
      theme_classic() +
      theme(legend.position = "none",
            axis.title = element_text(face = "bold",
                                      size = 15),
            axis.text = element_text(size = 12)) +
      ggtitle(titles)
    
  }
  
  if(gene_name) {
    require(ggrepel)
    volres <- volres +
      geom_text_repel(data = results,
                      aes(x= beta, y= -log10(Pval), label = gene.name),
                      min.segment.length = 0.05,
                      segment.alpha = 0.6,
                      max.overlaps = 10000)
  }
  
  FDR_lim <- 0
  if(min(results$FDR) < 0.05) {
    FDR_lim <- mean(c(results[max(which(results$FDR < 0.05)), "Pval"], results[min(which(results$FDR > 0.05)), "Pval"]))
  }
  
  if(FDR) {
    volres <- volres + 
      geom_point(data = results[which(results$FDR < 0.05),], aes(x=beta, y= -log10(Pval), color = change)) +
      geom_hline(yintercept = -log10(FDR_lim), linetype = 2)
  }
  
    volres <- volres + 
      geom_hline(yintercept = -log10(0.05), linetype = 2)
  
  return(volres)
}
