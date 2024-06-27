
########## Latent Class Mixed Model Class Plots所需函数 ##########

predY_plot <- function(latentclass, y, var.time, dat, bty = "l") {
  newdata <- dat
  colnames(newdata)[which(colnames(newdata) == y)] <- "variable_y"
  colnames(newdata)[which(colnames(newdata) == var.time)] <- "variable_x"
  
  plot(predictY(latentclass, var.time = "variable_x", newdata = newdata, bty = bty))
  
}


class_plot <- function(dat, leng, nclass = 4, xvar, yvar, xlab, ylab, xlim = NULL, ylim = NULL) {
  
  total_dat <- dat
  colnames(total_dat)[which(colnames(total_dat) == xvar)] <- "variable_x"
  colnames(total_dat)[which(colnames(total_dat) == yvar)] <- "variable_y"
  
  if(is.null(xlim)) {
    xlim_min <- floor(min(total_dat$variable_x) / 5) * 5
    xlim_max <- ceiling(max(total_dat$variable_x) / 5) * 5
    xlim <- c(xlim_min,xlim_max)
  }
  
  if(is.null(ylim)) {
    ylim_min <- floor(min(total_dat$variable_y) / 1000) * 1000
    ylim_max <- ceiling(max(total_dat$variable_y) / 1000) * 1000
    ylim <- c(ylim_min,ylim_max)
  }
  
  class_plots <- list()
  patch_class_plots <- list()
  all_class_plot <- ggplot()
  for (i in 1:nclass) {
    
    class_data <- total_dat %>%
      filter(class == i)
    
    assign(paste0("plot_class",i),
           ggplot()+
             geom_line(class_data, mapping = aes(x = variable_x, y = variable_y, group = as.character(pid)),alpha = 0.3, color = "gray40")+
             geom_smooth(data = class_data, mapping = aes(x = variable_x, y = variable_y),
                         color =  color_lcmm[i], fill= color_lcmm[i],
                         se = T,
                         alpha = 0.3,
                         lwd = 1.2)+
             scale_y_continuous(expand = c(0,0), limits = ylim, breaks = seq(ylim[1],ylim[2],by = 1000))+
             scale_x_continuous(expand = c(0,0), limits = xlim, breaks = seq(xlim[1],xlim[2],by = 10))+
             # geom_hline(yintercept = 100, lty = "dashed", color = "gray60", size = 1)+
             theme_classic()+
             theme(panel.border = element_rect(fill = NA), 
                   axis.line = element_line(),
                   axis.title = element_text(size = 14),
                   axis.text = element_text(size = 12))+
             theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = color_lcmm[i]),
                   legend.position = "top", 
                   legend.title=element_blank(), 
                   legend.text = element_text(size = 14, color = "gray50"))+
             labs(x = xlab, y = ylab))
    
    assign(paste0("plot_class",i,"_title"),
           get((paste0("plot_class",i)))+
             ggtitle(leng[i])+
             theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = color_lcmm[i])))
    
    all_class_plot <- all_class_plot + 
      geom_smooth(data = class_data, mapping = aes(x = variable_x, y = variable_y),
                  color =  color_lcmm[i], fill= color_lcmm[i],
                  se = T,
                  alpha = 0.3,
                  lwd = 1.2)
    
    class_plots[[i]] <- get(paste0("plot_class",i,"_title"))
    
    if(i == 1) {
      patch_class_plots <- get(paste0("plot_class",i,"_title"))
    } else {
      patch_class_plots <- patch_class_plots + get(paste0("plot_class",i,"_title"))
    }
  }
  
  all_class_plot <- all_class_plot + 
    scale_y_continuous(expand = c(0,0), limits = ylim, breaks = seq(ylim[1],ylim[2],by = 1000))+
    scale_x_continuous(expand = c(0,0), limits = xlim, breaks = seq(xlim[1],xlim[2],by = 10))+
    # geom_hline(yintercept = 100, lty = "dashed", color = "gray60", size = 1)+
    theme_classic()+
    theme(panel.border = element_rect(fill = NA), 
          axis.line = element_line(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))+
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = color_lcmm[i]),
          legend.position = "top", 
          legend.title=element_blank(), 
          legend.text = element_text(size = 14, color = "gray50"))+
    labs(x = xlab, y = ylab)
  
  class_plots[[nclass + 1]] <- patch_class_plots
  class_plots[[nclass + 2]] <- all_class_plot
  
  return(class_plots)
}










