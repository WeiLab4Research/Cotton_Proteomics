
########## Mendelian Randomization 所需函数 ##########

require(ggplot2)
require(ggsci)
#' Requires dev version of ggplot2
#' 
#' @param mr_results Output from \code{\link{mr}}.
#' @param dat Output from \code{\link{harmonise_data}}.
#' @export
#' @return List of plots
mr_scatter_plot <- function(mr_results, dat, xlab, ylab)
{
  # dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
  requireNamespace("ggplot2", quietly=TRUE)
  requireNamespace("plyr", quietly=TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE, key_glyph = "smooth") +
      ggsci::scale_color_nejm() +
      # scale_x_continuous(limits = c(min(xmin),max(xmax)),n.breaks = 5) +
      # scale_y_continuous(limits =c(min(ymin),max(ymax)), n.breaks = 5) +
      ggplot2::scale_colour_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF",
                                              "#6F99ADFF","#FFDC91FF","#EE4C97FF","black")) +
      ggplot2::labs(colour="MR Method", x = xlab, y = ylab) +
      ggplot2::theme_bw() + theme(legend.justification=c(0,0), legend.position=c(0.02,0.02), legend.direction="vertical", 
                                  legend.background = element_rect(color = "gray"),
                                  axis.line = element_line(colour = "black"),
                                  panel.border = element_blank(), panel.grid.major.x = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  text = element_text(size = 15)) + 
      ggplot2::guides(colour=ggplot2::guide_legend(ncol=1)) 
  })
  mrres
}
