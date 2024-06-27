# 2.1 p值函数
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    #      if (shapiro.test(y[g==levels(g)[1]])$p.value>=0.05){
    p <- round(t.test(y ~ g)$p.value,3)
    #      }else {
    #          p <- round(wilcox.test(y ~ g,alternative = "two.sided")$p.value,3)
    #      }
  } else if (is.ordered(x[[1]])){
    # for level variables, perform wilcoxn秩和检验
    p <- round(wilcox.test(as.numeric(y) ~ g,alternative = "two.sided")$p.value,3)
  }
  else {
    # For categorical variables, perform a chi-squared test of independence
    p <- tryCatch({round(chisq.test(table(y, g),correct=F)$p.value,3)},
                  warning = function(w){round(fisher.test(table(y, g))$p.value,3)},
                  finally = {round(chisq.test(table(y, g),correct=F)$p.value,3)})
    # Format the p-value, using an HTML entity for the less-than sign.
    # The initial empty string places the output on the line below the variable label.
    #c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
  }
  p<-format.pval(p, eps=0.001,digits=3)
}
# 2.2 统计量函数
statistic <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    #      if (shapiro.test(y[g==levels(g)[1]])$p.value>=0.05){
    statistic <- t.test(y ~ g)$statistic
    statistic <- paste0("t = ",round(statistic,3))
    #      }else {
    #          statistic <- wilcox.test(y ~ g,alternative = "two.sided")$statistic
    #          statistic <- paste0("W = ",statistic)
    #      }
  } else if (is.ordered(x[[1]])){
    # for level variables, perform wilcoxn秩和检验
    statistic <- wilcox.test(as.numeric(y) ~ g,alternative = "two.sided")$statistic
    statistic <- paste0("W = ",statistic)
  }
  else {
    # For categorical variables, perform a chi-squared test of independence
    statistic <- tryCatch({paste0("chisq = ",round(chisq.test(table(y, g),correct=F)$statistic,3))},warning = function(w){paste0("*")},finally = {paste0("chisq = ",round(chisq.test(table(y, g),correct=F)$statistic,3))})
  }}

# 2.3 格式函数
my.render.cont<-function(x) {
  with(stats.default(x), c("",
                           "Mean ± SD"=sprintf("%.4f ± %.4f", MEAN, SD),"Median (IQR)"=sprintf("%.4f (%.4f)", MEDIAN, IQR)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.1f %%)", FREQ, PCTnoNA))))
}
