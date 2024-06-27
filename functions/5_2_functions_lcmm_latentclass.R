
########## Latent Class Mixed Model latent_class分组所需函数 ##########

fit_latent_class <- function(lcmmdata, y, x, cov = NULL, 
                             random = 1,
                             link = "linear",
                             subject = "pid",
                             max_ng = 8,
                             maxiter = 2000,
                             nproc = 5) {
  temp_lcmmdata <<- lcmmdata
  colnames(temp_lcmmdata)[colnames(temp_lcmmdata) == y] <- "variable_y"
  colnames(temp_lcmmdata)[colnames(temp_lcmmdata) == x] <- "variable_x"
  if(is.null(cov)) {## 需要补充存在covariant情况
    m1 <<- lcmm(variable_y ~ variable_x,
               random = ~ 1,
               link = link,
               subject = subject,
               # classmb = ~ male81 + height,
               data = temp_lcmmdata,
               ng = 1,
               maxiter = maxiter,
               nproc = nproc)
    
    for (i in 2:max_ng) {
      cat("[",format(Sys.time(),"%H:%M:%S"),"]","........ i =", i, " \n")
      assign(paste0("m",i),
             lcmm(variable_y ~ variable_x,
                  random = ~ 1,## 需要调整
                  mixture = ~ variable_x,## 需要调整
                  link = link,
                  subject = subject,
                  classmb = ~ male81 + height,
                  data = temp_lcmmdata,
                  ng = i,
                  # 2-max类潜变量搜索最佳参数时以m1拟合的最佳参数为初始值
                  B = m1,
                  maxiter = maxiter,
                  nproc = nproc))
    }
  }
  
  sumtable <- summarytable(m1,m2,m3,m4,m5,m6,m7,m8,## 需要调整
                           which = c("npm","loglik","AIC","BIC", "SABIC","entropy", "%class","conv"))%>%
    as_tibble() %>%
    mutate(class = c(1:max_ng))
  
  latent_class_list <- list()
  latent_class_list[["sumtable"]] <- sumtable
  for (j in 1:max_ng) {
    latent_class_list[[paste0("m",j)]] <- get(paste0("m",j))
  }
  
  return(latent_class_list)
}
