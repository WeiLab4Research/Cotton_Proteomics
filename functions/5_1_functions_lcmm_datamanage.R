
########## Latent Class Mixed Model Data Management所需函数 ##########

latentclass_data <- function(dat, latentclass) {
  
  cotton_class <- latentclass$pprob[,c(1,2)]
  cotton_class_df <- left_join(lcmmdata, cotton_class)
  
  cotton_class_sum <- cotton_class %>%
    drop_na() %>%
    group_by(class) %>%
    summarise(n = n()) %>%
    mutate(freq = round(n/sum(n)*100,2))%>%
    as.data.frame()
  
  leng <- paste0(cotton_class_sum$class,
                 " (", "N = ", cotton_class_sum$n,
                 ", ", round_tidy(cotton_class_sum$freq,2), "%",
                 ")")
  
  latentclass_data <- list(cotton_class_df = cotton_class_df,
                           cotton_class_sum = cotton_class_sum,
                           leng = leng)
  return(latentclass_data)
}




