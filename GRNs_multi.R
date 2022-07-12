rm(list=ls())
args <- commandArgs()
# print(args[6])

source("mutipledatasets/muSuSiE-graph-multi.R")

# load data
load_data <- function(path_in) {
  dta <- read.table(file = path_in, 
                      header = T,
                      row.names = 1,
                      sep = ",",
                      as.is = TRUE)
  dta <- t(as.matrix(dta))
  return(dta)
}

files <- c("experiment/data/Curated/GSD/GSD-2000-1/ExpressionData.csv",
           "experiment/data/Curated/GSD/GSD-2000-2/ExpressionData.csv",
           "experiment/data/Curated/GSD/GSD-2000-3/ExpressionData.csv",
           "experiment/data/Curated/GSD/GSD-2000-4/ExpressionData.csv",
           "experiment/data/Curated/GSD/GSD-2000-5/ExpressionData.csv",
           "experiment/data/Curated/GSD/GSD-2000-6/ExpressionData.csv",
           "experiment/data/Curated/GSD/GSD-2000-7/ExpressionData.csv",
           "experiment/data/Curated/GSD/GSD-2000-8/ExpressionData.csv",
           "experiment/data/Curated/GSD/GSD-2000-9/ExpressionData.csv",
           "experiment/data/Curated/GSD/GSD-2000-10/ExpressionData.csv")
dta_list <- list()
for (i in seq_along(files)) {
  dta_list[[i]] <- load_data(files[i])
}

# muSuSiE
out <- joint_graph_multi(dta_list, 
                         scale_x = TRUE, 
                         intercept = TRUE,
                         sigma02_int = NULL, 
                         sigma2_int = NULL, 
                         prior_vec = NULL, #
                         itermax = 100, 
                         L_max = 10, #as.numeric(args[6]), #
                         tol = 1e-4, 
                         sigma0_low_bd = 1e-8,
                         residual_variance_lowerbound = NULL)
names(out)

# output
out_coeff <- list()
for (i in seq_along(files)) {
  out_coeff[[i]] <- out$A_list[[i]]
  out_coeff[[i]][abs(out_coeff[[i]]) < 0.01] <- 0.
  print(paste0("dropout:", sum(out_coeff[[i]]==0)/length(as.vector(out_coeff[[i]]))))
  print(paste0("non-zeros:", apply(out_coeff[[i]], 2, function(x) sum(x!=0))))
}

files <- c("coeff_1.txt",
           "coeff_2.txt",
           "coeff_3.txt",
           "coeff_4.txt",
           "coeff_5.txt",
           "coeff_6.txt",
           "coeff_7.txt",
           "coeff_8.txt",
           "coeff_9.txt",
           "coeff_10.txt")

save_coeff <- function(out, path_out) {
  write.table(out, 
              file = paste0("experiment/results/multi_coeff/", path_out),
              sep = ",",
              row.names = F, 
              col.names = F, 
              quote = F)
}
for (i in seq_along(files)) {
  save_coeff(out_coeff[[i]], files[i])
}




