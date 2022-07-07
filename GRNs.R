rm(list=ls())
args <- commandArgs()
# print(args[6])

source("twodatasets/muSuSiE-graph-two.R")
# source("pcNet.R")

dta_1 <- read.table(file = "./experiment/data/Curated/GSD/GSD-2000-1/ExpressionData.csv", 
                       header = T,
                       row.names = 1,
                       sep = ",",
                       as.is = TRUE)
dta_1 <- t(as.matrix(dta_1))
dta_2 <- read.table(file = "./experiment/data/Curated/GSD/GSD-2000-2/ExpressionData.csv", 
                    header = T,
                    row.names = 1,
                    sep = ",",
                    as.is = TRUE)
dta_2 <- t(as.matrix(dta_2))


out <- joint_graph_fun_two(dta_1, dta_2, 
                           scale_x = TRUE, 
                           intercept = TRUE,
                           sigma02_int = NULL, 
                           sigma2_int = NULL, 
                           prior_vec = NULL, #
                           itermax = 100, 
                           L_max = as.numeric(args[6]), #
                           tol = 1e-4, 
                           sigma0_low_bd = 1e-8,
                           residual_variance_lowerbound = NULL)
# names(out)
# print(sum(out$A_res_1==0)/length(as.vector(out$A_res_1))) # learned dropout (sparsity)
# print(sum(out$A_res_2==0)/length(as.vector(out$A_res_2)))
# apply(out$A_res_1, 2, function(x) sum(x!=0)) # effects for each col, change L_max?
# apply(out$A_res_2, 2, function(x) sum(x!=0))

# View(out$A_res_1)
out_1 <- out$A_res_1
out_2 <- out$A_res_2
out_1[abs(out_1) < 0.01] <- 0. # threshold based on prior
out_2[abs(out_2) < 0.01] <- 0.
# View(out_1)
# print(sum(out_1==0)/length(as.vector(out_1)))
# print(sum(out_2==0)/length(as.vector(out_2)))
# apply(out_1, 2, function(x) sum(x!=0))
# apply(out_2, 2, function(x) sum(x!=0))


write.table(out_1, file = "experiment/results/coeff_1.txt", sep = ",",
            row.names = F, col.names = F, quote = F)
write.table(out_2, file = "experiment/results/coeff_2.txt", sep = ",",
            row.names = F, col.names = F, quote = F)
# write.table(colnames(dta_1), file = "experiment/results/genes.txt", sep = "\t",
#             row.names = F, col.names = F, quote = F)

# # cov
# write.table(cov(dta_1), file = "experiment/results/coeff_1_cov.txt", sep = ",",
#             row.names = F, col.names = F, quote = F)
# write.table(cov(dta_2), file = "experiment/results/coeff_2_cov.txt", sep = ",",
#             row.names = F, col.names = F, quote = F)
# 
# #pcNet
# pcNet_1 <- pcNet(t(dta_1))
# pcNet_2 <- pcNet(t(dta_2))
# write.table(pcNet_1, file = "experiment/results/pcNet_1.txt", sep = ",",
#             row.names = F, col.names = F, quote = F)
# write.table(pcNet_2, file = "experiment/results/pcNet_2.txt", sep = ",",
#             row.names = F, col.names = F, quote = F)



# # prior
# p <- 19
# prior_vec <- c(1 / (2 * p^1.5), 1 / (p^2))
# prior_pi <- c(rep(prior_vec[1], 2 * p), rep(prior_vec[2], p))
# prior_pi <- c(prior_pi, 1 - sum(prior_pi))
# hist(prior_pi, breaks = p)
