## Create Table 2
## AG Schissler
## Created 27 Jul 2017

############################################################
## i. Set up environment

library(xtable)
source("~/Dropbox/Splice-n-of-1-pathways/Code/splice_functions.R")

count <- matrix(c(1, 4, 20, 0, 17, 24, 10, 2, 1, 5), nrow = 5, ncol = 2)
(rel1 <- c(1, 4, 20, 0, 17)/sum(c(1, 4, 20, 0, 17)))
(rel2 <- c(24, 10, 2, 1, 5)/sum(c(24, 10, 2, 1, 5)))
(rel <- matrix(c(rel1, rel2), nrow = 5, ncol = 2))
    
compute_hellinger(X = rel)
