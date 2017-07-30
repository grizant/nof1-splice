## Create Figure 1
## AG Schissler
## Created 27 Jul 2017

############################################################
## i. Set up environment

library(ggplot2)
library(cowplot)

setwd("~/Dropbox/splice-n-of-1-pathways/Figures")
load(file = "~/Dropbox/splice-n-of-1-pathways/Data/Figure1_dist_data.RData")
load(file = "~/Dropbox/splice-n-of-1-pathways/Data/Figure1_odds_ratio.RData")

(panel_b <- ggplot(data = dist_data, aes(x = gene_dist, fill = factor(call))) +
     geom_histogram(bins = 60) +
     xlab("Hellinger distance") + 
     ylab("Gene count") +
     background_grid(major = "xy", minor = "none") +
     ## theme_bw() +
     theme(legend.position = "none"))

odds_ratio <- data.frame(OR = odds_ratio)
fil_odds <- odds_ratio[!is.na(odds_ratio)]
tmp_locfdr <- locfdr::locfdr(zz = fil_odds, bre = ceiling(length(fil_odds)/8), df = 4,
                             pct = 0, pct0 = 1/4, nulltype = 2, plot = 1, mlests = c(1, sd(fil_odds)))

## could try to save as a ggplot object or recreate later
## https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/

(panel_c <- ggplot(data = odds_ratio, aes(x = OR)) +
     geom_histogram(bins = 30, color = "white") +
     xlab("Odds ratio across pathways") + 
     ylab("Pathway count") +
     background_grid(major = "xy", minor = "none") +
     geom_vline(xintercept = tmp_locfdr$z.2[2]))
     
     

(fig1 <- plot_grid(NULL, panel_b, panel_c, labels = c("A", "B", "C"), ncol = 3))

save_plot("Figure1_v1.pdf", fig1,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.1,
          base_width = 5
          )
