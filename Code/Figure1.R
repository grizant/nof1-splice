## Create Figure 1
## AG Schissler
## Created 27 Jul 2017

############################################################
## i. Set up environment

library(ggplot2)
library(cowplot)
library(png)
library(grid)

setwd("~/Dropbox/splice-n-of-1-pathways/Figures")
load(file = "~/Dropbox/splice-n-of-1-pathways/Data/Figure1_dist_data.RData")
load(file = "~/Dropbox/splice-n-of-1-pathways/Data/Figure1_odds_ratio.RData")
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/lusc_iso_kegg_data.RData") ## NEW ISOFORM 25 Jul 2017

tmp_data <- scores_list[[1]]
(tmp_pat <- names(scores_list)[1])


############################################################
## 1. Panel A

## https://stackoverflow.com/questions/9917049/inserting-an-image-to-ggplot2

img <- readPNG(system.file("img", "Rlogo.png", package="png"))
g <- rasterGrob(img, interpolate=TRUE)

(panel_a <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point())

############################################################
## 2. Panel B

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

############################################################
## 3. Panel C

############################################################
## 4. Panel D

(panel_d <- ggplot(data = odds_ratio, aes(x = OR)) +
     geom_histogram(bins = 30, color = "white") +
     xlab("Odds ratio across pathways") + 
     ylab("Pathway count") +
     background_grid(major = "xy", minor = "none") +
     geom_vline(xintercept = tmp_locfdr$z.2[2], color = "red"))

############################################################
## 5. Panel E

############################################################
## 6. Combine and save

(fig1 <- ggdraw() +
     ## draw_plot_label(c("A", "B", "C", "D", "E"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15)
     draw_plot(panel_a, 0, 0.5, 0.33, .5) +
     draw_plot(panel_b, 0.33, 0.5, 0.33, .5) +
     draw_plot(panel_b, 0.66, 0.5, 0.33, .5) +
     draw_plot(panel_d, 0, 0, .5, .5) +
     draw_plot(panel_b, 0.5, 0, .5, .5))

## (fig1 <- plot_grid(NULL, panel_b, panel_c, labels = c("A", "B", "C"), ncol = 3))
## (fig1 <- plot_grid(NULL, panel_b, panel_c, labels = c("A", "B", "C"), ncol = 3))

save_plot("Figure1_v2.pdf", fig1,
          ncol = 3, 
          nrow = 2,
          ## base_width = 5,
          base_aspect_ratio = 1.1
          )
