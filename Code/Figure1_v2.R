## Create Figure 1
## AG Schissler
## Created 27 Jul 2017

############################################################
## i. Set up environment

library(ggplot2)
library(cowplot)
library(xtable)

setwd("~/Dropbox/splice-n-of-1-pathways/Figures")
load(file = "~/Dropbox/splice-n-of-1-pathways/Data/Figure1_dist_data.RData")
load(file = "~/Dropbox/splice-n-of-1-pathways/Data/Figure1_odds_ratio.RData")
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/lusc_iso_kegg_data.RData") ## NEW ISOFORM 25 Jul 2017
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData") ## NEW ISOFORM 25 Jul 2017

tmp_data <- scores_list[[1]]
(tmp_pat <- names(scores_list)[1])


############################################################
## 1. Panel A: individual with Hellinger distance equation

## https://stackoverflow.com/questions/9917049/inserting-an-image-to-ggplot2

img <- readPNG(system.file("img", "Rlogo.png", package="png"))
g <- rasterGrob(img, interpolate=TRUE)

(panel_a <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point())

## panel_a will be removed by hand

############################################################
## 2. Panel B this will stay

(panel_b <- ggplot(data = dist_data, aes(x = gene_dist, fill = factor(call))) +
     geom_histogram(bins = 60) +
     xlab("Hellinger distance") + 
     ylab("Gene count") +
     background_grid(major = "xy", minor = "none") +
     ## theme_bw() +
     theme(legend.position = "none"))


############################################################
## 3. Panel C this will be a 2x2 table

## could try to save as a ggplot object or recreate later
## https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
str(dist_data)

## select random genes
set.seed(4)
index <- sample(1:nrow(dist_data), size = 50)

(tmp_genes <- rownames(dist_data)[index])
measured_genes <- tmp_genes[tmp_genes %in% rownames(dist_data)]
## genes outside pathway
not_genes <- rownames(dist_data)[!(rownames(dist_data) %in% measured_genes)]

(x11 <- sum(dist_data[measured_genes, "call"]))
(x21 <- sum(dist_data[measured_genes, "call"] == 0))
(x12 <- sum(dist_data[not_genes, "call"]))
(x22 <- sum(dist_data[not_genes, "call"] == 0))
(x11/x21)
(x12/x22)

(to_return <-  (x11/x21)/(x12/x22))

(panel_c <- data.frame(matrix(data = c(x11, x12, x21, x22), nrow = 2)))
names(panel_c) <- c("Alt spliced", "Not")
rownames(panel_c) <- c("Included in pathway", "Not")
panel_c

x <- xtable(panel_c)
align(x) <- xalign(x)
display(x) <- xdisplay(x)
print.xtable(x, include.rownames = T)


############################################################
## 4. Panel D


odds_ratio <- data.frame(OR = odds_ratio)
fil_odds <- odds_ratio[!is.na(odds_ratio)]

pdf("Fig1D.pdf", width = 10)
tmp_locfdr <- locfdr::locfdr(zz = fil_odds, bre = ceiling(length(fil_odds)/8), df = 4,
                             pct = 0, pct0 = 1/4, nulltype = 2, plot = 1, mlests = c(1, sd(fil_odds)))
dev.off()

(panel_d <- ggplot(data = odds_ratio, aes(x = OR)) +
     geom_histogram(bins = 30, color = "white") +
     xlab("Odds ratio across pathways") + 
     ylab("Pathway count") +
     background_grid(major = "xy", minor = "none") +
     geom_vline(xintercept = tmp_locfdr$z.2[2], color = "red"))

## replace panel D by with plot from locFDR


############################################################
## 5. Panel E

(panel_e <- rbind(tmp_data[1:3, c(8,1,3,4,5,6,7)], tmp_data[80,c(8,1,3,4,5,6,7)],
                  tmp_data[80,c(8,1,3,4,5,6,7)]))
## (panel_e[4,] <- t(data.frame(rep("...", 7))))

names(panel_e) <- c("Description", "Odds ratio", "FDR", "Num genes", "Num expressed", "critical", "Alt splice")
panel_e[,"Alt splice"] <- ifelse(panel_e[,"Alt splice"] == 1, "Yes", "No")

x <- xtable(panel_e)
align(x) <- xalign(x)
digits(x) <- xdigits(x, pad =F)
display(x) <- xdisplay(x)
print.xtable(x, include.rownames = F)

## 
## \begin{table}[ht]
## \centering
## \begin{tabular}{lrrrrrl}
##   \hline
## Pathway description & Odds ratio & FDR & Genes & Expressed & Critical OR & Alt splice \\ 
##   \hline
## Asthma & 4.05 & 0.00 &  30 &  15 & 1.93 & Yes \\ 
##   Allograft rejection & 2.90 & 0.00 &  37 &  26 & 1.93 & Yes \\ 
##   Taste transduction & 2.77 & 0.00 &  52 &  16 & 1.93 & Yes \\ 
##   \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\ 
##   Wnt signaling pathway & 1.21 & 1.00 & 150 & 101 & 1.93 & No \\ 
##   \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\ 
##    \hline
## \end{tabular}
## \end{table}
## 

############################################################
## 6. Combine and save

(fig1 <- ggdraw() +
     draw_plot(panel_a, 0, 0.5, 0.33, .5) +
     draw_plot(panel_b, 0.33, 0.5, 0.33, .5) +
     draw_plot(panel_b, 0.66, 0.5, 0.33, .5) +
     draw_plot(panel_d, 0, 0, .5, .5) +
     draw_plot(panel_b, 0.5, 0, .5, .5) +
     draw_plot_label(c("A", "B", "C", "D", "E"),
                     c(0, 0.33, 0.66, 0, 0.5),
                     c(1, 1, 1, 0.5, 0.5)))

## (fig1 <- plot_grid(NULL, panel_b, panel_c, labels = c("A", "B", "C"), ncol = 3))
## (fig1 <- plot_grid(NULL, panel_b, panel_c, labels = c("A", "B", "C"), ncol = 3))

save_plot("Figure1_v2.pdf", fig1,
          ncol = 3, 
          nrow = 2,
          ## base_width = 5,
          base_aspect_ratio = 1.1
          )
