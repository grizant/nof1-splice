## Create Figure 1
## AG Schissler
## Created 27 Jul 2017

############################################################
## i. Set up environment

library(ggplot2)
library(cowplot)

setwd("~/Dropbox/splice-n-of-1-pathways/Figures")

plot.mpg <- ggplot(mpg, aes(x = cty, y = hwy, colour = factor(cyl))) + 
  geom_point(size=2.5)
plot.mpg

fig3 <- ggdraw(plot.mpg) + 
    draw_plot_label("A", size = 14) + 
    draw_label("DRAFT!", angle = 35, size = 80, alpha = .2)

save_plot("Figure1.pdf", fig3,
          # each individual subplot should have an aspect ratio of 1.3
          )
