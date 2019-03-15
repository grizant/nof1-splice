## AG Schissler
## Explore N1PAS KEGG based simulation results
## 2 Mar 2019

#################################################################################
## 1. Retrieve simulation data
#################################################################################

setwd("~/Dropbox/Splice-n-of-1-pathways/Code/simulation")
system.time(sim_data <- readRDS("sim_data_n1pas_2mar2019.rds"))

rownames(sim_data) <- NULL

#################################################################################
## 2. Explore a bit
#################################################################################

## check the total amount of runs
(num_pat <- 7 + 19 + 51 + 52 + 58 + 59) ## 246 total
num_p <- 4
num_pi <- 5
reps <- 100
(num_rows_expected <- num_pat * num_p * num_pi * reps) == nrow(sim_data)

## misfits
sum(is.na(sim_data$target_captured)) ## 625 misfits
sum(is.na(sim_data$target_captured)) / nrow(sim_data) ## 0.00127

## remove the misfits to make coding easier later
sim_data <- sim_data[!is.na(sim_data$target_captured), ]

## capitalize tgca datasets
sim_data$dataset <- toupper(sim_data$dataset)
rownames(sim_data) <- NULL
head(sim_data)
table(sim_data$sim_p)

## explore asg background
str(unique(sim_data$prop_asg))

library("dplyr")

prop_asg_pat_data <- sim_data %>%
    filter(pi == 0) %>%
    group_by(patientID) %>%
    summarize(mean_prop_asg = mean(prop_asg))

summary(prop_asg_pat_data$mean_prop_asg)

sort(prop_asg_pat_data$mean_prop_asg)


#################################################################################
## 3. Develop false positive rate figure
#################################################################################

## DOES NOT MAKE SENSE TO EXPLORE THE GENE SET SIZE DIMENSION

## summarize at the patient level
summary(sim_data[sim_data$pi == 0, "false_positive_rate"])
## never more that 15.6%

fpr_pat_data <- sim_data %>%
    filter(pi == 0) %>%
    group_by(dataset, patientID) %>%
    summarize(mean_pat_fpr = mean(false_positive_rate))

fpr_pat_data

alpha <- 0.05

fpr_data <- fpr_pat_data %>%
    group_by(dataset) %>%
    summarize(mean_fpr = mean(mean_pat_fpr),
              sd_fpr = sd(mean_pat_fpr))

fpr_data <- fpr_data %>%
    mutate(lower = mean(mean_fpr) + qnorm(alpha / 2) * sd_fpr,
           upper = mean(mean_fpr) + qnorm(1 -( alpha / 2)) * sd_fpr)

fpr_data$upper[fpr_data$upper > 1] <- 1
fpr_data$lower[fpr_data$lower < 0] <- 0

fpr_data
head(fpr_data)
as.data.frame(fpr_data)


## figure settings
purple <- "#A078AA"
orange <- "#FFA040"
green <- "#72BF44"
grey <- "#777676"
yellow <- "#F9F080"

point_size <- 0.5
my_width <- 0.5
pd <- position_dodge(width = 7)

## p0 <- ggplot(data = fpr_data, aes(x = factor(dataset), y = mean_fpr, color = factor(sim_p)))
## p0 <- ggplot(data = fpr_data, aes(x = sim_p, y = mean_fpr, color = factor(dataset)))
## p0 <- ggplot(data = fpr_data, aes(x = dataset, y = mean_fpr, color = factor(dataset)))
p0 <- ggplot(data = fpr_data, aes(x = dataset, y = mean_fpr))
## p0 <- ggplot(data = fpr_data, aes(x = sim_p, y = mean_fpr, shape = factor(dataset)))
## p0 <- p0 + geom_hline(yintercept = 0.2) + geom_hline(yintercept = 0.05)
## p0 <- p0 + geom_hline(yintercept = 0.05, color = "red", linetype = "dashed")
## p1 <- p0 + geom_point(position = pd)
p1 <- p0 + geom_point()
## p2 <- p1 + geom_line(size = point_size, position = pd) + facet_wrap(facets = vars(dataset))
## p2 <- p1 + geom_line(size = point_size, position = pd)
p2 <- p1
p3 <- p2 + geom_errorbar(aes(ymin = lower,ymax = upper), width = my_width, size = point_size, position = pd)
p3 <- p2 + geom_errorbar(aes(ymin = lower,ymax = upper), width = my_width, size = point_size)
## p4 <- p3 + scale_color_manual(values=guide = guide_legend(title = "TCGA data set"))
## p4 <- p3 + guide_legend("TCGA data set")
p4 <- p3
p5 <- p4 + theme_bw(base_size = 12)
p6 <- p5 + labs(title="",x="TCGA data set", y="Average false-positive rates across patients")
p7 <- p6 + theme(legend.position = "bottom")
## p8 <- p3 + ylim(0, 0.0)
## p4 <- p3 + geom_hline(yintercept = 0.2)
## p8 <- p7 + ylim(0, 0.2)
## p8 <- p7 + guides(color = guide_legend(title="TCGA data set"))
p8 <- p7
p8

ggsave(plot = p8, height=8, width=7, dpi=200, filename = paste0("~/Dropbox/Splice-n-of-1-pathways/Figures/CurTrendsInTransBioinfor Figures/Figure6noGeneSetSize_gs.pdf"), useDingbats=FALSE)


#################################################################################
## 3. Develop power figure
#################################################################################

library("dplyr")

## number of simulations in power study
sum(sim_data$pi > 0)

## summarize at the patient level
power_pat_data <- sim_data %>%
    filter(pi > 0) %>%
    group_by(dataset, sim_p, pi, patientID) %>%
    summarize(power = sum(target_captured) / n())


str(power_pat_data)
tail(power_pat_data)

power_data <- power_pat_data %>%
    group_by(dataset, sim_p, pi) %>%
    summarize(mean_power = mean(power),
              sd_power = sd(power))

power_data <- power_data %>%
    group_by(dataset, sim_p, pi) %>%
    summarize(mean_power = mean_power,
              lower = mean(mean_power) + qnorm(alpha / 2) * sd_power,
              upper = mean(mean_power) + qnorm(1 -( alpha / 2)) * sd_power)

power_data$upper[power_data$upper > 1] <- 1
power_data$lower[power_data$lower < 0] <- 0

str(power_data)
tail(power_data)

as.data.frame(power_data[power_data$pi == 0.15,])
summary(power_data[power_data$pi == 0.15, "mean_power"])
summary(power_data[power_data$pi == 0.2, "mean_power"])


## develop figure
### lab colors:
purple <- "#A078AA"
orange <- "#FFA040"
green <- "#72BF44"
grey <- "#777676"
yellow <- "#F9F080"

point_size <- 1
my_width = 0.03
pd <- position_dodge(width = 0.03)

## use pi as x axis
p0 <- ggplot(data = power_data, aes(x = pi, y = mean_power, color = factor(sim_p)))
## p1 <- p0 + geom_point(size=point_size, col = "black") + geom_line(aes(linetype = factor(sim_p)), col = "black", size = point_size)
p1 <- p0 + geom_point(size=point_size, position = pd)
## p1 <- p1 + geom_line(aes(linetype = factor(sim_p)), size = point_size, position = pd)
p1 <- p1 + geom_line(size = point_size / 2, position = pd)
## p2 <- p1 + facet_wrap(facets = vars(dataset), nrow = 3)
p2 <- p1 + facet_wrap(facets = vars(dataset), nrow = 2)
p3 <- p2 + geom_errorbar(aes(ymin = lower, ymax = upper), width = my_width, size = point_size / 2, position = pd)
p4 <- p3 + scale_color_manual(values=c(purple, orange, green, grey), guide = guide_legend(title = "Expressed genes in pathway"))
## p4 <- p3 + scale_linetype_manual(values=c(4, 2, 3, 1), guide = guide_legend(title = "Pathway size"))
## p4 <- p4 + guides(color=guide_legend(), line = guide_legend())
p5 <- p4 + theme_bw(base_size = 12)
p6 <- p5 + labs(title="",x="Proportion of ASGs above background", y="Average power across patients")
## p7 <- p6 + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
## p7 <- p6 + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right")
p7 <- p6 + theme(legend.position = "bottom")
## p8 <- p7 + xlim(0.05, 0.20)
p8 <- p7
p8

ggsave(plot = p8 ,height=8, width=7, dpi=200, filename = paste0("~/Dropbox/Splice-n-of-1-pathways/Figures/CurTrendsInTransBioinfor Figures/Figure7.pdf"), useDingbats=FALSE)

