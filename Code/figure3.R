#################################################
#         Figure 3b, score-T2D
#################################################
t2dscore <- read.csv("~/iron/pool/metabolomics/final.T2DSCORE_STD.txt", sep="")
t2dscore <- t2dscore[504:1005, ]

options(bitmapType='cairo')
# plot(t2dscore$T2DSCORE_STD, t2dscore$Estimate)

t2dscore$T2DSCORE_STD <- as.numeric(t2dscore$T2DSCORE_STD)
t2dscore$Estimate <- as.numeric(t2dscore$Estimate)
t2dscore$Lower <- as.numeric(t2dscore$Lower)
t2dscore$Upper <- as.numeric(t2dscore$Upper)

library(ggplot2)
fig3b <- ggplot(data=t2dscore) +
  geom_line(aes(x=T2DSCORE_STD, y=Estimate), linewidth=1, linetype = "solid") +
  geom_line(aes(x=T2DSCORE_STD, y=Lower), linewidth=0.8, linetype = "dashed") +
  geom_line(aes(x=T2DSCORE_STD, y=Upper), linewidth=0.8, linetype = "dashed") +
  geom_ribbon(aes(ymin=Lower,ymax=Upper,x=T2DSCORE_STD), fill="#d7191c", alpha=0.15) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2), limits = c(-2.1,2.5)) +
  scale_y_continuous(breaks = c(-2,-1,0,1,2), limits = c(-2,2.4)) +
  xlab("Standardized T2D metabolomic score") +
  ylab("Log hazard ratio for T2D risk") +
  ggtitle("")+
  theme_bw() +
  theme(panel.border = element_rect(linewidth = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_line(colour = "black", linewidth = 1),
        axis.ticks.length=unit(.2, "cm"),
        plot.title = element_text(size = 20, hjust = 0, color="black"),
        axis.title = element_text(size = 16, color="black"),
        axis.text = element_text(size = 16, color="black")) +
  annotate(geom="text", x=-2, y=2.1, size=5.2, col="black", hjust=0,
           label="Hazard ratio per 1-SD score: \n2.12 (95% CI: 1.81, 2.48)") +
  annotate(geom="text", x=-0.2, y=-1.85, size=5.2, col="black", hjust=0,
           label="P for linearity <0.001")

#################################################
#         Figure 3c, heme-score
#################################################
heme_score=data.frame(quin=c("Q1 (Ref)", "Q2", "Q3", "Q4", "Q5"),
                      beta=c(0, 0.02860640, 0.08779815, 0.09686609, 0.16262381),
                      lci=c(0, -0.02590365, 0.03189386, 0.03853951, 0.09916916),
                      uci=c(0, 0.08311646, 0.14370243, 0.15519268, 0.22607846))

fig3c <- ggplot(aes(x=quin, y=beta, ymin=lci, ymax=uci), data=heme_score) +
  geom_hline(yintercept=0, linewidth=1, lty=2) +  # add a dotted line at x=1 after flip
  geom_errorbar(width=0.2, linewidth=1) +
  geom_point(size=5, shape=18, color="#d7191c") +
  ggtitle("")+
  xlab("Quintiles of heme iron intake") +
  scale_y_continuous(breaks=c(0,0.05,0.10,0.15,0.20,0.25), limits=c(-0.03, 0.26)) +
  ylab("Standardized T2D metabolomic score\n") +
  theme_bw() +  # use a white background
  theme(panel.border = element_rect(linewidth = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_line(colour = "black", linewidth = 1),
        axis.ticks.length = unit(.2, "cm"),
        plot.title = element_text(size = 20, hjust = 0, color="black"),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 16, vjust=-0.5, color="black"), 
        axis.text.y = element_text(size = 16, hjust=0, color="black")) +
  annotate(geom="text", x=0.7, y=0.23, size=5.2, col="black", hjust=0,
           label="Difference of standardized T2D score \nper 1 mg/d heme iron:\n0.16 (95% CI: 0.10, 0.21)") +
  annotate(geom="text", x=3.5, y=-0.02, size=5.2, col="black", hjust=0,
           label="P trend <0.001")

pdf("/udd/nhfwa/iron/archive/figures/figure3bc.pdf", 
    width = 9.6, height = 5.1, onefile = F) # Open a new pdf file
egg::ggarrange(fig3b, fig3c, nrow=1, widths = c(9.5,10))
dev.off() # Close the file

#################################################
#         Figure 3d, metabolite-T2D
#################################################
load("~/iron/pool/metabolomics/metscore_t2d_v2.RData")
gdata::keep(heme_metabolite4, sure=T)
library(readr)
library(tidyverse)

metname <- read_csv(file="/udd/nhfwa/iron/pool/metabolomics/metabolite_name.csv")
heme_metabolite4 <- merge(metname[, c("met", "biochemical_name")],
                          heme_metabolite4[, c(1,4:5)], by="met")

# Create dataset
data <- data.frame(met=heme_metabolite4$biochemical_name,
                   value=abs(heme_metabolite4$beta),
                   group=ifelse(heme_metabolite4$beta>0, "Positive", "Negative"))
data <- data %>% mutate(group=factor(group, levels = c("Positive", "Negative"))) %>% 
  arrange(group, value)
metorder <- data$met

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame(matrix(NA, empty_bar*nlevels(data$group), ncol(data)))
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)

# I substract 0.5 because the letter must have the angle of the center of the bars. 
# Not extreme right(1) or extreme left (0)
angle <- 90 - 360 * (label_data$id - 0.5) /number_of_bar     
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

colgp <- c("#D01C8B80","#4DAC2680")
names(colgp) <- c("Positive", "Negative")
# Make the plot
fig3d <- # Note that id is a factor. If x is numeric, there is some space between the first bar
  ggplot() +
  
  geom_bar(data=data, aes(x=as.factor(id), y=value, fill=group), stat="identity") +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0.4, xend = start, yend = 0.4),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.3, xend = start, yend = 0.3),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.2, xend = start, yend = 0.2),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.1, xend = start, yend = 0.1),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  
  geom_bar(data=data, aes(x=as.factor(id), y=value, fill=group), stat="identity") +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(0.1, 0.2, 0.3, 0.4), label = c("0.1", "0.2", "0.3", "0.4") , 
           color="grey", size=4 , angle=0, fontface="bold", hjust=1) +
  scale_fill_manual(values = colgp, name = "Association \nwith T2D") +
  ylim(-0.1, 0.7) +
  theme_minimal() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12, color="black"),
        legend.title = element_text(size = 12, color="black"),
        plot.margin = unit(c(-3,-1,-1,-1), "cm"),
        plot.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  coord_polar() + 
  geom_text(data=label_data, 
            aes(x=id, y=value + 0.15, label=met), 
            color="black", size=4, 
            angle= label_data$angle, inherit.aes = FALSE)

fig3d_legend <- # Note that id is a factor. If x is numeric, there is some space between the first bar
  ggplot() +
  
  geom_bar(data=data, aes(x=as.factor(id), y=value, fill=group), stat="identity") +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0.4, xend = start, yend = 0.4),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.3, xend = start, yend = 0.3),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.2, xend = start, yend = 0.2),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.1, xend = start, yend = 0.1),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  
  geom_bar(data=data, aes(x=as.factor(id), y=value, fill=group), stat="identity") +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(0.1, 0.2, 0.3, 0.4), label = c("0.1", "0.2", "0.3", "0.4") , 
           color="grey", size=4 , angle=0, fontface="bold", hjust=1) +
  scale_fill_manual(values = colgp, name = "Association \nwith T2D") +
  ylim(-0.1, 0.7) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12, color="black"),
        legend.title = element_text(size = 12, color="black"),
        plot.margin = unit(c(-3,-1,-1,-1), "cm"),
        plot.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  coord_polar() + 
  geom_text(data=label_data, 
            aes(x=id, y=value + 0.15, label=met), 
            color="black", size=4, 
            angle= label_data$angle, inherit.aes = FALSE)

pdf("/udd/nhfwa/iron/archive/figures/figure3d.pdf", 
    width = 6, height = 6, onefile = F) # Open a new pdf file
print(fig3d)
dev.off() # Close the file

pdf("/udd/nhfwa/iron/archive/figures/figure3d_legend.pdf", 
    width = 2, height = 1.5, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(fig3d_legend))
dev.off() # Close the file

#################################################
#         Figure 3e, heme-metabolite
#################################################
load("~/iron/pool/metabolomics/metscore_t2d_v2.RData")
gdata::keep(heme_metabolite4, metorder, sure=T)
library(readr)
library(tidyverse)

metname <- read_csv(file="/udd/nhfwa/iron/pool/metabolomics/metabolite_name.csv")
heme_metabolite4 <- merge(metname[, c("met", "biochemical_name")],
                          heme_metabolite4[, c(1,4:5)], by="met")

# Create dataset
data <- data.frame(met=heme_metabolite4$biochemical_name,
                   value=abs(heme_metabolite4$beta_heme),
                   group=ifelse(heme_metabolite4$beta_heme>0, "Positive", "Negative"))
data$met <- factor(data$met, levels = metorder)
data <- data %>% mutate(group=factor(group, levels = c("Positive", "Negative"))) %>% 
  arrange(group, met)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame(matrix(NA, empty_bar*nlevels(data$group), ncol(data)))
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)

# I substract 0.5 because the letter must have the angle of the center of the bars. 
# Not extreme right(1) or extreme left (0)
angle <- 90 - 360 * (label_data$id - 0.5) /number_of_bar     
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

colgp <- c("#E6610180", "#5E3C9980")
names(colgp) <- c("Positive", "Negative")
# Make the plot
fig3e <- # Note that id is a factor. If x is numeric, there is some space between the first bar
  ggplot() +
  
  geom_bar(data=data, aes(x=as.factor(id), y=value, fill=group), stat="identity") +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0.2, xend = start, yend = 0.2),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.15, xend = start, yend = 0.15),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.1, xend = start, yend = 0.1),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.05, xend = start, yend = 0.05),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  
  geom_bar(data=data, aes(x=as.factor(id), y=value, fill=group), stat="identity") +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(0.05, 0.1, 0.15, 0.2), label = c("0.05", "0.1", "0.15", "0.2") , 
           color="grey", size=4 , angle=0, fontface="bold", hjust=1) +
  scale_fill_manual(values = colgp, name = "Association \nwith heme iron") +
  ylim(-0.1, 0.7) +
  theme_minimal() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12, color="black"),
        legend.title = element_text(size = 12, color="black"),
        plot.margin = unit(c(-3,-1,-1,-1), "cm"),
        plot.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  coord_polar() + 
  geom_text(data=label_data, 
            aes(x=id, y=value + 0.15, label=met), 
            color="black", size=4, 
            angle= label_data$angle, inherit.aes = FALSE)

fig3e_legend <- # Note that id is a factor. If x is numeric, there is some space between the first bar
  ggplot() +
  
  geom_bar(data=data, aes(x=as.factor(id), y=value, fill=group), stat="identity") +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0.2, xend = start, yend = 0.2),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.15, xend = start, yend = 0.15),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.1, xend = start, yend = 0.1),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  geom_segment(data=grid_data, aes(x = end, y = 0.05, xend = start, yend = 0.05),
               colour = "grey", alpha=1, linewidth=0.7 , inherit.aes = FALSE) +
  
  geom_bar(data=data, aes(x=as.factor(id), y=value, fill=group), stat="identity") +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(0.05, 0.1, 0.15, 0.2), label = c("0.05", "0.1", "0.15", "0.2") , 
           color="grey", size=4 , angle=0, fontface="bold", hjust=1) +
  scale_fill_manual(values = colgp, name = "Association \nwith heme iron") +
  ylim(-0.1, 0.7) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12, color="black"),
        legend.title = element_text(size = 12, color="black"),
        plot.margin = unit(c(-3,-1,-1,-1), "cm"),
        plot.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  coord_polar() + 
  geom_text(data=label_data, 
            aes(x=id, y=value + 0.15, label=met), 
            color="black", size=4, 
            angle= label_data$angle, inherit.aes = FALSE)

pdf("/udd/nhfwa/iron/archive/figures/figure3e.pdf", 
    width = 6, height = 6, onefile = F) # Open a new pdf file
print(fig3e)
dev.off() # Close the file

pdf("/udd/nhfwa/iron/archive/figures/figure3e_legend.pdf", 
    width = 2, height = 1.5, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(fig3e_legend))
dev.off() # Close the file
