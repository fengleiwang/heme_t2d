library(tidyverse)

load("/udd/nhfwa/iron/pool/biomarker/biomarker.RData")

#### figure 2a
fig2a <- results %>% mutate(iron=case_when(type=="heme" ~ "Heme iron (per 1 mg/d)",
                                          type=="iron" ~ "Total iron (per 20 mg/d)",
                                          type=="nonheme" ~ "Nonheme iron (per 20 mg/d)",
                                          type=="iron_dt" ~ "Dietary iron (per 10 mg/d)",
                                          type=="iron_sp" ~ "Supplemental iron (per 20 mg/d)"),
                           group=case_when(bio=="adj_lna1c" ~ "HbA1c",
                                           bio=="adj_lnadipo" ~ "Adiponectin",
                                           bio=="adj_lnchol" ~ "Total cholesterol",
                                           bio=="adj_lncpep" ~ "C-peptide",
                                           bio=="adj_lncrp" ~ "CRP",
                                           bio=="adj_lnhdl" ~ "HDL-C",
                                           bio=="adj_lnldl" ~ "LDL-C",
                                           bio=="adj_lnleptin" ~ "Leptin",
                                           bio=="adj_lnratio" ~ "TAG/HDL-C ratio",
                                           bio=="adj_lntrig" ~ "TAG",
                                           bio=="adj_lnferritin" ~ "Ferritin",
                                           bio=="adj_lntrf" ~ "TfR/Ferritin ratio"))

fig2a$group <- factor(fig2a$group,
                        levels = c("C-peptide", "HbA1c", "Total cholesterol",
                                   "HDL-C", "LDL-C", "TAG", "TAG/HDL-C ratio",
                                   "CRP", "Adiponectin", "Leptin", "Ferritin", "TfR/Ferritin ratio"))
fig2a$iron <- factor(fig2a$iron,
                     levels=c("Total iron (per 20 mg/d)", "Heme iron (per 1 mg/d)", "Nonheme iron (per 20 mg/d)",
                              "Dietary iron (per 10 mg/d)", "Supplemental iron (per 20 mg/d)"))

ironcol <- c("#ffffbf", "#d7191c", "#abd9e9", "#fdae61", "#2c7bb6")
names(ironcol) <- c("Total iron (per 20 mg/d)", "Heme iron (per 1 mg/d)", "Nonheme iron (per 20 mg/d)",
                    "Dietary iron (per 10 mg/d)", "Supplemental iron (per 20 mg/d)")

library(ggplot2)
library(ggpubr)
library(grid)
p2a <- ggplot(aes(x=group, y=beta, ymin=lci, ymax=uci, fill=iron), data=fig2a) +
  geom_col(width = 0.9, color="black", position = "dodge") +
  geom_errorbar(width=0.4, size=1, position=position_dodge(width=0.9)) +
  ggtitle("a") +
  xlab("") +
  ylab("% Difference (95% CI)") +
  scale_y_continuous(breaks=c(-30,-15,0,15,30), limits=c(-32,38)) +
  scale_fill_manual(values = ironcol) +
  # coord_flip() +
  theme_bw() +  # use a white background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 2)) +
  theme(axis.title = element_text(size = 18), 
        axis.ticks = element_line(size=1.2),
        axis.ticks.length = unit(.2, "cm"),
        plot.title = element_text(size = 32, hjust = 0, color="black"),
        axis.text.x = element_text(size = 18, angle = 30, vjust=1, hjust=1,color="black"), 
        axis.text.y = element_text(size = 18, hjust=0, color="black"),
        strip.text.x = element_text(size = 18, angle = 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = "none")

p2a_legend <- ggplot(aes(x=group, y=beta, ymin=lci, ymax=uci, fill=iron), data=fig2a) +
  geom_col(width = 0.9, color="black", position = "dodge") +
  geom_errorbar(width=0.4, size=1, position=position_dodge(width=0.9)) +
  xlab("") +
  ylab("% Difference (95% CI)") +
  labs(fill="Types of iron intake") +
  scale_y_continuous(breaks=c(-30,-15,0,15,30), limits=c(-30,38)) +
  scale_fill_manual(values = ironcol) +
  # coord_flip() +
  theme_bw() +  # use a white background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 2)) +
  theme(axis.title = element_text(size = 18), 
        axis.ticks = element_line(size=1.2),
        axis.ticks.length = unit(.2, "cm"),
        axis.text.x = element_text(size = 18, angle = 30, vjust=1, hjust=1,color="black"), 
        axis.text.y = element_text(size = 18, hjust=0, color="black"),
        strip.text.x = element_text(size = 18, angle = 0),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "right")

pdf("/udd/nhfwa/iron/archive/figures/figure2a_legend.pdf", 
    width = 5.5, height = 5.5, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(p2a_legend))
dev.off() # Close the file

#### figure 2b
results1 <- fig2a %>% filter(type=="heme")
results2 <- results_con %>% filter(type=="heme")
results3 <- results_case %>% filter(type=="heme")
identical(results1$bio, results2$bio) # TRUE
identical(results1$bio, results3$bio) # TRUE

sum(results1$bon<0.05 & results2$bon<0.05) # 6

fig2b <- data.frame(bio=results1$group,
                    beta1=results1$beta,
                    beta2=results2$beta,
                    beta3=results3$beta,
                    group=ifelse(results1$bon<0.05, 
                                 "Significant (Bonferroni-adjusted P <0.05)", 
                                 "Nonsignificant"))
fig2b$group <- factor(fig2b$group, levels=c("Significant (Bonferroni-adjusted P <0.05)", 
                                            "Nonsignificant"))

dotcolor <- c("#d7191c", "grey")
names(dotcolor) <- c("Significant (Bonferroni-adjusted P <0.05)", 
                     "Nonsignificant")

library(ggrepel)
p2b=ggplot(aes(x=beta1, y=beta2, label=bio), data=fig2b) +
  geom_abline(slope=1, linetype=2, size=1.5) +
  geom_hline(yintercept = 0, linetype=2, size=1.5) +
  geom_vline(xintercept = 0, linetype=2, size=1.5) +
  geom_point(aes(color=group), size=3.5) +
  scale_y_continuous(breaks=c(-24, -12, 0, 12, 24), limits=c(-24, 27)) +
  scale_x_continuous(breaks=c(-24, -12, 0, 12, 24), limits=c(-24, 27)) +
  scale_color_manual(values=dotcolor) +
  labs(
    title = "b",
    x="% Difference\namong all participants",
    y="% Difference\namong controls") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=2.2)) +
  theme(plot.title = element_text(hjust = 0, size=32, color="black"),
        axis.title = element_text(size=18, color="black"),
        axis.text = element_text(size=18, color="black"),
        axis.ticks = element_line(size=1.2),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none") +
  geom_label_repel(size=5,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(aspect.ratio = 1)

p2b_legend=ggplot(aes(x=beta1, y=beta2, label=bio), data=fig2b) +
  geom_abline(slope=1, linetype=2, size=1.5) +
  geom_hline(yintercept = 0, linetype=2, size=1.5) +
  geom_vline(xintercept = 0, linetype=2, size=1.5) +
  geom_point(aes(color=group), size=3.5) +
  scale_y_continuous(breaks=c(-24, -12, 0, 12, 24), limits=c(-24, 27)) +
  scale_x_continuous(breaks=c(-24, -12, 0, 12, 24), limits=c(-24, 27)) +
  scale_color_manual(values=dotcolor) +
  labs(
    # title = "Associations between heme iron \nand biomakers among controls",
    x="Percentage of difference\namong all participants",
    y="Percentage of difference\namong controls",
    color="Significance of heme iron-biomarker \nassociations among all participants") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=2.2)) +
  theme(plot.title = element_text(hjust = 0.5, size=18, color="black"),
        axis.title = element_text(size=18, color="black"),
        axis.text = element_text(size=18, color="black"),
        axis.ticks = element_line(size=1.2),
        axis.ticks.length = unit(.2, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "right")

pdf("/udd/nhfwa/iron/archive/figures/figure2b_legend.pdf", 
    width = 5.5, height = 5.5, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(p2b_legend))
dev.off() # Close the file

blank <- ggplot()+theme_void()

library(patchwork)
pdf("/udd/nhfwa/iron/archive/figures/figure2.pdf", 
    width = 12, height = 12, onefile = F) # Open a new pdf file
p2a/(p2b+blank)
dev.off() # Close the file
