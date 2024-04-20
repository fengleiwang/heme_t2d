library(haven)

# readin spline data
t2d0 <- read.csv("~/iron/pool/final.HEME_C.txt", sep="")
for (i in 1:4){
  t2d0[,i] <- as.numeric(t2d0[,i])
}

summary(t2d0$HEME_C)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2361  0.7149  1.1975  1.1984  1.6802  2.1628

# readin distribution data
t2d_heme <- read_sas("iron/pool/t2d_heme.sas7bdat", NULL)
summary(t2d_heme$heme_c)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00515  0.83410  1.04087  1.07671  1.28224 13.78213 

library(ggplot2)
library(scales)
scaleFUN <- function(x) sprintf("%.1f", x)

p_t2d1=ggplot(data=t2d0) +
  geom_line(aes(x=HEME_C, y=Estimate), size=1, linetype = "solid") +
  geom_line(aes(x=HEME_C, y=Lower), size=0.8, linetype = "dashed") +
  geom_line(aes(x=HEME_C, y=Upper), size=0.8, linetype = "dashed") +
  geom_vline(xintercept=0.6590074, size=0.8, lty=3) +  # quintile 1
  geom_vline(xintercept=0.8773392, size=0.8, lty=3) +  # quintile 2
  geom_vline(xintercept=1.0393297, size=0.8, lty=3) +  # quintile 3
  geom_vline(xintercept=1.2213600, size=0.8, lty=3) +  # quintile 4
  geom_vline(xintercept=1.5409666, size=0.8, lty=3) +  # quintile 5
  geom_ribbon(aes(ymin=Lower,ymax=Upper,x=HEME_C), fill="#d7191c", alpha=0.15) +
  scale_x_continuous(breaks = c(0.5,1.0,1.5,2.0), limits = c(0.2,2.2)) +
  xlab("") +
  ylab("HR (95% CI) for type 2 diabetes") +
  ggtitle("")+
  theme_bw() +
  theme(panel.border = element_rect(size = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length=unit(.3, "cm"),
        plot.title = element_text(size = 36, hjust = 0, color="black"),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 20, color="black"),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20, color="black"))

p_t2d2 <- ggplot(t2d_heme, aes(x = heme_c)) + 
  geom_histogram(aes(y = ..density..), colour = 1, fill = "#d7191c", alpha=0.15) +
  geom_density() +
  scale_x_continuous(breaks = c(0.5,1.0,1.5,2.0), limits = c(0.2,2.2)) +
  xlab("Energy-adjusted heme iron\nintake (mg/d)") +
  ylab("Density") +
  theme_bw() +
  theme(panel.border = element_rect(size = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length=unit(.3, "cm"),
        plot.title = element_blank(),
        axis.title.x = element_text(size = 20, color="black"), 
        axis.title.y = element_text(size = 20, color="black"),
        axis.text.x = element_text(size = 20, color="black"), 
        axis.text.y = element_text(size = 20, color="black"))

####################################
#           Subgroup
####################################
library(readxl)
t2d <- read_excel("/udd/nhfwa/iron/meta_0905.xlsx", 
                  sheet = "t2d_sub", 
                  col_types = c("text", "text", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "text", "text"))
lab1=t2d$Lab
t2d$Lab=factor(t2d$Lab, levels=lab1)
t2d$Lab=factor(t2d$Lab, levels=rev(lab1))
t2d$factor=factor(t2d$factor, level=c("age", "sex", "bmi", "race",
                                      "meno", "alcohol", "phyt"))

plot1=ggplot(data=t2d, aes(x=Lab, y=HR, ymin=lci, ymax=uci)) +
  geom_point(size=5, shape=18) +
  geom_errorbar(width=0.4, size=0.9) +
  ggtitle("")+
  geom_hline(yintercept=1, size=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("") +
  scale_y_continuous(breaks=seq(1.0, 1.8, by=0.2), limits=c(0.9, 3)) +
  ylab("\nHR (95% CI)") +
  coord_flip() +
  theme_bw() +  # use a white background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  theme(axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.length.x=unit(.3, "cm"),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 36, hjust = 0, color="black"),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 20, vjust=-0.5, color="black"), 
        axis.text.y = element_text(size = 20, hjust=0, color="black")) +
  geom_text(aes(y = 1.9,label = HRci), size = 7.2, color="black", hjust=0) +
  geom_text(aes(y = 2.7,label = interp), size = 7.2, color="black", hjust=0)

library(patchwork)
pdf("/udd/nhfwa/iron/figure1ab.pdf", 
    width = 22, height = 10, onefile = F) # Open a new pdf file
((p_t2d1/p_t2d2 + plot_layout(heights = c(3,1))) | plot1) + plot_layout(widths = c(1,2))
dev.off() # Close the file

# pdf("/udd/nhfwa/iron/figure1a.pdf", 
#     width = 7, height = 10, onefile = F) # Open a new pdf file
# egg::ggarrange(p_t2d1, p_t2d2,
#                nrow=2, heights = c(3,1))  
# dev.off() # Close the file
# 
# pdf("/udd/nhfwa/iron/figure1b.pdf", 
#     width = 15, height = 10, onefile = F) # Open a new pdf file
# print(plot1)
# dev.off() # Close the file
