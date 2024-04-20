library(tidyverse)
load("~/iron/pool/metabolomics/metabolite_biomarker.RData")
bio_all %>% group_by(bio) %>% summarise(min=min(per),
                                        max=max(per),
                                        absmax=max(abs(per)))
library(readr)
metname <- read_csv(file="/udd/nhfwa/iron/pool/metabolomics/metabolite_name.csv")

####################################
#           heatmap 1
####################################
bio_all <- bio_all %>% mutate(type=case_when(bio=="adj_lna1c" ~ "HbA1c",
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
bio_all$type <- factor(bio_all$type, 
                       levels = c("C-peptide", "HbA1c", "Total cholesterol",
                                  "HDL-C", "LDL-C", "TAG", "TAG/HDL-C ratio",
                                  "CRP", "Adiponectin", "Leptin", "Ferritin", "TfR/Ferritin ratio"))

sum(bio_all$per>40) # 7
sum(bio_all$per<(-20)) # 3
bio_all$per[bio_all$per>40] <- 40
bio_all$per[bio_all$per<(-20)] <- (-20)

bio_all <- merge(metname[, c("met", "biochemical_name")],
                 bio_all[, -c(2:3)], by="met")
hh1 <- reshape(bio_all[, c("biochemical_name", "type", "per")] %>% arrange(type), 
               idvar = "biochemical_name", timevar = "type", direction = "wide")
rownames(hh1) <- hh1$biochemical_name
hh1 <- hh1[,-1]
colnames(hh1) <- gsub("per.", "", colnames(hh1))

hhp1 <- reshape(bio_all[, c("biochemical_name", "type", "bon")] %>% arrange(type),  
                idvar = "biochemical_name", timevar = "type", direction = "wide")
rownames(hhp1) <- hhp1$biochemical_name
hhp1 <- hhp1[,-1]
colnames(hhp1) <- gsub("bon.", "", colnames(hhp1))
identical(rownames(hh1), rownames(hhp1))

library(ComplexHeatmap)
library(circlize)
col_fun1 = colorRamp2(c(-20, 0, 40), c("#08519c", "white", "#bd0026"))

met_coef <- bio_all[bio_all$type=="CRP", c("biochemical_name", 'coef')]
colfunc_pos <- colorRampPalette(c("#8c510a", "white"))
colfunc_neg <- colorRampPalette(c("#01665e", "white"))
colfunc_pos(6)
colfunc_neg(5)
met_coef <- met_coef %>% mutate(colgp=case_when(coef<=(-0.15) ~ "(-0.20, -0.15]",
                                                coef>(-0.15) & coef<=(-0.10) ~ "(-0.15, -0.10]",
                                                coef>(-0.10) & coef<=(-0.05) ~ "(-0.10, -0.05]",
                                                coef>(-0.05) & coef<=(0) ~ "(-0.05, 0]",
                                                coef>(0) & coef<=(0.05) ~ "(0, 0.05]",
                                                coef>(0.05) & coef<=(0.10) ~ "(0.05, 0.10]",
                                                coef>(0.10) & coef<=(0.15) ~ "(0.10, 0.15]",
                                                coef>(0.15) & coef<=(0.20) ~ "(0.15, 0.20]",
                                                coef>(0.20) & coef<=(0.25) ~ "(0.20, 0.25)"))
met_coef$colgp <- factor(met_coef$colgp,
                         levels=c("(0.20, 0.25)",
                                  "(0.15, 0.20]",
                                  "(0.10, 0.15]",
                                  "(0.05, 0.10]",
                                  "(0, 0.05]",
                                  "(-0.05, 0]",
                                  "(-0.10, -0.05]",
                                  "(-0.15, -0.10]",
                                  "(-0.20, -0.15]"))

met_gp = rowAnnotation(coef = met_coef$colgp,
                       show_annotation_name = F,
                       col = list(coef = c("(0.20, 0.25)" = "#8C510A",
                                           "(0.15, 0.20]" = "#A3733B",
                                           "(0.10, 0.15]" = "#BA966C",
                                           "(0.05, 0.10]" = "#D1B99D",
                                           "(0, 0.05]" = "#E8DCCE",
                                           "(-0.05, 0]" = "#BFD8D6",
                                           "(-0.10, -0.05]" = "#80B2AE",
                                           "(-0.15, -0.10]" = "#408C86",
                                           "(-0.20, -0.15]" = "#01665E")),
                       annotation_legend_param = list(coef = list(
                         title = "Weight in the T2D \nmetabolomic score",
                         title_gp = gpar(col = "black", fontsize = 16),
                         labels_gp = gpar(col = "black", fontsize = 16),
                         title_position = "topleft"
                       )))

heat1 <- Heatmap(as.matrix(hh1), 
                 col = col_fun1,
                 row_dend_width = unit(2.5, "cm"),
                 row_dend_gp = gpar(lwd = 4),
                 show_row_names = T,
                 row_names_gp = gpar(fontsize = 18),
                 row_names_max_width = unit(25, "cm"),
                 column_labels = c("C-peptide", "HbA1c", "Total cholesterol",
                                   "HDL-C", "LDL-C", "TAG", "TAG/HDL-C ratio",
                                   "CRP", "Adiponectin", "Leptin", "Ferritin", "TfR/Ferritin ratio"),
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 18),
                 column_names_rot = 60,
                 column_title = NULL,
                 border = TRUE,
                 border_gp = gpar(col = "black", lwd = 4),
                 cluster_columns = FALSE,
                 cell_fun = function(j, i, x, y, w, h, fill){
                   if(hhp1[i, j] < 0.05) {
                     grid.text("*", x, y, vjust = 0.8,
                               gp = gpar(fontsize = 27))
                   }},
                 show_heatmap_legend = F,
                 right_annotation = met_gp)

lgd1 = Legend(title = "% Difference of biomarker \nper 1-SD metabolite", 
              col_fun = col_fun1,
              at=c(-20,0,20,40),
              labels = c("-20", "0", "20", "40"),
              legend_width = unit(2, "cm"),
              legend_height = unit(4, "cm"),
              grid_height = unit(0.7, "cm"),
              title_gp = gpar(col = "black", fontsize = 18),
              labels_gp = gpar(col = "black", fontsize = 18),
              direction = "vertical", 
              title_position = "topleft")

pdf("/udd/nhfwa/iron/archive/figures/metabolite_biomarker_heatmap_1.pdf", 
    width = 13, height = 21, onefile = F) # Open a new pdf file
draw(heat1, heatmap_legend_list = lgd1)
dev.off() # Close the file

####################################
#           heatmap 2
####################################
bio_all2 <- bio_all2 %>% mutate(type=case_when(bio=="adj_lna1c" ~ "HbA1c",
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
bio_all2$type <- factor(bio_all2$type, 
                        levels = c("C-peptide", "HbA1c", "Total cholesterol",
                                   "HDL-C", "LDL-C", "TAG", "TAG/HDL-C ratio",
                                   "CRP", "Adiponectin", "Leptin", "Ferritin", "TfR/Ferritin ratio"))

sum(bio_all2$per>40) # 10
sum(bio_all2$per<(-20)) # 5
bio_all2$per[bio_all2$per>40] <- 40
bio_all2$per[bio_all2$per<(-20)] <- (-20)

bio_all2 <- merge(metname[, c("met", "biochemical_name")],
                  bio_all2[, -c(2:3)], by="met")
hh2 <- reshape(bio_all2[, c("biochemical_name", "type", "per")] %>% arrange(type), 
               idvar = "biochemical_name", timevar = "type", direction = "wide")
rownames(hh2) <- hh2$biochemical_name
hh2 <- hh2[,-1]
colnames(hh2) <- gsub("per.", "", colnames(hh2))

hhp2 <- reshape(bio_all2[, c("biochemical_name", "type", "bon")] %>% arrange(type),  
                idvar = "biochemical_name", timevar = "type", direction = "wide")
rownames(hhp2) <- hhp2$biochemical_name
hhp2 <- hhp2[,-1]
colnames(hhp2) <- gsub("bon.", "", colnames(hhp2))
identical(rownames(hh2), rownames(hhp2))

heme_metabolite4 <- merge(metname[, c("met", "biochemical_name")],
                          heme_metabolite4[, c(1,4:5)], by="met")

met_beta_heme <- heme_metabolite4[, c("biochemical_name", 'beta_heme')]
met_beta_t2d <- heme_metabolite4[, c("biochemical_name", 'beta')]
summary(heme_metabolite4$beta_heme) # -0.2, 0.3
summary(heme_metabolite4$beta) # -0.4, 0.6

met_beta_heme <- met_beta_heme %>% 
  mutate(colgp=case_when(beta_heme<=0 ~ "Negative",
                         beta_heme>0 ~ "Positive"))
met_beta_heme$colgp <- factor(met_beta_heme$colgp,
                              levels=c("Negative", "Positive"))
met_beta_t2d <- met_beta_t2d %>% 
  mutate(colgp=case_when(beta<=0 ~ "Negative",
                         beta>0 ~ "Positive"))
met_beta_t2d$colgp <- factor(met_beta_t2d$colgp,
                             levels=c("Negative", "Positive"))

# annotation for medaition proportion
load("/udd/nhfwa/iron/pool/metabolomics/metscore_t2d_mediation.RData")
met_per <- merge(heme_metabolite4[, c("met", "biochemical_name")],
                 heme_metabolite_HR[, c("met", "per")])
met_per$per <- met_per$per/100

met_gp2 = rowAnnotation(beta_t2d = met_beta_t2d$colgp,
                        beta_heme = met_beta_heme$colgp,
                        porportion = anno_numeric(met_per$per, 
                                                  rg = c(0,0.2),
                                                  labels_gp = gpar(col = "black", fontsize = 14),
                                                  labels_format = function(x) paste0(sprintf("%.1f", x*100), "%"),
                                                  bg_gp = gpar(fill = "#bd002633", col = 0)),
                        show_annotation_name = F,
                        show_legend = F,
                        col = list(beta_heme = c("Positive" = "#E66101d0",
                                                 "Negative" = "#5E3C99d0"),
                                   beta_t2d = c("Positive" = "#D01C8Bd0",
                                                "Negative" = "#4DAC26d0")),
                        annotation_legend_param = list(
                          beta_heme = list(title = "Heme iron-metabolite\nassociation",
                                           title_gp = gpar(col = "black", fontsize = 16),
                                           labels_gp = gpar(col = "black", fontsize = 16),
                                           title_position = "topleft"),
                          beta_t2d = list(title = "Metabolite-T2D\nassociation",
                                          title_gp = gpar(col = "black", fontsize = 16),
                                          labels_gp = gpar(col = "black", fontsize = 16),
                                          title_position = "topleft")
                        ))

heat2 <- Heatmap(as.matrix(hh2), 
                 col = col_fun1,
                 row_dend_width = unit(2, "cm"),
                 row_dend_gp = gpar(lwd = 3),
                 show_row_names = T,
                 row_names_gp = gpar(fontsize = 16),
                 row_names_max_width = unit(25, "cm"),
                 column_labels = c("C-peptide", "HbA1c", "Total cholesterol",
                                   "HDL-C", "LDL-C", "TAG", "TAG/HDL-C ratio",
                                   "CRP", "Adiponectin", "Leptin", "Ferritin", "TfR/Ferritin ratio"),
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 16),
                 column_names_rot = 60,
                 column_title = NULL,
                 border = TRUE,
                 border_gp = gpar(col = "black", lwd = 3),
                 cluster_columns = FALSE,
                 cell_fun = function(j, i, x, y, w, h, fill){
                   if(hhp2[i, j] < 0.05) {
                     grid.text("*", x, y, vjust = 0.8,
                               gp = gpar(fontsize = 24))
                   }},
                 show_heatmap_legend = F,
                 right_annotation = met_gp2)

pdf("/udd/nhfwa/iron/archive/figures/metabolite_biomarker_heatmap_2.pdf", 
    width = 7.1, height = 6.9, onefile = F) # Open a new pdf file
draw(heat2)
dev.off() # Close the file
