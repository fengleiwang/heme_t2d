library(tidyverse)
library(dplyr)

nhs1 <- read.csv("/udd/nhfwa/iron/pool/metabolomics/data_nh_score.csv", header = T) # 1468
nhs2 <- read.csv("/udd/nhfwa/iron/pool/metabolomics/data_n2_score.csv", header = T) # 4962
hpfs <- read.csv("/udd/nhfwa/iron/pool/metabolomics/data_hp_score.csv", header = T) # 2739

nhs1$cohort="nh"
nhs2$cohort="n2"
hpfs$cohort="hp"

pool=bind_rows(nhs1, nhs2)
pool=bind_rows(pool, hpfs)
pool <- pool %>% filter(canbase==0 & hrtbase==0 & strbase==0 & dbbase==0)

#####################################################
#     metabolomics score - T2D associations
#####################################################
pool$hp <- ifelse(pool$cohort=="hp",1,0)
pool$n2 <- ifelse(pool$cohort=="n2",1,0)

pool$pmhc2 <- ifelse(pool$pmhc %in% c(2,5), 1, 0)
pool$pmhc3 <- ifelse(pool$pmhc==3, 1, 0)
pool$pmhc4 <- ifelse(pool$pmhc==4, 1, 0)
pool$pmhc9 <- ifelse(pool$pmhc==9, 1, 0)

var <- c("t2dscore", "ageyr", "white", "fast", "mv", "asp", "hbpbase", "cholbase",
         "htnrx", "hchtx", "pmhc", "smk2", "smk3", "alc2", "alc3", "alc4", "bmigp2", "bmigp3", 
         "bmigp4", "act2", "act3", "act4", "act5", "cvdfh", "dbfh",
         "calorn_avg", "ahei_noal", "cacoall", "cohort", "endpoint")
check_na <- data.frame(var=var,
                       na=colSums(is.na(pool[, var])))

pool <- pool[!is.na(pool$calorn_avg), ]
pool$ahei_noal[is.na(pool$ahei_noal)] <- median(pool$ahei_noal, na.rm = T)

pool$t2dscore_std <- scale(pool$t2dscore)

library(survival)
cox_t2d <- coxph(Surv(tdb2, db2case) ~  t2dscore_std + white + fast + mv + hbpbase + cholbase +
                   pmhc2 + pmhc3 + pmhc4 + smk2 + smk3 + alc2 + alc3 + alc4 + 
                   bmigp2 + bmigp3 + bmigp4 + act2 + act3 + act4 + act5 + dbfh +
                   calorn_avg + ahei_noal + strata(cacoall, cohort, endpoint, ageyr), # no need to include pmhc9 (hp)
                 data = pool)
options(scipen = 100)
summary(cox_t2d) # HR per SD: 2.12 (1.81, 2.48), n=851

summary(coxph(Surv(tdb2, db2case) ~  t2dscore_std + white + fast + mv + hbpbase + cholbase +
                   pmhc2 + pmhc3 + pmhc4 + smk2 + smk3 + alc2 + alc3 + alc4 + 
                   bmigp2 + bmigp3 + bmigp4 + act2 + act3 + act4 + act5 + dbfh +
                   calorn_avg + ahei_noal + fmtsvg_avg + strata(cacoall, cohort, endpoint, ageyr), # no need to include pmhc9 (hp)
                 data = pool))

cox_t2d_con <- coxph(Surv(tdb2, db2case) ~  t2dscore_std + white + fast + mv + hbpbase + cholbase +
                       pmhc2 + pmhc3 + pmhc4 + smk2 + smk3 + alc2 + alc3 + alc4 + 
                       bmigp2 + bmigp3 + bmigp4 + act2 + act3 + act4 + act5 + dbfh +
                       calorn_avg + ahei_noal + strata(cohort, endpoint, ageyr), # no need to include pmhc9 (hp)
                     data = pool[pool$cacoall==2, ])
summary(cox_t2d_con) # HR per SD: 2.41 (1.89, 3.06), n=447

cox_t2d_case <- coxph(Surv(tdb2, db2case) ~  t2dscore_std + white + fast + mv + hbpbase + cholbase +
                        pmhc2 + pmhc3 + pmhc4 + smk2 + smk3 + alc2 + alc3 + alc4 + 
                        bmigp2 + bmigp3 + bmigp4 + act2 + act3 + act4 + act5 + dbfh +
                        calorn_avg + ahei_noal + strata(cohort, endpoint, ageyr), # no need to include pmhc9 (hp)
                      data = pool[pool$cacoall==1, ])
summary(cox_t2d_case) # HR per SD: 2.03 (1.61, 2.57), n=404

# write.csv(pool[, c("tdb2", "db2case", "t2dscore_std", "t2dscore", "ageyr", "white", "fast", "mv", "hbpbase", "cholbase",
#                    "pmhc2", "pmhc3", "pmhc4", "smk2", "smk3", "alc2", "alc3", "alc4",
#                    "bmigp2", "bmigp3", "bmigp4", "act2", "act3", "act4", "act5", "dbfh",
#                    "calorn_avg", "ahei_noal", "cacoall", "cohort", "endpoint" )],
#           "/udd/nhfwa/iron/pool/metabolomics/metscore_t2d_spline.csv", na="", row.names = F)

###################################################
#       iron-metabolomics score associations
###################################################
results_heme <- data.frame(type="heme",
                           n=NA,
                           beta=NA,
                           lci=NA,
                           uci=NA,
                           p=NA)
results_nonheme <- data.frame(type="nonheme",
                              n=NA,
                              beta=NA,
                              lci=NA,
                              uci=NA,
                              p=NA)
results_iron <- data.frame(type="iron",
                           n=NA,
                           beta=NA,
                           lci=NA,
                           uci=NA,
                           p=NA)
results_iron_dt <- data.frame(type="iron_dt",
                              n=NA,
                              beta=NA,
                              lci=NA,
                              uci=NA,
                              p=NA)
results_iron_sp <- data.frame(type="iron_sp",
                              n=NA,
                              beta=NA,
                              lci=NA,
                              uci=NA,
                              p=NA)

pool$heme5 <- Hmisc::cut2(pool$hemea_avg, g=5)
pool$nonheme5 <- Hmisc::cut2(pool$nonhemea_avg, g=5)
pool$iron5 <- Hmisc::cut2(pool$irona_avg, g=5)
pool$iron_dt5 <- Hmisc::cut2(pool$irona_dt_avg, g=5)
pool$iron_sp5 <- Hmisc::cut2(pool$irona_sp_avg, g=5)

quantile(pool$hemea_avg, probs = 0.9) - quantile(pool$hemea_avg, probs = 0.1) # 0.88
quantile(pool$nonhemea_avg, probs = 0.9) - quantile(pool$nonhemea_avg, probs = 0.1) # 23.76
quantile(pool$irona_avg, probs = 0.9) - quantile(pool$irona_avg, probs = 0.1) # 23.74
quantile(pool$irona_dt_avg, probs = 0.9) - quantile(pool$irona_dt_avg, probs = 0.1) # 8.69 
quantile(pool$irona_sp_avg, probs = 0.9) - quantile(pool$irona_sp_avg, probs = 0.1) # 18.78

pool$nonhemea_avg <- pool$nonhemea_avg/20
pool$irona_avg <- pool$irona_avg/20
pool$irona_sp_avg <- pool$irona_sp_avg/20
pool$irona_dt_avg <- pool$irona_dt_avg/10

# heme iron
tempmod1 <- lm(t2dscore_std~hemea_avg+nonhemea_avg+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
                 hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
                 calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg, data=pool)
tempsum1 <- summary(tempmod1)

tempmod11 <- lm(t2dscore_std~heme5+nonheme5+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
                  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
                  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg, data=pool)
tempsum11 <- summary(tempmod11)

results_heme[1, "n"] <- dim(pool)[1]
results_heme[1, "beta"] <- tempsum1$coefficients[2,1]
results_heme[1, "lci"] <- confint(tempmod1)[2,1]
results_heme[1, "uci"] <- confint(tempmod1)[2,2]
results_heme[1, "p"] <- tempsum1$coefficients[2,4]

results_heme[2:5, "n"] <- rownames(tempsum11$coefficients)[2:5]
results_heme[2:5, "beta"] <- tempsum11$coefficients[2:5,1]
results_heme[2:5, "lci"] <- confint(tempmod11)[2:5,1]
results_heme[2:5, "uci"] <- confint(tempmod11)[2:5,2]
results_heme[2:5, "p"] <- tempsum11$coefficients[2:5,4]

# nonheme iron
tempmod2 <- lm(t2dscore_std~nonhemea_avg+hemea_avg+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
                 hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
                 calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg, data=pool)
tempsum2 <- summary(tempmod2)

tempmod22 <- lm(t2dscore_std~nonheme5+heme5+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
                  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
                  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg, data=pool)
tempsum22 <- summary(tempmod22)

results_nonheme[1, "n"] <- dim(pool)[1]
results_nonheme[1, "beta"] <- tempsum2$coefficients[2,1]
results_nonheme[1, "lci"] <- confint(tempmod2)[2,1]
results_nonheme[1, "uci"] <- confint(tempmod2)[2,2]
results_nonheme[1, "p"] <- tempsum2$coefficients[2,4]

results_nonheme[2:5, "n"] <- rownames(tempsum22$coefficients)[2:5]
results_nonheme[2:5, "beta"] <- tempsum22$coefficients[2:5,1]
results_nonheme[2:5, "lci"] <- confint(tempmod22)[2:5,1]
results_nonheme[2:5, "uci"] <- confint(tempmod22)[2:5,2]
results_nonheme[2:5, "p"] <- tempsum22$coefficients[2:5,4]

# total iron
tempmod3 <- lm(t2dscore_std~irona_avg+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
                 hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
                 calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg, data=pool)
tempsum3 <- summary(tempmod3)

tempmod33 <- lm(t2dscore_std~iron5+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
                  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
                  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg, data=pool)
tempsum33 <- summary(tempmod33)

results_iron[1, "n"] <- dim(pool)[1]
results_iron[1, "beta"] <- tempsum3$coefficients[2,1]
results_iron[1, "lci"] <- confint(tempmod3)[2,1]
results_iron[1, "uci"] <- confint(tempmod3)[2,2]
results_iron[1, "p"] <- tempsum3$coefficients[2,4]

results_iron[2:5, "n"] <- rownames(tempsum33$coefficients)[2:5]
results_iron[2:5, "beta"] <- tempsum33$coefficients[2:5,1]
results_iron[2:5, "lci"] <- confint(tempmod33)[2:5,1]
results_iron[2:5, "uci"] <- confint(tempmod33)[2:5,2]
results_iron[2:5, "p"] <- tempsum33$coefficients[2:5,4]

# dietary iron
tempmod4 <- lm(t2dscore_std~irona_dt_avg+irona_sp_avg+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
                 hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
                 calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg, data=pool)
tempsum4 <- summary(tempmod4)

tempmod44 <- lm(t2dscore_std~iron_dt5+iron_sp5+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
                  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
                  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg, data=pool)
tempsum44 <- summary(tempmod44)

results_iron_dt[1, "n"] <- dim(pool)[1]
results_iron_dt[1, "beta"] <- tempsum4$coefficients[2,1]
results_iron_dt[1, "lci"] <- confint(tempmod4)[2,1]
results_iron_dt[1, "uci"] <- confint(tempmod4)[2,2]
results_iron_dt[1, "p"] <- tempsum4$coefficients[2,4]

results_iron_dt[2:5, "n"] <- rownames(tempsum44$coefficients)[2:5]
results_iron_dt[2:5, "beta"] <- tempsum44$coefficients[2:5,1]
results_iron_dt[2:5, "lci"] <- confint(tempmod44)[2:5,1]
results_iron_dt[2:5, "uci"] <- confint(tempmod44)[2:5,2]
results_iron_dt[2:5, "p"] <- tempsum44$coefficients[2:5,4]

# supp iron
tempmod5 <- lm(t2dscore_std~irona_sp_avg+irona_dt_avg+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
                 hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
                 calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg, data=pool)
tempsum5 <- summary(tempmod5)

tempmod55 <- lm(t2dscore_std~iron_sp5+iron_dt5+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
                  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
                  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg, data=pool)
tempsum55 <- summary(tempmod55)

results_iron_sp[1, "n"] <- dim(pool)[1]
results_iron_sp[1, "beta"] <- tempsum5$coefficients[2,1]
results_iron_sp[1, "lci"] <- confint(tempmod5)[2,1]
results_iron_sp[1, "uci"] <- confint(tempmod5)[2,2]
results_iron_sp[1, "p"] <- tempsum5$coefficients[2,4]

results_iron_sp[2:5, "n"] <- rownames(tempsum55$coefficients)[2:5]
results_iron_sp[2:5, "beta"] <- tempsum55$coefficients[2:5,1]
results_iron_sp[2:5, "lci"] <- confint(tempmod55)[2:5,1]
results_iron_sp[2:5, "uci"] <- confint(tempmod55)[2:5,2]
results_iron_sp[2:5, "p"] <- tempsum55$coefficients[2:5,4]

results <- rbind(results_heme,
                 results_nonheme)
results <- rbind(results,
                 results_iron)
results <- rbind(results,
                 results_iron_dt)
results <- rbind(results,
                 results_iron_sp)

save.image(file="/udd/nhfwa/iron/pool/metabolomics/metscore_t2d_v1.RData")

#################################################################
#         metabolite-T2D, heme-metabolite associations
#################################################################
load("~/iron/pool/metabolomics/score_t2d/t2d_merge.RData")
gdata::keep(pool, met_coef, f, sure=T)

metabolite_t2d <- data.frame(met=met_coef$met,
                             name1=met_coef$metabolite_name,
                             name2=met_coef$biochemical_name,
                             beta=NA,
                             se=NA,
                             HR=NA,
                             p=NA)

heme_metabolite <- data.frame(met=met_coef$met,
                              name1=met_coef$metabolite_name,
                              name2=met_coef$biochemical_name,
                              beta=NA,
                              se=NA,
                              p=NA)

met_t2d <- list()
for (i in 1:length(met_coef$met)) {
  tempmet <- met_coef$met[i]
  pool[, tempmet] <- scale(pool[, tempmet])
  tempcox <- coxph(as.formula(paste0("Surv(tdb2, db2case) ~  ", tempmet, " + ageyr + white + fast + 
                                     mv + hbpbase + cholbase + pmhc2 + pmhc3 + pmhc4 + smk2 + smk3 + 
                                     alc2 + alc3 + alc4 + bmigp2 + bmigp3 + bmigp4 + act2 + act3 + 
                                     act4 + act5 + dbfh + calorn_avg + ahei_noal + 
                                     strata(cacoall, cohort, endpoint)")),
                   data = pool)
  met_t2d[[i]] <- tempcox
  tempsum <- summary(tempcox)
  hr_data <- tempsum$conf.int[1,c("exp(coef)", "lower .95", "upper .95")]
  metabolite_t2d[i, ]$beta <- tempsum$coefficients[1, 1]
  metabolite_t2d[i, ]$se <- tempsum$coefficients[1, 3]
  metabolite_t2d[i, ]$HR <- paste(format(round(hr_data[1], 2), nsmall=2), " (",
                               format(round(hr_data[2], 2), nsmall=2), ", ",
                               format(round(hr_data[3], 2), nsmall=2), ")", sep="")
  metabolite_t2d[i, ]$p <- tempsum$coefficients[1, 5]
  
  tempmod <- lm(paste0(tempmet, "~hemea_avg+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
                   hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
                   calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=pool)
  heme_metabolite[i, "beta"] <- summary(tempmod)$coefficients[2,1]
  heme_metabolite[i, "se"] <- summary(tempmod)$coefficients[2,2]
  heme_metabolite[i, "p"] <- summary(tempmod)$coefficients[2,4]
}

metabolite_t2d$fdr <- p.adjust(metabolite_t2d$p, method="fdr", length(metabolite_t2d$p))
metabolite_t2d$bon <- p.adjust(metabolite_t2d$p, method="bonferroni", length(metabolite_t2d$p))
sum(metabolite_t2d$bon<0.05) # 124

heme_metabolite2 <- heme_metabolite %>% 
  filter(met %in% metabolite_t2d[metabolite_t2d$bon<0.05, ]$met)
heme_metabolite2$fdr <- p.adjust(heme_metabolite2$p, method="fdr", length(heme_metabolite2$p))
heme_metabolite2$bon <- p.adjust(heme_metabolite2$p, method="bonferroni", length(heme_metabolite2$p))
sum(heme_metabolite2$bon<0.05) # 33

heme_metabolite3 <- heme_metabolite2[heme_metabolite2$bon<0.05, 1:4]
colnames(heme_metabolite3)[4] <- "beta_heme"
heme_metabolite3 <- merge(heme_metabolite3, metabolite_t2d[, c("met", "beta")], by="met")
heme_metabolite3$mediator <- ifelse(heme_metabolite3$beta_heme*heme_metabolite3$beta>0,
                                    1,
                                    0)
heme_metabolite4 <- heme_metabolite3[heme_metabolite3$mediator==1, ] # 17
save.image(file="/udd/nhfwa/iron/pool/metabolomics/metscore_t2d_v2.RData")

metname <- read_csv(file="/udd/nhfwa/iron/pool/metabolomics/metabolite_name.csv")
table_s5_1 <- merge(metname[, c("met", "biochemical_name")],
                    metabolite_t2d[, c(1,4:7,9)], by="met")
table_s5_2 <- merge(metname[, c("met", "biochemical_name")],
                    heme_metabolite2[, c(1,4,5,6,8)], by="met")

table_s5 <- merge(table_s5_1, table_s5_2[, -1], by="biochemical_name", all.x = T)
gdata::keep(table_s5, sure=T)

table_s5$beta.x <- format(round(table_s5$beta.x, 4), nsmall=4)
table_s5$se.x <- format(round(table_s5$se.x, 4), nsmall=4)
table_s5$beta.y <- format(round(table_s5$beta.y, 4), nsmall=4)
table_s5$se.y <- format(round(table_s5$se.y, 4), nsmall=4)

table_s5 <- table_s5 %>% 
  mutate(p.x=case_when(p.x<0.001 ~ "<0.001",
                       p.x>=0.001 & p.x<0.01 ~ format(round(p.x, 3), nsmall=3),
                       TRUE ~ format(round(p.x, 2), nsmall=2)),
         bon.x=case_when(bon.x==1 ~ "0.99",
                         bon.x<0.001 ~ "<0.001",
                         bon.x>=0.001 & bon.x<0.01 ~ format(round(bon.x, 3), nsmall=3),
                         TRUE ~ format(round(bon.x, 2), nsmall=2)),
         p.y=case_when(p.y<0.001 ~ "<0.001",
                       p.y>=0.001 & p.y<0.01 ~ format(round(p.y, 3), nsmall=3),
                       TRUE ~ format(round(p.y, 2), nsmall=2)),
         bon.y=case_when(bon.y==1 ~ "0.99",
                         bon.y<0.001 ~ "<0.001",
                         bon.y>=0.001 & bon.y<0.01 ~ format(round(bon.y, 3), nsmall=3),
                         TRUE ~ format(round(bon.y, 2), nsmall=2)))
write.csv(table_s5[, c(2,1,3:11)], file="/udd/nhfwa/iron/pool/metabolomics/table_s5.csv", row.names = F)

#################################################################
#         heme-T2D, adjusting for metabolites
#################################################################
met17 <- heme_metabolite4$met
hemecox <- coxph(Surv(tdb2, db2case) ~  hemea_avg + ageyr + white + fast + 
                   mv + hbpbase + cholbase + pmhc2 + pmhc3 + pmhc4 + smk2 + smk3 + 
                   alc2 + alc3 + alc4 + bmigp2 + bmigp3 + bmigp4 + act2 + act3 + 
                   act4 + act5 + dbfh + calorn_avg + magnn_avg + gln_avg + 
                   cerafn_avg + psn_avg + transn_avg +
                   strata(cacoall, cohort, endpoint), data = pool)
summary(hemecox) # 1.53 (1.23, 1.91)

heme_metabolite_HR <- data.frame(met=met17,
                                 name1=heme_metabolite4$name1,
                                 name2=heme_metabolite4$name2,
                                 HR=NA,
                                 hh=NA,
                                 lci=NA,
                                 uci=NA,
                                 p=NA)

for (i in 1:length(met17)) {
  tempmet <- met17[i]
  # pool[, tempmet] <- scale(pool[, tempmet])
  tempcox <- coxph(as.formula(paste0("Surv(tdb2, db2case) ~ hemea_avg + ", 
                   tempmet, " + ageyr + white + fast + 
                   mv + hbpbase + cholbase + pmhc2 + pmhc3 + pmhc4 + smk2 + smk3 + 
                   alc2 + alc3 + alc4 + bmigp2 + bmigp3 + bmigp4 + act2 + act3 + 
                   act4 + act5 + dbfh + calorn_avg + magnn_avg + gln_avg + 
                   cerafn_avg + psn_avg + transn_avg + 
                   strata(cacoall, cohort, endpoint)")),
                   data = pool)
  tempsum <- summary(tempcox)
  hr_data <- tempsum$conf.int[1,c("exp(coef)", "lower .95", "upper .95")]
  heme_metabolite_HR[i, ]$HR <- paste(format(round(hr_data[1], 2), nsmall=2), " (",
                                  format(round(hr_data[2], 2), nsmall=2), ", ",
                                  format(round(hr_data[3], 2), nsmall=2), ")", sep="")
  heme_metabolite_HR[i, ]$hh <- format(round(hr_data[1], 2), nsmall=2)
  heme_metabolite_HR[i, ]$lci <- format(round(hr_data[2], 2), nsmall=2)
  heme_metabolite_HR[i, ]$uci <- format(round(hr_data[3], 2), nsmall=2)
  heme_metabolite_HR[i, ]$p <- tempsum$coefficients[1, 5]
}  

write.csv(pool[, c("newid", "tdb2", "db2case", "hemea_avg", met17,
                   "white", "fast", "mv", "hbpbase", "cholbase",
                   "pmhc2", "pmhc3", "pmhc4", "smk2", "smk3", "alc2", "alc3", "alc4",
                   "bmigp2", "bmigp3", "bmigp4", "act2", "act3", "act4", "act5", "dbfh",
                   "calorn_avg", "magnn_avg", "gln_avg", "cerafn_avg", "psn_avg", "transn_avg",
                   "cacoall", "cohort", "endpoint", "ageyr")],
          "/udd/nhfwa/iron/pool/metabolomics/metscore_t2d_mediation.csv", na="", row.names = F)

save.image(file="/udd/nhfwa/iron/pool/metabolomics/metscore_t2d_v3.RData")

heme_metabolite_HR$per <- c(10.8, 9.9, 12.1, 17.0, 5.4, 5.1, 12.7, 11.9, 10.0,
                            11.5, 3.6, 8.3, 5.1, 8.5, 8.8, 4.5, 10.2)
metname <- read_csv(file="/udd/nhfwa/iron/pool/metabolomics/metabolite_name.csv")
heme_metabolite_HR <- merge(metname[, c("met", "biochemical_name")],
                            heme_metabolite_HR[, c(1,4:9)], by="met") %>% arrange(hh)
heme_metabolite_HR$biochemical_name <- paste0("  + ", heme_metabolite_HR$biochemical_name)
heme_t2d <- data.frame(met=NA,
                       biochemical_name="Heme iron",
                       HR="1.53 (1.23, 1.91)",
                       hh=1.53,
                       lci=1.23,
                       uci=1.91,
                       p=NA,
                       per=0)
heme_metabolite_HR <- rbind(heme_t2d, heme_metabolite_HR)

gdata::keep(heme_metabolite_HR, sure=T)
save.image(file="/udd/nhfwa/iron/pool/metabolomics/metscore_t2d_mediation.RData")

heme_metabolite_HR$biochemical_name <- factor(heme_metabolite_HR$biochemical_name, 
                                              levels=rev(heme_metabolite_HR$biochemical_name))

heme_metabolite_HR$hh <- as.numeric(heme_metabolite_HR$hh)
heme_metabolite_HR$lci <- as.numeric(heme_metabolite_HR$lci)
heme_metabolite_HR$uci <- as.numeric(heme_metabolite_HR$uci)

library(ggplot2)
options(bitmapType='cairo')

pdf("/udd/nhfwa/iron/archive/figures/figure_s2b.pdf", 
    width = 6, height = 7, onefile = F) # Open a new pdf file
ggplot(data=heme_metabolite_HR, aes(x=hh, y=biochemical_name, xmin=lci, xmax=uci)) +
  geom_point(size=5, shape=18) +
  geom_errorbar(width=0.3, size=1) +
  geom_vline(xintercept=1, size=0.9, lty=2) +  # add a dotted line at x=1 after flip
  geom_vline(xintercept=1.53, size=0.9, lty=2) +  # add a dotted line at x=1 after flip
  scale_x_continuous(breaks=seq(1, 1.9, by=0.3), limits=c(0.95, 3)) +
  xlab("\nHR (95% CI)") +
  theme_bw() +  # use a white background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(colour = "black", linewidth = 1),
        axis.ticks.length.x=unit(.2, "cm"),
        plot.title = element_text(size = 12, hjust = 0.5, color="black"),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, vjust=-1, color="black"), 
        axis.text.y = element_text(size = 12, hjust=0, color="black"))
dev.off() # Close the file
