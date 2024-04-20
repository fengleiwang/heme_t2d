library(tidyverse)
library(dplyr)

nhs1 <- read.csv("/udd/nhfwa/iron/pool/biomarker/data_nh_bio.csv", header = T)
nhs2 <- read.csv("/udd/nhfwa/iron/pool/biomarker/data_n2_bio.csv", header = T)
hpfs <- read.csv("/udd/nhfwa/iron/pool/biomarker/data_hp_bio.csv", header = T)

nhs1$cohort="nh"
nhs2$cohort="n2"
hpfs$cohort="hp"

pool=bind_rows(nhs1, nhs2)
pool=bind_rows(pool, hpfs)
pool <- pool %>% filter(canbase==0 & hrtbase==0 & strbase==0 & dbbase==0)

pool$adj_trf <- pool$adj_stfr/pool$adj_ferritin
pool$adj_lntrf <- log(pool$adj_trf)

##############################################
#       iron-biomarker associations
##############################################
pool$hp <- ifelse(pool$cohort=="hp",1,0)
pool$n2 <- ifelse(pool$cohort=="n2",1,0)

pool$pmhc2 <- ifelse(pool$pmhc %in% c(2,5), 1, 0)
pool$pmhc3 <- ifelse(pool$pmhc==3, 1, 0)
pool$pmhc4 <- ifelse(pool$pmhc==4, 1, 0)
pool$pmhc9 <- ifelse(pool$pmhc==9, 1, 0)

bio <- c("adj_lncpep", "adj_lna1c", "adj_lnchol",  "adj_lnhdl",  "adj_lnldl",  "adj_lntrig", "adj_lnratio",
         "adj_lncrp", "adj_lnadipo", "adj_lnleptin", "adj_lnferritin", "adj_lntrf") #, "adj_lntransf", "adj_lnstfr")

results_heme <- data.frame(type="heme",
                           bio=bio,
                           n=NA,
                           beta=NA,
                           lci=NA,
                           uci=NA,
                           p=NA)
results_nonheme <- data.frame(type="nonheme",
                              bio=bio,
                              n=NA,
                              beta=NA,
                              lci=NA,
                              uci=NA,
                              p=NA)
results_iron <- data.frame(type="iron",
                           bio=bio,
                           n=NA,
                           beta=NA,
                           lci=NA,
                           uci=NA,
                           p=NA)
results_iron_dt <- data.frame(type="iron_dt",
                           bio=bio,
                           n=NA,
                           beta=NA,
                           lci=NA,
                           uci=NA,
                           p=NA)
results_iron_sp <- data.frame(type="iron_sp",
                              bio=bio,
                              n=NA,
                              beta=NA,
                              lci=NA,
                              uci=NA,
                              p=NA)

quantile(pool$hemea_avg, probs = 0.9) - quantile(pool$hemea_avg, probs = 0.1) # 0.8972924
quantile(pool$nonhemea_avg, probs = 0.9) - quantile(pool$nonhemea_avg, probs = 0.1) # 23.91181
quantile(pool$irona_avg, probs = 0.9) - quantile(pool$irona_avg, probs = 0.1) # 23.8401
quantile(pool$irona_dt_avg, probs = 0.9) - quantile(pool$irona_dt_avg, probs = 0.1) # 8.635615 
quantile(pool$irona_sp_avg, probs = 0.9) - quantile(pool$irona_sp_avg, probs = 0.1) # 18.6159 

for (i in 1:12) {
  tempdata <- pool %>% filter(!is.na(hemea_avg))
  tempdata <- tempdata[!is.na(tempdata[,bio[i]]),]
  tempdata$nonhemea_avg <- tempdata$nonhemea_avg/20
  tempdata$irona_avg <- tempdata$irona_avg/20
  tempdata$irona_sp_avg <- tempdata$irona_sp_avg/20
  tempdata$irona_dt_avg <- tempdata$irona_dt_avg/10
  
  # heme iron
  tempmod1 <- lm(paste0(bio[i],"~hemea_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum1 <- summary(tempmod1)
  results_heme[i, "n"] <- dim(tempdata)[1]
  results_heme[i, "beta"] <- (exp(tempsum1$coefficients[2,1])-1)*100
  results_heme[i, "lci"] <- (exp(confint(tempmod1)[2,1])-1)*100
  results_heme[i, "uci"] <- (exp(confint(tempmod1)[2,2])-1)*100
  results_heme[i, "p"] <- tempsum1$coefficients[2,4]

  # nonheme iron
  tempmod2 <- lm(paste0(bio[i],"~nonhemea_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum2 <- summary(tempmod2)
  results_nonheme[i, "n"] <- dim(tempdata)[1]
  results_nonheme[i, "beta"] <- (exp(tempsum2$coefficients[2,1])-1)*100
  results_nonheme[i, "lci"] <- (exp(confint(tempmod2)[2,1])-1)*100
  results_nonheme[i, "uci"] <- (exp(confint(tempmod2)[2,2])-1)*100
  results_nonheme[i, "p"] <- tempsum2$coefficients[2,4]
  
  # total iron
  tempmod3 <- lm(paste0(bio[i],"~irona_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum3 <- summary(tempmod3)
  results_iron[i, "n"] <- dim(tempdata)[1]
  results_iron[i, "beta"] <- (exp(tempsum3$coefficients[2,1])-1)*100
  results_iron[i, "lci"] <- (exp(confint(tempmod3)[2,1])-1)*100
  results_iron[i, "uci"] <- (exp(confint(tempmod3)[2,2])-1)*100
  results_iron[i, "p"] <- tempsum3$coefficients[2,4]
  
  # dietary iron
  tempmod4 <- lm(paste0(bio[i],"~irona_dt_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum4 <- summary(tempmod4)
  results_iron_dt[i, "n"] <- dim(tempdata)[1]
  results_iron_dt[i, "beta"] <- (exp(tempsum4$coefficients[2,1])-1)*100
  results_iron_dt[i, "lci"] <- (exp(confint(tempmod4)[2,1])-1)*100
  results_iron_dt[i, "uci"] <- (exp(confint(tempmod4)[2,2])-1)*100
  results_iron_dt[i, "p"] <- tempsum4$coefficients[2,4]
  
  # supp iron
  tempmod5 <- lm(paste0(bio[i],"~irona_sp_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum5 <- summary(tempmod5)
  results_iron_sp[i, "n"] <- dim(tempdata)[1]
  results_iron_sp[i, "beta"] <- (exp(tempsum5$coefficients[2,1])-1)*100
  results_iron_sp[i, "lci"] <- (exp(confint(tempmod5)[2,1])-1)*100
  results_iron_sp[i, "uci"] <- (exp(confint(tempmod5)[2,2])-1)*100
  results_iron_sp[i, "p"] <- tempsum5$coefficients[2,4]
  
}

results <- rbind(results_heme,
                 results_nonheme)
results <- rbind(results,
                 results_iron)
results <- rbind(results,
                 results_iron_dt)
results <- rbind(results,
                 results_iron_sp)

results$bon <- results$p*12
results$per <- paste0(format(round(results$beta,digits=1),digits=1),
                      " (",
                      format(round(results$lci,digits=1),digits=1),
                      ", ",
                      format(round(results$uci,digits=1),digits=1),
                      ")")
sum(results$bon<0.05) # 10
View(results[results$bon<0.05,])

##############################################
#   iron-biomarker associations (control)
##############################################
pool2 <- pool[pool$caco==0,]

results_heme2 <- results_heme
results_nonheme2 <- results_nonheme
results_iron2 <- results_iron
results_iron_dt2 <- results_iron_dt
results_iron_sp2 <- results_iron_sp

for (i in 1:12) {
  tempdata <- pool2 %>% filter(!is.na(hemea_avg))
  tempdata <- tempdata[!is.na(tempdata[,bio[i]]),]
  tempdata$nonhemea_avg <- tempdata$nonhemea_avg/20
  tempdata$irona_avg <- tempdata$irona_avg/20
  tempdata$irona_sp_avg <- tempdata$irona_sp_avg/20
  tempdata$irona_dt_avg <- tempdata$irona_dt_avg/10
  
  # heme iron
  tempmod1 <- lm(paste0(bio[i],"~hemea_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum1 <- summary(tempmod1)
  results_heme2[i, "n"] <- dim(tempdata)[1]
  results_heme2[i, "beta"] <- (exp(tempsum1$coefficients[2,1])-1)*100
  results_heme2[i, "lci"] <- (exp(confint(tempmod1)[2,1])-1)*100
  results_heme2[i, "uci"] <- (exp(confint(tempmod1)[2,2])-1)*100
  results_heme2[i, "p"] <- tempsum1$coefficients[2,4]
  
  # nonheme iron
  tempmod2 <- lm(paste0(bio[i],"~nonhemea_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum2 <- summary(tempmod2)
  results_nonheme2[i, "n"] <- dim(tempdata)[1]
  results_nonheme2[i, "beta"] <- (exp(tempsum2$coefficients[2,1])-1)*100
  results_nonheme2[i, "lci"] <- (exp(confint(tempmod2)[2,1])-1)*100
  results_nonheme2[i, "uci"] <- (exp(confint(tempmod2)[2,2])-1)*100
  results_nonheme2[i, "p"] <- tempsum2$coefficients[2,4]
  
  # total iron
  tempmod3 <- lm(paste0(bio[i],"~irona_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum3 <- summary(tempmod3)
  results_iron2[i, "n"] <- dim(tempdata)[1]
  results_iron2[i, "beta"] <- (exp(tempsum3$coefficients[2,1])-1)*100
  results_iron2[i, "lci"] <- (exp(confint(tempmod3)[2,1])-1)*100
  results_iron2[i, "uci"] <- (exp(confint(tempmod3)[2,2])-1)*100
  results_iron2[i, "p"] <- tempsum3$coefficients[2,4]
  
  # dietary iron
  tempmod4 <- lm(paste0(bio[i],"~irona_dt_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum4 <- summary(tempmod4)
  results_iron_dt2[i, "n"] <- dim(tempdata)[1]
  results_iron_dt2[i, "beta"] <- (exp(tempsum4$coefficients[2,1])-1)*100
  results_iron_dt2[i, "lci"] <- (exp(confint(tempmod4)[2,1])-1)*100
  results_iron_dt2[i, "uci"] <- (exp(confint(tempmod4)[2,2])-1)*100
  results_iron_dt2[i, "p"] <- tempsum4$coefficients[2,4]
  
  # supp iron
  tempmod5 <- lm(paste0(bio[i],"~irona_sp_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum5 <- summary(tempmod5)
  results_iron_sp2[i, "n"] <- dim(tempdata)[1]
  results_iron_sp2[i, "beta"] <- (exp(tempsum5$coefficients[2,1])-1)*100
  results_iron_sp2[i, "lci"] <- (exp(confint(tempmod5)[2,1])-1)*100
  results_iron_sp2[i, "uci"] <- (exp(confint(tempmod5)[2,2])-1)*100
  results_iron_sp2[i, "p"] <- tempsum5$coefficients[2,4]
  
}

results_con <- rbind(results_heme2,
                     results_nonheme2)
results_con <- rbind(results_con,
                     results_iron2)
results_con <- rbind(results_con,
                     results_iron_dt2)
results_con <- rbind(results_con,
                     results_iron_sp2)

results_con$bon <- results_con$p*12
results_con$per <- paste0(format(round(results_con$beta,digits=1),digits=1),
                      " (",
                      format(round(results_con$lci,digits=1),digits=1),
                      ", ",
                      format(round(results_con$uci,digits=1),digits=1),
                      ")")
sum(results_con$bon<0.05) # 7
View(results_con[results_con$bon<0.05,])

##############################################
#   iron-biomarker associations (fast)
##############################################
pool3 <- pool[pool$fast==1,]

results_heme3 <- results_heme
results_nonheme3 <- results_nonheme
results_iron3 <- results_iron
results_iron_dt3 <- results_iron_dt
results_iron_sp3 <- results_iron_sp

for (i in 1:12) {
  tempdata <- pool3 %>% filter(!is.na(hemea_avg))
  tempdata <- tempdata[!is.na(tempdata[,bio[i]]),]
  tempdata$nonhemea_avg <- tempdata$nonhemea_avg/20
  tempdata$irona_avg <- tempdata$irona_avg/20
  tempdata$irona_sp_avg <- tempdata$irona_sp_avg/20
  tempdata$irona_dt_avg <- tempdata$irona_dt_avg/10
  
  # heme iron
  tempmod1 <- lm(paste0(bio[i],"~hemea_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum1 <- summary(tempmod1)
  results_heme3[i, "n"] <- dim(tempdata)[1]
  results_heme3[i, "beta"] <- (exp(tempsum1$coefficients[2,1])-1)*100
  results_heme3[i, "lci"] <- (exp(confint(tempmod1)[2,1])-1)*100
  results_heme3[i, "uci"] <- (exp(confint(tempmod1)[2,2])-1)*100
  results_heme3[i, "p"] <- tempsum1$coefficients[2,4]
  
  # nonheme iron
  tempmod2 <- lm(paste0(bio[i],"~nonhemea_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum2 <- summary(tempmod2)
  results_nonheme3[i, "n"] <- dim(tempdata)[1]
  results_nonheme3[i, "beta"] <- (exp(tempsum2$coefficients[2,1])-1)*100
  results_nonheme3[i, "lci"] <- (exp(confint(tempmod2)[2,1])-1)*100
  results_nonheme3[i, "uci"] <- (exp(confint(tempmod2)[2,2])-1)*100
  results_nonheme3[i, "p"] <- tempsum2$coefficients[2,4]
  
  # total iron
  tempmod3 <- lm(paste0(bio[i],"~irona_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum3 <- summary(tempmod3)
  results_iron3[i, "n"] <- dim(tempdata)[1]
  results_iron3[i, "beta"] <- (exp(tempsum3$coefficients[2,1])-1)*100
  results_iron3[i, "lci"] <- (exp(confint(tempmod3)[2,1])-1)*100
  results_iron3[i, "uci"] <- (exp(confint(tempmod3)[2,2])-1)*100
  results_iron3[i, "p"] <- tempsum3$coefficients[2,4]
 
  # dietary iron
  tempmod4 <- lm(paste0(bio[i],"~irona_dt_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum4 <- summary(tempmod4)
  results_iron_dt3[i, "n"] <- dim(tempdata)[1]
  results_iron_dt3[i, "beta"] <- (exp(tempsum4$coefficients[2,1])-1)*100
  results_iron_dt3[i, "lci"] <- (exp(confint(tempmod4)[2,1])-1)*100
  results_iron_dt3[i, "uci"] <- (exp(confint(tempmod4)[2,2])-1)*100
  results_iron_dt3[i, "p"] <- tempsum4$coefficients[2,4]
  
  # supp iron
  tempmod5 <- lm(paste0(bio[i],"~irona_sp_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum5 <- summary(tempmod5)
  results_iron_sp3[i, "n"] <- dim(tempdata)[1]
  results_iron_sp3[i, "beta"] <- (exp(tempsum5$coefficients[2,1])-1)*100
  results_iron_sp3[i, "lci"] <- (exp(confint(tempmod5)[2,1])-1)*100
  results_iron_sp3[i, "uci"] <- (exp(confint(tempmod5)[2,2])-1)*100
  results_iron_sp3[i, "p"] <- tempsum5$coefficients[2,4] 
}

results_fast <- rbind(results_heme3,
                      results_nonheme3)
results_fast <- rbind(results_fast,
                      results_iron3)
results_fast <- rbind(results_fast,
                      results_iron_dt3)
results_fast <- rbind(results_fast,
                      results_iron_sp3)

results_fast$bon <- results_fast$p*12
results_fast$per <- paste0(format(round(results_fast$beta,digits=1),digits=1),
                      " (",
                      format(round(results_fast$lci,digits=1),digits=1),
                      ", ",
                      format(round(results_fast$uci,digits=1),digits=1),
                      ")")
sum(results_fast$bon<0.05) # 8
View(results_fast[results_fast$bon<0.05,])

##############################################
#   iron-biomarker associations (healthy)
##############################################
pool4 <- pool[pool$cholbase==0 & pool$hbpbase==0,]

results_heme4 <- results_heme
results_nonheme4 <- results_nonheme
results_iron4 <- results_iron
results_iron_dt4 <- results_iron_dt
results_iron_sp4 <- results_iron_sp

for (i in 1:12) {
  tempdata <- pool4 %>% filter(!is.na(hemea_avg))
  tempdata <- tempdata[!is.na(tempdata[,bio[i]]),]
  tempdata$nonhemea_avg <- tempdata$nonhemea_avg/20
  tempdata$irona_avg <- tempdata$irona_avg/20
  tempdata$irona_sp_avg <- tempdata$irona_sp_avg/20
  tempdata$irona_dt_avg <- tempdata$irona_dt_avg/10
  
  # heme iron
  tempmod1 <- lm(paste0(bio[i],"~hemea_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum1 <- summary(tempmod1)
  results_heme4[i, "n"] <- dim(tempdata)[1]
  results_heme4[i, "beta"] <- (exp(tempsum1$coefficients[2,1])-1)*100
  results_heme4[i, "lci"] <- (exp(confint(tempmod1)[2,1])-1)*100
  results_heme4[i, "uci"] <- (exp(confint(tempmod1)[2,2])-1)*100
  results_heme4[i, "p"] <- tempsum1$coefficients[2,4]
  
  # nonheme iron
  tempmod2 <- lm(paste0(bio[i],"~nonhemea_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum2 <- summary(tempmod2)
  results_nonheme4[i, "n"] <- dim(tempdata)[1]
  results_nonheme4[i, "beta"] <- (exp(tempsum2$coefficients[2,1])-1)*100
  results_nonheme4[i, "lci"] <- (exp(confint(tempmod2)[2,1])-1)*100
  results_nonheme4[i, "uci"] <- (exp(confint(tempmod2)[2,2])-1)*100
  results_nonheme4[i, "p"] <- tempsum2$coefficients[2,4]
  
  # total iron
  tempmod3 <- lm(paste0(bio[i],"~irona_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum3 <- summary(tempmod3)
  results_iron4[i, "n"] <- dim(tempdata)[1]
  results_iron4[i, "beta"] <- (exp(tempsum3$coefficients[2,1])-1)*100
  results_iron4[i, "lci"] <- (exp(confint(tempmod3)[2,1])-1)*100
  results_iron4[i, "uci"] <- (exp(confint(tempmod3)[2,2])-1)*100
  results_iron4[i, "p"] <- tempsum3$coefficients[2,4]
  
  # dietary iron
  tempmod4 <- lm(paste0(bio[i],"~irona_dt_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum4 <- summary(tempmod4)
  results_iron_dt4[i, "n"] <- dim(tempdata)[1]
  results_iron_dt4[i, "beta"] <- (exp(tempsum4$coefficients[2,1])-1)*100
  results_iron_dt4[i, "lci"] <- (exp(confint(tempmod4)[2,1])-1)*100
  results_iron_dt4[i, "uci"] <- (exp(confint(tempmod4)[2,2])-1)*100
  results_iron_dt4[i, "p"] <- tempsum4$coefficients[2,4]
  
  # supp iron
  tempmod5 <- lm(paste0(bio[i],"~irona_sp_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum5 <- summary(tempmod5)
  results_iron_sp4[i, "n"] <- dim(tempdata)[1]
  results_iron_sp4[i, "beta"] <- (exp(tempsum5$coefficients[2,1])-1)*100
  results_iron_sp4[i, "lci"] <- (exp(confint(tempmod5)[2,1])-1)*100
  results_iron_sp4[i, "uci"] <- (exp(confint(tempmod5)[2,2])-1)*100
  results_iron_sp4[i, "p"] <- tempsum5$coefficients[2,4] 
}

results_health <- rbind(results_heme4,
                        results_nonheme4)
results_health <- rbind(results_health,
                        results_iron4)
results_health <- rbind(results_health,
                        results_iron_dt4)
results_health <- rbind(results_health,
                        results_iron_sp4)

results_health$bon <- results_health$p*12
results_health$per <- paste0(format(round(results_health$beta,digits=1),digits=1),
                           " (",
                           format(round(results_health$lci,digits=1),digits=1),
                           ", ",
                           format(round(results_health$uci,digits=1),digits=1),
                           ")")
sum(results_health$bon<0.05) # 11
View(results_health[results_health$bon<0.05,])

##############################################
#   iron-biomarker associations (nSES)
##############################################
results_heme5 <- results_heme
results_nonheme5 <- results_nonheme
results_iron5 <- results_iron
results_iron_dt5 <- results_iron_dt
results_iron_sp5 <- results_iron_sp

for (i in 1:12) {
  tempdata <- pool %>% filter(!is.na(hemea_avg))
  tempdata <- tempdata[!is.na(tempdata[,bio[i]]),]
  tempdata$nonhemea_avg <- tempdata$nonhemea_avg/20
  tempdata$irona_avg <- tempdata$irona_avg/20
  tempdata$irona_sp_avg <- tempdata$irona_sp_avg/20
  tempdata$irona_dt_avg <- tempdata$irona_dt_avg/10
  
  # heme iron
  tempmod1 <- lm(paste0(bio[i],"~hemea_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum1 <- summary(tempmod1)
  results_heme5[i, "n"] <- dim(tempdata)[1]
  results_heme5[i, "beta"] <- (exp(tempsum1$coefficients[2,1])-1)*100
  results_heme5[i, "lci"] <- (exp(confint(tempmod1)[2,1])-1)*100
  results_heme5[i, "uci"] <- (exp(confint(tempmod1)[2,2])-1)*100
  results_heme5[i, "p"] <- tempsum1$coefficients[2,4]
  
  # nonheme iron
  tempmod2 <- lm(paste0(bio[i],"~nonhemea_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum2 <- summary(tempmod2)
  results_nonheme5[i, "n"] <- dim(tempdata)[1]
  results_nonheme5[i, "beta"] <- (exp(tempsum2$coefficients[2,1])-1)*100
  results_nonheme5[i, "lci"] <- (exp(confint(tempmod2)[2,1])-1)*100
  results_nonheme5[i, "uci"] <- (exp(confint(tempmod2)[2,2])-1)*100
  results_nonheme5[i, "p"] <- tempsum2$coefficients[2,4]
  
  # total iron
  tempmod3 <- lm(paste0(bio[i],"~irona_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum3 <- summary(tempmod3)
  results_iron5[i, "n"] <- dim(tempdata)[1]
  results_iron5[i, "beta"] <- (exp(tempsum3$coefficients[2,1])-1)*100
  results_iron5[i, "lci"] <- (exp(confint(tempmod3)[2,1])-1)*100
  results_iron5[i, "uci"] <- (exp(confint(tempmod3)[2,2])-1)*100
  results_iron5[i, "p"] <- tempsum3$coefficients[2,4]
  
  # dietary iron
  tempmod4 <- lm(paste0(bio[i],"~irona_dt_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum4 <- summary(tempmod4)
  results_iron_dt5[i, "n"] <- dim(tempdata)[1]
  results_iron_dt5[i, "beta"] <- (exp(tempsum4$coefficients[2,1])-1)*100
  results_iron_dt5[i, "lci"] <- (exp(confint(tempmod4)[2,1])-1)*100
  results_iron_dt5[i, "uci"] <- (exp(confint(tempmod4)[2,2])-1)*100
  results_iron_dt5[i, "p"] <- tempsum4$coefficients[2,4]
  
  # supp iron
  tempmod5 <- lm(paste0(bio[i],"~irona_sp_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum5 <- summary(tempmod5)
  results_iron_sp5[i, "n"] <- dim(tempdata)[1]
  results_iron_sp5[i, "beta"] <- (exp(tempsum5$coefficients[2,1])-1)*100
  results_iron_sp5[i, "lci"] <- (exp(confint(tempmod5)[2,1])-1)*100
  results_iron_sp5[i, "uci"] <- (exp(confint(tempmod5)[2,2])-1)*100
  results_iron_sp5[i, "p"] <- tempsum5$coefficients[2,4] 
}

results_nSES <- rbind(results_heme5,
                      results_nonheme5)
results_nSES <- rbind(results_nSES,
                      results_iron5)
results_nSES <- rbind(results_nSES,
                      results_iron_dt5)
results_nSES <- rbind(results_nSES,
                      results_iron_sp5)

results_nSES$bon <- results_nSES$p*12
results_nSES$per <- paste0(format(round(results_nSES$beta,digits=1),digits=1),
                             " (",
                             format(round(results_nSES$lci,digits=1),digits=1),
                             ", ",
                             format(round(results_nSES$uci,digits=1),digits=1),
                             ")")
sum(results_nSES$bon<0.05) # 10

##############################################
#   iron-biomarker associations (case)
##############################################
pool5 <- pool[pool$caco==1,]

results_heme6 <- results_heme
results_nonheme6 <- results_nonheme
results_iron6 <- results_iron
results_iron_dt6 <- results_iron_dt
results_iron_sp6 <- results_iron_sp

for (i in 1:12) {
  tempdata <- pool5 %>% filter(!is.na(hemea_avg))
  tempdata <- tempdata[!is.na(tempdata[,bio[i]]),]
  tempdata$nonhemea_avg <- tempdata$nonhemea_avg/20
  tempdata$irona_avg <- tempdata$irona_avg/20
  tempdata$irona_sp_avg <- tempdata$irona_sp_avg/20
  tempdata$irona_dt_avg <- tempdata$irona_dt_avg/10
  
  # heme iron
  tempmod1 <- lm(paste0(bio[i],"~hemea_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum1 <- summary(tempmod1)
  results_heme6[i, "n"] <- dim(tempdata)[1]
  results_heme6[i, "beta"] <- (exp(tempsum1$coefficients[2,1])-1)*100
  results_heme6[i, "lci"] <- (exp(confint(tempmod1)[2,1])-1)*100
  results_heme6[i, "uci"] <- (exp(confint(tempmod1)[2,2])-1)*100
  results_heme6[i, "p"] <- tempsum1$coefficients[2,4]
  
  # nonheme iron
  tempmod2 <- lm(paste0(bio[i],"~nonhemea_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum2 <- summary(tempmod2)
  results_nonheme6[i, "n"] <- dim(tempdata)[1]
  results_nonheme6[i, "beta"] <- (exp(tempsum2$coefficients[2,1])-1)*100
  results_nonheme6[i, "lci"] <- (exp(confint(tempmod2)[2,1])-1)*100
  results_nonheme6[i, "uci"] <- (exp(confint(tempmod2)[2,2])-1)*100
  results_nonheme6[i, "p"] <- tempsum2$coefficients[2,4]
  
  # total iron
  tempmod3 <- lm(paste0(bio[i],"~irona_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum3 <- summary(tempmod3)
  results_iron6[i, "n"] <- dim(tempdata)[1]
  results_iron6[i, "beta"] <- (exp(tempsum3$coefficients[2,1])-1)*100
  results_iron6[i, "lci"] <- (exp(confint(tempmod3)[2,1])-1)*100
  results_iron6[i, "uci"] <- (exp(confint(tempmod3)[2,2])-1)*100
  results_iron6[i, "p"] <- tempsum3$coefficients[2,4]
  
  # dietary iron
  tempmod4 <- lm(paste0(bio[i],"~irona_dt_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum4 <- summary(tempmod4)
  results_iron_dt6[i, "n"] <- dim(tempdata)[1]
  results_iron_dt6[i, "beta"] <- (exp(tempsum4$coefficients[2,1])-1)*100
  results_iron_dt6[i, "lci"] <- (exp(confint(tempmod4)[2,1])-1)*100
  results_iron_dt6[i, "uci"] <- (exp(confint(tempmod4)[2,2])-1)*100
  results_iron_dt6[i, "p"] <- tempsum4$coefficients[2,4]
  
  # supp iron
  tempmod5 <- lm(paste0(bio[i],"~irona_sp_avg+age_blddraw+fast+caco+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
  tempsum5 <- summary(tempmod5)
  results_iron_sp6[i, "n"] <- dim(tempdata)[1]
  results_iron_sp6[i, "beta"] <- (exp(tempsum5$coefficients[2,1])-1)*100
  results_iron_sp6[i, "lci"] <- (exp(confint(tempmod5)[2,1])-1)*100
  results_iron_sp6[i, "uci"] <- (exp(confint(tempmod5)[2,2])-1)*100
  results_iron_sp6[i, "p"] <- tempsum5$coefficients[2,4]
  
}

results_case <- rbind(results_heme6,
                     results_nonheme6)
results_case <- rbind(results_case,
                     results_iron6)
results_case <- rbind(results_case,
                     results_iron_dt6)
results_case <- rbind(results_case,
                     results_iron_sp6)

results_case$bon <- results_case$p*12
results_case$per <- paste0(format(round(results_case$beta,digits=1),digits=1),
                      " (",
                      format(round(results_case$lci,digits=1),digits=1),
                      ", ",
                      format(round(results_case$uci,digits=1),digits=1),
                      ")")
sum(results_case$bon<0.05) # 6

save.image(file="/udd/nhfwa/iron/pool/biomarker/biomarker.RData")

library(openxlsx)
write.xlsx(list("results"=results[results$type=="heme",],
                "results_con"=results_con[results_con$type=="heme",],
                "results_case"=results_case[results_case$type=="heme",],
                "results_fast"=results_fast[results_fast$type=="heme",],
                "results_health"=results_health[results_health$type=="heme",],
                "results_nSES"=results_nSES[results_nSES$type=="heme",],
                "results_wthcrp"=results_wthcrp[results_wthcrp$type=="heme",]), 
           file = "/udd/nhfwa/iron/pool/biomarker/biomarker_results.xlsx")

