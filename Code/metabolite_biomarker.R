library(tidyverse)

load("/udd/nhfwa/iron/pool/metabolomics/score_t2d/t2d_merge.RData")
met_coef800 <- met_coef[met_coef$n>=800, ]
gdata::keep(met_coef, met_coef800, sure=T)

load("/udd/nhfwa/iron/pool/metabolomics/metscore_t2d_v2.RData")
pool_met <- pool
gdata::keep(pool_met, met_coef, met_coef800, heme_metabolite4, sure=T)
length(intersect(met_coef800$met, heme_metabolite4$met)) # 5
intersect(met_coef800$metabolite_name, heme_metabolite4$name1)
# [1] "urate/uric acid/Uric acid" "valine/Valine"             "C34:1 DAG/DG(34:1)"       
# [4] "C38:2 PE/PE(38:2)"         "C18:2 LPC/LPC(18:2)" 

load("~/iron/pool/biomarker/biomarker.RData")
pool_bio <- pool
gdata::keep(pool_bio, pool_met, met_coef, met_coef800, heme_metabolite4, bio, sure=T)

met51 <- met_coef800$met
met17 <- heme_metabolite4$met

#########################################################
#       T2D metabolites - biomarker associations 1
#########################################################
bio_results <- list()
for (i in 1:length(bio)){
  temp_results <- data.frame(bio=bio[i],
                             n=NA,
                             met=met51,
                             per=NA,
                             p=NA,
                             fdr=NA,
                             bon=NA)
  tempdata <- merge(pool_met, pool_bio[, c("newid", bio[i])], by="newid")
  tempdata <- tempdata[!is.na(tempdata[,bio[i]]),]
  temp_results$n <- dim(tempdata)[1]
  
  for (j in 1:length(met51)){
    tempdata[, met51[j]] <- scale(tempdata[, met51[j]])
    tempmod <- lm(paste0(bio[i],"~", met51[j], "+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
    tempsum <- summary(tempmod)
    temp_results[j, "per"] <- (exp(tempsum$coefficients[2,1])-1)*100
    temp_results[j, "p"] <- tempsum$coefficients[2,4]
  }
  
  temp_results$fdr <- p.adjust(temp_results$p, method="fdr", n=length(temp_results$p))
  temp_results$bon <- p.adjust(temp_results$p, method="bonferroni", n=length(temp_results$p))
  bio_results[[i]] <- temp_results
}

bio_all <- bind_rows(bio_results)
bio_all <- merge(met_coef800[, -4], bio_all, by="met")
table(bio_all[bio_all$bon<0.05, ]$bio)
# adj_lna1c    adj_lnadipo     adj_lnchol     adj_lncpep      adj_lncrp adj_lnferritin 
# 8             31             27             31             34             11 
# adj_lnhdl      adj_lnldl   adj_lnleptin    adj_lnratio     adj_lntrig 
# 21             12             22             28             34 

#########################################################
#       T2D metabolites - biomarker associations 2
#########################################################
bio_results2 <- list()
for (i in 1:length(bio)){
  temp_results <- data.frame(bio=bio[i],
                             n=NA,
                             met=met17,
                             per=NA,
                             p=NA,
                             fdr=NA,
                             bon=NA)
  tempdata <- merge(pool_met, pool_bio[, c("newid", bio[i])], by="newid")
  tempdata <- tempdata[!is.na(tempdata[,bio[i]]),]
  temp_results$n <- dim(tempdata)[1]
  
  for (j in 1:length(met17)){
    tempdata[, met17[j]] <- scale(tempdata[, met17[j]])
    tempmod <- lm(paste0(bio[i],"~", met17[j], "+ageyr+fast+cacoall+n2+white+pmhc2+pmhc3+pmhc4+pmhc9+
  hbpbase+cholbase+bmigp2+bmigp3+bmigp4+alc2+alc3+alc4+act2+act3+act4+act5+mv+asp+smk2+smk3+
  calorn_avg+magnn_avg+gln_avg+cerafn_avg+psn_avg+transn_avg"), data=tempdata)
    tempsum <- summary(tempmod)
    temp_results[j, "per"] <- (exp(tempsum$coefficients[2,1])-1)*100
    temp_results[j, "p"] <- tempsum$coefficients[2,4]
  }
  
  temp_results$fdr <- p.adjust(temp_results$p, method="fdr", n=length(temp_results$p))
  temp_results$bon <- p.adjust(temp_results$p, method="bonferroni", n=length(temp_results$p))
  bio_results2[[i]] <- temp_results
}

bio_all2 <- bind_rows(bio_results2)
bio_all2 <- merge(heme_metabolite4[, -4], bio_all2, by="met")
table(bio_all2[bio_all2$bon<0.05, ]$bio)
# adj_lna1c    adj_lnadipo     adj_lnchol     adj_lncpep      adj_lncrp adj_lnferritin 
# 2             15             12             16             13              6 
# adj_lnhdl      adj_lnldl   adj_lnleptin    adj_lnratio      adj_lntrf     adj_lntrig 
# 14              4             11             14              1             17 

save.image(file="/udd/nhfwa/iron/pool/metabolomics/metabolite_biomarker.RData")
write.csv(met_coef, file="/udd/nhfwa/iron/pool/metabolomics/metabolite_name.csv", row.names = F)
