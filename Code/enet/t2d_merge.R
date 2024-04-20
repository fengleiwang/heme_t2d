library(tidyverse)

t2dpath <- dir(path="~/iron/pool/metabolomics/score_t2d", pattern = ".RData")
load(paste0("~/iron/pool/metabolomics/score_t2d/", t2dpath[1]))

t2dscore <- bind_rows(score)
t2d_metcoef <- met_coef_db2
cindex <- record_db2[!is.na(record_db2$NoMetabs),]

for (i in 2:10) {
  load(paste0("~/iron/pool/metabolomics/score_t2d/", t2dpath[i]))
  temp <- bind_rows(score)
  t2dscore <- rbind(t2dscore, temp)
  t2d_metcoef <- cbind(t2d_metcoef, met_coef_db2[, 2:101])
  cindex <- rbind(cindex, record_db2[!is.na(record_db2$NoMetabs),])
  gdata::keep(t2dscore, t2d_metcoef, cindex, t2dpath, sure=T)
}
colnames(t2dscore)[2] <- "t2dscore"
mean(cindex$Cindex) # 0.7919

load("~/iron/pool/metabolomics/score_t2d/t2d_enet1.RData")
t2dscore <- merge(t2dscore, 
                  all[, c("newid", "fast", "endpoint", "ageyr", "ahei_noal", "dead", "tdead",
                          "dead_can", "dead_cvd", "db2case", "tdb2", "cvdcase", "tcvd",
                          "chdcase", "tchd", "strcase", "tstr", "cacoall",
                          met277)], 
                  by="newid")

gdata::keep(t2dscore, t2d_metcoef, cindex, f, sure=T)
write.csv(t2dscore, "/udd/nhfwa/iron/pool/metabolomics/score_t2d/t2dscore.csv", na="", row.names = F)

sum(rowSums(t2d_metcoef[, 2:1001]!=0)==1000) # 30
sum(rowSums(t2d_metcoef[, 2:1001]!=0)>950) # 41
sum(rowSums(t2d_metcoef[, 2:1001]!=0)>900) # 49
sum(rowSums(t2d_metcoef[, 2:1001]!=0)>800) # 51

met_coef <- data.frame(met=t2d_metcoef$HMDB,
                       n=rowSums(t2d_metcoef[, 2:1001]!=0))
met_coef$coef <- rowSums(t2d_metcoef[, 2:1001])/met_coef$n

f$met <- rownames(f)
met_coef <- merge(f[, c("met", "metabolite_name", "biochemical_name")], met_coef, by="met")

save.image("/udd/nhfwa/iron/pool/metabolomics/score_t2d/t2d_merge.RData")
