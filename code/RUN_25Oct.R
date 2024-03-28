gc()
library(readxl); library(writexl); library(ggplot2); library(scales); library(cowplot) #library(reshape2); library(plyr); library(splines); library(Hmisc); library(wesanderson); library(extrafont); #library(cowplot) # This last one is just a random color palett I found online as a test
library(ggseg); library(MASS); library(sfsmisc); library(rstanarm);  #library(ggseg3d)
library(mgcv)
library(tidymv)
library(plotrix)
library(table1)
library(tableone)
library(mclust)
library(data.table)
library(ggrepel)
library(ggalluvial)
library(dplyr)
library(mediation)
library(tidyr)
library(psych)
library(irr)
library(stargazer)
library(tibble)
library(stringr)
library(BSDA)
library(geiger)
library(ggpubr)
library(plyr)
library(gridExtra)
library(cowplot)
#clear memory

Apos_cutoff <- 18
Tpos_cutoff <- 1.3

source("./code/functions_25Oct.R")
source("./code/ABCDSclean_expanded.R")
scatterTheme<-theme(legend.key.size = unit(25, 'points'), legend.title = element_text(size=16), legend.text = element_text(size = 14), axis.title = element_text(size=18), axis.text = element_text(size=16))



#Create variable DFs ------
df_GFAP<-merge(demog_abcds, plasma_abcds[,c("subject_label", "event_code", "plasma_GFAP_value")], by=c("subject_label", "event_code"), all=FALSE)
df_GFAP<-df_GFAP[!is.na(df_GFAP$plasma_GFAP_value),]
df_pT217<-merge(demog_abcds, plasma_abcds[,c("subject_label", "event_code", "plasma_Ptau217_value")], by=c("subject_label", "event_code"), all=FALSE)
df_pT217<-df_pT217[!is.na(df_pT217$plasma_Ptau217_value),]
df_centiloid<-merge(demog_abcds, centiloid_abcds, by=c("subject_label", "event_code"), all=FALSE)
df_tauPET<-merge(demog_abcds, braak_abcds, by=c("subject_label", "event_code"), all=FALSE)
rm(AV45_abcds, centiloid_abcds, PIB_abcds, plasma_abcds, braak_abcds)



df_centiloid$Acat<-as.factor(ifelse(df_centiloid$ds_vs_control_flag=="control", "control",
                                    ifelse(df_centiloid$WUSTLcentiloid > Apos_cutoff, "A+", "A-")))
df_centiloid$Acat<-factor(df_centiloid$Acat, levels = c("control", "A-", "A+"))

df_tauPET$Tcat<-as.factor(ifelse(df_tauPET$ds_vs_control_flag=="control", "control",
                                      ifelse(df_tauPET$TempMeta < Tpos_cutoff, "T-", "T+")))
df_tauPET$Tcat<-factor(df_tauPET$Tcat, levels = c("control", "T-", "T+"))

df_temp<-subset(df_GFAP, ds_vs_control_flag == "DS")
df_temp<-df_temp[!is.na(df_temp$plasma_GFAP_value),]
BIC_GFAP = mclustBIC(df_temp$plasma_GFAP_value)
mdl_GFAP = Mclust(df_temp$plasma_GFAP_value, 3, modelNames = "V")
minGFAP<-min(df_temp$plasma_GFAP_value, na.rm = TRUE)
maxGFAP<-max(df_temp$plasma_GFAP_value, na.rm = TRUE)
newdata<-seq(minGFAP, maxGFAP, 0.01)
pred_GFAP<- predict.Mclust(mdl_GFAP, newdata)
high_cutoff<-max(newdata[pred_GFAP$classification==1])
df_GFAP$GFAPcat<-ifelse(df_GFAP$ds_vs_control_flag=="control", "control",
                        ifelse(df_GFAP$plasma_GFAP_value > high_cutoff, "high_GFAP", "low_GFAP"))
df_GFAP$GFAPcat<-factor(df_GFAP$GFAPcat, levels = c("control", "low_GFAP", "high_GFAP"))

rm(BIC_GFAP, mdl_GFAP, newdata, pred_GFAP, TPP, TPP_G, TPP_G_all, TPP_obj, TPP_Z, Staged_df_all, Staged_df_tau, sum_TPPs_all, sum_TPPs_tau)

colnames(df_GFAP)[colnames(df_GFAP)=="event_code"]<-"event_code_plasma"
colnames(df_GFAP)[colnames(df_GFAP)=="latency_in_days"]<-"latency_plasma"
colnames(df_pT217)[colnames(df_pT217)=="event_code"]<-"event_code_plasma"
colnames(df_pT217)[colnames(df_pT217)=="latency_in_days"]<-"latency_plasma"
colnames(df_centiloid)[colnames(df_centiloid)=="event_code"]<-"event_code_amy"
colnames(df_centiloid)[colnames(df_centiloid)=="latency_in_days"]<-"latency_amy"
colnames(df_tauPET)[colnames(df_tauPET)=="event_code"]<-"event_code_tau"
colnames(df_tauPET)[colnames(df_tauPET)=="latency_in_days"]<-"latency_tau"

#Only 1 visit per participant (will keep the first visit listed, usually baseline)
df_centiloid<-df_centiloid[!duplicated(df_centiloid$subject_label),]

source("./code/VARbyEYO_25Oct.R")
# pval_CentTauPET<-1 - sum(ifelse(SigDiffEYO_cent < SigDiffEYO_tauPET, 1, 0))/10000
# pval_CentpTau217<-1 - sum(ifelse(SigDiffEYO_cent < SigDiffEYO_pT217, 1, 0))/10000
# pval_CentGFAP<-1 - sum(ifelse(SigDiffEYO_cent < SigDiffEYO_GFAP, 1, 0))/10000
# pval_pTau217TauPET<-1 - sum(ifelse(SigDiffEYO_pT217 < SigDiffEYO_tauPET, 1, 0))/10000
# pval_pTau217GFAP<-1 - sum(ifelse(SigDiffEYO_pT217 < SigDiffEYO_GFAP, 1, 0))/10000
# pval_TauPETGFAP<-1 - sum(ifelse(SigDiffEYO_tauPET < SigDiffEYO_GFAP, 1, 0))/10000

 source("./code/TauMetaAnalysis_25July.R")

# #Demographics----
listVars<-c("age_at_visit", "EYO", "de_gender", "de_race", "karyotype", "apoe4", "allele_combo", "consensus")
catVars<-c("de_gender", "de_race", "karyotype", "apoe4", "allele_combo", "consensus")

demogVars<-c("subject_label","ds_vs_control_flag",  "age_at_visit", "EYO", "de_gender", "de_race", "HISPANIC", "karyotype", "apoe4", "allele_combo", "consensus")
df_all<-merge(df_GFAP, df_pT217, by=c(demogVars, "event_code_plasma"), all=TRUE)
df_all<-merge(df_all, df_centiloid, by.x = c(demogVars, "event_code_plasma"), by.y = c(demogVars, "event_code_amy"), all=TRUE)
df_all<-merge(df_all, df_tauPET, by.x = c(demogVars, "event_code_plasma"), by.y = c(demogVars, "event_code_tau"), all=TRUE)
DemogAll_CtlDS<-CreateTableOne(vars = c(listVars, "plasma_GFAP_value", "plasma_Ptau217_value", "WUSTLcentiloid", "TempMeta"), data = df_all, factorVars = catVars, strata = "ds_vs_control_flag")

DemogGFAP_CtlDS<-CreateTableOne(vars = c(listVars, "plasma_GFAP_value"), data = df_GFAP, factorVars = catVars, strata = "ds_vs_control_flag")
DemogpT217_CtlDS<-CreateTableOne(vars = c(listVars, "plasma_Ptau217_value"), data = df_pT217, factorVars = catVars, strata = "ds_vs_control_flag")
DemogAmy_CtlDS<-CreateTableOne(vars = c(listVars, "WUSTLcentiloid"), data = df_centiloid, factorVars = catVars, strata = "ds_vs_control_flag")
DemogTau_CtlDS<-CreateTableOne(vars = c(listVars, "TempMeta"), data = df_tauPET, factorVars = catVars, strata = "ds_vs_control_flag")

#################################################################################
#################################################################################
#BOXPLOT
#################################################################################
#################################################################################
df_all$Group <- as.factor(paste0(df_all$ds_vs_control_flag, df_all$Acat, df_all$Tcat))
levels(df_all$Group) <- c("Control", NA, NA, NA, NA, "A-/T-", NA, NA, "A+/T-", "A+/T+",
                          NA, NA, NA)
m_GFAP <- lm(plasma_GFAP_value ~ Group + de_gender + age_at_visit + latency_tau, data = df_all[!is.na(df_all$Group),] )
m_pTau <- lm(plasma_Ptau217_value ~ Group + de_gender + age_at_visit + latency_tau, data = df_all[!is.na(df_all$Group),] )

Tukey_GFAP <- data.frame(TukeyHSD(aov(m_GFAP), which = "Group")$Group)
Tukey_GFAP$Group <- row.names(Tukey_GFAP)
Tukey_GFAP$group1 <- substr(Tukey_GFAP$Group, start = 1, stop = 5)
Tukey_GFAP$group2 <- substr(Tukey_GFAP$Group, start = 7, stop = nchar(Tukey_GFAP$Group))
Tukey_GFAP$p.adj <- round(Tukey_GFAP$p.adj, 3)

Tukey_pTau <- data.frame(TukeyHSD(aov(m_pTau), which = "Group")$Group)
Tukey_pTau$Group <- row.names(Tukey_pTau)
Tukey_pTau$group1 <- substr(Tukey_pTau$Group, start = 1, stop = 5)
Tukey_pTau$group2 <- substr(Tukey_pTau$Group, start = 7, stop = nchar(Tukey_pTau$Group))
Tukey_pTau$p.adj <- round(Tukey_pTau$p.adj, 3)


tmp <- df_all[!is.na(df_all$Group), c("Group", "plasma_GFAP_value")]

p1 <- ggboxplot(tmp, x = "Group", y = "plasma_GFAP_value", colour = "Group", 
          fill = "Group", add = "jitter") + 
  stat_pvalue_manual(Tukey_GFAP, y.position = 400, step.increase = 0.05, 
                     label = "{scales::pvalue(p.adj)}") +
  ylab("Plasma GFAP (pg/mL)") + theme(legend.position = "none") +
  xlab("") + ggtitle("A.")


tmp <- df_all[!is.na(df_all$Group), c("Group", "plasma_Ptau217_value")]

p2 <- ggboxplot(tmp, x = "Group", y = "plasma_Ptau217_value", colour = "Group", 
          fill = "Group", add = "jitter") + 
  stat_pvalue_manual(Tukey_pTau, y.position = 1.35, step.increase = 0.05, 
                     label = "{scales::pvalue(p.adj)}") +
  ylab("Plasma pTau-217 (pg/mL)") + theme(legend.position = "none") + xlab("") +
  ggtitle("B.")

grid.arrange(p1, p2, nrow = 1)

#################################################################################
#################################################################################



#################################################################################
#################################################################################
#MEDIATION ANALYSIS
#################################################################################
#################################################################################
model.0 <- lm(tauopathy ~ WUSTLcentiloid + apoe4 + de_gender + latency_amy + latency_tau + latency_plasma.y, data = df_all[!is.na(df_all$tauopathy) & !is.na(df_all$WUSTLcentiloid) & !is.na(df_all$plasma_GFAP_value) & df_all$ds_vs_control_flag=="DS",])
model.1 <- lm(plasma_GFAP_value ~ WUSTLcentiloid + apoe4 + de_gender+ latency_amy + latency_tau + latency_plasma.y, data = df_all[!is.na(df_all$tauopathy) & !is.na(df_all$WUSTLcentiloid) & !is.na(df_all$plasma_GFAP_value)& df_all$ds_vs_control_flag=="DS",])
model.y <- lm(tauopathy ~ WUSTLcentiloid + plasma_GFAP_value + apoe4 + de_gender+ latency_amy + latency_tau + latency_plasma.y, data = df_all[!is.na(df_all$tauopathy) & !is.na(df_all$WUSTLcentiloid) & !is.na(df_all$plasma_GFAP_value)& df_all$ds_vs_control_flag=="DS",])

library(mediation)
results <- mediation::mediate(model.1, model.y, treat='WUSTLcentiloid', mediator='plasma_GFAP_value',
                   covariates = c("apoe4", "de_gender", "latency_amy", "latency_tau", "latency_plasma.y"),
                   boot=TRUE, sims=1000)


summary(results)


model.0_ptau <- lm(plasma_Ptau217_value ~ WUSTLcentiloid + apoe4 + de_gender+ latency_amy + latency_tau + latency_plasma.y, data = df_all[!is.na(df_all$tauopathy) & !is.na(df_all$WUSTLcentiloid) & !is.na(df_all$plasma_Ptau217_value)& df_all$ds_vs_control_flag=="DS",])
model.1_ptau <- lm(plasma_GFAP_value ~ WUSTLcentiloid + apoe4 + de_gender+ latency_amy + latency_tau + latency_plasma.y, data = df_all[!is.na(df_all$tauopathy) & !is.na(df_all$WUSTLcentiloid) & !is.na(df_all$plasma_Ptau217_value)& df_all$ds_vs_control_flag=="DS",])
model.y_ptau <- lm(plasma_Ptau217_value ~ WUSTLcentiloid + plasma_GFAP_value + apoe4 + de_gender+ latency_amy + latency_tau + latency_plasma.y, data = df_all[!is.na(df_all$tauopathy) & !is.na(df_all$WUSTLcentiloid) & !is.na(df_all$plasma_Ptau217_value)& df_all$ds_vs_control_flag=="DS",])

results_ptau <- mediation::mediate(model.1_ptau, model.y_ptau, treat='WUSTLcentiloid', mediator='plasma_GFAP_value',
                              covariates = c("apoe4", "de_gender", "latency_amy", "latency_tau", "latency_plasma.y"),
                              boot=TRUE, sims=1000)


summary(results_ptau)



#################################################################################
#################################################################################
