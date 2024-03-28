## Tmeta/Taumeta/TempMeta = Temporal Meta ROI pulled from Cho et al. (2016)

#Tau PET only analysis ------
dfGFAP_ATmeta<-merge(df_GFAP, df_tauPET[, c("subject_label", "event_code_tau", "latency_tau", "Tcat", "TempMeta")], by="subject_label", all=FALSE)
dfGFAP_ATmeta<-merge(dfGFAP_ATmeta, df_centiloid[,c("subject_label", "event_code_amy", "latency_amy", "Acat", "WUSTLcentiloid")], by="subject_label", all=FALSE)
dfGFAP_ATmeta<-dfGFAP_ATmeta[!duplicated(dfGFAP_ATmeta$subject_label),]
# dfGFAP_ATmeta$latency_CentTau<-abs(dfGFAP_ATmeta$latency_amy - dfGFAP_ATmeta$latency_tau)

## AT cat boxplot (excluding A-/T+) -----
dfGFAP_ATmeta$ATcat<-ifelse(dfGFAP_ATmeta$Acat=="control", "control",
                            ifelse(dfGFAP_ATmeta$Acat=="A-", ifelse(dfGFAP_ATmeta$Tcat=="T-", "A-/T-", "A-/T+"),
                                   ifelse(dfGFAP_ATmeta$Acat=="A+", ifelse(dfGFAP_ATmeta$Tcat=="T-", "A+/T-", "A+/T+"), NA)))
dfGFAP_ATmeta<-subset(dfGFAP_ATmeta, ATcat != "A-/T+")
dfGFAP_ATmeta$ATcat<-factor(dfGFAP_ATmeta$ATcat, levels = c("control", "A-/T-", "A+/T-", "A+/T+"))
Comp<-list(c("control", "A-/T-"), c("control", "A+/T-"), c("control", "A+/T+"), c("A-/T-", "A+/T-"), c("A-/T-", "A+/T+"), c("A+/T-", "A+/T+"))

ATcat_GFAP<-ggplot(dfGFAP_ATmeta, aes(x=ATcat, y=plasma_GFAP_value, color=ATcat, shape=ATcat)) + geom_boxplot() +
  ylab("Plasma GFAP (pg/mL)") + xlab("Amyloid/Tau Positivity") +
  geom_jitter(position=position_jitter(0.2)) + stat_compare_means(comparisons = Comp, size = 8) + 
  # scale_color_manual(values=c("#051405", "#1a661a", "#2eb82e", "#5cd65c")) +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size = 22), axis.title = element_text(size = 24))




## Cent~TauMeta by GFAPcat ------
AvTpet_GFAPcat<-ggplot(dfGFAP_ATmeta, aes(x=WUSTLcentiloid, y=TempMeta, color=GFAPcat)) + geom_point() + geom_smooth(method = "gam") +
  scatterTheme + ylab("Tau PET (Temporal-Meta SUVR)") + xlab("Amyloid PET (Centiloids)")

## Mediation Cent~TauMeta by GFAP ------
df<-merge(df_GFAP, df_tauPET[, c("subject_label", "event_code_tau", "latency_tau", "Tcat", "TempMeta")], by="subject_label", all=FALSE)
df<-merge(df, df_centiloid[,c("subject_label", "event_code_amy", "latency_amy", "Acat", "WUSTLcentiloid")], by="subject_label", all=FALSE)
df<-df[!duplicated(df$subject_label),]

medDF<-subset(df, ds_vs_control_flag=="DS")
medDF$latency_GFAPamy<-abs(medDF$latency_amy - medDF$latency_plasma)
medDF$latency_GFAPtau<-abs(medDF$latency_tau - medDF$latency_plasma)
medDF$latency_tauamy<-abs(medDF$latency_amy - medDF$latency_tau)

AmyGFAP_mdl<-lm(plasma_GFAP_value ~ WUSTLcentiloid + latency_GFAPtau + latency_tauamy + latency_GFAPamy + de_gender + apoe4, data = medDF)
TauGFAP_mdl<-lm(TempMeta ~ WUSTLcentiloid + plasma_GFAP_value + latency_GFAPtau + latency_tauamy + latency_GFAPamy + de_gender + apoe4, data=medDF)
TauGFAPAmy_mdl<-lm(TempMeta ~ WUSTLcentiloid*plasma_GFAP_value + latency_GFAPtau + latency_tauamy + latency_GFAPamy + de_gender + apoe4, data=medDF)
GFAPmed_TauMetaPET<-mediation::mediate(AmyGFAP_mdl, TauGFAP_mdl, treat='WUSTLcentiloid', mediator='plasma_GFAP_value', boot = TRUE)

summary(GFAPmed_TauMetaPET)


#pTau217 only analysis ------
dfGFAP_ApT217<-merge(df_GFAP, df_pT217[, c("subject_label", "event_code_plasma", "latency_plasma", "plasma_Ptau217_value")], by="subject_label", all=FALSE)
dfGFAP_ApT217<-merge(dfGFAP_ApT217, df_centiloid[,c("subject_label", "event_code_amy", "latency_amy", "Acat", "WUSTLcentiloid")], by="subject_label", all=FALSE)
dfGFAP_ApT217<-dfGFAP_ApT217[!duplicated(dfGFAP_ApT217$subject_label),]
## Cent~pT217 by GFAPcat ------
AvpT217_GFAPcat<-ggplot(dfGFAP_ApT217, aes(x=WUSTLcentiloid, y=plasma_Ptau217_value, color=GFAPcat)) + geom_point() + geom_smooth(method = "gam") +
  scatterTheme + ylab("pTau-217 (pg/mL)") + xlab("Amyloid PET (Centiloids)")
# jpeg("./output/CentpT217_GFAPcat.jpeg", width = 800, height = 600)
# print(AvpT217_GFAPcat)
# dev.off()

## Mediation Cent~pT217 by GFAP ------
medDF<-subset(dfGFAP_ApT217, ds_vs_control_flag=="DS")
medDF$latency_GFAPamy<-abs(medDF$latency_amy - medDF$latency_plasma.x)

AmyGFAP_mdl_pt<-lm(plasma_GFAP_value ~ WUSTLcentiloid + latency_GFAPamy + de_gender + apoe4, data = medDF)
TauGFAP_mdl_pt<-lm(plasma_Ptau217_value ~ WUSTLcentiloid + plasma_GFAP_value + latency_GFAPamy + de_gender + apoe4, data=medDF)
TauGFAPAmy_mdl_pt<-lm(plasma_Ptau217_value ~ WUSTLcentiloid*plasma_GFAP_value, data=medDF)
GFAPmed_pT217<-mediation::mediate(AmyGFAP_mdl_pt, TauGFAP_mdl_pt, treat='WUSTLcentiloid', mediator='plasma_GFAP_value', boot = TRUE)

summary(GFAPmed_pT217)


library(psych)
library(lavaan)
library(ggplot2)
library(readxl)
library(semPlot)


mediation_model <- '
  # Direct effects
  plasma_GFAP_value ~ a * WUSTLcentiloid
  plasma_Ptau217_value ~ c * WUSTLcentiloid + b * plasma_GFAP_value

  # Indirect effect (a * b)
  indirect := a * b

  # Total effect (c + indirect)
  total := c + indirect
'
# Estimate the mediation model
mediation_results <- sem(mediation_model, data = medDF)

# Summarize the results
summary(mediation_results, standardized = TRUE, fit.measures = TRUE)
semPaths(mediation_results, whatLabels = "est", style = "lisrel", intercepts = FALSE)

summary(AmyGFAP_mdl_pt)
