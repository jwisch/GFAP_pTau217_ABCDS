source("./code/functions.R")
# #~ EYO scatterplots-----
##Centiloid ~ EYO-----
z=1.96
df <-expand.grid("EYO" = seq(from = -25, to = 10, by = 0.1))
df$EYO<-round(df$EYO, digits = 1)
itr<-10000
SigDiffEYO<-matrix(0, itr, 1)
MaxCentiloid<-matrix(0, itr, 1)
AOC<-matrix(0, itr, 1)
amyPETdata_ctrl<-subset(df_centiloid, ds_vs_control_flag=="control")
amyPETdata_DS<-subset(df_centiloid, ds_vs_control_flag=="DS")


# 
# SigDiffEYO <- list()
# for (i in 1:itr){
#   print(i)
#   SigDiffEYO[[i]]<-get_takeOffPoint(amyPETdata_ctrl, amyPETdata_DS, 1.96, "WUSTLcentiloid")
# }
# 
# SigDiffEYO_GFAP <- list()
# for (i in 1:itr){
#   print(i)
#   SigDiffEYO_GFAP[[i]]<-get_takeOffPoint(df_GFAP[df_GFAP$ds_vs_control_flag == "control",],
#                                  df_GFAP[df_GFAP$ds_vs_control_flag == "DS",], 1.96, "plasma_GFAP_value")
# }
# 
# SigDiffEYO_tauPET <- list()
# for (i in 1:itr){
#   print(i)
#   SigDiffEYO_tauPET[[i]]<-get_takeOffPoint(df_tauPET[df_tauPET$ds_vs_control_flag == "control",],
#                                  df_tauPET[df_tauPET$ds_vs_control_flag == "DS",], 1.96, "TempMeta")
# }
# 
# SigDiffEYO_pT217 <- list()
# for (i in 1:itr){
#   print(i)
#   SigDiffEYO_pT217[[i]]<-get_takeOffPoint(df_pT217[df_pT217$ds_vs_control_flag == "control",],
#                                         df_pT217[df_pT217$ds_vs_control_flag == "DS",], 1.96, "plasma_Ptau217_value")
# }
# # 
# mean(unlist(SigDiffEYO))
# mean(unlist(SigDiffEYO_GFAP))
# mean(unlist(SigDiffEYO_tauPET))
# mean(unlist(SigDiffEYO_pT217))
# 
# 
# SigDiffEYO <- unlist(SigDiffEYO)
# SigDiffEYO_GFAP <- unlist(SigDiffEYO_GFAP)
# SigDiffEYO_tauPET <- unlist(SigDiffEYO_tauPET)
# SigDiffEYO_pT217 <- unlist(SigDiffEYO_pT217)



# 
# saveRDS((SigDiffEYO), "./data/generated_WUSTLcentiloidMedianTakeoffpoints_10k.RDS")
# saveRDS((SigDiffEYO_GFAP), "./data/generated_GFAPMedianTakeoffpoints_10k.RDS")
# saveRDS((SigDiffEYO_tauPET), "./data/generated_tauPETMedianTakeoffpoints_10k.RDS")
# saveRDS((SigDiffEYO_pT217), "./data/generated_pT217MedianTakeoffpoints_10k.RDS")

SigDiffEYO <- readRDS("./data/generated_WUSTLcentiloidMedianTakeoffpoints_10k.RDS")
SigDiffEYO_GFAP <- readRDS("./data/generated_GFAPMedianTakeoffpoints_10k.RDS")
SigDiffEYO_tauPET <- readRDS("./data/generated_tauPETMedianTakeoffpoints_10k.RDS")
SigDiffEYO_pT217 <- readRDS("./data/generated_pT217MedianTakeoffpoints_10k.RDS")

ks.test(SigDiffEYO, SigDiffEYO_GFAP)
ks.test(SigDiffEYO, SigDiffEYO_pT217)
ks.test(SigDiffEYO_GFAP, SigDiffEYO_pT217)
ks.test(SigDiffEYO, SigDiffEYO_tauPET)

tmp <- data.frame("val" = c((SigDiffEYO_GFAP),(SigDiffEYO_tauPET)),
                  "measure" = c(rep("GFAP", 10000), rep("tau", 10000) ) )

ggplot(tmp, aes(x = val, group = measure, fill = measure)) + geom_histogram(alpha = 0.3, position = "identity") +
  theme_bw() + 
  geom_vline(xintercept = mean(tmp[tmp$measure == "GFAP", "val"]), linetype = "dashed") +
  geom_vline(xintercept = median(tmp[tmp$measure == "GFAP", "val"]), linetype = "solid") +
  geom_vline(xintercept = mean(tmp[tmp$measure == "tau", "val"]), linetype = "dashed", colour = "blue") +
  geom_vline(xintercept = median(tmp[tmp$measure == "tau", "val"]), linetype = "solid", colour = "blue") 
  

AllSigDiffEYO_centiloid<-SigDiffEYO
SigDiffEYO_centiloid<-round(mean(SigDiffEYO), 1)

scatterEYO_centiloid<-
  ggplot(df_centiloid, aes(x=(age_at_visit), y=WUSTLcentiloid, color=ds_vs_control_flag)) + geom_point() + geom_smooth(method="gam") +
  scale_colour_manual(values = c( "#F8766D", "#00BFC4"), labels = c("Control", "Down Syndrome")) +
  scatterTheme + scale_x_continuous(breaks=c(30, (52.5 + round(SigDiffEYO_centiloid, 1)), 52.5, 60),
                                    sec.axis = sec_axis(~ . - 52.5, breaks = c(-22.5, round(SigDiffEYO_centiloid, 1), 0, 7.5), 
                                                        name = "Estimated Years to Symptom Onset"), 
                                    name = "Age", limits = c(24.8, 72.5)) +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(colour = c('black', "#00BFC4", 'black',  'black',  'black'), 
                                   size=18), legend.title = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 22),
        plot.title = element_text(size = 24, face = "bold")) + geom_vline(xintercept = 52.5, linetype = "dashed") + 
  ylab("Amyloid PET (Centiloids)") 
# jpeg("./output/CentEYO_DSvCtrl_annot.jpeg", width = 800, height = 600)
# print(scatterEYO_centiloid)
# dev.off()

## Cent~EYO by GFAPcat -------
df<-merge(df_GFAP, df_centiloid[, c("subject_label", "event_code_amy", "latency_amy", "Acat", "WUSTLcentiloid")], by="subject_label", all=FALSE)
scatterEYO_centiloid_byGFAP<-  ggplot(df, aes(x=(age_at_visit-52.5), y=WUSTLcentiloid)) + geom_point(aes(color=GFAPcat)) + geom_smooth(method="gam", aes(linetype=ds_vs_control_flag)) +
  xlim(-30,30) +
  scatterTheme + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(colour = c('black', "#00BFC4", 'black', 'black'), size=22), legend.title = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold")) +  geom_vline(xintercept = 0, linetype = "dashed") + ylab("Amyloid PET (Centiloids)") + xlab("EYO (AAO = 52.5 yrs)")
# jpeg("./output/CentEYO_byGFAP.jpeg", width = 800, height = 600)
# print(scatterEYO_centiloid_byGFAP)
# dev.off()

##GFAP ~ EYO-----

AllSigDiffEYO_GFAP<-SigDiffEYO_GFAP
SigDiffEYO_GFAP<-round(mean(SigDiffEYO_GFAP), 1)

scatterEYO_GFAP<-
  ggplot(df_GFAP, aes(x=(age_at_visit), y=plasma_GFAP_value, color=ds_vs_control_flag)) + geom_point() + geom_smooth(method="gam") +
  scale_colour_manual(values = c("#F8766D", "#00BFC4"), labels = c("Control", "Down Syndrome")) +
  scatterTheme + 
  scale_x_continuous(breaks=c(30, 52.5 + SigDiffEYO_GFAP, 52.5, 60),
                         sec.axis = sec_axis(~ . - 52.5, breaks = c(-22.5, SigDiffEYO_GFAP, 0, 7.5), 
                                             name = "Estimated Years to Symptom Onset"), 
                     name = "Age", limits = c(24.8, 72.5)) +
    theme(legend.position = "bottom", 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          axis.text.x = element_text(colour = c('black', "#00BFC4", 'black',  'black',  'black'), 
                                     size=18), legend.title = element_blank(),
          axis.text.y = element_text(size = 18),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 22),
          plot.title = element_text(size = 24, face = "bold")) + geom_vline(xintercept = 52.5, linetype = "dashed") + 
    ylab("Plasma GFAP (pg/mL)")
# jpeg("./output/GFAPEYO_DSvCtrl_annot.jpeg", width = 800, height = 600)
# print(scatterEYO_GFAP)
# dev.off()

## TauPET ~ EYO-----

AllSigDiffEYO_tauPET<-SigDiffEYO_tauPET
SigDiffEYO_tauPET<-round(mean(SigDiffEYO_tauPET), 1)

scatterEYO_tauPET<-
  ggplot(df_tauPET, aes(x=(age_at_visit), y=TempMeta, color=ds_vs_control_flag)) + geom_point() + geom_smooth(method="gam") +
  scale_colour_manual(values = c( "#F8766D", "#00BFC4"), labels = c("Control", "Down Syndrome")) +
  scale_x_continuous(breaks=c(30, 52.5 + SigDiffEYO_tauPET, 52.5, 60),
                     sec.axis = sec_axis(~ . - 52.5, breaks = c(-22.5, SigDiffEYO_tauPET, 0, 7.5), 
                                         name = "Estimated Years to Symptom Onset"), 
                     name = "Age", limits = c(24.8, 72.5)) +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(colour = c('black', '#00BFC4','black',  'black',  'black'), 
                                   size=18), legend.title = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 22),
        plot.title = element_text(size = 24, face = "bold")) + geom_vline(xintercept = 52.5, linetype = "dashed") +  
  ylab("Tau PET (Temporal-Meta SUVR)") + xlab("Age") 

# jpeg("./output/TauPETEYO_DSvCtrl_annot.jpeg", width = 800, height = 600)
# print(scatterEYO_tauPET)
# dev.off()

## Tau-meta~EYO by GFAPcat -------
df<-merge(df_GFAP, df_tauPET[, c("subject_label", "event_code_tau", "latency_tau", "Tcat", "TempMeta")], by="subject_label", all=FALSE)
scatterEYO_MetatauPET<-ggplot(df, aes(x=(age_at_visit-52.5), y=TempMeta)) + geom_point(aes(color=GFAPcat)) + geom_smooth(method="gam", aes(linetype=ds_vs_control_flag)) +
  xlim(-30,30) +
  scatterTheme + geom_vline(xintercept = 0, linetype = "dashed") + ylab("Tau PET (Temporal-Meta SUVR)") + xlab("EYO (AAO = 52.5 yrs)") + theme(legend.title = element_blank())
# jpeg("./output/MetaTauPETEYO_byGFAP.jpeg", width = 800, height = 600)
# print(scatterEYO_MetatauPET)
# dev.off()

jpeg(file = "./ScatterPlot_EYO.jpg", width = 2000, height = 600)
lemon::grid_arrange_shared_legend(scatterEYO_centiloid, scatterEYO_GFAP, scatterEYO_tauPET)
dev.off()
##pTau217 ~ EYO-----
# z=1.96
# df <-expand.grid("EYO" = seq(from = -25, to = 10, by = 0.1))
# df$EYO<-round(df$EYO, digits = 1)
# itr<-10000
# SigDiffEYO<-matrix(0, itr, 1)
# MaxCentiloid<-matrix(0, itr, 1)
# AOC<-matrix(0, itr, 1)
# Taudata_ctrl<-subset(df_pT217, ds_vs_control_flag=="control")
# Taudata_DS<-subset(df_pT217, ds_vs_control_flag=="DS")
# 
# for (i in 1:itr){
#   sampleDF<-Taudata_ctrl[sample(nrow(Taudata_ctrl), replace = TRUE),]
#   sampleDF<-rbind(sampleDF, Taudata_DS[sample(nrow(Taudata_DS), replace = TRUE),])
#   GAM_EYO<-gam(plasma_Ptau217_value~ds_vs_control_flag + s(EYO, bs="cs", k=4, by=ds_vs_control_flag), data=sampleDF)
#   predictGAM<-predict_gam(GAM_EYO, values=list(EYO=df$EYO))
#   
#   predictctrl<-subset(predictGAM, ds_vs_control_flag=="control")
#   predictDS<-subset(predictGAM, ds_vs_control_flag=="DS")
#   
#   MaxCentiloid[i,1]<-max(predictDS$fit)
#   AOC[i,1]<-geiger:::.area.between.curves(predictDS$EYO, predictctrl$fit, predictDS$fit, xrange = c(-25,0))
#   
#   smoothDiff_DSvCtrl<-get_difference(GAM_EYO, list(ds_vs_control_flag=c("DS", "control")), cond = list(EYO=df$EYO), rm.ranef = FALSE, se = TRUE, f = z)
#   smoothDiff_DSvCtrl$CI_upper<-smoothDiff_DSvCtrl$difference + smoothDiff_DSvCtrl$CI
#   smoothDiff_DSvCtrl$CI_lower<-smoothDiff_DSvCtrl$difference - smoothDiff_DSvCtrl$CI
#   smoothDiff_DSvCtrl$SigDiff_DSvCtrl<-as.factor(ifelse(smoothDiff_DSvCtrl$CI_upper<0, "Significant", ifelse(smoothDiff_DSvCtrl$CI_lower>0, "Significant", "NotSig")))
#   SigEYO_DSvCtrl<-subset(smoothDiff_DSvCtrl, difference>0)
#   SigEYO_DSvCtrl<-subset(SigEYO_DSvCtrl, SigDiff_DSvCtrl=="Significant")
#   SigEYO_DSvCtrl<-SigEYO_DSvCtrl[!(duplicated(SigEYO_DSvCtrl$SigDiff_DSvCtrl)), "EYO"]
#   
#   SigDiffEYO[i,1]<-SigEYO_DSvCtrl
# }

AllSigDiffEYO_pT217<-SigDiffEYO_pT217
SigDiffEYO_pT217<-round(mean(SigDiffEYO_pT217), 1)

scatterEYO_pT217<-
ggplot(df_pT217, aes(x=(age_at_visit), y=plasma_Ptau217_value, color=ds_vs_control_flag)) + geom_point() + geom_smooth(method="gam") +
  scale_colour_manual(values = c("#F8766D", "#00BFC4"), labels = c("Control", "Down Syndrome")) +
  scale_x_continuous(breaks=c(30, 52.5 + SigDiffEYO_pT217, 52.5, 60),
                     sec.axis = sec_axis(~ . - 52.5, breaks = c(-22.5, SigDiffEYO_pT217, 0, 7.5), 
                                         name = "Estimated Years to Symptom Onset"), 
                     name = "Age", limits = c(24.8, 72.5)) +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(colour = c('black', '#00BFC4','black',  'black',  'black'), 
                                   size=18), legend.title = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 22),
        plot.title = element_text(size = 24, face = "bold")) + geom_vline(xintercept = 52.5, linetype = "dashed") + ylab("pTau-217 (pg/mL)") + xlab("EYO (AAO = 52.5 yrs)")
# jpeg("./output/pT217EYO_DSvCtrl_annot.jpeg", width = 800, height = 600)
# print(scatterEYO_pT217)
# dev.off()


lemon::grid_arrange_shared_legend(scatterEYO_centiloid + ggtitle("A."), 
                                  scatterEYO_pT217 + ggtitle("B."),
                                  scatterEYO_GFAP + ggtitle("C."),
                                  scatterEYO_tauPET + ggtitle("D."), 
                                   nrow = 2, ncol = 2)


## pTau217~EYO by GFAPcat -------
df<-merge(df_GFAP, df_pT217[, c("subject_label", "plasma_Ptau217_value")], by="subject_label", all=FALSE)
scatterEYO_pT217_byGFAP<-ggplot(df, aes(x=(age_at_visit-52.5), y=plasma_Ptau217_value)) + geom_point(aes(color=GFAPcat)) + geom_smooth(method="gam", aes(linetype=ds_vs_control_flag)) +
  xlim(-30,30) +
  scatterTheme + geom_vline(xintercept = 0, linetype = "dashed") + ylab("pTau-217 (pg/mL)") + xlab("EYO (AAO = 52.5 yrs)")
# jpeg("./output/pT217EYO_byGFAP.jpeg", width = 800, height = 600)
# print(scatterEYO_pT217_byGFAP)
# dev.off()

