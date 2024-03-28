library(readxl)
demogRAW_abcds<-read.csv("./data/DEMOGRAPH_abcdsDFAnnaSummer_26July.csv", na.strings = c("", "NA"))
geneRAW_abcds<-read.csv("./data/GENETICS_abcdsDFAnnaSummer.csv", na.strings = c("", "NA"))
clinRAW_abcds<-read.csv("./data/CLINICAL_ds_abcdsDFAnnaSummer.csv", na.strings = c("", "NA"))
AV45RAW_abcds<-read.csv("./data/amyAV45_abcdsDFAnnaSummer.csv", na.strings = c("", "NA"))
PIBRAW_abcds<-read.csv("./data/amyPIB_abcdsDFAnnaSummer.csv", na.strings = c("", "NA"))
FDGRAW_abcds<-read.csv("./data/FDG_abcdsDFAnnaSummer.csv", na.strings = c("", "NA"))
tauRAW_abcds<-read.csv("./data/tauAV1451_abcdsDFAnnaSummer.csv", na.strings = c("", "NA"))
bioRAW_abcds<-read.csv("./data/BIOSPECIMEN_abcdsDFAnnaSummer.csv", na.strings = c("", "NA"))
CentlRAW_abcds<-read.csv("./data/amyWUSTLcentiloid_abcdsDFAnnaSummer.csv", na.strings = c("", "NA"))


demog_abcds<-demogRAW_abcds[,c("subject_label", "event_code", "ds_vs_control_flag", "age_at_visit", "latency_in_days", "de_gender", "de_ethnicity", "de_race")]

gene_abcds<-geneRAW_abcds[!duplicated(geneRAW_abcds$subject_label),]
demog_abcds<-merge(demog_abcds, gene_abcds[,c("subject_label", "karyotype", "allele_combo")], by="subject_label", all=FALSE)
demog_abcds$karyotype<-ifelse(demog_abcds$ds_vs_control_flag=="Control", "control", demog_abcds$karyotype)
## 1=Trisomy 21;2=Mosaicism;3=Translocation;4=Other
demog_abcds$karyotype<-factor(demog_abcds$karyotype, levels = c(1,2,3,4,"control"), labels = c("Trisomy 21", "Mosaicism", "Translocation", "Other", "control"))
demog_abcds$allele_combo<-as.factor(demog_abcds$allele_combo)
demog_abcds$apoe4<-as.factor(ifelse(demog_abcds$allele_combo == "E2/E4"|demog_abcds$allele_combo == "E4/E3"|demog_abcds$allele_combo == "E4/E4", 1, 0))
demog_abcds<-merge(demog_abcds, clinRAW_abcds[,c("subject_label", "event_code", "consensus")], by=c("subject_label", "event_code"), all=TRUE)
# rm(demogRAW_abcds, geneRAW_abcds, clinRAW_abcds)

demog_abcds$event_code<-as.factor(demog_abcds$event_code)
demog_abcds$ds_vs_control_flag<-factor(demog_abcds$ds_vs_control_flag, levels = c("Control", "DS"), labels = c("control", "DS"))
## 1=Male;2=Female
demog_abcds$de_gender<-factor(demog_abcds$de_gender, levels = c(1, 2), labels = c("male", "female"))
## Hispanic/Latino: 1=yes; 0=no
demog_abcds$de_ethnicity<-factor(demog_abcds$de_ethnicity, levels = c(0,1), labels = c("No", "Yes"))
colnames(demog_abcds)[colnames(demog_abcds)=="de_ethnicity"]<-"HISPANIC"
## 1=White ;2=Black/African-American/African/Caribbean/Black British;3=American Indian/Alaskan Native;4=Asian/Asian British;5=Native Hawaiian/Other Pacific Islander;6=Unknown or Not Reported
demog_abcds$de_race<-as.factor(demog_abcds$de_race)
levels(demog_abcds$de_race)[!(levels(demog_abcds$de_race)%in%c(1,2,3,4,5,6))]<-"multi/other"
demog_abcds$de_race<-factor(demog_abcds$de_race, levels = c(1,2,3,4,5,"multi/other",6), labels = c("White", "Black or African American", "American Indian or Alaska Native", "Asian", "Native Hawaiian or Other Pacific Islander", "Multi/Other", "Unknown"))
## 0=No MCI and no Dementia;1=MCI;2=Dementia;3=Unable to Determine
demog_abcds$consensus<-ifelse(demog_abcds$ds_vs_control_flag=="control", "control", demog_abcds$consensus)
demog_abcds$consensus<-factor(demog_abcds$consensus, levels = c("control", 0, 1, 2, 3), labels = c("control", "asym", "MCI", "dem", "noConsensus"))
demog_abcds$cog<-factor(demog_abcds$consensus, levels = c("control", "asym", "MCI", "dem", "noConsensus"), labels = c("control", "aDS", "sDS", "sDS", "ncDS"))
demog_abcds$EYO<-demog_abcds$age_at_visit - 52.5
demog_abcds$latency_in_days<-ifelse(demog_abcds$event_code == "bl", 0, demog_abcds$latency_in_days)

plasma_abcds<-bioRAW_abcds[,c("subject_label", "event_code", "plasma_GFAP_value", "plasma_NT1tau_value", "plasma_pTau181_value", "plasma_Ptau217_value", "plasma_Tau_value", "plasma_Ab40_value")]
plasma_abcds<-plasma_abcds[!with(plasma_abcds, is.na(plasma_GFAP_value) & is.na(plasma_NT1tau_value) & is.na(plasma_pTau181_value) & is.na(plasma_Ptau217_value) & is.na(plasma_Tau_value) & is.na(plasma_Ab40_value)),]

centiloid_abcds<-CentlRAW_abcds[,c("subject_label" ,"event_code", "WUSTLcentiloid")]
centiloid_abcds<-centiloid_abcds[!is.na(centiloid_abcds$WUSTLcentiloid),]


BrianOrigData =read_xlsx("./data/data_for_anna.xlsx", col_types="guess")
ROIs<-colnames(BrianOrigData)[19:60] #only need ROI names
ROIs_tau<-ROIs[ROIs!="CHORPLEX"]
ROIs_tau<-ROIs_tau[ROIs_tau!="VENTRALDC"]
# ROIs_tau<-ROIs_tau[ROIs_tau!="ENTORHINAL"]
PUPregions_rsfSUVR<-colnames(AV45RAW_abcds)[startsWith(colnames(AV45RAW_abcds), "PET_fSUVR_rsf_TOT_")]
PUPregions_rsfSUVR<-PUPregions_rsfSUVR[!startsWith(PUPregions_rsfSUVR, "PET_fSUVR_rsf_TOT_WM")]
PUPregions_rsfSUVR<-PUPregions_rsfSUVR[!startsWith(PUPregions_rsfSUVR, "PET_fSUVR_rsf_TOT_CBLL_")]
# PUPregions_rsfSUVR<-PUPregions_rsfSUVR[!startsWith(PUPregions_rsfSUVR, "PET_fSUVR_rsf_TOT_CORTMEAN")]

TempMeta_ROIs<-c("PET_fSUVR_rsf_TOT_CTX_ENTORHINAL",
              "PET_fSUVR_rsf_TOT_CTX_INFERTMP",
              "PET_fSUVR_rsf_TOT_CTX_MIDTMP",
              "PET_fSUVR_rsf_TOT_CTX_FUSIFORM",
              "PET_fSUVR_rsf_TOT_CTX_PARAHPCMPL",
              "PET_fSUVR_rsf_TOT_AMYGDALA")
tauopath_ROIs<-c("PET_fSUVR_rsf_TOT_AMYGDALA", "PET_fSUVR_rsf_TOT_CTX_ENTORHINAL", "PET_fSUVR_rsf_TOT_CTX_INFERTMP", "PET_fSUVR_rsf_TOT_CTX_LATOCC")
Braak_ROIS<-unique(c(TempMeta_ROIs, tauopath_ROIs))

AV45_abcds<-AV45RAW_abcds[,c("subject_label", "event_code", PUPregions_rsfSUVR)]
colnames(AV45_abcds)<-gsub("PET_fSUVR_rsf_TOT_CTX_","",colnames(AV45_abcds))
colnames(AV45_abcds)<-gsub("PET_fSUVR_rsf_TOT_","",colnames(AV45_abcds))
AV45_abcds<-AV45_abcds[,c("subject_label", "event_code", ROIs, "CORTMEAN")]

PIB_abcds<-PIBRAW_abcds[,c("subject_label", "event_code", PUPregions_rsfSUVR)]
colnames(PIB_abcds)<-gsub("PET_fSUVR_rsf_TOT_CTX_","",colnames(PIB_abcds))
colnames(PIB_abcds)<-gsub("PET_fSUVR_rsf_TOT_","",colnames(PIB_abcds))
PIB_abcds<-PIB_abcds[,c("subject_label", "event_code", ROIs, "CORTMEAN")]
PIB_abcds<-PIB_abcds[!is.na(PIB_abcds$SSTSBANK),]

tau_abcds<-tauRAW_abcds[,c("subject_label", "event_code", PUPregions_rsfSUVR)]
colnames(tau_abcds)<-gsub("PET_fSUVR_rsf_TOT_CTX_","",colnames(tau_abcds))
colnames(tau_abcds)<-gsub("PET_fSUVR_rsf_TOT_","",colnames(tau_abcds))
tau_abcds<-tau_abcds[,c("subject_label", "event_code", ROIs_tau)]
tau_abcds<-tau_abcds[!is.na(tau_abcds$SSTSBANK),]

braak_abcds<-tauRAW_abcds[,c("subject_label", "event_code", Braak_ROIS)]
braak_abcds<-braak_abcds[!is.na(braak_abcds$PET_fSUVR_rsf_TOT_CTX_ENTORHINAL),]
braak_abcds$TempMeta<-rowMeans(braak_abcds[,TempMeta_ROIs], na.rm=FALSE)
braak_abcds$tauopathy<-rowMeans(braak_abcds[,tauopath_ROIs], na.rm=FALSE)


rm(AV45RAW_abcds,PIBRAW_abcds, tauRAW_abcds, FDGRAW_abcds, BrianOrigData)
rm(bioRAW_abcds, CentlRAW_abcds, clinRAW_abcds, demogRAW_abcds, geneRAW_abcds, gene_abcds)
