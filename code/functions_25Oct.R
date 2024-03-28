get_takeOffPoint <- function(DF_ctrl, DF_DS, z = 1.96, resp){
  sampleDF <- rbind(DF_ctrl[sample(nrow(DF_ctrl), replace = TRUE),],
                    DF_DS[sample(nrow(DF_DS), replace = TRUE),])
  
  GAM_EYO<-gam(get(resp)~ds_vs_control_flag + s(EYO, bs="cs", k=4, by=ds_vs_control_flag), data=sampleDF)
  predictGAM<-predict_gam(GAM_EYO, values=list(EYO=df$EYO))
  
  predictctrl<-predictGAM[predictGAM$ds_vs_control_flag == "control",]
  predictDS<-predictGAM[predictGAM$ds_vs_control_flag == "DS",]
  
  smoothDiff_DSvCtrl<-get_difference(GAM_EYO, list(ds_vs_control_flag=c("DS", "control")), cond = list(EYO=df$EYO), rm.ranef = FALSE, se = TRUE, f = z)
  smoothDiff_DSvCtrl$CI_upper<-smoothDiff_DSvCtrl$difference + smoothDiff_DSvCtrl$CI
  smoothDiff_DSvCtrl$CI_lower<-smoothDiff_DSvCtrl$difference - smoothDiff_DSvCtrl$CI
  smoothDiff_DSvCtrl$SigDiff_DSvCtrl<-as.factor(ifelse(smoothDiff_DSvCtrl$CI_upper<0, "Significant", ifelse(smoothDiff_DSvCtrl$CI_lower>0, "Significant", "NotSig")))
  SigEYO_DSvCtrl<-subset(smoothDiff_DSvCtrl, difference>0)
  SigEYO_DSvCtrl<-subset(SigEYO_DSvCtrl, SigDiff_DSvCtrl=="Significant")
  SigEYO_DSvCtrl<-SigEYO_DSvCtrl[!(duplicated(SigEYO_DSvCtrl$SigDiff_DSvCtrl)), "EYO"]
  
  return(SigEYO_DSvCtrl)
  
}



get_TPP_scores <- function(df, region){
  
  fit <- Mclust(df[[region]], G = 1:2, model = "V")
  
  if(fit$G == 2){
    #calculate group means
    g1_mean <- mean(df[[region]][fit$classification == 1])
    g2_mean <- mean(df[[region]][fit$classification == 2])
    
    tau_data_g <- fit$classification
    tau_data_tpp <- fit$z[,2]
    
    #reclassify the very low SUVRs to be G == 2
    tau_data_g[df[[region]] < g1_mean] = 1
    tau_data_tpp[df[[region]] < g1_mean] = min(fit$z[,2])
    
    #recalculate group means and sd's
    g1_mean = mean(df[[region]][tau_data_g==1])
    g2_mean = mean(df[[region]][tau_data_g==2])
    g1_sd = sd(df[[region]][tau_data_g==1])
    g2_sd = sd(df[[region]][tau_data_g==2])
    
    
    #save z transformed SUVRs
    tau_data_z = (df[[region]] - g1_mean) / g1_sd
    
  }
  
  #if optimal G is 1, save as tpp = 0 or near 0 with random noise
  if(fit$G == 1){
    tau_data_tpp = rnorm(length(df), mean = 0.01, sd = .001) #add some random noise to the tpp values
    tau_data_g = fit$classification
    
    #calculate group mean
    g1_mean = mean(df[[region]][fit$classification==1])
    g1_sd = sd(df[[region]][fit$classification==1])
    #save z transformed SUVRs
    tau_data_z = (df[[region]] - g1_mean) / g1_sd
    
  }
  
  
  return(list(tau_data_g, tau_data_z, tau_data_tpp))}


get_bound_dataframe <- function(TPP, region_names){
  TPP <- dplyr::bind_cols(TPP)
  colnames(TPP) <- region_names
  return(TPP)
}


get_Staged_by_Region <- function(TPP_G_Group, df, GROUP){
  
  sum_TPPs_DS <- colSums(TPP_G_Group)
  sum_TPPs_DS <- data.frame("Region" = names(sum_TPPs_DS),
                            "Count" = as.numeric(sum_TPPs_DS))
  sum_TPPs_DS$Percent <- sum_TPPs_DS$Count / length(df[df$ds_vs_control_flag == GROUP,]$subject_label)
  sum_TPPs_DS <- sum_TPPs_DS[with(sum_TPPs_DS, order(-Count)),]
  sum_TPPs_DS$running_total <- cumsum(sum_TPPs_DS$Count)
  equal_frequency <- c(sum(sum_TPPs_DS$Count) / 4,
                       2 * sum(sum_TPPs_DS$Count) / 4,
                       3 * sum(sum_TPPs_DS$Count) / 4,
                       sum(sum_TPPs_DS$Count))
  sum_TPPs_DS$Stage_DS <- as.factor(ifelse(sum_TPPs_DS$running_total <= equal_frequency[1], 1,
                                           ifelse(sum_TPPs_DS$running_total > equal_frequency[1] & sum_TPPs_DS$running_total <= equal_frequency[2], 2,
                                                  ifelse(sum_TPPs_DS$running_total > equal_frequency[2] & sum_TPPs_DS$running_total <= equal_frequency[3], 3, 4))))
  #Need to check and see if two regions with equal percentage are being arbitrarily split between stages
  #If they are, bump them both into the smaller stage
  min_values <- data.table(sum_TPPs_DS)[, min(Percent), by = Stage_DS]
  sum_TPPs_DS$Stage_DS <- as.factor(ifelse(sum_TPPs_DS$Percent >= min_values$V1[1], 1,
                                           ifelse(sum_TPPs_DS$Percent >= min_values$V1[2] & sum_TPPs_DS$Percent < min_values$V1[1], 2,
                                                  ifelse(sum_TPPs_DS$Percent >= min_values$V1[3] & sum_TPPs_DS$Percent < min_values$V1[2], 3, 4))))
  sum_TPPs_DS$Stage_DS <- ifelse(sum_TPPs_DS$Count == 0, 0, sum_TPPs_DS$Stage_DS)  
  names(sum_TPPs_DS)[5] <- paste0("Stage_", GROUP)
  
  return(sum_TPPs_DS)}


get_Staged_by_ID <- function(sum_TPPs_DS, TPP_G_DS, Stage_Name, df){
  
  Stage_1_Regions <- as.character(sum_TPPs_DS[sum_TPPs_DS[,Stage_Name] == 1, "Region"])
  Stage_2_Regions <- as.character(sum_TPPs_DS[sum_TPPs_DS[,Stage_Name] == 2, "Region"])
  Stage_3_Regions <- as.character(sum_TPPs_DS[sum_TPPs_DS[,Stage_Name] == 3, "Region"])
  Stage_4_Regions <- as.character(sum_TPPs_DS[sum_TPPs_DS[,Stage_Name] == 4, "Region"])
  
  Stage_1_check <- ifelse(rowSums(TPP_G_DS[,names(TPP_G_DS) %in% Stage_1_Regions, drop = FALSE]) / length(Stage_1_Regions) >= 0.5, 1, 0)
  Stage_2_check <- ifelse(rowSums(TPP_G_DS[,names(TPP_G_DS) %in% Stage_2_Regions, drop = FALSE]) / length(Stage_2_Regions) >= 0.5, 1, 0)
  Stage_3_check <- ifelse(rowSums(TPP_G_DS[,names(TPP_G_DS) %in% Stage_3_Regions, drop = FALSE]) / length(Stage_3_Regions) >= 0.5, 1, 0)
  Stage_4_check <- ifelse(rowSums(TPP_G_DS[,names(TPP_G_DS) %in% Stage_4_Regions, drop = FALSE]) / length(Stage_4_Regions) >= 0.5, 1, 0)
  
  Staged_df <- data.frame(cbind(df[, c("subject_label", "event_code", "ds_vs_control_flag", "age_at_visit", "apoe4", "de_gender", "de_race", "consensus")],
                                Stage_1_check, Stage_2_check, Stage_3_check, Stage_4_check))
  
  Staged_df$Stage <-ifelse(rowSums(Staged_df[, c("Stage_1_check", "Stage_2_check", "Stage_3_check",
                                                 "Stage_4_check")]) == 0, "Stage 0",
                           ifelse(Staged_df$Stage_1_check == 1 & rowSums(Staged_df[, c("Stage_2_check", "Stage_3_check",
                                                                                       "Stage_4_check")]) == 0, "Stage 1",
                                  ifelse(Staged_df$Stage_1_check == 1 &Staged_df$Stage_2_check == 1 &
                                           rowSums(Staged_df[, c( "Stage_3_check","Stage_4_check")]) == 0, "Stage 2",
                                         ifelse(Staged_df$Stage_1_check == 1 &Staged_df$Stage_2_check == 1 &
                                                  Staged_df$Stage_3_check == 1 & Staged_df$Stage_4_check == 0, "Stage 3",
                                                ifelse(rowSums(Staged_df[, c("Stage_1_check", "Stage_2_check", "Stage_3_check",
                                                                             "Stage_4_check")]) == 4, "Stage 4", "Discordant")
                                         ))))
  
  Staged_df$Stage_pooled <- ifelse(Staged_df$Stage == "Stage 2" | Staged_df$Stage == "Stage 3", "Stage 2/3",
                                   Staged_df$Stage)
  return(Staged_df)}



