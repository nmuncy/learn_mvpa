

### --- Notes
#
# Target = Face, 1 = Correct
# Foil = Scene, 2 = Correct


# Set up
data_dir <- "/Users/nmuncy/Projects/learn_mvpa/behAnalysis/"
study_list <- c("Study_BE", "Study_CE", "Study_FP")

# get files, subj list
all_list <- list.files(data_dir)
subj_list <- vector(mode = "list", length = length(all_list)/length(study_list))
for(i in 1:length(all_list)){
  subj_list[i] <- gsub("_.*","",all_list[i])
}
subj_list <- unique(subj_list)

# make data frame
for(study in study_list){
  # study <- study_list[1]

  # Make master df
  df_master <- as.data.frame(matrix(NA,nrow = length(subj_list), ncol = 2))
  colnames(df_master) <- c("Subj", "d.prime")
  m_count <- 1
  
  for(subj in subj_list){
    # subj <- subj_list[1]
    df_master$Subj[m_count] <- subj
    
    # read data, calc z-scores
    df_subj_face <- read.csv(paste0(data_dir, subj, "_", study, "_Face.txt"), header = T)
    hit_num <- grep(1, df_subj_face$Pred)
    hit_prop <- length(hit_num)/dim(df_subj_face)[1]
    hit_z <- round(qnorm(hit_prop),3)
    
    df_subj_scene <- read.csv(paste0(data_dir, subj, "_", study, "_Scene.txt"), header = T)
    fa_num <- grep(2, df_subj_face$Pred)
    fa_prop <- length(fa_num)/dim(df_subj_scene)[1]
    fa_z <- round(qnorm(fa_prop),3)
    
    # fill d-prime
    df_master$d.prime[m_count] <- hit_z-fa_z
    m_count <- m_count+1
  }
  
  write.csv(df_master, paste0(data_dir, "Stat_dprime_", gsub(".*_","",study), ".txt"), row.names = F, quote = F)
}




