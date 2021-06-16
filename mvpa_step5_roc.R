library("lme4")


### --- Notes
#
# Target = Face, 1 = Correct
# Foil = Scene, 2 = Correct


# Set up
data_dir <- "/Users/nmuncy/Projects/learn_mvpa/mvpa_pred/"
study_list <- c("Study_single")
cat_list <- c("all", "Face", "Scene")

# get files, subj list
all_list <- list.files(data_dir)
subj_list <- vector(
  mode = "list", 
  length = length(all_list)/length(study_list)
)
for(i in 1:length(all_list)){
  subj_list[i] <- gsub("_.*","",all_list[i])
}
subj_list <- unique(subj_list)


# logistic regression for all & e/category
# for(cat in ccat_list){
#   
# }
cat <- cat_list[1]
df_cat <- as.data.frame(matrix(NA, nrow=1, ncol=5))
colnames(df_cat) <- c("Subj", "Volume", "True", "Pred", "Acc")

for(subj in subj_list){
  
  df_subj <- read.csv(
    paste0(data_dir, 
     "/", 
     subj, 
     "_", 
     study_list[1], 
     "_", 
     cat, 
     ".txt"
   )
  )
  df_subj$Subj <- rep(subj, dim(df_subj)[1])
  df_subj$Volume <- 1:dim(df_subj)[1]
  
  df_subj$Acc <- NA
  for(row in 1:dim(df_subj)[1]){
    if(df_subj[row,]$Pred == df_subj[row,]$True){
      df_subj[row,]$Acc <- 1
    }else{
      df_subj[row,]$Acc <- 0
    }
  }
  
  df_subj <- df_subj[,c(3,4,1,2,5)]
  df_cat <- rbind(df_cat, df_subj)
}
df_cat <- df_cat[-1,]
row.names(df_cat) <- 1:dim(df_cat)[1]

fit <- glmer(Acc ~ Volume + (1|Subj), data=df_cat, family=binomial(link="logit"))

# Unequal number of volumes, should probably do a bayesian model



# do d-prime for e/subj
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




