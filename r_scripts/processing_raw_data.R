if (!require("readxl")) install.packages("readxl")
if (!require("tidyverse")) install.packages("tidyverse")

library(tidyverse)
library(readxl)

# One-hot encoding된 lesion anatomy column (LM, pLAD, mLAD, ...) => collapsing
collapse_lesion_anatomy <- function(df) {
  return (df %>% 
            mutate(lesion_anatomy = case_when(
              LM == 1 ~ "LM",
              pLAD == 1 ~ "pLAD",
              mLAD == 1 ~ "mLAD",
              dLAD == 1 ~ "dLAD",
              Dg == 1 ~ "Dg",
              septal_br == 1 ~ "septal_br",
              pLCX == 1 ~ "pLCX",
              dLCX == 1 ~ "dLCX",
              OM == 1 ~ "OM",
              RI == 1 ~ "RI",
              pRCA == 1 ~ "pRCA",
              mRCA == 1 ~ "mRCA",
              dRCA == 1 ~ "dRCA",
              PDA == 1 ~ "PDA",
              PLB == 1 ~ "PLB"),
              lesion_anatomy = as.factor(lesion_anatomy)
              )
          )
}

simplify_lesion_anatomy <- function(df) {
  return (df %>%
    mutate(lesion_anatomy = case_when(
      # 1. Diagonal Branch 통합
      lesion_anatomy %in% c("D1", "D2", "D3") ~ "Dg",
      
      # 2. Septal Branch 통합
      lesion_anatomy %in% c("S1", "S2", "S3") ~ "septal_br",
      
      # 3. Obtuse Marginal Branch 통합
      lesion_anatomy %in% c("OM1", "OM2", "OM3") ~ "OM",
      
      # [중요] 위 조건에 해당하지 않는 나머지는 원래 값 유지
      TRUE ~ lesion_anatomy)
      )
  )
}

# Restore the lesion level data from the saved excel file
get_lesion_level_data <- function(lesion_data_file="excel_files/df_backup_251220.xlsx") {
  df <- read_excel(lesion_data_file)
  
  df <- collapse_lesion_anatomy(df)
  
  df <- df %>%
    mutate(across(c("cath_no", "pt_id", "cath_no"), as.integer))
  
  return(df)
}

# DB에서 PCI data를 불러들임
# 이유는 기존 data에는 DP vs BP data만 있으므로 revascularization data를 추가로 확보하기 위해서
get_db_pci_data <- function(pci_data_file="excel_files/outcome_new_251220.xlsx") {
  df <- read_excel(pci_data_file)
  
  df <- simplify_lesion_anatomy(df)
  
  # column selection and renaming
  df <- df %>%
    select(pt_id, cath_no, lesion_anatomy, ISR_lesion, is_ACS, ACS_type, 
           # outcomes
           ISR, ST, TLR, TLR_date, TVR, TVR_date) %>%
    rename(ISR_outcome=ISR, ST_outcome=ST)
  
  # datatype conversion
  df <- df %>%
    mutate(across(c(pt_id, cath_no), as.integer),
           across(c(lesion_anatomy, ISR_lesion, 
                    ISR_outcome, ST_outcome, TLR, TVR), as.factor),
           across(c(TLR_date, TVR_date), as.Date))
  
  return(df)
}

# DB에서 CAG 정보를 불러오기, PCI를 시행하지 않았던 경우도 확인하기 위해서
get_db_cag_data <- function(cag_data_file="excel_files/cag_data.xlsx") {
  df <- read_excel(cag_data_file)
  
  df <- df %>%
    mutate(across(c("cath_no", "pt_id"), as.integer),
           across(c("CAG_date"), as.Date),
           across(c("is_ACS", "ACS_type", "num_CAOD", "L_LAD", "L_LM", "L_LCX", "L_RCA",
                    "is_PCI_done", "rec_CABG", "IVUS_used", "OCT_used", "VH_used", "TVC_used",
                    "MDCT_prev_done"), as.factor))
  
  return(df)
}

get_join_data_with_pci <- function(df, pci_data) {
  df.joined <- df %>% 
    left_join(pci_data, by = c("cath_no", "lesion_anatomy")) %>%
    # removing duplicated columns and renaming renamed columns in the join process
    select(-c("TLR.x", "TLR_date.x", "TV_MI", "TV_MI_date", "pt_id.y")) %>%
    rename(pt_id=pt_id.x, TLR=TLR.y, TLR_date=TLR_date.y) %>%
    mutate(across(c("CAG_date", "death_date"), as.Date))
  
  return(df.joined)
}

df.lesion.all <- get_join_data_with_pci(get_lesion_level_data(), get_db_pci_data())

df.cag.all <- get_db_cag_data()

################################################################################
# Study flow diagram 작성용

all_pcis <- get_db_pci_data() # 32,243
df.cag.all

count_unique_id(all_pcis$pt_id) # 16,190

count_unique_id(df.lesion.all$pt_id) # 9,586

sum(df.pts$BP == 0) # DP
sum(df.pts$BP == 1) # BP

