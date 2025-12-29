# 1. 패키지 로드 (없으면 install.packages("gtsummary") 실행 필요)
library(dplyr)
library(gtsummary)
library(dplyr)
library(tidyr)
library(dplyr)
library(tidyr)

add_procedural_data <- function(data, df_lesion) {
  
  # 1. 시술 부위(Anatomy) 매핑 정의
  target_LAD <- c("pLAD", "mLAD", "dLAD", "Dg", "septal_br") 
  target_LCX <- c("pLCX", "dLCX", "OM", "RI")
  target_RCA <- c("pRCA", "mRCA", "dRCA", "PLB", "PDA") 
  target_LM  <- c("LM")
  
  # 2. Index Procedure에 해당하는 병변만 추출 및 요약
  proc_summary <- df_lesion %>%
    # 환자의 Index procedure 날짜와 일치하는 스텐트만 필터링
    inner_join(data %>% select(pt_id, CAG_date), by = c("pt_id", "CAG_date")) %>%
    
    group_by(pt_id) %>%
    summarise(
      # A. Procedural Characteristics (Continuous)
      num_stents = n(),
      total_stent_length = sum(stent_length, na.rm = TRUE),
      avg_stent_diameter = mean(stent_diameter, na.rm = TRUE),
      
      # B. Target Vessel (0 or 1)
      is_LM = max(ifelse(lesion_anatomy %in% target_LM, 1, 0), na.rm = TRUE),
      is_LAD = max(ifelse(lesion_anatomy %in% target_LAD, 1, 0), na.rm = TRUE),
      is_LCX = max(ifelse(lesion_anatomy %in% target_LCX, 1, 0), na.rm = TRUE),
      is_RCA = max(ifelse(lesion_anatomy %in% target_RCA, 1, 0), na.rm = TRUE),
      
      # C. Lesion Complexity (0 or 1)
      is_bifurcation = max(ifelse(!is.na(bifurcation) & bifurcation == 1, 1, 0), na.rm = TRUE),
      is_complex_lesion = max(ifelse(ACC_type %in% c("B2", "C"), 1, 0), na.rm = TRUE),
      is_CTO = max(ifelse(!is.na(CTO) & CTO == 1, 1, 0), na.rm = TRUE),
      is_ISR = max(ifelse(!is.na(ISR_lesion) & ISR_lesion == 1, 1, 0), na.rm = TRUE)
      
    ) %>%
    mutate(
      # D. Multivessel PCI 정의 (0 or 1)
      vessel_count = is_LAD + is_LCX + is_RCA,
      is_multivessel = ifelse(vessel_count >= 2, 1, 0)
    ) %>%
    ungroup()
  
  # 3. 원본 데이터에 붙이기 (Left Join)
  data_final <- left_join(data, proc_summary, by = "pt_id") %>%
    mutate(
      # 매칭 안 된 환자(NA)를 0으로 처리 (이미 1인 값은 그대로 1 유지)
      num_stents = ifelse(is.na(num_stents), 0, num_stents),
      total_stent_length = ifelse(is.na(total_stent_length), 0, total_stent_length),
      avg_stent_diameter = ifelse(is.na(avg_stent_diameter), 0, avg_stent_diameter),
      
      # Binary 변수들: NA -> 0으로 변환 (Factor 변환 제거)
      is_LM = ifelse(is.na(is_LM), 0, is_LM),
      is_LAD = ifelse(is.na(is_LAD), 0, is_LAD),
      is_LCX = ifelse(is.na(is_LCX), 0, is_LCX),
      is_RCA = ifelse(is.na(is_RCA), 0, is_RCA),
      
      is_bifurcation = ifelse(is.na(is_bifurcation), 0, is_bifurcation),
      is_complex_lesion = ifelse(is.na(is_complex_lesion), 0, is_complex_lesion),
      is_CTO = ifelse(is.na(is_CTO), 0, is_CTO),
      is_ISR = ifelse(is.na(is_ISR), 0, is_ISR),
      is_multivessel = ifelse(is.na(is_multivessel), 0, is_multivessel)
    ) %>%
    mutate(
      num_stents = as.numeric(as.character(num_stents))
    )
  
  return(data_final)
}

# 2. Baseline Table 생성 함수 정의
get_baseline_table <- function(data, group_var = NULL) {
  # 추가 column 생성
  data <- data %>%
    mutate(num_CAOD = as.numeric(as.character(num_CAOD)))
  
  data <- add_procedural_data(data, df.lesion.all)
  
  data <- data %>%
    mutate(LVEFlt40 = factor(ifelse(LVEF < 40, 1, 0)),
           # 진단명을 하나의 컬럼으로 정리
           clinical_presentation = case_when(
             ACS_type == "STEMI" ~ "STEMI",
             ACS_type == "NSTEMI" ~ "NSTEMI",
             (ACS_type == "UA" | ACS_type == "Unstable Angina") ~ "Unstable Angina",
             # ACS가 아니면(또는 Stable CAD라고 명시되어 있다면) Stable CAD로 분류
             # (데이터 상황에 맞춰 조건 수정 필요)
             is.na(ACS_type) | ACS_type == "" | is_ACS == 0 ~ "Stable CAD",
             TRUE ~ NA_character_ 
           ),
           
           # 순서 지정 (Stable을 맨 위로 올릴지, STEMI를 위로 올릴지 결정)
           clinical_presentation = factor(clinical_presentation, 
                                          levels = c("STEMI", "NSTEMI", "Unstable Angina", "Stable CAD")),
           
           # Study Period (기존 유지)
           study_period = case_when(
             procedure_year >= 2010 & procedure_year <= 2013 ~ "2010-2013",
             procedure_year >= 2014 & procedure_year <= 2017 ~ "2014-2017",
             procedure_year >= 2018 & procedure_year <= 2021 ~ "2018-2021",
             TRUE ~ NA_character_),
           
           num_CAOD = case_when(
             is.na(num_CAOD) | num_CAOD == 1 ~ "1 vessel",  # NA이거나 1이면 '1 vessel'
             num_CAOD == 2 ~ "2 vessels",                   # 2면 '2 vessels'
             num_CAOD >= 3 ~ "3 vessels"                    # 3 이상이면 '3 vessels'
           )) %>%
    mutate(num_CAOD = factor(num_CAOD, levels = c("1 vessel", "2 vessels", "3 vessels")),
           target_vessel = case_when(
             is_LM  == 1 ~ "LM",
             is_LAD == 1 ~ "LAD",
             is_LCX == 1 ~ "LCX",
             is_RCA == 1 ~ "RCA"
           ),
           target_vessel = factor(target_vessel, 
                                  levels = c("LM", "LAD", "LCX", "RCA"), 
                                  labels = c("Left main", "Left anterior descending", "Left circumflex", "Right coronary artery")))
  
  # 약제 정보 추가
  temp_drug <- df.lesion.all %>% select(pt_id, antiplatelet, anticoagulants, `_statin`, beta_blocker, ARB_ACEi) %>% distinct()
  
  data <- left_join(data, temp_drug, by="pt_id") %>%
    mutate(ARB_ACEi = as.factor(ifelse(ARB_ACEi >= 1, 1, 0)))
  
  # 논문에 들어갈 변수 리스트 (계속 추가하시면 됩니다)
  # 괄호 안의 column name과 data의 column name이 일치해야 합니다.
  vars_to_include <- c("age", "male", "BMI",       # Demographics
                       
                       "HTN", "DM", "DL", "smoking", "CKD", # Risk Factors
                       
                       "prior_MI", "prior_PCI", "prior_CABG", "prior_stroke",
                       "LVEF", "LVEFlt40", # Cardiac history
                       
                       "clinical_presentation",
                       
                       "TC", "LDL", "HDL", "TG", # lab
                       
                       "antiplatelet", "anticoagulants", "_statin", "beta_blocker", "ARB_ACEi",
                       
                       "num_CAOD",
                       
                       "is_LM", "target_vessel", "is_complex_lesion", "is_bifurcation", "is_CTO", "is_ISR",
                       
                       "is_multivessel", "num_stents", "total_stent_length", "avg_stent_diameter",
                       
                       "study_period") 
  
  # 테이블 생성
  tbl <- data %>%
    select(all_of(c(vars_to_include, group_var))) %>% # 필요한 변수만 선택
    tbl_summary(
      by = all_of(group_var), # 그룹 변수 (예: 'treatment') 없으면 NULL
      
      # 1. 통계량 표시 형식 지정 (연속형: Mean ± SD, 범주형: N (%))
      statistic = list(
        all_continuous() ~ "{mean} \u00B1 {sd}", 
        all_categorical() ~ "{n} ({p}%)"
      ),
      
      # 'Yes'인 경우만 출력하도록 설정 (이분형 변수 깔끔하게)
      type = list(
        c(male, HTN, DM, DL, smoking, CKD,
          prior_MI, prior_PCI, prior_CABG, prior_stroke, LVEFlt40,
          antiplatelet, anticoagulants, `_statin`, beta_blocker, ARB_ACEi,
          is_LM, is_complex_lesion, is_bifurcation, is_CTO, is_ISR, 
          is_multivessel
          ) ~ "dichotomous",
        num_stents ~ "continuous"
      ),
      value = list(
        c(male, HTN, DM, DL, smoking, CKD,
          prior_MI, prior_PCI, prior_CABG, prior_stroke, LVEFlt40,
          antiplatelet, anticoagulants, `_statin`, beta_blocker, ARB_ACEi,
          is_LM, is_complex_lesion, is_bifurcation, is_CTO, is_ISR, 
          is_multivessel) ~ 1
      ),
      
      # 2. 결측치 처리 (논문에서는 보통 결측치 줄을 따로 표시하지 않음)
      missing = "no",
      
      # 3. 변수 라벨링 (논문에 출력될 이름)
      label = list(
        age ~ "Age (years)",
        male ~ "Male sex",
        BMI ~ "Body Mass Index (kg/m²)",
        HTN ~ "Hypertension",
        DM ~ "Diabetes Mellitus",
        DL ~ "Dyslipidemia",
        smoking ~ "Current Smoking",
        CKD ~ "Chronic Kidney Disease",
        prior_MI ~ "Previous myocardial infarction",
        prior_PCI ~ "Previous PCI",
        prior_CABG ~ "Previous CABG",
        prior_stroke ~ "Previous stroke",
        LVEF ~ "LVEF",
        LVEFlt40 ~ "LVEF <40%",
        clinical_presentation ~ "Clinical Presentation",
        TC ~ "Total cholesterol",
        LDL ~ "LDL cholesterol",
        HDL ~ "HDL cholesterol",
        TG ~ "Triglycerides",
        antiplatelet ~ "Antiplatelet", 
        anticoagulants ~ "Anticoagulant", 
        `_statin` ~ "Statin", 
        beta_blocker ~ "Beta-blocker", 
        ARB_ACEi ~ "ACE inhibitor/ARB",
        num_CAOD ~ "Number of diseased vessels",
        is_LM ~ "Left main involvement",
        target_vessel ~ "Target vessel",
        is_complex_lesion ~ "Type B2/C lesion", 
        is_bifurcation ~ "Bifurcation lesion", 
        is_CTO ~ "Chronic total occlusion", 
        is_ISR ~ "In-stent restenosis", 
        is_multivessel ~ "Multivessel PCI",
        num_stents ~ "Number of stents per patient",
        total_stent_length ~ "Total stent length",
        avg_stent_diameter ~ "Average stent diameter",
        study_period ~ "Study Period"
      )
    ) %>%
    
    # 4. P-value 추가
    add_p(pvalue_fun = ~style_pvalue(.x, digits = 3)) %>%
    
    # 5. 전체 컬럼 추가
    add_overall() %>%
    
    # 6. 중요: PSM 논문용 SMD(Standardized Mean Difference) 추가
    add_stat_label(location = "row") %>%
    bold_labels()
  
  # 그룹 변수가 있을 때만 SMD를 계산해서 붙임
  if (!is.null(group_var)) {
    tbl <- tbl %>% add_difference(test = everything() ~ "smd") 
  }
  
  return(tbl)
}

table1 <- get_baseline_table(df.pts, group_var = "BP")

# 결과 출력
table1

df.pts %>% left_join(df.lesion.all, by=c("pt_id", "cath_no"))

df.lesion.all$CTO
