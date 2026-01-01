library(MatchIt)
library(dplyr)
library(tableone)
library(cobalt)
library(skimr)

# -----------------------------------------------------------------------------
# 함수 1: 데이터 전처리 (Era 생성 및 Treatment 변수 정리)
# -----------------------------------------------------------------------------
prepare_psm_data <- function(data, treat_col = "BP", year_col = "procedure_year") {
  data %>%
    # 1. Treatment 변수 Factor 변환 (0, 1 보장)
    mutate(
      !!sym(treat_col) := as.integer(as.character(!!sym(treat_col))),
      !!sym(treat_col) := factor(!!sym(treat_col), levels = c(0, 1))
    ) %>%
    filter(!is.na(!!sym(treat_col)))
  # era는 이미 생성되었으므로 사용하지 않음
}

# -----------------------------------------------------------------------------
# 함수 2: PSM 실행 및 매칭된 데이터 반환
# -----------------------------------------------------------------------------
run_psm <- function(data, treat_col, match_vars, keep_vars, 
                    ratio = 1, caliper = 0.05, seed = 1234) {
  
  # 1. 데이터 필터링 (공변량에 NA가 있는 행만 제거)
  # match_vars와 keep_vars를 합친 후 중복 제거
  target_cols <- unique(c(match_vars, keep_vars, treat_col))
  
  data_clean <- data %>%
    select(all_of(target_cols)) %>%
    filter(if_all(all_of(match_vars), ~ !is.na(.)))
  
  message(paste("Original N:", nrow(data)))
  message(paste("N for PSM (NA removed):", nrow(data_clean)))
  
  # 2. Formula 생성
  ps_formula <- as.formula(paste(treat_col, "~", paste(match_vars, collapse = " + ")))
  
  # 3. MatchIt 실행
  set.seed(seed)
  m.out <- matchit(
    formula = ps_formula,
    data = data_clean,
    method = "nearest", 
    distance = "glm",
    ratio = ratio,
    caliper = caliper,
    std.caliper = TRUE
  )
  
  # 4. 결과 요약 출력
  print(summary(m.out))
  
  # 5. 매칭된 데이터 추출
  matched_data <- match.data(m.out)
  
  return(list(model = m.out, matched_data = matched_data))
}

# -----------------------------------------------------------------------------
# 함수 3: 밸런스 확인 (Table 1 & Love Plot)
# -----------------------------------------------------------------------------
check_balance <- function(psm_result, vars_to_check, cat_vars) {
  
  m.out <- psm_result$model
  matched_data <- psm_result$matched_data
  
  library(cobalt)
  
  # 1. 변수명 매핑 리스트 생성 (순서는 상관없음)
  new_labels <- c(
    "distance" = "Propensity Score",
    "age" = "Age (years)",
    "male" = "Male Sex",
    "BMI" = "BMI",
    "smoking" = "Current Smoking",
    "HTN" = "HTN",
    "DM" = "DM",
    "CKD" = "CKD",
    "DL" = "Dyslipidemia",
    "LVEF" = "LVEF (%)",
    "is_ACS" = "ACS",
    "prior_PCI" = "Previous PCI",
    "prior_MI" = "Previous MI",
    "prior_stroke" = "Previous Stroke",
    
    # 범주형 변수가 원-핫 인코딩된 경우 구체적인 이름 지정 필요
    "num_CAOD_1" = "1-Vessel Disease",
    "num_CAOD_2" = "2-Vessel Disease",
    "num_CAOD_3" = "3-Vessel Disease",
    
    "is_LM" = "Left Main Disease",
    "is_bifurcation" = "Bifurcation Lesion",
    "is_complex_lesion" = "Complex Lesion",
    "is_CTO" = "CTO",
    "is_ISR" = "ISR lesion",
    
    "era_2010-2013" = "Study Period: 2010-2013",
    "era_2014-2016" = "Study Period: 2014-2016",
    "era_2017-2021" = "Study Period: 2017-2021"
  )
  
  # 2. love.plot 그리기
  print(love.plot(
    x = psm_result$model, # MatchIt 결과 객체 등
    stat = "mean.diffs",
    abs = TRUE,
    threshold = 0.1,
    
    # ★ 여기서 이름 변경 적용
    var.names = new_labels, 
    
    # (선택) Propensity Score(distance)를 그래프에서 빼고 싶다면
    # drop.distance = TRUE 
  ))
  
  
  # 2. Matched Table 1 생성
  table1 <- CreateTableOne(vars = vars_to_check, 
                           strata = "BP", 
                           data = matched_data, 
                           factorVars = cat_vars,
                           test = TRUE)
  
  # 결과 출력
  print(table1, smd = TRUE, showAllLevels = FALSE)
  
  return(table1)
}

# ==============================================================================
# 1. 설정 (Configuration)
# ==============================================================================
# 원본 데이터 로드
df.toPSM <- df.pts # patient-level data

df.toPSM <- add_procedural_data(df.toPSM, df.lesion.all)

# 매칭에 사용할 변수 (NA 있으면 안됨)
match_vars_list <- c("age", "male", "BMI", "smoking",
                     "HTN", "DM", "CKD", "DL",
                     "LVEF", "is_ACS", "prior_PCI", "prior_MI", "prior_stroke", 
                     #"prior_CABG", 
                     "num_CAOD", "is_LM", "is_bifurcation", "is_complex_lesion", "is_CTO", "is_ISR", 
                     "era")

# 결과 분석에 필요한 변수 (NA 있어도 됨)
keep_vars_list <- c("pt_id", "procedure_year", "CAG_date", "fu_days", "ACS_type",
                    "TC", "LDL", "HDL", "TG", "prior_CABG",
                    "is_LAD", "is_LCX", "is_RCA", "is_multivessel", 
                    "num_stents", "total_stent_length", "avg_stent_diameter", 
                    "last_fu_date",
                    "TVR", "TVR_date",
                    "TLR", "TLR_date", 
                    "TVMI", "TVMI_date",
                    "cardiac_death", "cardiac_death_date", 
                    "any_death", "death_date", 
                    "next_MI", "next_MI_date",
                    "composite_event", "composite_date",
                    "ST_outcome", "ST_date")

# Table 1에 표시할 변수 (보통 match_vars와 비슷하지만 더 많을 수도 있음)
table1_vars <- match_vars_list
table1_cat_vars <- c("male", "smoking", "HTN", "DM", "CKD", "DL", 
                     "is_ACS", 
                     "prior_PCI", "prior_MI", "prior_stroke",
                     "num_CAOD", 
                     "is_bifurcation", "is_complex_lesion", "is_CTO", "is_ISR",
                     "era")

# ==============================================================================
# 2. 실행 (Execution)
# ==============================================================================

# A. 전처리 (Era 생성 등)
df.ready <- prepare_psm_data(df.toPSM, treat_col = "BP", year_col = "procedure_year")

# B. PSM 실행
psm_result <- run_psm(
  data = df.ready,
  treat_col = "BP",
  match_vars = match_vars_list,
  keep_vars = keep_vars_list,
  caliper = 0.04,
  ratio = 1
)

# C. 결과 확인 (SMD Plot & Table 1)
tbl1_result <- check_balance(
  psm_result = psm_result, 
  vars_to_check = table1_vars, 
  cat_vars = table1_cat_vars
)

# D. 추가 데이터 확인 (환자 수 등)
matched_data <- psm_result$matched_data
print(paste("BP=0 Unique Patients:", count_unique_id(matched_data[matched_data$BP == 0,]$pt_id)))
print(paste("BP=1 Unique Patients:", count_unique_id(matched_data[matched_data$BP == 1,]$pt_id)))

matched_data %>% count(era)

################################################################################
# Creating time-to-outcome variables
################################################################################
# 추가 분석을 위해서 time_to_event varaible 생성
matched_data <- matched_data %>%
  mutate(
    # --------------------------------------------------------------------------
    # 1. TLR (Target Lesion Revascularization)
    # --------------------------------------------------------------------------
    # 날짜 차이 계산
    diff_TLR = as.numeric(TLR_date - CAG_date),
    
    # 시간 변수 결정 (이벤트=1이면 발생일까지, 아니면 추적기간까지)
    time_to_TLR = if_else(
      TLR == 1 & !is.na(diff_TLR), 
      diff_TLR, 
      as.numeric(fu_days)
    ),
    
    # 0 이하 보정 (Safety buffer)
    time_to_TLR = if_else(time_to_TLR <= 0, 0.1, time_to_TLR),
    
    # --------------------------------------------------------------------------
    # 2. TVR (Target Vessel Revascularization)
    # --------------------------------------------------------------------------
    # 날짜 차이 계산
    diff_TVR = as.numeric(TVR_date - CAG_date),
    
    # 시간 변수 결정
    time_to_TVR = if_else(
      TVR == 1 & !is.na(diff_TVR), 
      diff_TVR, 
      as.numeric(fu_days)
    ),
    
    # 0 이하 보정
    time_to_TVR = if_else(time_to_TVR <= 0, 0.1, time_to_TVR)
    
  ) %>%
  # 계산에 사용한 임시 날짜 차이 컬럼 제거
  select(-diff_TLR, -diff_TVR)

# 결과 확인 (제대로 생성되었는지 앞부분만 출력)
matched_data %>% 
  select(pt_id, TLR, time_to_TLR, TVR, time_to_TVR, fu_days) %>% 
  head()



