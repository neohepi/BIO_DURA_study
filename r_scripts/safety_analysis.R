# 1. matched_data에 'death_30d' 변수를 먼저 생성하여 업데이트(저장)합니다.
matched_data <- matched_data %>%
  mutate(
    time_to_death = as.numeric(death_date - CAG_date),
    # 사망했고(1), 사망 시점이 30일 이내인 경우만 1, 나머지는 0
    # (NA 처리 안전장치: time_to_death가 NA여도 0이 되도록)
    death_30d = if_else(any_death == 1 & !is.na(time_to_death) & time_to_death <= 30, 1, 0)
  )

# 2. 결과 요약 (생성된 변수 활용)
result_30d_mortality <- matched_data %>%
  group_by(BP) %>%
  summarize(
    n = n(),
    n_deaths = sum(death_30d, na.rm = TRUE),
    rate = n_deaths / n * 100
  )
print(result_30d_mortality)

# 3. Statistical Test (이제 matched_data 안에 death_30d가 있으므로 작동합니다)
table_30d <- table(matched_data$BP, matched_data$death_30d)

print(table_30d)

# Chi-square test (혹은 Fisher's exact test)
chisq_result <- chisq.test(table_30d)
print(chisq_result)

# 만약 빈도가 적어 Fisher test가 필요하다면:
# fisher.test(table_30d)

# Odds ratio
library(epitools)
oddsratio(table_30d)

################################################################################
# Competing risk analysis with death
install.packages("cmprsk" ,"tidycmprsk")
library(cmprsk)
library(tidycmprsk)

library(survival)
library(tidycmprsk)
library(dplyr)

# 1. Competing Risk용 데이터 전처리
matched_data_cr <- matched_data %>%
  mutate(
    # (A) 사망까지의 시간 계산 (사망 안 했으면 fu_days)
    time_to_death_clean = if_else(any_death == 1, as.numeric(death_date - CAG_date), as.numeric(fu_days)),
    time_to_death_clean = if_else(time_to_death_clean <= 0, 0.1, time_to_death_clean), # 0 보정
    
    # (B) 365일 이내 발생 여부 확인 (1년 시점 Landmark)
    # MI가 1년 내 발생했으면 그 날짜, 아니면 NA
    t_mi_1yr = if_else(next_MI == 1 & time_to_MI <= 365, time_to_MI, NA_real_),
    # 사망이 1년 내 발생했으면 그 날짜, 아니면 NA
    t_death_1yr = if_else(any_death == 1 & time_to_death_clean <= 365, time_to_death_clean, NA_real_),
    
    # (C) Event Type 결정 (Time to FIRST Event)
    # MI와 사망이 둘 다 기록된 경우, 날짜가 더 빠른 것을 우선시함
    final_event = case_when(
      # 1. 둘 다 1년 내 발생 -> 더 빨리 온 것이 진짜 이벤트
      !is.na(t_mi_1yr) & !is.na(t_death_1yr) & t_mi_1yr <= t_death_1yr ~ "MI",
      !is.na(t_mi_1yr) & !is.na(t_death_1yr) & t_death_1yr < t_mi_1yr ~ "Death",
      
      # 2. 하나만 발생
      !is.na(t_mi_1yr) ~ "MI",
      !is.na(t_death_1yr) ~ "Death",
      
      # 3. 아무 일도 없음 (Censored)
      TRUE ~ "Censored"
    ),
    
    # (D) Time 결정 (결정된 이벤트의 시간 사용)
    final_time = case_when(
      final_event == "MI" ~ t_mi_1yr,
      final_event == "Death" ~ t_death_1yr,
      # Censored인 경우: 1년(365일) 혹은 그 전 추적종료일 중 빠른 날
      TRUE ~ pmin(fu_days, 365, na.rm = TRUE)
    ),
    
    # [핵심] Factor 변환 (반드시 Censored가 첫 번째 레벨이어야 함)
    event_factor = factor(final_event, levels = c("Censored", "MI", "Death"))
  )

# 2. Cumulative Incidence Function 실행
cif <- cuminc(Surv(final_time, event_factor) ~ BP, 
              data = matched_data_cr)

# 3. 결과 출력 (Table)
print(cif)

# 1. 시각화 전용 패키지 설치 및 로드
if(!require("ggsurvfit")) install.packages("ggsurvfit")
library(ggsurvfit)

# 2. 그래프 그리기 (옵션 조금 더 추가해서 예쁘게)
cif %>%
  ggcuminc(outcome = "MI") +             # 관심 Outcome 지정
  labs(
    title = "Cumulative Incidence of MI (1-year Landmark)",
    subtitle = "Competing Risk Analysis (Death as competing event)",
    x = "Days since PCI",
    y = "Cumulative Incidence Probability"
  ) +
  add_confidence_interval() +            # 95% 신뢰구간 (음영)
  add_risktable(                         # Risk Table 추가
    risktable_stats = "n.risk"           # 남은 환자 수만 표시
  ) +
  add_pvalue(caption = "Gray's Test") +  # P-value 자동 추가 (중요!)
  scale_ggsurvfit()                      # 깔끔한 테마 적용

################################################################################
# ST analysis
# Stent Thrombosis Analysis with Landmark
library(dplyr)
library(lubridate)

library(dplyr)
library(lubridate)

create_ST_variables_corrected <- function(data,
                                          index_date_col = "CAG_date",
                                          st_date_col = "ST_date", 
                                          last_fu_col = "last_FU_date",
                                          st_indicator_col = "ST_outcome", 
                                          max_followup_days = 1825) {  # 5 years
  
  result_data <- data %>%
    mutate(
      # 날짜 형식 변환
      index_date = as.Date(.data[[index_date_col]]),
      st_date_event = as.Date(.data[[st_date_col]]), # st_date 이름 충돌 방지
      last_fu_date = as.Date(.data[[last_fu_col]]),
      
      # 1. Raw Time 계산 (제한 없는 실제 일수)
      raw_days_to_event = case_when(
        # ST 발생 시: ST 날짜 기준
        !is.na(.data[[st_indicator_col]]) & .data[[st_indicator_col]] == 1 & !is.na(st_date_event) ~ 
          as.numeric(st_date_event - index_date),
        
        # ST 미발생 시: 마지막 추적일 기준
        TRUE ~ as.numeric(last_fu_date - index_date)
      ),
      
      # 음수 값 보정 (데이터 오류 방지)
      raw_days_to_event = pmax(raw_days_to_event, 0, na.rm = TRUE),
      
      # 2. Time Truncation (최대 기간으로 자르기)
      time_to_ST = pmin(raw_days_to_event, max_followup_days),
      
      # 3. Outcome Censoring (가장 중요!)
      # 원래 이벤트가 있었더라도(1), 기간(1825일)을 넘겨서 발생했다면 0으로 변경
      ST_outcome_censored = case_when(
        !is.na(.data[[st_indicator_col]]) & .data[[st_indicator_col]] == 1 & raw_days_to_event <= max_followup_days ~ 1,
        TRUE ~ 0
      ),
      
      # 기존 변수명과의 호환성을 위해 ST_outcome으로 저장 (필요시 변경 가능)
      ST_outcome = ST_outcome_censored,
      
      # 4. ST Timing 재분류 (Outcome이 0이 되면 Timing도 수정 필요)
      ST_timing = case_when(
        ST_outcome == 1 & time_to_ST <= 1 ~ "Acute (≤24h)",
        ST_outcome == 1 & time_to_ST <= 30 ~ "Subacute (1-30d)",
        ST_outcome == 1 & time_to_ST <= 365 ~ "Late (30d-1y)",
        ST_outcome == 1 & time_to_ST > 365 ~ "Very Late (>1y)",
        TRUE ~ "No ST" # Censored 된 경우 포함
      )
    )
  
  # 요약 출력
  cat("\n=== ST Variable Creation (Censored at", max_followup_days, "days) ===\n")
  cat("Total patients:", nrow(result_data), "\n")
  
  # Censoring 전후 비교 (데이터 검증용)
  original_events <- sum(data[[st_indicator_col]], na.rm=TRUE)
  new_events <- sum(result_data$ST_outcome, na.rm=TRUE)
  
  cat("Original ST events:", original_events, "\n")
  cat("Final ST events (within 5 years):", new_events, "\n")
  cat("-> Events censored (occurred > 5 years):", original_events - new_events, "\n\n")
  
  cat("ST rate (5-year):", round(mean(result_data$ST_outcome, na.rm = TRUE) * 100, 2), "%\n")
  
  return(result_data)
}

# 사용 예시
matched_data_ST <- create_ST_variables_corrected(
  data = matched_data,
  index_date_col = "CAG_date",
  st_date_col = "ST_date",
  last_fu_col = "last_fu_date",
  st_indicator_col = "ST_outcome",
  max_followup_days = 1825
)


# ============================================================================
# Most flexible version: Multiple ST definitions
# ============================================================================

create_comprehensive_ST_variables <- function(data,
                                              index_date_col = "CAG_date",
                                              last_fu_col = "last_FU_date",
                                              definite_ST_col = "definite_ST",
                                              probable_ST_col = "probable_ST",
                                              definite_ST_date_col = NULL,
                                              probable_ST_date_col = NULL,
                                              max_followup_days = 1825) {
  
  result_data <- data %>%
    mutate(
      index_date = as.Date(.data[[index_date_col]]),
      last_fu_date = as.Date(.data[[last_fu_col]])
    )
  
  # Definite ST
  if (!is.null(definite_ST_date_col) && definite_ST_date_col %in% names(data)) {
    result_data <- result_data %>%
      mutate(
        definite_ST_date = as.Date(.data[[definite_ST_date_col]]),
        time_to_definite_ST = case_when(
          !is.na(.data[[definite_ST_col]]) & .data[[definite_ST_col]] == 1 & !is.na(definite_ST_date) ~
            as.numeric(definite_ST_date - index_date),
          TRUE ~ as.numeric(last_fu_date - index_date)
        ),
        time_to_definite_ST = pmin(pmax(time_to_definite_ST, 0), max_followup_days),
        definite_ST_outcome = ifelse(!is.na(.data[[definite_ST_col]]) & 
                                       .data[[definite_ST_col]] == 1, 1, 0)
      )
  } else {
    result_data <- result_data %>%
      mutate(
        time_to_definite_ST = as.numeric(last_fu_date - index_date),
        time_to_definite_ST = pmin(pmax(time_to_definite_ST, 0), max_followup_days),
        definite_ST_outcome = ifelse(!is.na(.data[[definite_ST_col]]) & 
                                       .data[[definite_ST_col]] == 1, 1, 0)
      )
  }
  
  # Definite or Probable ST (composite)
  result_data <- result_data %>%
    mutate(
      # Composite ST outcome
      ST_composite = case_when(
        (!is.na(.data[[definite_ST_col]]) & .data[[definite_ST_col]] == 1) |
          (!is.na(.data[[probable_ST_col]]) & .data[[probable_ST_col]] == 1) ~ 1,
        TRUE ~ 0
      ),
      
      # Time to composite ST
      time_to_ST = case_when(
        ST_composite == 1 & !is.na(definite_ST_date) ~ 
          as.numeric(definite_ST_date - index_date),
        ST_composite == 1 ~ as.numeric(last_fu_date - index_date),
        TRUE ~ as.numeric(last_fu_date - index_date)
      ),
      
      time_to_ST = pmin(pmax(time_to_ST, 0), max_followup_days),
      ST_outcome = ST_composite
    )
  
  # Summary
  cat("\n=== Comprehensive ST Analysis ===\n")
  cat("Total patients:", nrow(result_data), "\n\n")
  
  cat("Definite ST:\n")
  cat("  Events:", sum(result_data$definite_ST_outcome, na.rm = TRUE), "\n")
  cat("  Rate:", round(mean(result_data$definite_ST_outcome, na.rm = TRUE) * 100, 2), "%\n\n")
  
  cat("Definite or Probable ST:\n")
  cat("  Events:", sum(result_data$ST_outcome, na.rm = TRUE), "\n")
  cat("  Rate:", round(mean(result_data$ST_outcome, na.rm = TRUE) * 100, 2), "%\n\n")
  
  cat("Follow-up duration:\n")
  cat("  Median:", median(result_data$time_to_ST, na.rm = TRUE), "days\n")
  cat("  Mean:", round(mean(result_data$time_to_ST, na.rm = TRUE), 1), "days\n")
  cat("  Max:", max(result_data$time_to_ST, na.rm = TRUE), "days\n\n")
  
  return(result_data)
}

# ST 발생 날짜가 별도로 기록된 경우
matched_data <- create_ST_variables(
  data = matched_data,
  index_date_col = "CAG_date",
  st_date_col = "ST_date",          # ST 발생일
  last_fu_col = "last_fu_date",
  st_indicator_col = "ST_outcome",          # ST 발생 여부 (0/1)
  max_followup_days = 1825
)

# 결과 확인
head(matched_data %>% select(pt_id, CAG_date, ST_date, last_fu_date, 
                             time_to_ST, ST_outcome, ST_timing))

################################################################################
library(dplyr)
library(lubridate)

create_ST_variables_corrected <- function(data,
                                          index_date_col = "CAG_date",
                                          st_date_col = "ST_date", 
                                          last_fu_col = "last_FU_date",
                                          st_indicator_col = "ST_outcome", 
                                          max_followup_days = 1825) {  # 5 years
  
  result_data <- data %>%
    mutate(
      # 날짜 형식 변환
      index_date = as.Date(.data[[index_date_col]]),
      st_date_event = as.Date(.data[[st_date_col]]), # st_date 이름 충돌 방지
      last_fu_date = as.Date(.data[[last_fu_col]]),
      
      # 1. Raw Time 계산 (제한 없는 실제 일수)
      raw_days_to_event = case_when(
        # ST 발생 시: ST 날짜 기준
        !is.na(.data[[st_indicator_col]]) & .data[[st_indicator_col]] == 1 & !is.na(st_date_event) ~ 
          as.numeric(st_date_event - index_date),
        
        # ST 미발생 시: 마지막 추적일 기준
        TRUE ~ as.numeric(last_fu_date - index_date)
      ),
      
      # 음수 값 보정 (데이터 오류 방지)
      raw_days_to_event = pmax(raw_days_to_event, 0, na.rm = TRUE),
      
      # 2. Time Truncation (최대 기간으로 자르기)
      time_to_ST = pmin(raw_days_to_event, max_followup_days),
      
      
      # 3. Outcome Censoring (가장 중요!)
      # 원래 이벤트가 있었더라도(1), 기간(1825일)을 넘겨서 발생했다면 0으로 변경
      ST_outcome_censored = case_when(
        !is.na(.data[[st_indicator_col]]) & .data[[st_indicator_col]] == 1 & raw_days_to_event <= max_followup_days ~ 1,
        TRUE ~ 0
      ),
      
      # 기존 변수명과의 호환성을 위해 ST_outcome으로 저장 (필요시 변경 가능)
      ST_outcome = ST_outcome_censored,
      
      # 4. ST Timing 재분류 (Outcome이 0이 되면 Timing도 수정 필요)
      ST_timing = case_when(
        ST_outcome == 1 & time_to_ST <= 1 ~ "Acute (≤24h)",
        ST_outcome == 1 & time_to_ST <= 30 ~ "Subacute (1-30d)",
        ST_outcome == 1 & time_to_ST <= 365 ~ "Late (30d-1y)",
        ST_outcome == 1 & time_to_ST > 365 ~ "Very Late (>1y)",
        TRUE ~ "No ST" # Censored 된 경우 포함
      )
    )
  
  # 요약 출력
  cat("\n=== ST Variable Creation (Censored at", max_followup_days, "days) ===\n")
  cat("Total patients:", nrow(result_data), "\n")
  
  # Censoring 전후 비교 (데이터 검증용)
  original_events <- sum(data[[st_indicator_col]], na.rm=TRUE)
  new_events <- sum(result_data$ST_outcome, na.rm=TRUE)
  
  cat("Original ST events:", original_events, "\n")
  cat("Final ST events (within 5 years):", new_events, "\n")
  cat("-> Events censored (occurred > 5 years):", original_events - new_events, "\n\n")
  
  cat("ST rate (5-year):", round(mean(result_data$ST_outcome, na.rm = TRUE) * 100, 2), "%\n")
  
  return(result_data)
}

# 사용 예시
matched_data <- create_ST_variables_corrected(
  data = matched_data,
  index_date_col = "CAG_date",
  st_date_col = "ST_date",
  last_fu_col = "last_fu_date",
  st_indicator_col = "ST_outcome",
  max_followup_days = 1825
)


