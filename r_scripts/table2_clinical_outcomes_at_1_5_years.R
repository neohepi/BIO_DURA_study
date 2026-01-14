################################################################################
# Outcome Table Generation
# 1-Year and 5-Year Cumulative Incidence with HR
################################################################################

library(survival)
library(dplyr)

################################################################################
# 1. 분석 함수 정의
################################################################################

# Cumulative Incidence 및 95% CI 계산 함수
get_cumulative_incidence <- function(fit, time_point, group_idx) {
  # survfit summary에서 특정 시점의 survival 추출
  surv_summary <- summary(fit, times = time_point, extend = TRUE)
  
  # Cumulative incidence = 1 - survival
  surv <- surv_summary$surv[group_idx]
  lower <- surv_summary$lower[group_idx]
  upper <- surv_summary$upper[group_idx]
  
  ci <- 1 - surv
  ci_lower <- 1 - upper  # 역전: survival의 upper가 CI의 lower
  
  ci_upper <- 1 - lower
  
  # Percentage로 변환
  return(list(
    ci = ci * 100,
    ci_lower = ci_lower * 100,
    ci_upper = ci_upper * 100
  ))
}

# HR 및 95% CI 계산 함수 (특정 기간)
get_hr_at_timepoint <- function(data, time_var, status_var, time_cap) {
  # Time capping
  analysis_data <- data %>%
    mutate(
      time_capped = pmin(!!sym(time_var), time_cap),
      status_capped = if_else(!!sym(time_var) > time_cap, 0, !!sym(status_var))
    )
  
  # Cox model
  cox_model <- tryCatch({
    coxph(Surv(time_capped, status_capped) ~ BP, data = analysis_data)
  }, error = function(e) NULL)
  
  if (is.null(cox_model)) {
    return(list(hr = NA, ci_lower = NA, ci_upper = NA, p = NA))
  }
  
  summ <- summary(cox_model)
  hr <- summ$conf.int[1, "exp(coef)"]
  ci_lower <- summ$conf.int[1, "lower .95"]
  ci_upper <- summ$conf.int[1, "upper .95"]
  p <- summ$coefficients[1, "Pr(>|z|)"]
  
  return(list(hr = hr, ci_lower = ci_lower, ci_upper = ci_upper, p = p))
}

# 단일 Outcome 분석 함수
analyze_single_outcome <- function(data, time_var, status_var, outcome_name) {
  
  # 데이터 준비 (5년 cap)
  analysis_data <- data %>%
    mutate(
      time_5yr = pmin(!!sym(time_var), 1825),
      status_5yr = if_else(!!sym(time_var) > 1825, 0, as.numeric(!!sym(status_var)))
    )
  
  # Survfit (5년 기준)
  fit <- survfit(Surv(time_5yr, status_5yr) ~ BP, data = analysis_data)
  
  # 1년 Cumulative Incidence
  ci_1yr_dp <- get_cumulative_incidence(fit, 365, 1)  # BP=0 (DP-DES)
  ci_1yr_bp <- get_cumulative_incidence(fit, 365, 2)  # BP=1 (BP-DES)
  
  # 5년 Cumulative Incidence
  ci_5yr_dp <- get_cumulative_incidence(fit, 1825, 1)
  ci_5yr_bp <- get_cumulative_incidence(fit, 1825, 2)
  
  # 1년 HR
  hr_1yr <- get_hr_at_timepoint(data, time_var, status_var, 365)
  
  # 5년 HR
  hr_5yr <- get_hr_at_timepoint(data, time_var, status_var, 1825)
  
  # 결과 반환
  return(data.frame(
    Outcome = outcome_name,
    
    # 1-Year
    `1yr_DP_CI` = sprintf("%.1f%% (%.1f–%.1f%%)", 
                          ci_1yr_dp$ci, ci_1yr_dp$ci_lower, ci_1yr_dp$ci_upper),
    `1yr_BP_CI` = sprintf("%.1f%% (%.1f–%.1f%%)", 
                          ci_1yr_bp$ci, ci_1yr_bp$ci_lower, ci_1yr_bp$ci_upper),
    `1yr_HR` = sprintf("%.2f (%.2f–%.2f)", hr_1yr$hr, hr_1yr$ci_lower, hr_1yr$ci_upper),
    `1yr_P` = sprintf("%.3f", hr_1yr$p),
    
    # 5-Year
    `5yr_DP_CI` = sprintf("%.1f%% (%.1f–%.1f%%)", 
                          ci_5yr_dp$ci, ci_5yr_dp$ci_lower, ci_5yr_dp$ci_upper),
    `5yr_BP_CI` = sprintf("%.1f%% (%.1f–%.1f%%)", 
                          ci_5yr_bp$ci, ci_5yr_bp$ci_lower, ci_5yr_bp$ci_upper),
    `5yr_HR` = sprintf("%.2f (%.2f–%.2f)", hr_5yr$hr, hr_5yr$ci_lower, hr_5yr$ci_upper),
    `5yr_P` = sprintf("%.3f", hr_5yr$p),
    
    stringsAsFactors = FALSE
  ))
}


################################################################################
# 2. Time 변수 생성 (fu_days_correct 사용)
################################################################################

# fu_days_correct가 없으면 생성
if (!"fu_days_correct" %in% names(matched_data)) {
  matched_data <- matched_data %>%
    mutate(fu_days_correct = as.numeric(last_fu_date - CAG_date))
}

# 각 outcome별 time 변수 생성
matched_data <- matched_data %>%
  mutate(
    # Composite (TLF)
    time_to_composite = case_when(
      composite_event == 1 ~ as.numeric(composite_date - CAG_date),
      TRUE ~ fu_days_correct
    ),
    
    # Cardiac death
    time_to_cardiac_death = case_when(
      cardiac_death == 1 ~ as.numeric(cardiac_death_date - CAG_date),
      TRUE ~ fu_days_correct
    ),
    
    # TVMI
    time_to_TVMI = case_when(
      TVMI == 1 ~ as.numeric(TVMI_date - CAG_date),
      TRUE ~ fu_days_correct
    ),
    
    # TLR
    time_to_TLR = case_when(
      TLR == 1 ~ as.numeric(TLR_date - CAG_date),
      TRUE ~ fu_days_correct
    ),
    
    # Any death
    time_to_any_death = case_when(
      any_death == 1 ~ as.numeric(death_date - CAG_date),
      TRUE ~ fu_days_correct
    ),
    
    # ST
    time_to_ST = case_when(
      ST_outcome == 1 ~ as.numeric(ST_date - CAG_date),
      TRUE ~ fu_days_correct
    )
  )

# NA/음수 처리
matched_data <- matched_data %>%
  mutate(across(starts_with("time_to_"), ~ if_else(is.na(.) | . <= 0, 0.1, .)))

# Event 변수를 numeric으로 변환
matched_data <- matched_data %>%
  mutate(
    composite_event = as.numeric(as.character(composite_event)),
    cardiac_death = as.numeric(as.character(cardiac_death)),
    TVMI = as.numeric(as.character(TVMI)),
    TLR = as.numeric(as.character(TLR)),
    any_death = as.numeric(as.character(any_death)),
    ST_outcome = as.numeric(as.character(ST_outcome)),
    BP = as.numeric(as.character(BP))
  )


################################################################################
# 3. Outcome Table 생성
################################################################################

cat("\n", rep("=", 80), "\n")
cat("GENERATING OUTCOME TABLE\n")
cat(rep("=", 80), "\n\n")

# 각 Outcome 분석
result_tlf <- analyze_single_outcome(
  matched_data, "time_to_composite", "composite_event", "Target lesion failure"
)

result_cd <- analyze_single_outcome(
  matched_data, "time_to_cardiac_death", "cardiac_death", "  Cardiac death"
)

result_tvmi <- analyze_single_outcome(
  matched_data, "time_to_TVMI", "TVMI", "  Target vessel MI"
)

result_tlr <- analyze_single_outcome(
  matched_data, "time_to_TLR", "TLR", "  Clinically driven TLR"
)

result_death <- analyze_single_outcome(
  matched_data, "time_to_any_death", "any_death", "All-cause death"
)

result_st <- analyze_single_outcome(
  matched_data, "time_to_ST", "ST_outcome", "Stent thrombosis"
)

# 결과 합치기
outcome_table <- bind_rows(
  result_tlf,
  result_cd,
  result_tvmi,
  result_tlr,
  result_death,
  result_st
)

# 열 이름 정리
colnames(outcome_table) <- c(
  "Outcome",
  "1-Year DP-DES, % (95% CI)",
  "1-Year BP-DES, % (95% CI)", 
  "1-Year HR (95% CI)",
  "P-value (1yr)",
  "5-Year DP-DES, % (95% CI)",
  "5-Year BP-DES, % (95% CI)",
  "5-Year HR (95% CI)",
  "P-value (5yr)"
)

# 출력
cat("\n=== OUTCOME TABLE ===\n\n")
print(outcome_table, row.names = FALSE)

# CSV 저장
write.csv(outcome_table, "Table_Outcomes_1yr_5yr.csv", row.names = FALSE)


################################################################################
# 4. 깔끔한 형식으로 출력
################################################################################

cat("\n\n", rep("=", 100), "\n")
cat("FORMATTED OUTCOME TABLE\n")
cat(rep("=", 100), "\n\n")

# 헤더 출력
cat(sprintf("%-25s | %-20s | %-20s | %-18s | %-8s | %-20s | %-20s | %-18s | %-8s\n",
            "Outcome", 
            "1-Yr DP-DES", "1-Yr BP-DES", "1-Yr HR", "P",
            "5-Yr DP-DES", "5-Yr BP-DES", "5-Yr HR", "P"))
cat(rep("-", 160), "\n")

# 각 행 출력
for (i in 1:nrow(outcome_table)) {
  cat(sprintf("%-25s | %-20s | %-20s | %-18s | %-8s | %-20s | %-20s | %-18s | %-8s\n",
              outcome_table$Outcome[i],
              outcome_table$`1-Year DP-DES, % (95% CI)`[i],
              outcome_table$`1-Year BP-DES, % (95% CI)`[i],
              outcome_table$`1-Year HR (95% CI)`[i],
              outcome_table$`P-value (1yr)`[i],
              outcome_table$`5-Year DP-DES, % (95% CI)`[i],
              outcome_table$`5-Year BP-DES, % (95% CI)`[i],
              outcome_table$`5-Year HR (95% CI)`[i],
              outcome_table$`P-value (5yr)`[i]))
}

cat(rep("=", 160), "\n")


################################################################################
# 5. Event Counts 추가 (Optional)
################################################################################

cat("\n\n=== EVENT COUNTS ===\n\n")

event_counts <- data.frame(
  Outcome = c("Target lesion failure", "Cardiac death", "Target vessel MI", 
              "Clinically driven TLR", "All-cause death", "Stent thrombosis"),
  DP_Events = c(
    sum(matched_data$composite_event == 1 & matched_data$BP == 0 & matched_data$time_to_composite <= 1825, na.rm = TRUE),
    sum(matched_data$cardiac_death == 1 & matched_data$BP == 0 & matched_data$time_to_cardiac_death <= 1825, na.rm = TRUE),
    sum(matched_data$TVMI == 1 & matched_data$BP == 0 & matched_data$time_to_TVMI <= 1825, na.rm = TRUE),
    sum(matched_data$TLR == 1 & matched_data$BP == 0 & matched_data$time_to_TLR <= 1825, na.rm = TRUE),
    sum(matched_data$any_death == 1 & matched_data$BP == 0 & matched_data$time_to_any_death <= 1825, na.rm = TRUE),
    sum(matched_data$ST_outcome == 1 & matched_data$BP == 0 & matched_data$time_to_ST <= 1825, na.rm = TRUE)
  ),
  BP_Events = c(
    sum(matched_data$composite_event == 1 & matched_data$BP == 1 & matched_data$time_to_composite <= 1825, na.rm = TRUE),
    sum(matched_data$cardiac_death == 1 & matched_data$BP == 1 & matched_data$time_to_cardiac_death <= 1825, na.rm = TRUE),
    sum(matched_data$TVMI == 1 & matched_data$BP == 1 & matched_data$time_to_TVMI <= 1825, na.rm = TRUE),
    sum(matched_data$TLR == 1 & matched_data$BP == 1 & matched_data$time_to_TLR <= 1825, na.rm = TRUE),
    sum(matched_data$any_death == 1 & matched_data$BP == 1 & matched_data$time_to_any_death <= 1825, na.rm = TRUE),
    sum(matched_data$ST_outcome == 1 & matched_data$BP == 1 & matched_data$time_to_ST <= 1825, na.rm = TRUE)
  ),
  DP_N = sum(matched_data$BP == 0),
  BP_N = sum(matched_data$BP == 1)
)

event_counts <- event_counts %>%
  mutate(
    DP_Rate = sprintf("%d/%d (%.1f%%)", DP_Events, DP_N, DP_Events/DP_N*100),
    BP_Rate = sprintf("%d/%d (%.1f%%)", BP_Events, BP_N, BP_Events/BP_N*100)
  )

print(event_counts %>% select(Outcome, DP_Rate, BP_Rate))

write.csv(event_counts, "Table_Event_Counts.csv", row.names = FALSE)

cat("\n", rep("=", 80), "\n")
cat("Tables saved: Table_Outcomes_1yr_5yr.csv, Table_Event_Counts.csv\n")
cat(rep("=", 80), "\n")