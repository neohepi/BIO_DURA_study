################################################################################
# Supplementary Table S2: Landmark Analysis (1-5 Years)
# 기존 survival_analysis_corrected.R 실행 후 이 코드 실행
################################################################################

library(survival)
library(dplyr)

cat("\n", rep("=", 70), "\n")
cat("SUPPLEMENTARY TABLE S2: LANDMARK ANALYSIS (1-5 YEARS)\n")
cat(rep("=", 70), "\n\n")

################################################################################
# 1. Time 변수 확인/생성
################################################################################

# fu_days_correct 생성 (없으면)
if (!"fu_days_correct" %in% names(matched_data)) {
  matched_data <- matched_data %>%
    mutate(fu_days_correct = as.numeric(last_fu_date - CAG_date))
}

# 각 outcome별 time 변수 생성/확인
matched_data <- matched_data %>%
  mutate(
    # Composite (TLF)
    time_to_composite = case_when(
      composite_event == 1 ~ as.numeric(composite_date - CAG_date),
      TRUE ~ fu_days_correct
    ),
    time_to_composite = pmin(time_to_composite, 1825),
    time_to_composite = if_else(is.na(time_to_composite) | time_to_composite <= 0, 0.1, time_to_composite),
    
    # Cardiac death
    time_to_cardiac_death = case_when(
      cardiac_death == 1 ~ as.numeric(cardiac_death_date - CAG_date),
      TRUE ~ fu_days_correct
    ),
    time_to_cardiac_death = pmin(time_to_cardiac_death, 1825),
    time_to_cardiac_death = if_else(is.na(time_to_cardiac_death) | time_to_cardiac_death <= 0, 0.1, time_to_cardiac_death),
    
    # TVMI
    time_to_TVMI = case_when(
      TVMI == 1 ~ as.numeric(TVMI_date - CAG_date),
      TRUE ~ fu_days_correct
    ),
    time_to_TVMI = pmin(time_to_TVMI, 1825),
    time_to_TVMI = if_else(is.na(time_to_TVMI) | time_to_TVMI <= 0, 0.1, time_to_TVMI),
    
    # TLR
    time_to_TLR = case_when(
      TLR == 1 ~ as.numeric(TLR_date - CAG_date),
      TRUE ~ fu_days_correct
    ),
    time_to_TLR = pmin(time_to_TLR, 1825),
    time_to_TLR = if_else(is.na(time_to_TLR) | time_to_TLR <= 0, 0.1, time_to_TLR),
    
    # Any death
    time_to_any_death = case_when(
      any_death == 1 ~ as.numeric(death_date - CAG_date),
      TRUE ~ fu_days_correct
    ),
    time_to_any_death = pmin(time_to_any_death, 1825),
    time_to_any_death = if_else(is.na(time_to_any_death) | time_to_any_death <= 0, 0.1, time_to_any_death)
  )

# Event 변수를 numeric으로
matched_data <- matched_data %>%
  mutate(
    composite_event = as.numeric(as.character(composite_event)),
    cardiac_death = as.numeric(as.character(cardiac_death)),
    TVMI = as.numeric(as.character(TVMI)),
    TLR = as.numeric(as.character(TLR)),
    any_death = as.numeric(as.character(any_death)),
    BP = as.numeric(as.character(BP))
  )


################################################################################
# 2. Landmark Analysis 함수 (1년 ~ 5년)
################################################################################

analyze_landmark_period <- function(data, time_var, event_var, 
                                    landmark_start = 365, landmark_end = 1825) {
  
  # 1년 시점에 event-free인 환자만 포함
  # time > landmark_start 이거나 event가 0인 경우
  data_landmark <- data %>%
    filter(!!sym(time_var) > landmark_start | !!sym(event_var) == 0)
  
  # 1년 이후 ~ 5년 사이에 발생한 event만 계산
  data_landmark <- data_landmark %>%
    mutate(
      # 1년 이후부터의 시간
      time_lm = pmin(!!sym(time_var), landmark_end) - landmark_start,
      time_lm = if_else(time_lm <= 0, 0.1, time_lm),
      
      # 1년 이후 5년 이내에 발생한 event
      event_lm = if_else(
        !!sym(event_var) == 1 & 
          !!sym(time_var) > landmark_start & 
          !!sym(time_var) <= landmark_end, 
        1, 0
      )
    )
  
  # Survfit for cumulative incidence at 4 years (1460 days from landmark)
  fit <- survfit(Surv(time_lm, event_lm) ~ BP, data = data_landmark)
  
  # 4년 시점 (1825 - 365 = 1460일)의 cumulative incidence
  analysis_time <- landmark_end - landmark_start
  surv_summary <- summary(fit, times = analysis_time, extend = TRUE)
  
  # DP-DES (BP=0) 
  if (length(surv_summary$surv) >= 1) {
    ci_dp <- (1 - surv_summary$surv[1]) * 100
    ci_dp_lower <- (1 - surv_summary$upper[1]) * 100
    ci_dp_upper <- (1 - surv_summary$lower[1]) * 100
  } else {
    ci_dp <- ci_dp_lower <- ci_dp_upper <- NA
  }
  
  # BP-DES (BP=1)
  if (length(surv_summary$surv) >= 2) {
    ci_bp <- (1 - surv_summary$surv[2]) * 100
    ci_bp_lower <- (1 - surv_summary$upper[2]) * 100
    ci_bp_upper <- (1 - surv_summary$lower[2]) * 100
  } else {
    ci_bp <- ci_bp_lower <- ci_bp_upper <- NA
  }
  
  # Cox model for HR
  cox_model <- tryCatch({
    coxph(Surv(time_lm, event_lm) ~ BP, data = data_landmark)
  }, error = function(e) NULL)
  
  if (!is.null(cox_model)) {
    summ <- summary(cox_model)
    hr <- summ$conf.int[1, "exp(coef)"]
    hr_lower <- summ$conf.int[1, "lower .95"]
    hr_upper <- summ$conf.int[1, "upper .95"]
    p_value <- summ$coefficients[1, "Pr(>|z|)"]
  } else {
    hr <- hr_lower <- hr_upper <- p_value <- NA
  }
  
  # Event counts
  events_dp <- sum(data_landmark$event_lm[data_landmark$BP == 0], na.rm = TRUE)
  events_bp <- sum(data_landmark$event_lm[data_landmark$BP == 1], na.rm = TRUE)
  n_dp <- sum(data_landmark$BP == 0)
  n_bp <- sum(data_landmark$BP == 1)
  
  return(list(
    ci_dp = sprintf("%.1f%% (%.1f-%.1f%%)", ci_dp, ci_dp_lower, ci_dp_upper),
    ci_bp = sprintf("%.1f%% (%.1f-%.1f%%)", ci_bp, ci_bp_lower, ci_bp_upper),
    hr = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
    p_value = sprintf("%.3f", p_value),
    events_dp = events_dp,
    events_bp = events_bp,
    n_dp = n_dp,
    n_bp = n_bp
  ))
}


################################################################################
# 3. 각 Outcome에 대해 Landmark 분석 수행
################################################################################

cat("Analyzing landmark period (365 days to 1825 days)...\n\n")

landmark_tlf <- analyze_landmark_period(matched_data, "time_to_composite", "composite_event")
landmark_cd <- analyze_landmark_period(matched_data, "time_to_cardiac_death", "cardiac_death")
landmark_tvmi <- analyze_landmark_period(matched_data, "time_to_TVMI", "TVMI")
landmark_tlr <- analyze_landmark_period(matched_data, "time_to_TLR", "TLR")
landmark_death <- analyze_landmark_period(matched_data, "time_to_any_death", "any_death")


################################################################################
# 4. 결과 테이블 생성
################################################################################

landmark_table <- data.frame(
  Outcome = c("TLF", "Cardiac death", "TV-MI", "Clinically driven TLR", "All-cause death"),
  DP_DES = c(landmark_tlf$ci_dp, landmark_cd$ci_dp, landmark_tvmi$ci_dp, 
             landmark_tlr$ci_dp, landmark_death$ci_dp),
  BP_DES = c(landmark_tlf$ci_bp, landmark_cd$ci_bp, landmark_tvmi$ci_bp,
             landmark_tlr$ci_bp, landmark_death$ci_bp),
  HR = c(landmark_tlf$hr, landmark_cd$hr, landmark_tvmi$hr,
         landmark_tlr$hr, landmark_death$hr),
  P_value = c(landmark_tlf$p_value, landmark_cd$p_value, landmark_tvmi$p_value,
              landmark_tlr$p_value, landmark_death$p_value),
  stringsAsFactors = FALSE
)


################################################################################
# 5. 출력
################################################################################

cat("\n")
cat("Supplementary Table S2. Landmark Analysis of Clinical Outcomes Between 1 and 5 Years\n")
cat(rep("=", 100), "\n")
cat(sprintf("%-25s | %-20s | %-20s | %-18s | %-8s\n",
            "Outcome", "DP-DES, % (95% CI)", "BP-DES, % (95% CI)", "HR (95% CI)", "P-value"))
cat(rep("-", 100), "\n")

for (i in 1:nrow(landmark_table)) {
  cat(sprintf("%-25s | %-20s | %-20s | %-18s | %-8s\n",
              landmark_table$Outcome[i],
              landmark_table$DP_DES[i],
              landmark_table$BP_DES[i],
              landmark_table$HR[i],
              landmark_table$P_value[i]))
}

cat(rep("=", 100), "\n")
cat("Analysis includes patients event-free at 1 year (landmark analysis)\n")
cat("HR = Hazard Ratio for BP-DES vs DP-DES (reference)\n")
cat(rep("=", 100), "\n")


################################################################################
# 6. Event Counts 추가 정보
################################################################################

cat("\n\n=== EVENT COUNTS (Landmark Period: 1-5 Years) ===\n\n")

event_counts <- data.frame(
  Outcome = c("TLF", "Cardiac death", "TV-MI", "Clinically driven TLR", "All-cause death"),
  DP_Events = c(landmark_tlf$events_dp, landmark_cd$events_dp, landmark_tvmi$events_dp,
                landmark_tlr$events_dp, landmark_death$events_dp),
  DP_N = c(landmark_tlf$n_dp, landmark_cd$n_dp, landmark_tvmi$n_dp,
           landmark_tlr$n_dp, landmark_death$n_dp),
  BP_Events = c(landmark_tlf$events_bp, landmark_cd$events_bp, landmark_tvmi$events_bp,
                landmark_tlr$events_bp, landmark_death$events_bp),
  BP_N = c(landmark_tlf$n_bp, landmark_cd$n_bp, landmark_tvmi$n_bp,
           landmark_tlr$n_bp, landmark_death$n_bp)
)

event_counts <- event_counts %>%
  mutate(
    DP_Rate = sprintf("%d/%d", DP_Events, DP_N),
    BP_Rate = sprintf("%d/%d", BP_Events, BP_N)
  )

cat(sprintf("%-25s | %-15s | %-15s\n", "Outcome", "DP-DES (n/N)", "BP-DES (n/N)"))
cat(rep("-", 60), "\n")
for (i in 1:nrow(event_counts)) {
  cat(sprintf("%-25s | %-15s | %-15s\n",
              event_counts$Outcome[i],
              event_counts$DP_Rate[i],
              event_counts$BP_Rate[i]))
}


################################################################################
# 7. CSV 저장
################################################################################

colnames(landmark_table) <- c("Outcome", "DP-DES, % (95% CI)", "BP-DES, % (95% CI)", 
                              "HR (95% CI)", "P-value")
write.csv(landmark_table, "Table_S2_Landmark_Analysis.csv", row.names = FALSE)

cat("\n\n✓ Table saved: Table_S2_Landmark_Analysis.csv\n")
cat(rep("=", 70), "\n")

