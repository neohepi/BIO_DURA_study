################################################################################
# Table 3. Stent Thrombosis at 1 Year and 5 Years
################################################################################

library(survival)
library(dplyr)
library(knitr)

################################################################################
# 1. 데이터 준비
################################################################################

# fu_days_correct 생성 (없으면)
if (!"fu_days_correct" %in% names(matched_data)) {
  matched_data <- matched_data %>%
    mutate(fu_days_correct = as.numeric(last_fu_date - CAG_date))
}

# time_to_ST 생성/수정
matched_data <- matched_data %>%
  mutate(
    time_to_ST = case_when(
      ST_outcome == 1 ~ as.numeric(ST_date - CAG_date),
      TRUE ~ fu_days_correct
    ),
    time_to_ST = if_else(is.na(time_to_ST) | time_to_ST <= 0, 0.1, time_to_ST),
    
    # 5년 capping
    time_to_ST_5yr = pmin(time_to_ST, 1825),
    ST_outcome_5yr = if_else(time_to_ST > 1825, 0, as.numeric(ST_outcome))
  )

# ST timing 분류
matched_data <- matched_data %>%
  mutate(
    ST_timing = case_when(
      ST_outcome == 1 & time_to_ST <= 1 ~ "Acute (≤24 hours)",
      ST_outcome == 1 & time_to_ST > 1 & time_to_ST <= 30 ~ "Subacute (>24h-30d)",
      ST_outcome == 1 & time_to_ST > 30 & time_to_ST <= 365 ~ "Late (30d-1yr)",
      ST_outcome == 1 & time_to_ST > 365 ~ "Very late (>1yr)",
      TRUE ~ NA_character_
    )
  )

# BP를 numeric으로
matched_data <- matched_data %>%
  mutate(BP = as.numeric(as.character(BP)))


################################################################################
# 2. ST Timing Distribution
################################################################################

cat("\n", rep("=", 70), "\n")
cat("TABLE 3. STENT THROMBOSIS AT 1 YEAR AND 5 YEARS\n")
cat(rep("=", 70), "\n\n")

# Timing별 이벤트 수
timing_counts <- matched_data %>%
  filter(ST_outcome == 1) %>%
  group_by(BP, ST_timing) %>%
  summarise(n = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = BP, values_from = n, values_fill = 0) %>%
  rename(DP_DES = `0`, BP_DES = `1`)

# 순서 정렬
timing_order <- c("Acute (≤24 hours)", "Subacute (>24h-30d)", 
                  "Late (30d-1yr)", "Very late (>1yr)")

timing_counts <- timing_counts %>%
  mutate(ST_timing = factor(ST_timing, levels = timing_order)) %>%
  arrange(ST_timing)


################################################################################
# 3. Cox Regression 분석
################################################################################

# N counts
n_dp <- sum(matched_data$BP == 0)
n_bp <- sum(matched_data$BP == 1)

# --- Overall (0-5 years) ---
cox_overall <- coxph(Surv(time_to_ST_5yr, ST_outcome_5yr) ~ BP, data = matched_data)
summ_overall <- summary(cox_overall)

hr_overall <- summ_overall$conf.int[1, "exp(coef)"]
ci_lower_overall <- summ_overall$conf.int[1, "lower .95"]
ci_upper_overall <- summ_overall$conf.int[1, "upper .95"]
p_overall <- summ_overall$coefficients[1, "Pr(>|z|)"]

events_dp_overall <- sum(matched_data$ST_outcome_5yr == 1 & matched_data$BP == 0)
events_bp_overall <- sum(matched_data$ST_outcome_5yr == 1 & matched_data$BP == 1)

# --- 1 Year (0-365 days) ---
data_1yr <- matched_data %>%
  mutate(
    time_1yr = pmin(time_to_ST, 365),
    event_1yr = if_else(time_to_ST > 365, 0, as.numeric(ST_outcome))
  )

cox_1yr <- coxph(Surv(time_1yr, event_1yr) ~ BP, data = data_1yr)
summ_1yr <- summary(cox_1yr)

hr_1yr <- summ_1yr$conf.int[1, "exp(coef)"]
ci_lower_1yr <- summ_1yr$conf.int[1, "lower .95"]
ci_upper_1yr <- summ_1yr$conf.int[1, "upper .95"]
p_1yr <- summ_1yr$coefficients[1, "Pr(>|z|)"]

events_dp_1yr <- sum(data_1yr$event_1yr == 1 & data_1yr$BP == 0)
events_bp_1yr <- sum(data_1yr$event_1yr == 1 & data_1yr$BP == 1)

# --- Landmark (365 days - 5 years) ---
data_landmark <- matched_data %>%
  filter(time_to_ST > 365 | ST_outcome == 0) %>%
  mutate(
    time_landmark = pmin(time_to_ST, 1825) - 365,
    time_landmark = if_else(time_landmark < 0, 0.1, time_landmark),
    event_landmark = if_else(ST_outcome == 1 & time_to_ST > 365 & time_to_ST <= 1825, 1, 0)
  )

cox_landmark <- coxph(Surv(time_landmark, event_landmark) ~ BP, data = data_landmark)
summ_landmark <- summary(cox_landmark)

hr_landmark <- summ_landmark$conf.int[1, "exp(coef)"]
ci_lower_landmark <- summ_landmark$conf.int[1, "lower .95"]
ci_upper_landmark <- summ_landmark$conf.int[1, "upper .95"]
p_landmark <- summ_landmark$coefficients[1, "Pr(>|z|)"]

events_dp_landmark <- sum(data_landmark$event_landmark == 1 & data_landmark$BP == 0)
events_bp_landmark <- sum(data_landmark$event_landmark == 1 & data_landmark$BP == 1)


################################################################################
# 4. 표 생성
################################################################################

# 결과 테이블 생성
result_table <- data.frame(
  Period = c(
    "Overall (0-5 years)",
    "  Acute (≤24 hours)",
    "  Subacute (>24h-30d)",
    "  Late (30d-1yr)",
    "  Very late (>1yr)",
    "1 year (0-365 days)",
    "Landmark (365d-5yr)*"
  ),
  DP_DES = c(
    sprintf("%d (%.1f%%)", events_dp_overall, events_dp_overall/n_dp*100),
    as.character(timing_counts$DP_DES[timing_counts$ST_timing == "Acute (≤24 hours)"]),
    as.character(timing_counts$DP_DES[timing_counts$ST_timing == "Subacute (>24h-30d)"]),
    as.character(timing_counts$DP_DES[timing_counts$ST_timing == "Late (30d-1yr)"]),
    as.character(timing_counts$DP_DES[timing_counts$ST_timing == "Very late (>1yr)"]),
    sprintf("%d (%.1f%%)", events_dp_1yr, events_dp_1yr/n_dp*100),
    sprintf("%d (%.1f%%)", events_dp_landmark, events_dp_landmark/n_dp*100)
  ),
  BP_DES = c(
    sprintf("%d (%.1f%%)", events_bp_overall, events_bp_overall/n_bp*100),
    as.character(timing_counts$BP_DES[timing_counts$ST_timing == "Acute (≤24 hours)"]),
    as.character(timing_counts$BP_DES[timing_counts$ST_timing == "Subacute (>24h-30d)"]),
    as.character(timing_counts$BP_DES[timing_counts$ST_timing == "Late (30d-1yr)"]),
    as.character(timing_counts$BP_DES[timing_counts$ST_timing == "Very late (>1yr)"]),
    sprintf("%d (%.1f%%)", events_bp_1yr, events_bp_1yr/n_bp*100),
    sprintf("%d (%.1f%%)", events_bp_landmark, events_bp_landmark/n_bp*100)
  ),
  HR_95CI = c(
    sprintf("%.2f (%.2f-%.2f)", hr_overall, ci_lower_overall, ci_upper_overall),
    "", "", "", "",
    sprintf("%.2f (%.2f-%.2f)", hr_1yr, ci_lower_1yr, ci_upper_1yr),
    sprintf("%.2f (%.2f-%.2f)", hr_landmark, ci_lower_landmark, ci_upper_landmark)
  ),
  P_value = c(
    sprintf("%.3f", p_overall),
    "", "", "", "",
    sprintf("%.3f", p_1yr),
    sprintf("%.3f", p_landmark)
  ),
  stringsAsFactors = FALSE
)

# NA를 0으로 대체
result_table$DP_DES <- gsub("NA", "0", result_table$DP_DES)
result_table$BP_DES <- gsub("NA", "0", result_table$BP_DES)


################################################################################
# 5. 출력
################################################################################

cat("\nTable 3. Stent Thrombosis at 1 Year and 5 Years\n")
cat(rep("-", 80), "\n")

# 깔끔한 형식으로 출력
cat(sprintf("%-28s | %-15s | %-15s | %-18s | %-8s\n",
            "Period", "DP-DES, n (%)", "BP-DES, n (%)", "HR (95% CI)", "P-value"))
cat(rep("-", 80), "\n")

for (i in 1:nrow(result_table)) {
  cat(sprintf("%-28s | %-15s | %-15s | %-18s | %-8s\n",
              result_table$Period[i],
              result_table$DP_DES[i],
              result_table$BP_DES[i],
              result_table$HR_95CI[i],
              result_table$P_value[i]))
}

cat(rep("-", 80), "\n")
cat("* Landmark analysis includes patients event-free at 1 year\n")
cat("ST defined as definite or probable per ARC criteria\n")
cat(rep("=", 80), "\n")


################################################################################
# 6. Markdown/Publication 형식 출력
################################################################################

cat("\n\n=== PUBLICATION FORMAT ===\n\n")

print(kable(result_table,
            col.names = c("Period", "DP-DES, n (%)", "BP-DES, n (%)", 
                          "HR (95% CI)", "P-value"),
            align = c('l', 'c', 'c', 'c', 'c'),
            format = "markdown"))


################################################################################
# 7. CSV 저장
################################################################################

write.csv(result_table, "Table3_Stent_Thrombosis.csv", row.names = FALSE)
cat("\n✓ Table saved: Table3_Stent_Thrombosis.csv\n")


################################################################################
# 8. 추가 분석: Timing별 상세 분포
################################################################################

cat("\n\n=== TIMING DISTRIBUTION DETAILS ===\n\n")

timing_detail <- matched_data %>%
  filter(ST_outcome == 1) %>%
  mutate(
    Group = if_else(BP == 0, "DP-DES", "BP-DES"),
    Days_to_ST = time_to_ST
  ) %>%
  group_by(Group, ST_timing) %>%
  summarise(
    N = n(),
    Mean_days = round(mean(Days_to_ST), 0),
    Median_days = median(Days_to_ST),
    .groups = 'drop'
  ) %>%
  arrange(Group, factor(ST_timing, levels = timing_order))

print(kable(timing_detail,
            col.names = c("Group", "Timing", "N", "Mean (days)", "Median (days)"),
            align = c('l', 'l', 'c', 'c', 'c')))


################################################################################
# 9. Incidence Rate (per 1000 patient-years)
################################################################################

cat("\n\n=== INCIDENCE RATES ===\n\n")

incidence <- matched_data %>%
  group_by(Group = if_else(BP == 0, "DP-DES", "BP-DES")) %>%
  summarise(
    N = n(),
    Events = sum(ST_outcome_5yr),
    Person_Years = sum(pmin(time_to_ST, 1825)) / 365.25,
    IR_per_1000py = round(Events / Person_Years * 1000, 2),
    .groups = 'drop'
  )

print(kable(incidence,
            col.names = c("Group", "N", "Events", "Person-Years", "IR per 1000 PY"),
            align = c('l', 'c', 'c', 'c', 'c')))

cat("\n", rep("=", 70), "\n")