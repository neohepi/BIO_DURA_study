################################################################################
# Follow-up Rate Analysis - CORRECTED VERSION
################################################################################

library(dplyr)
library(survival)

################################################################################
# 1. Follow-up 변수 수정
################################################################################

cat("\n", rep("=", 70), "\n")
cat("FOLLOW-UP ANALYSIS (CORRECTED)\n")
cat(rep("=", 70), "\n\n")

# fu_days 재계산 (last_fu_date - CAG_date)
matched_data <- matched_data %>%
  mutate(
    fu_days_correct = as.numeric(last_fu_date - CAG_date),
    fu_years_correct = fu_days_correct / 365.25
  )

# 기본 통계 확인
cat("=== Corrected Follow-up Statistics ===\n")
cat("Min (days):", min(matched_data$fu_days_correct, na.rm = TRUE), "\n")
cat("Max (days):", max(matched_data$fu_days_correct, na.rm = TRUE), "\n")
cat("Median (days):", round(median(matched_data$fu_days_correct, na.rm = TRUE), 1), "\n")
cat("Median (years):", round(median(matched_data$fu_years_correct, na.rm = TRUE), 2), "\n")
cat("IQR (years):", 
    round(quantile(matched_data$fu_years_correct, 0.25, na.rm = TRUE), 2), "-",
    round(quantile(matched_data$fu_years_correct, 0.75, na.rm = TRUE), 2), "\n")


################################################################################
# 2. Median Follow-up (Reverse Kaplan-Meier Method)
################################################################################

cat("\n=== Median Follow-up (Reverse KM Method) ===\n")

# Reverse KM: censoring을 event로, death를 censoring으로
fu_surv_reverse <- survfit(
  Surv(fu_years_correct, 1 - any_death) ~ 1, 
  data = matched_data
)

print(fu_surv_reverse)

# 그룹별
fu_surv_by_group <- survfit(
  Surv(fu_years_correct, 1 - any_death) ~ BP, 
  data = matched_data
)
cat("\nBy Group:\n")
print(fu_surv_by_group)


################################################################################
# 3. Number at Risk at Specific Timepoints
################################################################################

cat("\n=== Number at Risk ===\n")

time_points_days <- c(365, 1825)  # 1 year, 5 years
time_labels <- c("1 year", "5 years")

n_total <- nrow(matched_data)
n_bp <- sum(matched_data$BP == 1)
n_dp <- sum(matched_data$BP == 0)

for (i in 1:length(time_points_days)) {
  t <- time_points_days[i]
  
  # fu_days_correct >= t인 환자
  at_risk_total <- sum(matched_data$fu_days_correct >= t, na.rm = TRUE)
  at_risk_bp <- sum(matched_data$fu_days_correct >= t & matched_data$BP == 1, na.rm = TRUE)
  at_risk_dp <- sum(matched_data$fu_days_correct >= t & matched_data$BP == 0, na.rm = TRUE)
  
  cat("\n", time_labels[i], ":\n")
  cat("  Total:", at_risk_total, "/", n_total, 
      sprintf("(%.1f%%)", at_risk_total/n_total*100), "\n")
  cat("  BP-DES:", at_risk_bp, "/", n_bp, 
      sprintf("(%.1f%%)", at_risk_bp/n_bp*100), "\n")
  cat("  DP-DES:", at_risk_dp, "/", n_dp, 
      sprintf("(%.1f%%)", at_risk_dp/n_dp*100), "\n")
}


################################################################################
# 4. Follow-up Complete Rate
################################################################################

cat("\n\n=== Follow-up Complete Rate ===\n")
cat("(이벤트 발생 또는 목표 기간까지 추적 완료된 비율)\n\n")

matched_data <- matched_data %>%
  mutate(
    # 1년 follow-up 완료: 1년 이상 추적 OR 1년 내 이벤트/사망
    fu_complete_1yr_correct = case_when(
      fu_days_correct >= 365 ~ 1,
      any_death == 1 ~ 1,
      composite_event == 1 ~ 1,
      TRUE ~ 0
    ),
    
    # 5년 follow-up 완료
    fu_complete_5yr_correct = case_when(
      fu_days_correct >= 1825 ~ 1,
      any_death == 1 & fu_days_correct < 1825 ~ 1,
      composite_event == 1 & time_to_composite < 1825 ~ 1,
      TRUE ~ 0
    )
  )

# 1 Year
fu_1yr_total <- mean(matched_data$fu_complete_1yr_correct, na.rm = TRUE) * 100
fu_1yr_bp <- mean(matched_data$fu_complete_1yr_correct[matched_data$BP == 1], na.rm = TRUE) * 100
fu_1yr_dp <- mean(matched_data$fu_complete_1yr_correct[matched_data$BP == 0], na.rm = TRUE) * 100

cat("At 1 Year:\n")
cat(sprintf("  Total: %.1f%%\n", fu_1yr_total))
cat(sprintf("  BP-DES: %.1f%%\n", fu_1yr_bp))
cat(sprintf("  DP-DES: %.1f%%\n", fu_1yr_dp))

# 5 Years
fu_5yr_total <- mean(matched_data$fu_complete_5yr_correct, na.rm = TRUE) * 100
fu_5yr_bp <- mean(matched_data$fu_complete_5yr_correct[matched_data$BP == 1], na.rm = TRUE) * 100
fu_5yr_dp <- mean(matched_data$fu_complete_5yr_correct[matched_data$BP == 0], na.rm = TRUE) * 100

cat("\nAt 5 Years:\n")
cat(sprintf("  Total: %.1f%%\n", fu_5yr_total))
cat(sprintf("  BP-DES: %.1f%%\n", fu_5yr_bp))
cat(sprintf("  DP-DES: %.1f%%\n", fu_5yr_dp))

# 차이 검정
fu_5yr_test <- chisq.test(table(matched_data$BP, matched_data$fu_complete_5yr_correct))
cat(sprintf("\n5-year follow-up difference: P = %.3f\n", fu_5yr_test$p.value))


################################################################################
# 5. Lost to Follow-up 분석
################################################################################

cat("\n\n=== Lost to Follow-up ===\n")

# Study end date 확인
study_end <- max(matched_data$last_fu_date, na.rm = TRUE)
cat("Study end date:", as.character(study_end), "\n")

# 추적 소실: 이벤트 없이 연구 종료 전에 follow-up 끝난 경우
# Administrative censoring threshold: study_end - 180일
admin_threshold <- as.numeric(study_end - as.Date("2024-01-01"))  # 또는 적절한 날짜

matched_data <- matched_data %>%
  mutate(
    # 연구 종료 전 추적 소실 (5년 이전)
    lost_before_5yr = if_else(
      any_death == 0 & composite_event == 0 & fu_days_correct < 1825,
      1, 0
    )
  )

ltfu_total <- sum(matched_data$lost_before_5yr, na.rm = TRUE)
ltfu_bp <- sum(matched_data$lost_before_5yr[matched_data$BP == 1], na.rm = TRUE)
ltfu_dp <- sum(matched_data$lost_before_5yr[matched_data$BP == 0], na.rm = TRUE)

cat("\nLost to follow-up before 5 years (no event):\n")
cat(sprintf("  Total: %d/%d (%.1f%%)\n", ltfu_total, n_total, ltfu_total/n_total*100))
cat(sprintf("  BP-DES: %d/%d (%.1f%%)\n", ltfu_bp, n_bp, ltfu_bp/n_bp*100))
cat(sprintf("  DP-DES: %d/%d (%.1f%%)\n", ltfu_dp, n_dp, ltfu_dp/n_dp*100))


################################################################################
# 6. Summary Table
################################################################################

cat("\n\n=== SUMMARY TABLE ===\n\n")

summary_table <- data.frame(
  Metric = c(
    "Total patients, n",
    "Median follow-up, years (IQR)",
    "",
    "Follow-up at 1 year",
    "  Number at risk, n (%)",
    "  Follow-up complete, %",
    "",
    "Follow-up at 5 years", 
    "  Number at risk, n (%)",
    "  Follow-up complete, %",
    "",
    "Lost to follow-up before 5 years, n (%)"
  ),
  Total = c(
    n_total,
    sprintf("%.1f (%.1f-%.1f)", 
            median(matched_data$fu_years_correct, na.rm = TRUE),
            quantile(matched_data$fu_years_correct, 0.25, na.rm = TRUE),
            quantile(matched_data$fu_years_correct, 0.75, na.rm = TRUE)),
    "",
    "",
    sprintf("%d (%.1f)", sum(matched_data$fu_days_correct >= 365, na.rm = TRUE),
            sum(matched_data$fu_days_correct >= 365, na.rm = TRUE)/n_total*100),
    sprintf("%.1f", fu_1yr_total),
    "",
    "",
    sprintf("%d (%.1f)", sum(matched_data$fu_days_correct >= 1825, na.rm = TRUE),
            sum(matched_data$fu_days_correct >= 1825, na.rm = TRUE)/n_total*100),
    sprintf("%.1f", fu_5yr_total),
    "",
    sprintf("%d (%.1f)", ltfu_total, ltfu_total/n_total*100)
  ),
  BP_DES = c(
    n_bp,
    sprintf("%.1f (%.1f-%.1f)",
            median(matched_data$fu_years_correct[matched_data$BP == 1], na.rm = TRUE),
            quantile(matched_data$fu_years_correct[matched_data$BP == 1], 0.25, na.rm = TRUE),
            quantile(matched_data$fu_years_correct[matched_data$BP == 1], 0.75, na.rm = TRUE)),
    "",
    "",
    sprintf("%d (%.1f)", sum(matched_data$fu_days_correct >= 365 & matched_data$BP == 1, na.rm = TRUE),
            sum(matched_data$fu_days_correct >= 365 & matched_data$BP == 1, na.rm = TRUE)/n_bp*100),
    sprintf("%.1f", fu_1yr_bp),
    "",
    "",
    sprintf("%d (%.1f)", sum(matched_data$fu_days_correct >= 1825 & matched_data$BP == 1, na.rm = TRUE),
            sum(matched_data$fu_days_correct >= 1825 & matched_data$BP == 1, na.rm = TRUE)/n_bp*100),
    sprintf("%.1f", fu_5yr_bp),
    "",
    sprintf("%d (%.1f)", ltfu_bp, ltfu_bp/n_bp*100)
  ),
  DP_DES = c(
    n_dp,
    sprintf("%.1f (%.1f-%.1f)",
            median(matched_data$fu_years_correct[matched_data$BP == 0], na.rm = TRUE),
            quantile(matched_data$fu_years_correct[matched_data$BP == 0], 0.25, na.rm = TRUE),
            quantile(matched_data$fu_years_correct[matched_data$BP == 0], 0.75, na.rm = TRUE)),
    "",
    "",
    sprintf("%d (%.1f)", sum(matched_data$fu_days_correct >= 365 & matched_data$BP == 0, na.rm = TRUE),
            sum(matched_data$fu_days_correct >= 365 & matched_data$BP == 0, na.rm = TRUE)/n_bp*100),
    sprintf("%.1f", fu_1yr_dp),
    "",
    "",
    sprintf("%d (%.1f)", sum(matched_data$fu_days_correct >= 1825 & matched_data$BP == 0, na.rm = TRUE),
            sum(matched_data$fu_days_correct >= 1825 & matched_data$BP == 0, na.rm = TRUE)/n_dp*100),
    sprintf("%.1f", fu_5yr_dp),
    "",
    sprintf("%d (%.1f)", ltfu_dp, ltfu_dp/n_dp*100)
  )
)

print(summary_table)

write.csv(summary_table, "Follow_up_Summary_Corrected.csv", row.names = FALSE)


################################################################################
# 7. 논문 문구 제안
################################################################################

cat("\n\n=== SUGGESTED MANUSCRIPT TEXT ===\n\n")

median_fu <- round(median(matched_data$fu_years_correct, na.rm = TRUE), 1)
iqr_low <- round(quantile(matched_data$fu_years_correct, 0.25, na.rm = TRUE), 1)
iqr_high <- round(quantile(matched_data$fu_years_correct, 0.75, na.rm = TRUE), 1)

at_risk_1yr_pct <- round(sum(matched_data$fu_days_correct >= 365, na.rm = TRUE)/n_total*100, 1)
at_risk_5yr_pct <- round(sum(matched_data$fu_days_correct >= 1825, na.rm = TRUE)/n_total*100, 1)

cat(sprintf(
  "Median follow-up duration was %.1f years (interquartile range %.1f–%.1f years) 
in the matched cohort. At 1 year, %.1f%% of patients remained at risk, 
and at 5 years, %.1f%% remained at risk.",
  median_fu, iqr_low, iqr_high, at_risk_1yr_pct, at_risk_5yr_pct
))

cat("\n\n")
cat(rep("=", 70), "\n")


# 실제 추적 기간 < 5년인데, time_to_composite = 1825인 환자 수
mismatch <- matched_data %>%
  filter(fu_days_correct < 1825 & time_to_composite >= 1825 & composite_event == 0)

cat("Mismatch cases:", nrow(mismatch), "\n")
cat("Percentage:", round(nrow(mismatch)/nrow(matched_data)*100, 1), "%\n")
