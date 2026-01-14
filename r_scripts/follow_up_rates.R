################################################################################
# Follow-up Rate Analysis
################################################################################

library(dplyr)
library(survival)

################################################################################
# 1. 기본 Follow-up 통계
################################################################################

cat("\n", rep("=", 70), "\n")
cat("FOLLOW-UP ANALYSIS\n")
cat(rep("=", 70), "\n\n")

# Median follow-up (Reverse Kaplan-Meier method)
# 사망을 censoring으로, censoring을 event로 처리
matched_data <- matched_data %>%
  mutate(
    # Follow-up time in years
    fu_years = fu_days / 365.25,
    
    # Reverse coding for median follow-up calculation
    # 실제 follow-up 완료 = censored (event 없이 마지막까지 추적)
    fu_complete = if_else(any_death == 0 & composite_event == 0, 1, 0)
  )

# Median follow-up (reverse KM)
fu_surv <- survfit(Surv(fu_years, 1 - any_death) ~ 1, data = matched_data)
cat("=== Median Follow-up Duration (Reverse KM Method) ===\n")
print(fu_surv)

# Simple median
cat("\n=== Simple Follow-up Statistics ===\n")
cat("Median follow-up (years):", round(median(matched_data$fu_years), 2), "\n")
cat("IQR:", round(quantile(matched_data$fu_years, 0.25), 2), "-", 
    round(quantile(matched_data$fu_years, 0.75), 2), "\n")
cat("Mean follow-up (years):", round(mean(matched_data$fu_years), 2), "\n")
cat("Range:", round(min(matched_data$fu_years), 2), "-", 
    round(max(matched_data$fu_years), 2), "\n")


################################################################################
# 2. Number at Risk vs Follow-up Complete Rate 비교
################################################################################

cat("\n\n=== Number at Risk Analysis ===\n")
cat("(환자가 해당 시점까지 '살아있고 추적 가능한' 수)\n\n")

# Number at risk at specific time points
time_points <- c(365, 365*5)  # 1 year, 5 years in days
time_labels <- c("1 year", "5 years")

for (i in 1:length(time_points)) {
  t <- time_points[i]
  
  # 해당 시점에서 아직 추적 중인 환자 수
  # fu_days >= t 이면 해당 시점에 at risk
  at_risk_total <- sum(matched_data$fu_days >= t)
  at_risk_bp <- sum(matched_data$fu_days >= t & matched_data$BP == 1)
  at_risk_dp <- sum(matched_data$fu_days >= t & matched_data$BP == 0)
  
  n_total <- nrow(matched_data)
  n_bp <- sum(matched_data$BP == 1)
  n_dp <- sum(matched_data$BP == 0)
  
  cat(time_labels[i], ":\n")
  cat("  Total at risk:", at_risk_total, "/", n_total, 
      sprintf("(%.1f%%)", at_risk_total/n_total*100), "\n")
  cat("  BP-DES at risk:", at_risk_bp, "/", n_bp, 
      sprintf("(%.1f%%)", at_risk_bp/n_bp*100), "\n")
  cat("  DP-DES at risk:", at_risk_dp, "/", n_dp, 
      sprintf("(%.1f%%)", at_risk_dp/n_dp*100), "\n\n")
}


################################################################################
# 3. Follow-up Complete Rate (Administrative Censoring 제외)
################################################################################

cat("\n=== Follow-up Complete Rate ===\n")
cat("(이벤트 발생 또는 연구 종료까지 추적 완료된 비율)\n")
cat("(중도 탈락/추적 소실 제외)\n\n")

# 연구 종료일 확인 (가장 긴 follow-up 기준)
max_fu <- max(matched_data$fu_days)
cat("Maximum follow-up in dataset:", round(max_fu/365.25, 2), "years\n\n")

# Follow-up complete 정의:
# 1) 이벤트 발생 (사망 또는 composite event)
# 2) 또는 administrative censoring (연구 종료까지 추적)

# Administrative censoring threshold (예: max_fu - 30일 이내면 완료로 간주)
admin_censor_threshold <- max_fu - 30

matched_data <- matched_data %>%
  mutate(
    # 1년 시점 follow-up 완료
    fu_complete_1yr = case_when(
      fu_days >= 365 ~ 1,  # 1년 이상 추적됨
      any_death == 1 & fu_days < 365 ~ 1,  # 1년 내 사망 (이벤트로 추적 완료)
      composite_event == 1 & time_to_composite < 365 ~ 1,  # 1년 내 이벤트
      TRUE ~ 0  # 1년 전 추적 소실
    ),
    
    # 5년 시점 follow-up 완료
    fu_complete_5yr = case_when(
      fu_days >= 365*5 ~ 1,  # 5년 이상 추적됨
      any_death == 1 & fu_days < 365*5 ~ 1,  # 5년 내 사망
      composite_event == 1 & time_to_composite < 365*5 ~ 1,  # 5년 내 이벤트
      TRUE ~ 0  # 5년 전 추적 소실
    )
  )

# 1 Year
cat("At 1 Year:\n")
fu_1yr_total <- mean(matched_data$fu_complete_1yr) * 100
fu_1yr_bp <- mean(matched_data$fu_complete_1yr[matched_data$BP == 1]) * 100
fu_1yr_dp <- mean(matched_data$fu_complete_1yr[matched_data$BP == 0]) * 100

cat(sprintf("  Total: %.1f%%\n", fu_1yr_total))
cat(sprintf("  BP-DES: %.1f%%\n", fu_1yr_bp))
cat(sprintf("  DP-DES: %.1f%%\n", fu_1yr_dp))

# 5 Years
cat("\nAt 5 Years:\n")
fu_5yr_total <- mean(matched_data$fu_complete_5yr) * 100
fu_5yr_bp <- mean(matched_data$fu_complete_5yr[matched_data$BP == 1]) * 100
fu_5yr_dp <- mean(matched_data$fu_complete_5yr[matched_data$BP == 0]) * 100

cat(sprintf("  Total: %.1f%%\n", fu_5yr_total))
cat(sprintf("  BP-DES: %.1f%%\n", fu_5yr_bp))
cat(sprintf("  DP-DES: %.1f%%\n", fu_5yr_dp))


################################################################################
# 4. Lost to Follow-up 분석
################################################################################

cat("\n\n=== Lost to Follow-up Analysis ===\n")

# Lost to follow-up = 이벤트 없이 중도 탈락
matched_data <- matched_data %>%
  mutate(
    # 추적 소실 여부 (이벤트 없이 max_fu보다 일찍 종료)
    lost_to_fu = if_else(
      any_death == 0 & composite_event == 0 & fu_days < admin_censor_threshold,
      1, 0
    )
  )

ltfu_total <- sum(matched_data$lost_to_fu)
ltfu_bp <- sum(matched_data$lost_to_fu[matched_data$BP == 1])
ltfu_dp <- sum(matched_data$lost_to_fu[matched_data$BP == 0])

cat("Lost to follow-up (before study end):\n")
cat(sprintf("  Total: %d/%d (%.1f%%)\n", 
            ltfu_total, nrow(matched_data), ltfu_total/nrow(matched_data)*100))
cat(sprintf("  BP-DES: %d/%d (%.1f%%)\n", 
            ltfu_bp, sum(matched_data$BP == 1), ltfu_bp/sum(matched_data$BP == 1)*100))
cat(sprintf("  DP-DES: %d/%d (%.1f%%)\n", 
            ltfu_dp, sum(matched_data$BP == 0), ltfu_dp/sum(matched_data$BP == 0)*100))

# Chi-square test for difference
ltfu_table <- table(matched_data$BP, matched_data$lost_to_fu)
ltfu_test <- chisq.test(ltfu_table)
cat(sprintf("\nDifference between groups: P = %.3f\n", ltfu_test$p.value))


################################################################################
# 5. KM Curve의 Number at Risk 재현
################################################################################

cat("\n\n=== Reproducing KM Number at Risk ===\n")
cat("(그래프 하단의 숫자와 비교)\n\n")

# Survfit 객체에서 number at risk 추출
# composite event 기준
surv_fit <- survfit(Surv(fu_days, composite_event) ~ BP, data = matched_data)

# 특정 시점의 number at risk
summary_times <- c(0, 365, 730, 1095, 1460, 1825)  # 0, 1, 2, 3, 4, 5 years in days
surv_summary <- summary(surv_fit, times = summary_times)

cat("Time (days) | DP-DES at risk | BP-DES at risk\n")
cat(rep("-", 50), "\n")

# strata 순서 확인
strata_levels <- levels(factor(matched_data$BP))

for (i in 1:length(summary_times)) {
  idx_dp <- i  # BP=0
  idx_bp <- i + length(summary_times)  # BP=1
  
  cat(sprintf("%11d | %14d | %14d\n", 
              summary_times[i],
              surv_summary$n.risk[idx_dp],
              surv_summary$n.risk[idx_bp]))
}


################################################################################
# 6. 상세 Follow-up Table
################################################################################

cat("\n\n=== Detailed Follow-up Summary Table ===\n\n")

fu_summary <- data.frame(
  Metric = c(
    "Total patients",
    "Median follow-up (years)",
    "IQR (years)",
    "",
    "Number at risk at 1 year",
    "Number at risk at 5 years",
    "",
    "Follow-up complete at 1 year (%)",
    "Follow-up complete at 5 years (%)",
    "",
    "Lost to follow-up (%)"
  ),
  Total = c(
    nrow(matched_data),
    round(median(matched_data$fu_years), 2),
    paste0(round(quantile(matched_data$fu_years, 0.25), 2), "-", 
           round(quantile(matched_data$fu_years, 0.75), 2)),
    "",
    sum(matched_data$fu_days >= 365),
    sum(matched_data$fu_days >= 365*5),
    "",
    sprintf("%.1f", fu_1yr_total),
    sprintf("%.1f", fu_5yr_total),
    "",
    sprintf("%.1f", ltfu_total/nrow(matched_data)*100)
  ),
  BP_DES = c(
    sum(matched_data$BP == 1),
    round(median(matched_data$fu_years[matched_data$BP == 1]), 2),
    paste0(round(quantile(matched_data$fu_years[matched_data$BP == 1], 0.25), 2), "-",
           round(quantile(matched_data$fu_years[matched_data$BP == 1], 0.75), 2)),
    "",
    sum(matched_data$fu_days >= 365 & matched_data$BP == 1),
    sum(matched_data$fu_days >= 365*5 & matched_data$BP == 1),
    "",
    sprintf("%.1f", fu_1yr_bp),
    sprintf("%.1f", fu_5yr_bp),
    "",
    sprintf("%.1f", ltfu_bp/sum(matched_data$BP == 1)*100)
  ),
  DP_DES = c(
    sum(matched_data$BP == 0),
    round(median(matched_data$fu_years[matched_data$BP == 0]), 2),
    paste0(round(quantile(matched_data$fu_years[matched_data$BP == 0], 0.25), 2), "-",
           round(quantile(matched_data$fu_years[matched_data$BP == 0], 0.75), 2)),
    "",
    sum(matched_data$fu_days >= 365 & matched_data$BP == 0),
    sum(matched_data$fu_days >= 365*5 & matched_data$BP == 0),
    "",
    sprintf("%.1f", fu_1yr_dp),
    sprintf("%.1f", fu_5yr_dp),
    "",
    sprintf("%.1f", ltfu_dp/sum(matched_data$BP == 0)*100)
  )
)

print(fu_summary)

# Export
write.csv(fu_summary, "Follow_up_Summary.csv", row.names = FALSE)

cat("\n", rep("=", 70), "\n")
cat("NOTE:\n")
cat("- 'Number at risk' = 해당 시점에 아직 이벤트 없이 추적 중인 환자 수\n")
cat("- 'Follow-up complete' = 이벤트 발생 또는 연구 종료까지 추적된 비율\n")
cat("- KM 그래프의 Number at risk는 '이벤트 없이 추적 중'인 수를 보여줌\n")
cat("- Follow-up complete rate는 '추적 소실' 없이 완료된 비율\n")
cat(rep("=", 70), "\n")

