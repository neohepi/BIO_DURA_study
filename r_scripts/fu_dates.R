################################################################################
# Follow-up 변수 진단
################################################################################

library(dplyr)

cat("\n", rep("=", 70), "\n")
cat("FOLLOW-UP VARIABLE DIAGNOSIS\n")
cat(rep("=", 70), "\n\n")

################################################################################
# 1. fu_days 기본 통계
################################################################################

cat("=== fu_days Basic Statistics ===\n")
cat("Min:", min(matched_data$fu_days, na.rm = TRUE), "\n")
cat("Max:", max(matched_data$fu_days, na.rm = TRUE), "\n")
cat("Median:", median(matched_data$fu_days, na.rm = TRUE), "\n")
cat("Mean:", mean(matched_data$fu_days, na.rm = TRUE), "\n")
cat("Q1:", quantile(matched_data$fu_days, 0.25, na.rm = TRUE), "\n")
cat("Q3:", quantile(matched_data$fu_days, 0.75, na.rm = TRUE), "\n")

cat("\n=== fu_days Distribution ===\n")
cat("< 365:", sum(matched_data$fu_days < 365, na.rm = TRUE), "\n")
cat("365-730:", sum(matched_data$fu_days >= 365 & matched_data$fu_days < 730, na.rm = TRUE), "\n")
cat("730-1825:", sum(matched_data$fu_days >= 730 & matched_data$fu_days < 1825, na.rm = TRUE), "\n")
cat(">= 1825 (5 years):", sum(matched_data$fu_days >= 1825, na.rm = TRUE), "\n")

################################################################################
# 2. 단위 확인
################################################################################

cat("\n=== Unit Check ===\n")
# 만약 fu_days가 실제로 days면 max는 보통 3000-5000 정도
# 만약 fu_days가 years면 max는 10-15 정도
# 결과를 보면 55년이므로, 아마 다른 문제가 있음

# Sample values
cat("\nFirst 20 fu_days values:\n")
print(head(matched_data$fu_days, 20))

cat("\nLast 20 fu_days values:\n")
print(tail(matched_data$fu_days, 20))

################################################################################
# 3. 날짜 변수 확인
################################################################################

cat("\n=== Date Variables Check ===\n")

# CAG_date 확인
if ("CAG_date" %in% names(matched_data)) {
  cat("\nCAG_date:\n")
  cat("  Class:", class(matched_data$CAG_date), "\n")
  cat("  Min:", as.character(min(matched_data$CAG_date, na.rm = TRUE)), "\n")
  cat("  Max:", as.character(max(matched_data$CAG_date, na.rm = TRUE)), "\n")
  cat("  Sample:\n")
  print(head(matched_data$CAG_date, 5))
}

# last_fu_date 또는 유사 변수 확인
date_vars <- names(matched_data)[grepl("date|Date|fu|FU", names(matched_data), ignore.case = TRUE)]
cat("\nDate-related variables in dataset:\n")
print(date_vars)

################################################################################
# 4. time_to_composite와 비교
################################################################################

cat("\n=== time_to_composite Check ===\n")
if ("time_to_composite" %in% names(matched_data)) {
  cat("Min:", min(matched_data$time_to_composite, na.rm = TRUE), "\n")
  cat("Max:", max(matched_data$time_to_composite, na.rm = TRUE), "\n")
  cat("Median:", median(matched_data$time_to_composite, na.rm = TRUE), "\n")
  
  # time_to_composite와 fu_days 비교
  cat("\ntime_to_composite vs fu_days comparison:\n")
  cat("Correlation:", cor(matched_data$time_to_composite, matched_data$fu_days, use = "complete.obs"), "\n")
  
  # time_to_composite가 days 단위인지 확인
  cat("\ntime_to_composite distribution:\n")
  cat("< 365:", sum(matched_data$time_to_composite < 365, na.rm = TRUE), "\n")
  cat("365-1825:", sum(matched_data$time_to_composite >= 365 & matched_data$time_to_composite < 1825, na.rm = TRUE), "\n")
  cat(">= 1825:", sum(matched_data$time_to_composite >= 1825, na.rm = TRUE), "\n")
}

################################################################################
# 5. 실제 follow-up 계산 (날짜 변수가 있다면)
################################################################################

cat("\n=== Recalculating Follow-up ===\n")

# 마지막 추적 날짜 변수 찾기
if ("last_fu_date" %in% names(matched_data)) {
  fu_date_var <- "last_fu_date"
} else if ("fu_date" %in% names(matched_data)) {
  fu_date_var <- "fu_date"
} else if ("last_contact_date" %in% names(matched_data)) {
  fu_date_var <- "last_contact_date"
} else {
  fu_date_var <- NULL
  cat("No follow-up date variable found. Available variables:\n")
  print(names(matched_data)[1:50])
}

if (!is.null(fu_date_var) && "CAG_date" %in% names(matched_data)) {
  cat("\nRecalculating follow-up from dates...\n")
  matched_data <- matched_data %>%
    mutate(
      fu_days_recalc = as.numeric(!!sym(fu_date_var) - CAG_date)
    )
  
  cat("Recalculated fu_days:\n")
  cat("  Min:", min(matched_data$fu_days_recalc, na.rm = TRUE), "\n")
  cat("  Max:", max(matched_data$fu_days_recalc, na.rm = TRUE), "\n")
  cat("  Median:", median(matched_data$fu_days_recalc, na.rm = TRUE), "\n")
}

################################################################################
# 6. KM 그래프의 Number at Risk 수동 계산
################################################################################

cat("\n=== Manual Number at Risk Calculation ===\n")
cat("Using time_to_composite variable:\n\n")

if ("time_to_composite" %in% names(matched_data)) {
  time_points <- c(0, 365, 730, 1095, 1460, 1825)
  
  cat("Time (days) | Total at risk | BP-DES | DP-DES\n")
  cat(rep("-", 55), "\n")
  
  for (t in time_points) {
    # composite event가 발생하지 않았고, time_to_composite >= t인 환자
    at_risk_total <- sum(matched_data$time_to_composite >= t, na.rm = TRUE)
    at_risk_bp <- sum(matched_data$time_to_composite >= t & matched_data$BP == 1, na.rm = TRUE)
    at_risk_dp <- sum(matched_data$time_to_composite >= t & matched_data$BP == 0, na.rm = TRUE)
    
    cat(sprintf("%11d | %13d | %6d | %6d\n", t, at_risk_total, at_risk_bp, at_risk_dp))
  }
}

################################################################################
# 7. 이벤트 발생 시점 분포
################################################################################

cat("\n=== Event Timing Distribution ===\n")

if ("composite_event" %in% names(matched_data) && "time_to_composite" %in% names(matched_data)) {
  events_only <- matched_data %>% filter(composite_event == 1)
  
  cat("Total events:", nrow(events_only), "\n")
  cat("Events within 1 year:", sum(events_only$time_to_composite <= 365), "\n")
  cat("Events 1-5 years:", sum(events_only$time_to_composite > 365 & events_only$time_to_composite <= 1825), "\n")
  cat("Events after 5 years:", sum(events_only$time_to_composite > 1825), "\n")
  
  cat("\nEvent time distribution (days):\n")
  print(summary(events_only$time_to_composite))
}

cat("\n", rep("=", 70), "\n")