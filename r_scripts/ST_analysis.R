library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)
library(knitr)

# ============================================================================
# Complete ST Analysis Function
# ============================================================================

analyze_ST_comprehensive <- function(data,
                                     time_var = "time_to_ST",
                                     event_var = "ST_outcome",
                                     group_var = "BP",
                                     timing_var = "ST_timing",
                                     landmark_day = 365) {
  
  cat("\n╔═══════════════════════════════════════════════════════════╗\n")
  cat("║        COMPREHENSIVE STENT THROMBOSIS ANALYSIS           ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n\n")
  
  # Prepare data
  analysis_data <- data %>%
    mutate(
      time = .data[[time_var]],
      event = .data[[event_var]],
      group = .data[[group_var]],
      timing = .data[[timing_var]]
    )
  
  # ========================================================================
  # 1. OVERALL SUMMARY
  # ========================================================================
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("1. OVERALL SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  overall_summary <- analysis_data %>%
    group_by(group) %>%
    summarise(
      N = n(),
      ST_events = sum(event),
      ST_rate = round(mean(event) * 100, 2),
      Median_FU_days = median(time),
      Mean_FU_days = round(mean(time), 1),
      .groups = 'drop'
    ) %>%
    mutate(
      Group = ifelse(group == 0, "DP-DES", "BP-DES")
    ) %>%
    select(Group, N, ST_events, ST_rate, Median_FU_days, Mean_FU_days)
  
  print(kable(overall_summary, 
              col.names = c("Group", "N", "ST Events", "Rate (%)", "Median FU (d)", "Mean FU (d)"),
              align = 'c'))
  cat("\n")
  
  # ========================================================================
  # 2. ST TIMING DISTRIBUTION
  # ========================================================================
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("2. STENT THROMBOSIS TIMING DISTRIBUTION\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  timing_summary <- analysis_data %>%
    filter(event == 1) %>%
    group_by(Group = ifelse(group == 0, "DP-DES", "BP-DES"), timing) %>%
    summarise(N = n(), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = timing, values_from = N, values_fill = 0) %>%
    mutate(Total = rowSums(select(., -Group)))
  
  print(kable(timing_summary, align = 'c'))
  cat("\n")
  
  # ========================================================================
  # 3. COX REGRESSION - OVERALL AND LANDMARK
  # ========================================================================
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("3. COX PROPORTIONAL HAZARDS ANALYSIS\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # Overall (0-5y)
  fit_overall <- coxph(Surv(time, event) ~ group, data = analysis_data)
  hr_overall <- exp(coef(fit_overall))
  ci_overall <- exp(confint(fit_overall))
  p_overall <- summary(fit_overall)$coefficients[,"Pr(>|z|)"]
  
  # Early landmark (0 to landmark_day)
  data_early <- analysis_data %>%
    mutate(
      event_early = ifelse(time <= landmark_day & event == 1, 1, 0),
      time_early = pmin(time, landmark_day)
    )
  
  fit_early <- coxph(Surv(time_early, event_early) ~ group, data = data_early)
  hr_early <- exp(coef(fit_early))
  ci_early <- exp(confint(fit_early))
  p_early <- summary(fit_early)$coefficients[,"Pr(>|z|)"]
  
  # Late landmark (after landmark_day)
  data_late <- analysis_data %>%
    filter(time > landmark_day | (time <= landmark_day & event == 0)) %>%
    mutate(
      time_late = time - landmark_day,
      event_late = ifelse(time > landmark_day & event == 1, 1, 0)
    ) %>%
    filter(time_late >= 0)
  
  fit_late <- coxph(Surv(time_late, event_late) ~ group, data = data_late)
  hr_late <- exp(coef(fit_late))
  ci_late <- exp(confint(fit_late))
  p_late <- summary(fit_late)$coefficients[,"Pr(>|z|)"]
  
  # Event counts
  events_overall_dp <- sum(analysis_data$event[analysis_data$group == 0])
  events_overall_bp <- sum(analysis_data$event[analysis_data$group == 1])
  events_early_dp <- sum(data_early$event_early[data_early$group == 0])
  events_early_bp <- sum(data_early$event_early[data_early$group == 1])
  events_late_dp <- sum(data_late$event_late[data_late$group == 0])
  events_late_bp <- sum(data_late$event_late[data_late$group == 1])
  
  # Cox table
  cox_table <- data.frame(
    Period = c("Overall (0-5y)", 
               sprintf("Landmark (0-%dd)", landmark_day),
               sprintf("Landmark (%dd-5y)", landmark_day)),
    DP_DES = sprintf("%d (%.1f%%)", 
                     c(events_overall_dp, events_early_dp, events_late_dp),
                     c(events_overall_dp, events_early_dp, events_late_dp) / 
                       c(sum(analysis_data$group == 0), 
                         sum(data_early$group == 0),
                         sum(data_late$group == 0)) * 100),
    BP_DES = sprintf("%d (%.1f%%)", 
                     c(events_overall_bp, events_early_bp, events_late_bp),
                     c(events_overall_bp, events_early_bp, events_late_bp) / 
                       c(sum(analysis_data$group == 1),
                         sum(data_early$group == 1),
                         sum(data_late$group == 1)) * 100),
    HR = sprintf("%.2f", c(hr_overall, hr_early, hr_late)),
    CI_95 = sprintf("%.2f-%.2f", 
                    c(ci_overall[1], ci_early[1], ci_late[1]),
                    c(ci_overall[2], ci_early[2], ci_late[2])),
    P_value = sprintf("%.3f", c(p_overall, p_early, p_late))
  )
  
  print(kable(cox_table, 
              col.names = c("Period", "DP-DES", "BP-DES", "HR", "95% CI", "P-value"),
              align = 'lccccr'))
  cat("\n")
  
  # ========================================================================
  # 4. INCIDENCE RATES
  # ========================================================================
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("4. INCIDENCE RATES (per 100 patient-years)\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  incidence_data <- analysis_data %>%
    group_by(Group = ifelse(group == 0, "DP-DES", "BP-DES")) %>%
    summarise(
      Events = sum(event),
      Person_Years = sum(time) / 365.25,
      IR_per_100py = round(Events / (Person_Years / 100), 2),
      .groups = 'drop'
    )
  
  print(kable(incidence_data,
              col.names = c("Group", "Events", "Person-Years", "IR per 100 PY"),
              align = 'lccc'))
  cat("\n")
  
  # ========================================================================
  # 5. KAPLAN-MEIER CURVES
  # ========================================================================
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("5. GENERATING KAPLAN-MEIER CURVES\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # KM fit
  fit_km <- survfit(Surv(time, event) ~ group, data = analysis_data)
  
  # Get cumulative incidence at specific timepoints
  km_summary <- summary(fit_km, times = c(30, 365, 1825))
  
  timepoints_table <- data.frame(
    Timepoint = c("30 days", "1 year", "5 years"),
    DP_DES = sprintf("%.2f%%", (1 - km_summary$surv[km_summary$strata == "group=0"]) * 100),
    BP_DES = sprintf("%.2f%%", (1 - km_summary$surv[km_summary$strata == "group=1"]) * 100)
  )
  
  print(kable(timepoints_table,
              col.names = c("Timepoint", "DP-DES", "BP-DES"),
              align = 'lcc'))
  cat("\n")
  
  # Main KM plot
  p_km <- ggsurvplot(
    fit_km,
    data = analysis_data,
    fun = "event",
    palette = c("#E64B35", "#4DBBD5"),
    legend.title = "Stent Type",
    legend.labs = c("DP-DES", "BP-DES"),
    xlab = "Time (days)",
    ylab = "Cumulative Incidence (%)",
    xlim = c(0, 1825),
    break.time.by = 365,
    pval = TRUE,
    pval.coord = c(200, 0.008),
    pval.size = 4,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.fontsize = 3.5,
    tables.theme = theme_cleantable(),
    ggtheme = theme_classic(base_size = 12),
    font.main = c(14, "bold"),
    font.x = c(12, "plain"),
    font.y = c(12, "plain")
  )
  
  # Add landmark line
  p_km$plot <- p_km$plot + 
    geom_vline(xintercept = landmark_day, linetype = "dashed", 
               color = "black", linewidth = 0.5) +
    annotate("text", x = landmark_day + 100, y = 0.009,
             label = sprintf("Landmark\n(%d days)", landmark_day),
             hjust = 0, size = 3.5, fontface = "bold") +
    scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                       limits = c(0, 0.01)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = c(0.15, 0.85),
      legend.background = element_rect(fill = "white", color = "black")
    ) +
    labs(title = "Stent Thrombosis: Cumulative Incidence")
  
  # ========================================================================
  # 6. TIMING-SPECIFIC ANALYSIS
  # ========================================================================
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("6. TIMING-SPECIFIC STENT THROMBOSIS RATES\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  timing_detailed <- analysis_data %>%
    filter(event == 1) %>%
    group_by(Group = ifelse(group == 0, "DP-DES", "BP-DES"), timing) %>%
    summarise(N = n(), .groups = 'drop') %>%
    group_by(Group) %>%
    mutate(
      Total = sum(N),
      Percentage = round(N / Total * 100, 1)
    ) %>%
    arrange(Group, timing)
  
  print(kable(timing_detailed,
              col.names = c("Group", "Timing", "N", "Total", "%"),
              align = 'llccc'))
  cat("\n")
  
  # Timing bar plot
  p_timing <- ggplot(timing_detailed, aes(x = Group, y = N, fill = timing)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = N), position = position_dodge(width = 0.7),
              vjust = -0.5, size = 3.5) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Stent Thrombosis by Timing",
      x = "",
      y = "Number of Events",
      fill = "ST Timing"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  # ========================================================================
  # 7. STATISTICAL INTERPRETATION
  # ========================================================================
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("7. STATISTICAL INTERPRETATION\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # Power calculation (post-hoc)
  total_events <- sum(analysis_data$event)
  expected_hr <- hr_overall
  
  cat(sprintf("Total ST events: %d\n", total_events))
  cat(sprintf("Overall HR: %.2f (95%% CI %.2f-%.2f)\n", 
              hr_overall, ci_overall[1], ci_overall[2]))
  cat(sprintf("P-value: %.3f\n\n", p_overall))
  
  if (p_overall > 0.05) {
    cat("INTERPRETATION:\n")
    cat("• No statistically significant difference detected (p > 0.05)\n")
    cat(sprintf("• Study underpowered for ST (n=%d events)\n", total_events))
    cat("• ST rate very low (<1%) in both groups - reassuring safety signal\n")
    cat("• Wide confidence interval indicates uncertainty\n")
  } else {
    cat("INTERPRETATION:\n")
    cat("• Statistically significant difference detected (p < 0.05)\n")
    cat(sprintf("• BP-DES associated with %.0f%% relative risk reduction\n", 
                (1 - hr_overall) * 100))
  }
  
  cat("\n")
  
  # ========================================================================
  # 8. RETURN RESULTS
  # ========================================================================
  
  results <- list(
    overall_summary = overall_summary,
    timing_summary = timing_summary,
    cox_table = cox_table,
    incidence_rates = incidence_data,
    timepoints = timepoints_table,
    timing_detailed = timing_detailed,
    km_plot = p_km,
    timing_plot = p_timing,
    fit_overall = fit_overall,
    fit_early = fit_early,
    fit_late = fit_late,
    fit_km = fit_km
  )
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("✓ Analysis completed successfully!\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  return(invisible(results))
}


# ============================================================================
# Publication-Ready Table Generator
# ============================================================================

create_ST_publication_table <- function(results) {
  
  cat("\n╔═══════════════════════════════════════════════════════════╗\n")
  cat("║          PUBLICATION-READY TABLE                         ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n\n")
  
  cat("Table X. Stent Thrombosis at 5 Years\n")
  cat("─────────────────────────────────────────────────────────────\n\n")
  
  print(kable(results$cox_table,
              format = "markdown",
              col.names = c("Period", "DP-DES n (%)", "BP-DES n (%)", 
                            "Hazard Ratio", "95% CI", "P Value")))
  
  cat("\n\nDefinitions:\n")
  cat("ST = Stent thrombosis (definite or probable by ARC criteria)\n")
  cat("Values are n (%) unless otherwise specified\n")
  cat("HR = Hazard ratio from Cox proportional hazards model\n")
  cat("\n")
}


# ============================================================================
# Save All Outputs
# ============================================================================

save_ST_analysis <- function(results, output_prefix = "ST_analysis") {
  
  # Save KM plot
  ggsave(
    paste0(output_prefix, "_KM_curve.png"),
    plot = results$km_plot$plot,
    width = 10, height = 7, dpi = 300
  )
  
  ggsave(
    paste0(output_prefix, "_KM_curve.tiff"),
    plot = results$km_plot$plot,
    width = 10, height = 7, dpi = 300,
    compression = "lzw"
  )
  
  # Save timing plot
  ggsave(
    paste0(output_prefix, "_timing.png"),
    plot = results$timing_plot,
    width = 8, height = 6, dpi = 300
  )
  
  # Save tables to CSV
  write.csv(results$cox_table, 
            paste0(output_prefix, "_cox_results.csv"),
            row.names = FALSE)
  
  write.csv(results$overall_summary,
            paste0(output_prefix, "_summary.csv"),
            row.names = FALSE)
  
  cat("\n✓ All outputs saved with prefix:", output_prefix, "\n")
  cat("  - KM curves (PNG & TIFF)\n")
  cat("  - Timing plot (PNG)\n")
  cat("  - Results tables (CSV)\n\n")
}

# ============================================================================
# COMPLETE WORKFLOW
# ============================================================================

# Step 1: Create ST variables (이미 완료)
# matched_data <- create_ST_variables(
#   data = matched_data,
#   index_date_col = "CAG_date",
#   st_date_col = "ST_date",
#   last_fu_col = "last_fu_date",
#   st_indicator_col = "ST",
#   max_followup_days = 1825
# )

# Step 2: Run comprehensive analysis
results_ST <- analyze_ST_comprehensive(
  data = matched_data,
  time_var = "time_to_ST",
  event_var = "ST_outcome",
  group_var = "BP",
  timing_var = "ST_timing",
  landmark_day = 365
)

# Step 3: View KM plot
print(results_ST$km_plot)

# Step 4: View timing plot
print(results_ST$timing_plot)

# Step 5: Create publication table
create_ST_publication_table(results_ST)

# Step 6: Save all outputs
save_ST_analysis(results_ST, output_prefix = "ESC_ST_analysis")

# Step 7: Access individual results
results_ST$cox_table
results_ST$overall_summary
results_ST$incidence_rates
results_ST$timepoints

################################################################################
# Very Late ST (>1 year) detailed analysis
very_late_st <- matched_data %>%
  filter(ST_outcome == 1, time_to_ST > 365) %>%
  group_by(BP) %>%
  summarise(
    N = n(),
    Group = ifelse(BP[1] == 0, "DP-DES", "BP-DES"),
    Mean_time = round(mean(time_to_ST), 0),
    Median_time = median(time_to_ST),
    Min_time = min(time_to_ST),
    Max_time = max(time_to_ST)
  ) %>%
  select(Group, N, Mean_time, Median_time, Min_time, Max_time)

print(kable(very_late_st,
            col.names = c("Group", "N", "Mean (d)", "Median (d)", "Min (d)", "Max (d)"),
            caption = "Very Late ST (>365 days) Characteristics"))

################################################################################
# Median ST-free survival
fit_km <- survfit(Surv(time_to_ST, ST_outcome) ~ BP, data = matched_data)

cat("\nMedian ST-free survival:\n")
print(fit_km)

################################################################################
# ST by ACS vs stable CAD
subgroup_st <- matched_data %>%
  group_by(is_ACS, Group = ifelse(BP == 0, "DP-DES", "BP-DES")) %>%
  summarise(
    N = n(),
    ST_events = sum(ST_outcome),
    ST_rate = round(mean(ST_outcome) * 100, 2),
    .groups = 'drop'
  ) %>%
  mutate(Presentation = ifelse(is_ACS == 1, "ACS", "Stable CAD"))

print(kable(subgroup_st))





################################################################################
library(survival)
library(survminer)
library(scales)

result_ST <- analyze_landmark_survival(
  data = matched_data,
  group_var = "BP",
  date_cag_var = "CAG_date",
  date_event_var = "ST_date",   # <--- 여기만 변경
  event_status_var = "ST_outcome",      # <--- 여기만 변경
  fu_days_var = "fu_days",
  cap_years = 5,
  landmark_days = 0,
  yrange = c(0, 0.01),
  #  plot_title = "Primary Composite"
)

print(result_ST$cox_table)
print(result_ST$plot)


# 1. 패키지 로드 (없으면 install.packages("survminer"))
library(survminer)
library(survival)

# 2. Survival Object 생성 (기존에 fit 하신 게 있다면 그걸 쓰시면 됩니다)
# 예시: fit <- survfit(Surv(time, status) ~ group, data = data)

# 3. 그래프 그리기 (ggsurvplot 사용)
p <- ggsurvplot(
  result_ST$fit,                    # 교수님의 survfit 객체 이름
  data = matched_data,          # 사용한 데이터프레임
  
  # 디자인 설정 (jskm 스타일과 비슷하게)
  risk.table = TRUE,      # 리스크 테이블 표시
  pval = TRUE,            # P-value 표시
  conf.int = FALSE,       # 신뢰구간 표시 여부 (선택)
  xlim = c(0, 1820),      # 5년(1825일)까지 표시 (필요에 따라 수정)
  ylim = c(0, 0.02),
  break.time.by = 365,    # 1년 단위로 눈금 끊기
  
  # 색상 및 스타일
  palette = c("red", "blue"),    # DP(Red), BP(Blue) - 순서 확인 필요
  legend.labs = c("DP-DES", "BP-DES"), 
  legend.title = "Strata",
  
  # 축 라벨
  xlab = "Time-to-event (days)",
  ylab = "Cumulative incidence (%)",
  fun = "event"           # 누적 발생률(Cumulative Incidence)로 변환
)

# 4. 여기에 365일 점선 추가 (핵심!)
# ggsurvplot은 결과가 리스트($plot, $table)로 나오므로 $plot에 추가해야 합니다.
p$plot <- p$plot + 
  geom_vline(xintercept = 365, linetype = "dashed", color = "black", linewidth = 0.5) +
  annotate("text", x = 365, y = 0.02, label = "1 Year", hjust = -0.2)

# 5. 최종 출력
print(p)

library(dplyr)

# 1. ST가 발생한 환자만 추출하여 내림차순 정렬
check_tail_events <- matched_data %>%
  # (1) ST 이벤트가 있는 사람만 필터링
  filter(ST_outcome == 1) %>%
  
  # (2) 시간 계산 (이미 되어 있다면 생략 가능)
  # Index Date(CAG_date)와 ST 발생일(ST_date) 차이 계산
  mutate(
    days_to_ST = as.numeric(as.Date(ST_date) - as.Date(CAG_date))
  ) %>%
  
  # (3) 내림차순 정렬 (가장 늦게 생긴 이벤트부터)
  arrange(desc(days_to_ST)) %>%
  
  # (4) 확인하기 편하게 주요 컬럼만 선택
  select(pt_id, BP, CAG_date, ST_date, days_to_ST)

# 2. 상위 20개 출력 (그래프 끝부분 확인)
print(head(check_tail_events, 20))

# -----------------------------------------------------------
# [추가 분석] 5년(1825일) 근처에 이벤트가 몇 건이나 몰려있는지 확인
# -----------------------------------------------------------
cat("\n=== Distribution of Events near 5 Years (1800~1830 days) ===\n")

check_tail_events %>%
  filter(days_to_ST >= 1800) %>%  # 1800일 이후 발생한 건만
  count(days_to_ST) %>%           # 날짜별 발생 건수 카운트
  arrange(desc(days_to_ST)) %>%
  print()


################################################################################
# time_to_death
library(dplyr)
library(lubridate)

# 예: index_date가 Date 타입이어야 함
matched_data <- matched_data %>%
  mutate(
    index_date   = as.Date(index_date),
    death_date   = as.Date(death_date),
    time_to_death = ifelse(any_death == 1,
                           as.numeric(death_date - index_date),
                           NA_real_)
  )

admin_censor <- 1825  # 5 years

dat_cr <- matched_data %>%
  mutate(
    # NA 처리용
    t_st    = ifelse(is.na(time_to_ST),    Inf, time_to_ST),
    t_death = ifelse(is.na(time_to_death), Inf, time_to_death),
    
    # 실제 분석 시간: ST/사망/행정censor 중 가장 빠른 것
    ftime = pmin(t_st, t_death, admin_censor),
    
    # 같은 날 ST와 death가 같이 있으면(매우 드뭄) 우선순위 규칙 필요
    # 여기서는 ST를 우선(원하면 death 우선으로 바꿔도 됨)
    fstatus = case_when(
      ST_outcome == 1 & t_st <= t_death & t_st <= admin_censor ~ 1L,
      any_death  == 1 & t_death <  t_st & t_death <= admin_censor ~ 2L,
      TRUE ~ 0L
    )
  )

library(cmprsk)

# BP가 0/1이면 factor로
ci <- with(dat_cr, cuminc(ftime, fstatus, group = BP))

print(ci)     # Gray test p-value 포함
plot(ci, xlab="Days", ylab="Cumulative incidence", lty=1)

library(survival)

dat_cr <- dat_cr %>% mutate(pair_id = subclass)

# 예: pair_id가 없고 MatchIt subclass가 있다면
# dat_cr <- dat_cr %>% mutate(pair_id = subclass)

dat_cr <- dat_cr %>%
  mutate(
    fstatus = case_when(
      fstatus == 0 ~ "censor",
      fstatus == 1 ~ "ST",
      fstatus == 2 ~ "death",
      TRUE ~ NA_character_
    ),
    fstatus = factor(fstatus, levels = c("censor", "ST", "death")),
    BP = factor(BP, levels = c(0,1), labels = c("DP-DES","BP-DES"))
  )

table(dat_cr$fstatus, useNA = "ifany")   # 반드시 NA 없어야 함
summary(dat_cr$ftime)                   # ftime도 NA 없어야 함


fgdat <- finegray(Surv(ftime, fstatus) ~ BP + pair_id,
                  data = dat_cr,
                  etype = "ST")


fg_fit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ BP + cluster(pair_id),
                weights = fgwt, data = fgdat)

summary(fg_fit)
exp(coef(fg_fit))                 # sHR
exp(confint(fg_fit))              # 95% CI


library(survival)
library(survminer)
library(ggplot2)
library(scales)

fit <- survfit(Surv(time_to_ST, ST_outcome) ~ BP, data = matched_data)

p <- ggsurvplot(
  fit, data = matched_data,
  fun = "event",              # 1 - S(t) : 누적발생
#  conf.int = TRUE,
  risk.table = FALSE,         # ST는 risk table 빼는 게 보통 더 깔끔
  xlab = "Days",
  ylab = "Cumulative incidence of stent thrombosis"
)

p$plot <- p$plot +
  geom_vline(xintercept = 365, linetype = "dashed") +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    limits = c(0, 0.015)      # 0~1.5% (필요시 0.02로)
  ) +
  labs(color = "Stent type")

print(p)

