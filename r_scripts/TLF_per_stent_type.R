################################################################################
#
#  Individual Stent-specific Cox Regression Analysis
#  
#  목적: 어떤 스텐트가 실제로 outcome에 영향을 미쳤는지 확인
#  Reference: Xience (가장 많이 사용된 DP-DES)
#
################################################################################

library(dplyr)
library(tidyr)
library(survival)
library(ggplot2)
library(forestplot)
library(broom)

################################################################################
# 1. 데이터 준비
################################################################################

cat("\n", rep("=", 80), "\n")
cat("INDIVIDUAL STENT-SPECIFIC COX REGRESSION ANALYSIS\n")
cat(rep("=", 80), "\n\n")

# stent_names 변수 확인
if (!"stent_names" %in% names(matched_data)) {
  stop("stent_names variable not found in matched_data")
}

# Era 변수 확인/생성
if (!"era" %in% names(matched_data)) {
  matched_data <- matched_data %>%
    mutate(
      procedure_year = as.numeric(format(CAG_date, "%Y")),
      era = case_when(
        procedure_year >= 2010 & procedure_year <= 2013 ~ "2010-2013",
        procedure_year >= 2014 & procedure_year <= 2016 ~ "2014-2016",
        procedure_year >= 2017 & procedure_year <= 2021 ~ "2017-2021",
        TRUE ~ NA_character_
      ),
      era = factor(era, levels = c("2010-2013", "2014-2016", "2017-2021"))
    )
}

# 스텐트 분포 확인
cat("=== Overall Stent Distribution ===\n")
stent_dist <- matched_data %>%
  group_by(stent_names) %>%
  summarise(
    N = n(),
    TLF = sum(composite_event, na.rm = TRUE),
    TLF_rate = round(TLF / N * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(N))

print(stent_dist, n = 30)


################################################################################
# 2. 주요 스텐트 식별 및 그룹화
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 1: STENT CLASSIFICATION\n")
cat(rep("=", 80), "\n")

# 주요 단일 스텐트만 추출 (combination 제외, N >= 30)
major_stents <- stent_dist %>%
  filter(N >= 30) %>%
  filter(!grepl(",", stent_names)) %>%  # combination 제외
  pull(stent_names)

cat("\nMajor stents (N >= 30, single stent):\n")
print(major_stents)

# 스텐트 타입 분류 (확장)
bp_stents <- c("SYNERGY", "Orsiro", "Nobori", "BioMatrix", "GENOSS DES", 
               "BioMime", "DESyne X2", "Ultimaster", "SYNERGY XD", "Firehawk")
dp_stents <- c("XIENCE Prime", "XIENCE Xpedition", "XIENCE Alpine", "XIENCE Sierra",
               "XIENCE", "Resolute Integrity", "Resolute Onyx", "PROMUS Element",
               "Promus PREMIER", "Promus", "Endeavor", "DESyne", "D+Storm")

# 분석용 데이터 생성 (단일 스텐트만)
analysis_data <- matched_data %>%
  filter(!grepl(",", stent_names)) %>%  # combination 제외
  mutate(
    stent_category = case_when(
      stent_names %in% bp_stents ~ "BP-DES",
      stent_names %in% dp_stents ~ "DP-DES",
      TRUE ~ "Other"
    )
  )

cat("\n=== Stent Distribution After Filtering ===\n")
cat("Total patients with single stent:", nrow(analysis_data), "\n")
cat("Excluded (combination stents):", nrow(matched_data) - nrow(analysis_data), "\n\n")

# 스텐트별 분포 재확인
stent_summary <- analysis_data %>%
  group_by(stent_names, stent_category) %>%
  summarise(
    N = n(),
    TLF = sum(composite_event, na.rm = TRUE),
    TLF_rate = round(TLF / N * 100, 1),
    Deaths = sum(any_death, na.rm = TRUE),
    Death_rate = round(Deaths / N * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(stent_category, desc(N))

cat("Stent Summary:\n")
print(stent_summary, n = 30)


################################################################################
# 3. Time-to-event 변수 확인/생성
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 2: TIME-TO-EVENT PREPARATION\n")
cat(rep("=", 80), "\n")

# 사용 가능한 time 변수 확인
cat("Available time variables:\n")
time_vars <- grep("time|fu_days|days", names(analysis_data), value = TRUE, ignore.case = TRUE)
print(time_vars)

# Time 변수 생성
analysis_data <- analysis_data %>%
  mutate(
    # TLF time 계산
    time_to_tlf = case_when(
      composite_event == 1 & !is.na(composite_date) ~ 
        pmin(as.numeric(composite_date - CAG_date), 1825),
      !is.na(fu_days_correct) ~ pmin(fu_days_correct, 1825),
      !is.na(fu_days) ~ pmin(fu_days, 1825),
      TRUE ~ pmin(as.numeric(last_fu_date - CAG_date), 1825)
    ),
    # Event (5년 내)
    event_tlf = case_when(
      composite_event == 1 & time_to_tlf <= 1825 ~ 1,
      TRUE ~ 0
    )
  )

cat("Time variable summary:\n")
cat("Median follow-up:", round(median(analysis_data$time_to_tlf, na.rm = TRUE)/365.25, 2), "years\n")
cat("Total TLF events:", sum(analysis_data$event_tlf, na.rm = TRUE), "\n")


################################################################################
# 4. Reference 스텐트 설정 (Xience 계열)
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 3: SETTING REFERENCE STENT\n")
cat(rep("=", 80), "\n")

# Xience 계열 합치기
analysis_data <- analysis_data %>%
  mutate(
    stent_grouped = case_when(
      grepl("XIENCE", stent_names) ~ "XIENCE (all)",
      grepl("Resolute", stent_names) ~ "Resolute (all)",
      grepl("PROMUS", stent_names) ~ "PROMUS (all)",
      TRUE ~ stent_names
    )
  )

# 그룹별 분포
cat("\nGrouped Stent Distribution:\n")
grouped_dist <- analysis_data %>%
  group_by(stent_grouped, stent_category) %>%
  summarise(
    N = n(),
    TLF = sum(event_tlf, na.rm = TRUE),
    TLF_rate = round(TLF / N * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(N))

print(grouped_dist, n = 20)

# Reference 설정: XIENCE (가장 많이 사용된 DP-DES)
# N >= 50인 스텐트만 분석에 포함
stents_for_analysis <- grouped_dist %>%
  filter(N >= 50) %>%
  pull(stent_grouped)

cat("\nStents included in analysis (N >= 50):\n")
print(stents_for_analysis)

# Factor 설정 (XIENCE를 reference로)
analysis_data <- analysis_data %>%
  filter(stent_grouped %in% stents_for_analysis) %>%
  mutate(
    stent_factor = factor(stent_grouped),
    stent_factor = relevel(stent_factor, ref = "XIENCE (all)")
  )

cat("\nReference stent: XIENCE (all)\n")
cat("Patients in final analysis:", nrow(analysis_data), "\n")


################################################################################
# 5. Univariable Cox Regression (각 스텐트 vs XIENCE)
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 4: UNIVARIABLE COX REGRESSION\n")
cat(rep("=", 80), "\n")

# Univariable Cox model
cox_uni <- coxph(Surv(time_to_tlf, event_tlf) ~ stent_factor, data = analysis_data)

cat("\n=== Univariable Cox Model (Reference: XIENCE) ===\n")
print(summary(cox_uni))

# 결과 추출
uni_results <- tidy(cox_uni, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(
    stent = gsub("stent_factor", "", term),
    HR = round(estimate, 2),
    CI_lower = round(conf.low, 2),
    CI_upper = round(conf.high, 2),
    P_value = round(p.value, 4),
    HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, CI_lower, CI_upper)
  ) %>%
  select(stent, HR, CI_lower, CI_upper, HR_CI, P_value)

cat("\n=== Univariable Results ===\n")
print(uni_results)


################################################################################
# 6. Multivariable Cox Regression (adjusted)
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 5: MULTIVARIABLE COX REGRESSION\n")
cat(rep("=", 80), "\n")

# Adjustment 변수
# age, male, DM, CKD, ACS, multivessel, bifurcation, total_stent_length

# 변수 준비
analysis_data <- analysis_data %>%
  mutate(
    male_n = as.numeric(as.character(male)),
    DM_n = as.numeric(as.character(DM)),
    CKD_n = as.numeric(as.character(CKD)),
    ACS_n = as.numeric(as.character(is_ACS)),
    multivessel_n = as.numeric(is_multivessel),
    bifurcation_n = as.numeric(is_bifurcation)
  )

# Multivariable Cox model
cox_multi <- coxph(
  Surv(time_to_tlf, event_tlf) ~ stent_factor + age + male_n + DM_n + CKD_n + 
    ACS_n + multivessel_n + bifurcation_n + total_stent_length,
  data = analysis_data
)

cat("\n=== Multivariable Cox Model ===\n")
print(summary(cox_multi))

# 스텐트 관련 결과만 추출
multi_results <- tidy(cox_multi, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(grepl("stent_factor", term)) %>%
  mutate(
    stent = gsub("stent_factor", "", term),
    HR = round(estimate, 2),
    CI_lower = round(conf.low, 2),
    CI_upper = round(conf.high, 2),
    P_value = round(p.value, 4),
    HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, CI_lower, CI_upper)
  ) %>%
  select(stent, HR, CI_lower, CI_upper, HR_CI, P_value)

cat("\n=== Multivariable Results (Stent Effects) ===\n")
print(multi_results)


################################################################################
# 7. 결과 합치기 및 Forest Plot 데이터 준비
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 6: COMBINED RESULTS AND FOREST PLOT\n")
cat(rep("=", 80), "\n")

# Univariable과 Multivariable 결과 합치기
combined_results <- uni_results %>%
  rename(
    HR_uni = HR,
    CI_lower_uni = CI_lower,
    CI_upper_uni = CI_upper,
    HR_CI_uni = HR_CI,
    P_uni = P_value
  ) %>%
  left_join(
    multi_results %>%
      rename(
        HR_adj = HR,
        CI_lower_adj = CI_lower,
        CI_upper_adj = CI_upper,
        HR_CI_adj = HR_CI,
        P_adj = P_value
      ),
    by = "stent"
  )

# 스텐트 카테고리 추가
combined_results <- combined_results %>%
  mutate(
    category = case_when(
      stent %in% c("SYNERGY", "Orsiro", "Nobori", "BioMatrix", "GENOSS DES", 
                   "BioMime", "DESyne X2", "Ultimaster") ~ "BP-DES",
      TRUE ~ "DP-DES"
    )
  )

# N과 Event 수 추가
stent_counts <- analysis_data %>%
  group_by(stent_grouped) %>%
  summarise(
    N = n(),
    Events = sum(event_tlf, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(stent = stent_grouped)

# Reference (XIENCE) 정보 추가
xience_info <- analysis_data %>%
  filter(stent_grouped == "XIENCE (all)") %>%
  summarise(
    N = n(),
    Events = sum(event_tlf, na.rm = TRUE)
  )

reference_row <- data.frame(
  stent = "XIENCE (all)",
  HR_uni = 1.00,
  CI_lower_uni = NA,
  CI_upper_uni = NA,
  HR_CI_uni = "1.00 (Reference)",
  P_uni = NA,
  HR_adj = 1.00,
  CI_lower_adj = NA,
  CI_upper_adj = NA,
  HR_CI_adj = "1.00 (Reference)",
  P_adj = NA,
  category = "DP-DES"
)

final_results <- bind_rows(reference_row, combined_results) %>%
  left_join(stent_counts, by = "stent") %>%
  mutate(
    N = ifelse(stent == "XIENCE (all)", xience_info$N, N),
    Events = ifelse(stent == "XIENCE (all)", xience_info$Events, Events)
  ) %>%
  arrange(category, desc(N))

cat("\n=== Final Combined Results ===\n")
print(final_results %>% select(stent, category, N, Events, HR_CI_uni, P_uni, HR_CI_adj, P_adj))

# CSV 저장
write.csv(final_results, "Table_Individual_Stent_Cox_Results.csv", row.names = FALSE)
cat("\nSaved: Table_Individual_Stent_Cox_Results.csv\n")


################################################################################
# 8. Forest Plot 생성 (Adjusted HR)
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 7: FOREST PLOT\n")
cat(rep("=", 80), "\n")

# Forest plot 데이터 준비
forest_data <- final_results %>%
  filter(!is.na(HR_adj) | stent == "XIENCE (all)") %>%
  arrange(category, desc(N)) %>%
  mutate(
    y_order = row_number(),
    label = sprintf("%s (N=%d, %d events)", stent, N, Events),
    HR_display = ifelse(stent == "XIENCE (all)", 1.00, HR_adj),
    CI_lower_display = ifelse(stent == "XIENCE (all)", 1.00, CI_lower_adj),
    CI_upper_display = ifelse(stent == "XIENCE (all)", 1.00, CI_upper_adj)
  )

# ggplot forest plot
p_forest <- ggplot(forest_data, aes(x = HR_display, y = reorder(label, -y_order))) +
  # Reference line
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  # Confidence intervals
  geom_errorbarh(
    aes(xmin = CI_lower_display, xmax = CI_upper_display),
    height = 0.2,
    color = ifelse(forest_data$category == "BP-DES", "#D73027", "#4575B4")
  ) +
  # Point estimates
  geom_point(
    aes(color = category, shape = category),
    size = 3
  ) +
  # Category separation
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  # Scales
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
    limits = c(0.4, 3.5)
  ) +
  scale_color_manual(values = c("BP-DES" = "#D73027", "DP-DES" = "#4575B4")) +
  scale_shape_manual(values = c("BP-DES" = 16, "DP-DES" = 15)) +
  # Labels
  labs(
    title = "Adjusted Hazard Ratios for TLF by Individual Stent",
    subtitle = "Reference: XIENCE (all types)",
    x = "Hazard Ratio (95% CI)",
    y = "",
    caption = "Adjusted for age, sex, DM, CKD, ACS, multivessel, bifurcation, total stent length"
  ) +
  # Theme
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 9),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

print(p_forest)

ggsave("Figure_Individual_Stent_Forest_Plot.png", p_forest, 
       width = 10, height = 8, dpi = 300)
ggsave("Figure_Individual_Stent_Forest_Plot.pdf", p_forest, 
       width = 10, height = 8)
cat("Saved: Figure_Individual_Stent_Forest_Plot.png/pdf\n")


################################################################################
# 9. Era별 스텐트 분석
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 8: ERA-SPECIFIC STENT ANALYSIS\n")
cat(rep("=", 80), "\n")

# 각 Era별 스텐트 성적
era_stent_results <- list()

for (era_name in c("2010-2013", "2014-2016", "2017-2021")) {
  
  cat("\n", rep("-", 60), "\n")
  cat("Era:", era_name, "\n")
  cat(rep("-", 60), "\n")
  
  era_data <- analysis_data %>% filter(era == era_name)
  
  # 해당 era에서 스텐트 분포
  era_stent_dist <- era_data %>%
    group_by(stent_grouped, stent_category) %>%
    summarise(
      N = n(),
      Events = sum(event_tlf, na.rm = TRUE),
      TLF_rate = round(Events / N * 100, 1),
      .groups = "drop"
    ) %>%
    arrange(desc(N))
  
  cat("\nStent distribution in", era_name, ":\n")
  print(era_stent_dist)
  
  # Cox model (해당 era 내에서)
  # N >= 30인 스텐트만
  era_major_stents <- era_stent_dist %>%
    filter(N >= 30) %>%
    pull(stent_grouped)
  
  if (length(era_major_stents) >= 2 && "XIENCE (all)" %in% era_major_stents) {
    
    era_analysis <- era_data %>%
      filter(stent_grouped %in% era_major_stents) %>%
      mutate(
        stent_factor = factor(stent_grouped),
        stent_factor = relevel(stent_factor, ref = "XIENCE (all)")
      )
    
    # Univariable
    cox_era <- tryCatch({
      coxph(Surv(time_to_tlf, event_tlf) ~ stent_factor, data = era_analysis)
    }, error = function(e) NULL)
    
    if (!is.null(cox_era)) {
      era_results <- tidy(cox_era, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(grepl("stent_factor", term)) %>%
        mutate(
          stent = gsub("stent_factor", "", term),
          HR = round(estimate, 2),
          CI_lower = round(conf.low, 2),
          CI_upper = round(conf.high, 2),
          P_value = round(p.value, 3),
          era = era_name
        ) %>%
        select(era, stent, HR, CI_lower, CI_upper, P_value)
      
      cat("\nCox Results (Reference: XIENCE):\n")
      print(era_results)
      
      era_stent_results[[era_name]] <- era_results
    }
  }
}

# Era별 결과 합치기
if (length(era_stent_results) > 0) {
  all_era_results <- bind_rows(era_stent_results)
  write.csv(all_era_results, "Table_Era_Specific_Stent_Cox.csv", row.names = FALSE)
  cat("\nSaved: Table_Era_Specific_Stent_Cox.csv\n")
}


################################################################################
# 10. BP-DES 내 스텐트 비교 (Orsiro를 reference로)
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 9: WITHIN BP-DES COMPARISON (Reference: Orsiro)\n")
cat(rep("=", 80), "\n")

bp_analysis <- analysis_data %>%
  filter(stent_category == "BP-DES") %>%
  filter(stent_grouped %in% stents_for_analysis)

# Orsiro를 reference로
bp_analysis <- bp_analysis %>%
  mutate(
    stent_factor = factor(stent_grouped),
    stent_factor = relevel(stent_factor, ref = "Orsiro")
  )

cat("\nBP-DES Stent Distribution:\n")
bp_dist <- bp_analysis %>%
  group_by(stent_grouped) %>%
  summarise(
    N = n(),
    Events = sum(event_tlf, na.rm = TRUE),
    TLF_rate = round(Events / N * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(N))
print(bp_dist)

# Cox model within BP-DES
cox_bp <- tryCatch({
  coxph(Surv(time_to_tlf, event_tlf) ~ stent_factor, data = bp_analysis)
}, error = function(e) NULL)

if (!is.null(cox_bp)) {
  cat("\n=== BP-DES Internal Comparison (Reference: Orsiro) ===\n")
  print(summary(cox_bp))
  
  bp_results <- tidy(cox_bp, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(grepl("stent_factor", term)) %>%
    mutate(
      stent = gsub("stent_factor", "", term),
      HR = round(estimate, 2),
      CI_lower = round(conf.low, 2),
      CI_upper = round(conf.high, 2),
      P_value = round(p.value, 3)
    ) %>%
    select(stent, HR, CI_lower, CI_upper, P_value)
  
  cat("\nBP-DES Internal Comparison Results:\n")
  print(bp_results)
}


################################################################################
# 11. DP-DES 내 스텐트 비교
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 10: WITHIN DP-DES COMPARISON (Reference: XIENCE)\n")
cat(rep("=", 80), "\n")

dp_analysis <- analysis_data %>%
  filter(stent_category == "DP-DES") %>%
  filter(stent_grouped %in% stents_for_analysis)

# XIENCE를 reference로
dp_analysis <- dp_analysis %>%
  mutate(
    stent_factor = factor(stent_grouped),
    stent_factor = relevel(stent_factor, ref = "XIENCE (all)")
  )

cat("\nDP-DES Stent Distribution:\n")
dp_dist <- dp_analysis %>%
  group_by(stent_grouped) %>%
  summarise(
    N = n(),
    Events = sum(event_tlf, na.rm = TRUE),
    TLF_rate = round(Events / N * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(N))
print(dp_dist)

# Cox model within DP-DES
cox_dp <- tryCatch({
  coxph(Surv(time_to_tlf, event_tlf) ~ stent_factor, data = dp_analysis)
}, error = function(e) NULL)

if (!is.null(cox_dp)) {
  cat("\n=== DP-DES Internal Comparison (Reference: XIENCE) ===\n")
  print(summary(cox_dp))
}


################################################################################
# 12. Summary Table
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 11: SUMMARY\n")
cat(rep("=", 80), "\n")

# 최종 요약 테이블
summary_table <- final_results %>%
  select(
    Stent = stent,
    Category = category,
    N,
    Events,
    `TLF Rate (%)` = Events,
    `Unadjusted HR (95% CI)` = HR_CI_uni,
    `P (unadj)` = P_uni,
    `Adjusted HR (95% CI)` = HR_CI_adj,
    `P (adj)` = P_adj
  ) %>%
  mutate(
    `TLF Rate (%)` = round(Events / N * 100, 1)
  )

cat("\n=== FINAL SUMMARY TABLE ===\n")
print(summary_table)

write.csv(summary_table, "Table_Stent_Comparison_Summary.csv", row.names = FALSE)
cat("\nSaved: Table_Stent_Comparison_Summary.csv\n")


################################################################################
# 13. Key Findings
################################################################################

cat("\n", rep("=", 80), "\n")
cat("KEY FINDINGS\n")
cat(rep("=", 80), "\n\n")

# 유의하게 나쁜 스텐트 (HR > 1, P < 0.1)
worse_stents <- final_results %>%
  filter(HR_adj > 1 & P_adj < 0.1)

if (nrow(worse_stents) > 0) {
  cat("Stents with WORSE outcomes vs XIENCE (P < 0.1):\n")
  print(worse_stents %>% select(stent, category, N, HR_CI_adj, P_adj))
}

# 유의하게 좋은 스텐트 (HR < 1, P < 0.1)
better_stents <- final_results %>%
  filter(HR_adj < 1 & P_adj < 0.1)

if (nrow(better_stents) > 0) {
  cat("\nStents with BETTER outcomes vs XIENCE (P < 0.1):\n")
  print(better_stents %>% select(stent, category, N, HR_CI_adj, P_adj))
}

# TLF rate 높은 순
cat("\nStents ranked by TLF rate (descending):\n")
print(
  final_results %>%
    mutate(TLF_rate = round(Events / N * 100, 1)) %>%
    arrange(desc(TLF_rate)) %>%
    select(stent, category, N, TLF_rate, HR_CI_adj)
)


cat("\n", rep("=", 80), "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n")