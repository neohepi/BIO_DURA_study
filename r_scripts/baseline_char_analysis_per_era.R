################################################################################
#
#  Era-specific Baseline Characteristics Comparison
#  
#  목적: 2014-2016 기간 BP-DES 성적이 나빴던 이유가 환자 특성 때문인지 확인
#
################################################################################

library(dplyr)
library(tidyr)
library(tableone)
library(knitr)
library(ggplot2)

################################################################################
# 1. 데이터 준비
################################################################################

cat("\n", rep("=", 80), "\n")
cat("ERA-SPECIFIC BASELINE CHARACTERISTICS ANALYSIS\n")
cat(rep("=", 80), "\n\n")

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

# Stent group 생성
matched_data <- matched_data %>%
  mutate(
    stent_group = factor(ifelse(BP == 1, "BP-DES", "DP-DES"), 
                         levels = c("DP-DES", "BP-DES"))
  )

# 분석용 변수 확인
cat("=== Data Overview ===\n")
cat("Total patients:", nrow(matched_data), "\n\n")

cat("Era distribution:\n")
print(table(matched_data$era, matched_data$stent_group))


################################################################################
# 2. Era별 전체 환자 특성 비교 (BP-DES + DP-DES combined)
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 1: OVERALL PATIENT CHARACTERISTICS BY ERA\n")
cat(rep("=", 80), "\n")

# 분석 변수 목록
vars_baseline <- c("age", "male", "BMI", "HTN", "DM", "DL", "smoking", 
                   "CKD", "prior_MI", "prior_PCI", "prior_CABG", "prior_stroke",
                   "LVEF", "is_ACS")

vars_lesion <- c("is_LM", "is_LAD", "is_LCX", "is_RCA", "is_multivessel",
                 "is_bifurcation", "is_CTO", "is_ISR", "is_complex_lesion",
                 "num_stents", "total_stent_length", "avg_stent_diameter")

vars_all <- c(vars_baseline, vars_lesion)

# 사용 가능한 변수만 선택
vars_available <- vars_all[vars_all %in% names(matched_data)]
cat("\nAvailable variables for analysis:", length(vars_available), "\n")

# Factor 변수 지정
factor_vars <- c("male", "HTN", "DM", "DL", "smoking", "CKD", 
                 "prior_MI", "prior_PCI", "prior_CABG", "prior_stroke",
                 "is_ACS", "is_LM", "is_LAD", "is_LCX", "is_RCA", 
                 "is_multivessel", "is_bifurcation", "is_CTO", "is_ISR",
                 "is_complex_lesion")
factor_vars <- factor_vars[factor_vars %in% vars_available]

# TableOne 생성 - Era별 전체 비교
tab_era_overall <- CreateTableOne(
  vars = vars_available,
  strata = "era",
  data = matched_data,
  factorVars = factor_vars,
  addOverall = TRUE
)

cat("\n=== Overall Patient Characteristics by Era ===\n")
print(tab_era_overall, showAllLevels = TRUE, formatOptions = list(big.mark = ","))


################################################################################
# 3. Era × Stent Type 별 비교 (핵심 분석)
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 2: ERA × STENT TYPE STRATIFIED ANALYSIS\n")
cat(rep("=", 80), "\n")

# Era와 Stent type 조합 변수 생성
matched_data <- matched_data %>%
  mutate(
    era_stent = paste(era, stent_group, sep = "_"),
    era_stent = factor(era_stent, levels = c(
      "2010-2013_DP-DES", "2010-2013_BP-DES",
      "2014-2016_DP-DES", "2014-2016_BP-DES",
      "2017-2021_DP-DES", "2017-2021_BP-DES"
    ))
  )

# 각 Era 내에서 BP-DES vs DP-DES 비교
for (era_name in c("2010-2013", "2014-2016", "2017-2021")) {
  
  cat("\n", rep("-", 60), "\n")
  cat("Era:", era_name, "\n")
  cat(rep("-", 60), "\n")
  
  era_data <- matched_data %>% filter(era == era_name)
  
  tab_era <- CreateTableOne(
    vars = vars_available,
    strata = "stent_group",
    data = era_data,
    factorVars = factor_vars
  )
  
  print(tab_era, showAllLevels = FALSE, formatOptions = list(big.mark = ","))
}


################################################################################
# 4. 핵심 복잡도 지표 Era별 트렌드 분석
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 3: COMPLEXITY INDICATORS TREND BY ERA\n")
cat(rep("=", 80), "\n")

# 복잡도 지표 요약
complexity_summary <- matched_data %>%
  group_by(era, stent_group) %>%
  summarise(
    N = n(),
    
    # Demographics
    Age_mean = round(mean(age, na.rm = TRUE), 1),
    Male_pct = round(mean(as.numeric(as.character(male)), na.rm = TRUE) * 100, 1),
    DM_pct = round(mean(as.numeric(as.character(DM)), na.rm = TRUE) * 100, 1),
    CKD_pct = round(mean(as.numeric(as.character(CKD)), na.rm = TRUE) * 100, 1),
    
    # Clinical presentation
    ACS_pct = round(mean(as.numeric(as.character(is_ACS)), na.rm = TRUE) * 100, 1),
    
    # Lesion complexity
    LM_pct = round(mean(is_LM, na.rm = TRUE) * 100, 1),
    Multivessel_pct = round(mean(is_multivessel, na.rm = TRUE) * 100, 1),
    Bifurcation_pct = round(mean(is_bifurcation, na.rm = TRUE) * 100, 1),
    CTO_pct = round(mean(is_CTO, na.rm = TRUE) * 100, 1),
    ISR_pct = round(mean(is_ISR, na.rm = TRUE) * 100, 1),
    
    # Procedural
    Stents_mean = round(mean(num_stents, na.rm = TRUE), 2),
    Length_mean = round(mean(total_stent_length, na.rm = TRUE), 1),
    
    .groups = "drop"
  )

cat("\n=== Complexity Indicators by Era and Stent Type ===\n")
print(complexity_summary, n = 20)

# CSV 저장
write.csv(complexity_summary, "Table_Era_Complexity_Summary.csv", row.names = FALSE)
cat("\nSaved: Table_Era_Complexity_Summary.csv\n")


################################################################################
# 5. 복잡도 지표 시각화
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 4: VISUALIZATION OF COMPLEXITY TRENDS\n")
cat(rep("=", 80), "\n")

# Long format으로 변환
complexity_long <- complexity_summary %>%
  select(era, stent_group, N, ACS_pct, LM_pct, Multivessel_pct, 
         Bifurcation_pct, CTO_pct, ISR_pct, DM_pct, CKD_pct) %>%
  pivot_longer(
    cols = c(ACS_pct, LM_pct, Multivessel_pct, Bifurcation_pct, 
             CTO_pct, ISR_pct, DM_pct, CKD_pct),
    names_to = "Variable",
    values_to = "Percentage"
  ) %>%
  mutate(
    Variable = gsub("_pct", "", Variable),
    Variable = factor(Variable, levels = c("ACS", "DM", "CKD", "LM", 
                                           "Multivessel", "Bifurcation", 
                                           "CTO", "ISR"))
  )

# Bar plot
p_complexity <- ggplot(complexity_long, aes(x = era, y = Percentage, fill = stent_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("DP-DES" = "#4575B4", "BP-DES" = "#D73027"),
                    name = "Stent Type") +
  labs(
    title = "Patient and Lesion Complexity by Era and Stent Type",
    x = "Era",
    y = "Percentage (%)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_complexity)

ggsave("Figure_Era_Complexity_Trends.png", p_complexity, 
       width = 12, height = 8, dpi = 300)
cat("Saved: Figure_Era_Complexity_Trends.png\n")


################################################################################
# 6. 통계적 검정: Era별 복잡도 차이
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 5: STATISTICAL TESTING FOR ERA DIFFERENCES\n")
cat(rep("=", 80), "\n")

# Era 간 복잡도 차이 검정 (전체 환자)
cat("\n=== Trend Test Across Eras (All Patients) ===\n")

complexity_vars <- c("is_ACS", "DM", "CKD", "is_LM", "is_multivessel", 
                     "is_bifurcation", "is_CTO", "is_ISR", "num_stents", 
                     "total_stent_length")

trend_results <- data.frame(
  Variable = character(),
  Era_2010_2013 = character(),
  Era_2014_2016 = character(),
  Era_2017_2021 = character(),
  P_trend = numeric(),
  stringsAsFactors = FALSE
)

for (var in complexity_vars) {
  if (!var %in% names(matched_data)) next
  
  # Era를 numeric으로 변환하여 trend test
  matched_data$era_num <- as.numeric(matched_data$era)
  
  # 변수 타입에 따른 처리
  if (var %in% c("num_stents", "total_stent_length", "age", "BMI", "LVEF")) {
    # Continuous variable - linear regression for trend
    formula_trend <- as.formula(paste(var, "~ era_num"))
    fit <- lm(formula_trend, data = matched_data)
    p_trend <- summary(fit)$coefficients["era_num", "Pr(>|t|)"]
    
    # Era별 mean (SD)
    era_stats <- matched_data %>%
      group_by(era) %>%
      summarise(
        stat = sprintf("%.1f (%.1f)", 
                       mean(get(var), na.rm = TRUE),
                       sd(get(var), na.rm = TRUE)),
        .groups = "drop"
      )
  } else {
    # Categorical variable - Cochran-Armitage trend test approximation
    var_numeric <- as.numeric(as.character(matched_data[[var]]))
    fit <- glm(var_numeric ~ era_num, data = matched_data, family = binomial)
    p_trend <- summary(fit)$coefficients["era_num", "Pr(>|z|)"]
    
    # Era별 n (%)
    era_stats <- matched_data %>%
      group_by(era) %>%
      summarise(
        n_event = sum(as.numeric(as.character(get(var))), na.rm = TRUE),
        n_total = n(),
        stat = sprintf("%d (%.1f%%)", n_event, n_event/n_total * 100),
        .groups = "drop"
      )
  }
  
  trend_results <- rbind(trend_results, data.frame(
    Variable = var,
    Era_2010_2013 = era_stats$stat[1],
    Era_2014_2016 = era_stats$stat[2],
    Era_2017_2021 = era_stats$stat[3],
    P_trend = round(p_trend, 4),
    stringsAsFactors = FALSE
  ))
}

cat("\nTrend Analysis Results:\n")
print(trend_results)

write.csv(trend_results, "Table_Era_Trend_Analysis.csv", row.names = FALSE)
cat("\nSaved: Table_Era_Trend_Analysis.csv\n")


################################################################################
# 7. BP-DES 내에서 Era별 비교
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 6: BP-DES ONLY - ERA COMPARISON\n")
cat(rep("=", 80), "\n")

bp_only <- matched_data %>% filter(BP == 1)

cat("\n=== BP-DES Patient Characteristics by Era ===\n")

tab_bp_era <- CreateTableOne(
  vars = vars_available,
  strata = "era",
  data = bp_only,
  factorVars = factor_vars
)

print(tab_bp_era, showAllLevels = FALSE, formatOptions = list(big.mark = ","))

# BP-DES 내 Era별 복잡도 요약
bp_complexity <- bp_only %>%
  group_by(era) %>%
  summarise(
    N = n(),
    ACS_pct = round(mean(as.numeric(as.character(is_ACS)), na.rm = TRUE) * 100, 1),
    DM_pct = round(mean(as.numeric(as.character(DM)), na.rm = TRUE) * 100, 1),
    CKD_pct = round(mean(as.numeric(as.character(CKD)), na.rm = TRUE) * 100, 1),
    LM_pct = round(mean(is_LM, na.rm = TRUE) * 100, 1),
    Multivessel_pct = round(mean(is_multivessel, na.rm = TRUE) * 100, 1),
    Bifurcation_pct = round(mean(is_bifurcation, na.rm = TRUE) * 100, 1),
    CTO_pct = round(mean(is_CTO, na.rm = TRUE) * 100, 1),
    ISR_pct = round(mean(is_ISR, na.rm = TRUE) * 100, 1),
    Stents_mean = round(mean(num_stents, na.rm = TRUE), 2),
    Length_mean = round(mean(total_stent_length, na.rm = TRUE), 1),
    .groups = "drop"
  )

cat("\n=== BP-DES Complexity Summary by Era ===\n")
print(bp_complexity)


################################################################################
# 8. DP-DES 내에서 Era별 비교
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 7: DP-DES ONLY - ERA COMPARISON\n")
cat(rep("=", 80), "\n")

dp_only <- matched_data %>% filter(BP == 0)

cat("\n=== DP-DES Patient Characteristics by Era ===\n")

tab_dp_era <- CreateTableOne(
  vars = vars_available,
  strata = "era",
  data = dp_only,
  factorVars = factor_vars
)

print(tab_dp_era, showAllLevels = FALSE, formatOptions = list(big.mark = ","))


################################################################################
# 9. 2014-2016 기간 심층 분석
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 8: DEEP DIVE INTO 2014-2016 PERIOD\n")
cat(rep("=", 80), "\n")

period2 <- matched_data %>% filter(era == "2014-2016")

cat("\n=== 2014-2016 Period: BP-DES vs DP-DES ===\n")
cat("Total patients in 2014-2016:", nrow(period2), "\n")
cat("BP-DES:", sum(period2$BP == 1), "\n")
cat("DP-DES:", sum(period2$BP == 0), "\n\n")

# 2014-2016 상세 비교
tab_period2 <- CreateTableOne(
  vars = vars_available,
  strata = "stent_group",
  data = period2,
  factorVars = factor_vars
)

print(tab_period2, showAllLevels = FALSE, formatOptions = list(big.mark = ","))

# 2014-2016 내 BP-DES 스텐트 분포
if ("stent_names" %in% names(period2)) {
  cat("\n=== Stent Distribution in 2014-2016 (BP-DES) ===\n")
  bp_period2 <- period2 %>% filter(BP == 1)
  stent_dist <- bp_period2 %>%
    group_by(stent_names) %>%
    summarise(
      N = n(),
      TLF_events = sum(composite_event, na.rm = TRUE),
      TLF_rate = round(TLF_events / N * 100, 1),
      .groups = "drop"
    ) %>%
    arrange(desc(N))
  
  print(stent_dist)
}


################################################################################
# 10. Outcome과 복잡도 관계 분석
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 9: COMPLEXITY AND OUTCOME RELATIONSHIP\n")
cat(rep("=", 80), "\n")

# 복잡도 점수 생성 (간단한 버전)
matched_data <- matched_data %>%
  mutate(
    complexity_score = as.numeric(is_LM) + 
      as.numeric(is_multivessel) + 
      as.numeric(is_bifurcation) + 
      as.numeric(is_CTO) + 
      as.numeric(is_ISR) +
      as.numeric(as.character(DM)) +
      as.numeric(as.character(CKD)),
    complexity_group = case_when(
      complexity_score <= 1 ~ "Low (0-1)",
      complexity_score <= 3 ~ "Moderate (2-3)",
      complexity_score >= 4 ~ "High (≥4)"
    ),
    complexity_group = factor(complexity_group, 
                              levels = c("Low (0-1)", "Moderate (2-3)", "High (≥4)"))
  )

# Era별 복잡도 분포
cat("\n=== Complexity Score Distribution by Era ===\n")
complexity_dist <- matched_data %>%
  group_by(era, complexity_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(era) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

print(complexity_dist)

# 복잡도별 TLF rate
cat("\n=== TLF Rate by Complexity Group and Stent Type ===\n")
tlf_by_complexity <- matched_data %>%
  group_by(complexity_group, stent_group) %>%
  summarise(
    N = n(),
    TLF = sum(composite_event, na.rm = TRUE),
    TLF_rate = round(TLF / N * 100, 1),
    .groups = "drop"
  )

print(tlf_by_complexity)

# 복잡도별 TLF Forest plot 데이터
cat("\n=== TLF by Complexity Group, Era, and Stent Type ===\n")
tlf_complex_era <- matched_data %>%
  group_by(era, complexity_group, stent_group) %>%
  summarise(
    N = n(),
    TLF = sum(composite_event, na.rm = TRUE),
    TLF_rate = round(TLF / N * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(era, complexity_group, stent_group)

print(tlf_complex_era, n = 30)


################################################################################
# 11. Summary Table 생성
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 10: SUMMARY\n")
cat(rep("=", 80), "\n")

cat("\n=== KEY FINDINGS ===\n\n")

# Era별 주요 지표 비교
summary_table <- matched_data %>%
  group_by(era, stent_group) %>%
  summarise(
    N = n(),
    
    # Patient complexity
    ACS = sprintf("%.1f%%", mean(as.numeric(as.character(is_ACS)), na.rm = TRUE) * 100),
    DM = sprintf("%.1f%%", mean(as.numeric(as.character(DM)), na.rm = TRUE) * 100),
    CKD = sprintf("%.1f%%", mean(as.numeric(as.character(CKD)), na.rm = TRUE) * 100),
    
    # Lesion complexity  
    Multivessel = sprintf("%.1f%%", mean(is_multivessel, na.rm = TRUE) * 100),
    Bifurcation = sprintf("%.1f%%", mean(is_bifurcation, na.rm = TRUE) * 100),
    CTO = sprintf("%.1f%%", mean(is_CTO, na.rm = TRUE) * 100),
    
    # Outcome
    TLF_rate = sprintf("%.1f%%", sum(composite_event, na.rm = TRUE) / n() * 100),
    
    .groups = "drop"
  )

cat("Summary Table:\n")
print(summary_table, n = 20)

write.csv(summary_table, "Table_Era_Summary_Final.csv", row.names = FALSE)
cat("\nSaved: Table_Era_Summary_Final.csv\n")


################################################################################
# 12. Heatmap 시각화
################################################################################

# Heatmap을 위한 데이터 준비
heatmap_data <- matched_data %>%
  group_by(era, stent_group) %>%
  summarise(
    ACS = mean(as.numeric(as.character(is_ACS)), na.rm = TRUE) * 100,
    DM = mean(as.numeric(as.character(DM)), na.rm = TRUE) * 100,
    CKD = mean(as.numeric(as.character(CKD)), na.rm = TRUE) * 100,
    LM = mean(is_LM, na.rm = TRUE) * 100,
    Multivessel = mean(is_multivessel, na.rm = TRUE) * 100,
    Bifurcation = mean(is_bifurcation, na.rm = TRUE) * 100,
    CTO = mean(is_CTO, na.rm = TRUE) * 100,
    ISR = mean(is_ISR, na.rm = TRUE) * 100,
    TLF = sum(composite_event, na.rm = TRUE) / n() * 100,
    .groups = "drop"
  ) %>%
  mutate(era_stent = paste(era, stent_group, sep = "\n"))

# Long format
heatmap_long <- heatmap_data %>%
  select(era_stent, ACS, DM, CKD, LM, Multivessel, Bifurcation, CTO, ISR, TLF) %>%
  pivot_longer(cols = -era_stent, names_to = "Variable", values_to = "Value") %>%
  mutate(
    Variable = factor(Variable, levels = c("ACS", "DM", "CKD", "LM", 
                                           "Multivessel", "Bifurcation", 
                                           "CTO", "ISR", "TLF"))
  )

p_heatmap <- ggplot(heatmap_long, aes(x = Variable, y = era_stent, fill = Value)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.1f", Value)), color = "black", size = 3) +
  scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
                       midpoint = 15, name = "Percentage (%)") +
  labs(
    title = "Patient and Lesion Complexity Heatmap by Era and Stent Type",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )

print(p_heatmap)

ggsave("Figure_Era_Complexity_Heatmap.png", p_heatmap, 
       width = 10, height = 6, dpi = 300)
cat("Saved: Figure_Era_Complexity_Heatmap.png\n")


cat("\n", rep("=", 80), "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n")