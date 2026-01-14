################################################################################
#
#  2014-2016 Orsiro Deep Dive Analysis
#  
#  목적: 왜 2014-2016 기간에 Orsiro가 유의하게 나빴는지 분석
#
################################################################################

library(dplyr)
library(tidyr)
library(survival)
library(ggplot2)
library(survminer)
library(broom)

################################################################################
# 1. 데이터 준비
################################################################################

cat("\n", rep("=", 80), "\n")
cat("2014-2016 ORSIRO DEEP DIVE ANALYSIS\n")
cat(rep("=", 80), "\n\n")

# Era 변수 확인
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

# 단일 스텐트만 필터링
analysis_data <- matched_data %>%
  filter(!grepl(",", stent_names))

# Time-to-event 변수 생성
analysis_data <- analysis_data %>%
  mutate(
    time_to_tlf = case_when(
      composite_event == 1 & !is.na(composite_date) ~ 
        pmin(as.numeric(composite_date - CAG_date), 1825),
      !is.na(fu_days_correct) ~ pmin(fu_days_correct, 1825),
      !is.na(fu_days) ~ pmin(fu_days, 1825),
      TRUE ~ pmin(as.numeric(last_fu_date - CAG_date), 1825)
    ),
    event_tlf = case_when(
      composite_event == 1 & time_to_tlf <= 1825 ~ 1,
      TRUE ~ 0
    )
  )

# Numeric 변수 생성
analysis_data <- analysis_data %>%
  mutate(
    male_n = as.numeric(as.character(male)),
    DM_n = as.numeric(as.character(DM)),
    CKD_n = as.numeric(as.character(CKD)),
    ACS_n = as.numeric(as.character(is_ACS)),
    multivessel_n = as.numeric(is_multivessel),
    bifurcation_n = as.numeric(is_bifurcation),
    complex_n = as.numeric(is_complex_lesion)
  )


################################################################################
# 2. Orsiro 사용 Era별 분포
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 1: ORSIRO USAGE BY ERA\n")
cat(rep("=", 80), "\n")

orsiro_era <- analysis_data %>%
  filter(stent_names == "Orsiro") %>%
  group_by(era) %>%
  summarise(
    N = n(),
    TLF = sum(event_tlf, na.rm = TRUE),
    TLF_rate = round(TLF / N * 100, 1),
    .groups = "drop"
  )

cat("\nOrsiro Usage by Era:\n")
print(orsiro_era)

# XIENCE 비교
xience_era <- analysis_data %>%
  filter(grepl("XIENCE", stent_names)) %>%
  group_by(era) %>%
  summarise(
    N = n(),
    TLF = sum(event_tlf, na.rm = TRUE),
    TLF_rate = round(TLF / N * 100, 1),
    .groups = "drop"
  )

cat("\nXIENCE Usage by Era:\n")
print(xience_era)


################################################################################
# 3. 2014-2016 Orsiro vs XIENCE 환자 특성 비교
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 2: 2014-2016 ORSIRO vs XIENCE BASELINE COMPARISON\n")
cat(rep("=", 80), "\n")

# 2014-2016 데이터
period2 <- analysis_data %>%
  filter(era == "2014-2016") %>%
  filter(stent_names == "Orsiro" | grepl("XIENCE", stent_names)) %>%
  mutate(
    stent_group = ifelse(stent_names == "Orsiro", "Orsiro", "XIENCE")
  )

cat("\n2014-2016 Orsiro vs XIENCE:\n")
cat("Orsiro N:", sum(period2$stent_group == "Orsiro"), "\n")
cat("XIENCE N:", sum(period2$stent_group == "XIENCE"), "\n")

# 기본 특성 비교
baseline_compare <- period2 %>%
  group_by(stent_group) %>%
  summarise(
    N = n(),
    Age = round(mean(age, na.rm = TRUE), 1),
    Male_pct = round(mean(male_n, na.rm = TRUE) * 100, 1),
    DM_pct = round(mean(DM_n, na.rm = TRUE) * 100, 1),
    CKD_pct = round(mean(CKD_n, na.rm = TRUE) * 100, 1),
    ACS_pct = round(mean(ACS_n, na.rm = TRUE) * 100, 1),
    Multivessel_pct = round(mean(multivessel_n, na.rm = TRUE) * 100, 1),
    Bifurcation_pct = round(mean(bifurcation_n, na.rm = TRUE) * 100, 1),
    LM_pct = round(mean(is_LM, na.rm = TRUE) * 100, 1),
    CTO_pct = round(mean(is_CTO, na.rm = TRUE) * 100, 1),
    ISR_pct = round(mean(is_ISR, na.rm = TRUE) * 100, 1),
    Complex_pct = round(mean(complex_n, na.rm = TRUE) * 100, 1),
    Num_stents = round(mean(num_stents, na.rm = TRUE), 2),
    Total_length = round(mean(total_stent_length, na.rm = TRUE), 1),
    Avg_diameter = round(mean(avg_stent_diameter, na.rm = TRUE), 2),
    TLF_rate = round(sum(event_tlf) / n() * 100, 1),
    .groups = "drop"
  )

cat("\n=== Baseline Characteristics (2014-2016) ===\n")
print(t(baseline_compare))

# 통계적 검정
cat("\n=== Statistical Tests ===\n")

# Continuous variables
cat("\nAge:", t.test(age ~ stent_group, data = period2)$p.value %>% round(3))
cat("\nNum stents:", t.test(num_stents ~ stent_group, data = period2)$p.value %>% round(3))
cat("\nTotal length:", t.test(total_stent_length ~ stent_group, data = period2)$p.value %>% round(3))
cat("\nAvg diameter:", t.test(avg_stent_diameter ~ stent_group, data = period2)$p.value %>% round(3))

# Categorical variables
cat("\n\nDM:", chisq.test(table(period2$stent_group, period2$DM_n))$p.value %>% round(3))
cat("\nACS:", chisq.test(table(period2$stent_group, period2$ACS_n))$p.value %>% round(3))
cat("\nMultivessel:", chisq.test(table(period2$stent_group, period2$multivessel_n))$p.value %>% round(3))
cat("\nBifurcation:", chisq.test(table(period2$stent_group, period2$bifurcation_n))$p.value %>% round(3))
cat("\n")


################################################################################
# 4. Small Vessel Analysis
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 3: SMALL VESSEL ANALYSIS\n")
cat(rep("=", 80), "\n")

# Small vessel 정의: avg_stent_diameter <= 2.75mm
period2 <- period2 %>%
  mutate(
    small_vessel = ifelse(avg_stent_diameter <= 2.75, 1, 0),
    vessel_size = ifelse(avg_stent_diameter <= 2.75, "Small (≤2.75mm)", "Large (>2.75mm)")
  )

cat("\n=== Vessel Size Distribution ===\n")
vessel_dist <- period2 %>%
  group_by(stent_group, vessel_size) %>%
  summarise(
    N = n(),
    TLF = sum(event_tlf),
    TLF_rate = round(TLF / N * 100, 1),
    .groups = "drop"
  )
print(vessel_dist)

# Small vessel에서 Orsiro vs XIENCE
cat("\n=== Small Vessel Cox Regression ===\n")
small_vessel_data <- period2 %>% filter(small_vessel == 1)
cat("Small vessel patients: Orsiro =", sum(small_vessel_data$stent_group == "Orsiro"),
    ", XIENCE =", sum(small_vessel_data$stent_group == "XIENCE"), "\n")

if (nrow(small_vessel_data) >= 30) {
  cox_small <- coxph(Surv(time_to_tlf, event_tlf) ~ stent_group, data = small_vessel_data)
  cat("\nSmall Vessel Results:\n")
  print(summary(cox_small))
}

# Large vessel에서 Orsiro vs XIENCE
cat("\n=== Large Vessel Cox Regression ===\n")
large_vessel_data <- period2 %>% filter(small_vessel == 0)
cat("Large vessel patients: Orsiro =", sum(large_vessel_data$stent_group == "Orsiro"),
    ", XIENCE =", sum(large_vessel_data$stent_group == "XIENCE"), "\n")

if (nrow(large_vessel_data) >= 30) {
  cox_large <- coxph(Surv(time_to_tlf, event_tlf) ~ stent_group, data = large_vessel_data)
  cat("\nLarge Vessel Results:\n")
  print(summary(cox_large))
}

# Interaction test
cat("\n=== Interaction Test (Stent × Vessel Size) ===\n")
cox_interaction <- coxph(Surv(time_to_tlf, event_tlf) ~ stent_group * vessel_size, 
                         data = period2)
print(summary(cox_interaction))


################################################################################
# 5. ACS vs Stable Angina Analysis
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 4: ACS vs STABLE ANGINA ANALYSIS\n")
cat(rep("=", 80), "\n")

# ACS 분포
cat("\n=== Clinical Presentation Distribution ===\n")
acs_dist <- period2 %>%
  group_by(stent_group, ACS_n) %>%
  summarise(
    N = n(),
    TLF = sum(event_tlf),
    TLF_rate = round(TLF / N * 100, 1),
    .groups = "drop"
  ) %>%
  mutate(Presentation = ifelse(ACS_n == 1, "ACS", "Stable"))
print(acs_dist)

# ACS에서 Orsiro vs XIENCE
cat("\n=== ACS Patients Cox Regression ===\n")
acs_data <- period2 %>% filter(ACS_n == 1)
cat("ACS patients: Orsiro =", sum(acs_data$stent_group == "Orsiro"),
    ", XIENCE =", sum(acs_data$stent_group == "XIENCE"), "\n")

if (nrow(acs_data) >= 30 && length(unique(acs_data$stent_group)) == 2) {
  cox_acs <- coxph(Surv(time_to_tlf, event_tlf) ~ stent_group, data = acs_data)
  acs_result <- tidy(cox_acs, conf.int = TRUE, exponentiate = TRUE)
  cat("\nACS Results: HR =", round(acs_result$estimate, 2), 
      "(", round(acs_result$conf.low, 2), "-", round(acs_result$conf.high, 2), ")",
      "P =", round(acs_result$p.value, 3), "\n")
}

# Stable에서 Orsiro vs XIENCE
cat("\n=== Stable Angina Patients Cox Regression ===\n")
stable_data <- period2 %>% filter(ACS_n == 0)
cat("Stable patients: Orsiro =", sum(stable_data$stent_group == "Orsiro"),
    ", XIENCE =", sum(stable_data$stent_group == "XIENCE"), "\n")

if (nrow(stable_data) >= 30 && length(unique(stable_data$stent_group)) == 2) {
  cox_stable <- coxph(Surv(time_to_tlf, event_tlf) ~ stent_group, data = stable_data)
  stable_result <- tidy(cox_stable, conf.int = TRUE, exponentiate = TRUE)
  cat("\nStable Results: HR =", round(stable_result$estimate, 2), 
      "(", round(stable_result$conf.low, 2), "-", round(stable_result$conf.high, 2), ")",
      "P =", round(stable_result$p.value, 3), "\n")
}


################################################################################
# 6. Complex Lesion Analysis
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 5: COMPLEX LESION ANALYSIS\n")
cat(rep("=", 80), "\n")

# Complex lesion 분포
cat("\n=== Complex Lesion Distribution ===\n")
complex_dist <- period2 %>%
  group_by(stent_group, complex_n) %>%
  summarise(
    N = n(),
    TLF = sum(event_tlf),
    TLF_rate = round(TLF / N * 100, 1),
    .groups = "drop"
  ) %>%
  mutate(Lesion = ifelse(complex_n == 1, "Complex (B2/C)", "Simple (A/B1)"))
print(complex_dist)

# Complex에서 Orsiro vs XIENCE
cat("\n=== Complex Lesion Cox Regression ===\n")
complex_data <- period2 %>% filter(complex_n == 1)
cat("Complex lesion patients: Orsiro =", sum(complex_data$stent_group == "Orsiro"),
    ", XIENCE =", sum(complex_data$stent_group == "XIENCE"), "\n")

if (nrow(complex_data) >= 30 && length(unique(complex_data$stent_group)) == 2) {
  cox_complex <- coxph(Surv(time_to_tlf, event_tlf) ~ stent_group, data = complex_data)
  complex_result <- tidy(cox_complex, conf.int = TRUE, exponentiate = TRUE)
  cat("\nComplex Lesion Results: HR =", round(complex_result$estimate, 2), 
      "(", round(complex_result$conf.low, 2), "-", round(complex_result$conf.high, 2), ")",
      "P =", round(complex_result$p.value, 3), "\n")
}

# Simple에서 Orsiro vs XIENCE
cat("\n=== Simple Lesion Cox Regression ===\n")
simple_data <- period2 %>% filter(complex_n == 0)
cat("Simple lesion patients: Orsiro =", sum(simple_data$stent_group == "Orsiro"),
    ", XIENCE =", sum(simple_data$stent_group == "XIENCE"), "\n")

if (nrow(simple_data) >= 30 && length(unique(simple_data$stent_group)) == 2) {
  cox_simple <- coxph(Surv(time_to_tlf, event_tlf) ~ stent_group, data = simple_data)
  simple_result <- tidy(cox_simple, conf.int = TRUE, exponentiate = TRUE)
  cat("\nSimple Lesion Results: HR =", round(simple_result$estimate, 2), 
      "(", round(simple_result$conf.low, 2), "-", round(simple_result$conf.high, 2), ")",
      "P =", round(simple_result$p.value, 3), "\n")
}


################################################################################
# 7. Early vs Late Event Analysis (Learning Curve)
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 6: EARLY vs LATE EVENT ANALYSIS (LEARNING CURVE)\n")
cat(rep("=", 80), "\n")

# 2014-2016 내에서 연도별 분석
period2 <- period2 %>%
  mutate(
    procedure_year = as.numeric(format(CAG_date, "%Y"))
  )

yearly_orsiro <- period2 %>%
  filter(stent_group == "Orsiro") %>%
  group_by(procedure_year) %>%
  summarise(
    N = n(),
    TLF = sum(event_tlf),
    TLF_rate = round(TLF / N * 100, 1),
    .groups = "drop"
  )

cat("\n=== Orsiro TLF Rate by Year (2014-2016) ===\n")
print(yearly_orsiro)

yearly_xience <- period2 %>%
  filter(stent_group == "XIENCE") %>%
  group_by(procedure_year) %>%
  summarise(
    N = n(),
    TLF = sum(event_tlf),
    TLF_rate = round(TLF / N * 100, 1),
    .groups = "drop"
  )

cat("\n=== XIENCE TLF Rate by Year (2014-2016) ===\n")
print(yearly_xience)

# Early event (30일 이내) vs Late event 분석
period2 <- period2 %>%
  mutate(
    early_event = ifelse(event_tlf == 1 & time_to_tlf <= 30, 1, 0),
    late_event = ifelse(event_tlf == 1 & time_to_tlf > 30, 1, 0),
    very_late_event = ifelse(event_tlf == 1 & time_to_tlf > 365, 1, 0)
  )

cat("\n=== Event Timing Analysis ===\n")
event_timing <- period2 %>%
  group_by(stent_group) %>%
  summarise(
    N = n(),
    Total_TLF = sum(event_tlf),
    Early_30d = sum(early_event),
    Early_pct = round(Early_30d / Total_TLF * 100, 1),
    Late_30d_1yr = sum(late_event) - sum(very_late_event),
    Late_pct = round((Late_30d_1yr) / Total_TLF * 100, 1),
    Very_late_1yr = sum(very_late_event),
    Very_late_pct = round(Very_late_1yr / Total_TLF * 100, 1),
    .groups = "drop"
  )
print(event_timing)


################################################################################
# 8. TLF Component Analysis
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 7: TLF COMPONENT ANALYSIS\n")
cat(rep("=", 80), "\n")

# TLF 구성요소 분석
component_analysis <- period2 %>%
  group_by(stent_group) %>%
  summarise(
    N = n(),
    TLF = sum(event_tlf),
    TLF_rate = round(TLF / N * 100, 1),
    
    # Components (if available)
    CD = sum(cardiac_death, na.rm = TRUE),
    CD_rate = round(CD / N * 100, 1),
    
    MI = sum(MI, na.rm = TRUE),
    MI_rate = round(MI / N * 100, 1),
    
    TLR = sum(TLR, na.rm = TRUE),
    TLR_rate = round(TLR / N * 100, 1),
    
    ST = sum(any_ST, na.rm = TRUE),
    ST_rate = round(ST / N * 100, 1),
    
    .groups = "drop"
  )

cat("\n=== TLF Component Analysis (2014-2016) ===\n")
print(component_analysis)


################################################################################
# 9. Era 간 Orsiro 비교 (2014-2016 vs 2017-2021)
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 8: ORSIRO COMPARISON ACROSS ERAS\n")
cat(rep("=", 80), "\n")

# Orsiro만 추출
orsiro_all <- analysis_data %>%
  filter(stent_names == "Orsiro") %>%
  filter(era %in% c("2014-2016", "2017-2021"))

cat("\nOrsiro patients by era:\n")
cat("2014-2016:", sum(orsiro_all$era == "2014-2016"), "\n")
cat("2017-2021:", sum(orsiro_all$era == "2017-2021"), "\n")

# Era별 환자 특성 비교
cat("\n=== Orsiro Patient Characteristics by Era ===\n")
orsiro_compare <- orsiro_all %>%
  group_by(era) %>%
  summarise(
    N = n(),
    Age = round(mean(age, na.rm = TRUE), 1),
    Male_pct = round(mean(male_n, na.rm = TRUE) * 100, 1),
    DM_pct = round(mean(DM_n, na.rm = TRUE) * 100, 1),
    CKD_pct = round(mean(CKD_n, na.rm = TRUE) * 100, 1),
    ACS_pct = round(mean(ACS_n, na.rm = TRUE) * 100, 1),
    Multivessel_pct = round(mean(multivessel_n, na.rm = TRUE) * 100, 1),
    Bifurcation_pct = round(mean(bifurcation_n, na.rm = TRUE) * 100, 1),
    Complex_pct = round(mean(complex_n, na.rm = TRUE) * 100, 1),
    Avg_diameter = round(mean(avg_stent_diameter, na.rm = TRUE), 2),
    Total_length = round(mean(total_stent_length, na.rm = TRUE), 1),
    Small_vessel_pct = round(mean(avg_stent_diameter <= 2.75, na.rm = TRUE) * 100, 1),
    TLF_rate = round(sum(event_tlf) / n() * 100, 1),
    .groups = "drop"
  )
print(t(orsiro_compare))

# Cox regression: 2014-2016 vs 2017-2021
cat("\n=== Orsiro: 2014-2016 vs 2017-2021 Cox Regression ===\n")
orsiro_all$era_factor <- factor(orsiro_all$era, levels = c("2017-2021", "2014-2016"))

cox_orsiro_era <- coxph(Surv(time_to_tlf, event_tlf) ~ era_factor, data = orsiro_all)
print(summary(cox_orsiro_era))

# Adjusted
cox_orsiro_era_adj <- coxph(
  Surv(time_to_tlf, event_tlf) ~ era_factor + age + male_n + DM_n + CKD_n + 
    ACS_n + multivessel_n + total_stent_length,
  data = orsiro_all
)
cat("\n=== Adjusted Model ===\n")
print(summary(cox_orsiro_era_adj))


################################################################################
# 10. Kaplan-Meier Plot
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 9: KAPLAN-MEIER PLOTS\n")
cat(rep("=", 80), "\n")

# 2014-2016 Orsiro vs XIENCE KM plot
fit_2014_2016 <- survfit(Surv(time_to_tlf, event_tlf) ~ stent_group, data = period2)

p_km <- ggsurvplot(
  fit_2014_2016,
  data = period2,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlim = c(0, 1825),
  break.time.by = 365,
  xlab = "Days from Index PCI",
  ylab = "TLF-free Survival",
  title = "2014-2016: Orsiro vs XIENCE",
  legend.title = "Stent",
  legend.labs = c("Orsiro", "XIENCE"),
  palette = c("#D73027", "#4575B4"),
  risk.table.height = 0.25,
  ggtheme = theme_bw()
)

print(p_km)
ggsave("Figure_2014-2016_Orsiro_vs_XIENCE_KM.png", p_km$plot, 
       width = 10, height = 8, dpi = 300)

# Orsiro Era 비교 KM plot
fit_orsiro_era <- survfit(Surv(time_to_tlf, event_tlf) ~ era, data = orsiro_all)

p_km_era <- ggsurvplot(
  fit_orsiro_era,
  data = orsiro_all,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlim = c(0, 1825),
  break.time.by = 365,
  xlab = "Days from Index PCI",
  ylab = "TLF-free Survival",
  title = "Orsiro: 2014-2016 vs 2017-2021",
  legend.title = "Era",
  palette = c("#FDB863", "#5AB4AC"),
  risk.table.height = 0.25,
  ggtheme = theme_bw()
)

print(p_km_era)
ggsave("Figure_Orsiro_Era_Comparison_KM.png", p_km_era$plot, 
       width = 10, height = 8, dpi = 300)


################################################################################
# 11. Subgroup Forest Plot
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SECTION 10: SUBGROUP ANALYSIS FOREST PLOT\n")
cat(rep("=", 80), "\n")

# Subgroup 분석 함수
run_subgroup_cox <- function(data, subgroup_var, subgroup_name) {
  results <- data.frame()
  
  for (level in unique(data[[subgroup_var]])) {
    subset_data <- data[data[[subgroup_var]] == level, ]
    
    if (nrow(subset_data) >= 20 && 
        length(unique(subset_data$stent_group)) == 2 &&
        sum(subset_data$event_tlf) >= 5) {
      
      cox_fit <- tryCatch({
        coxph(Surv(time_to_tlf, event_tlf) ~ stent_group, data = subset_data)
      }, error = function(e) NULL)
      
      if (!is.null(cox_fit)) {
        res <- tidy(cox_fit, conf.int = TRUE, exponentiate = TRUE)
        
        n_orsiro <- sum(subset_data$stent_group == "Orsiro")
        n_xience <- sum(subset_data$stent_group == "XIENCE")
        events_orsiro <- sum(subset_data$stent_group == "Orsiro" & subset_data$event_tlf == 1)
        events_xience <- sum(subset_data$stent_group == "XIENCE" & subset_data$event_tlf == 1)
        
        results <- rbind(results, data.frame(
          Subgroup = subgroup_name,
          Level = as.character(level),
          N_Orsiro = n_orsiro,
          Events_Orsiro = events_orsiro,
          N_XIENCE = n_xience,
          Events_XIENCE = events_xience,
          HR = round(res$estimate, 2),
          CI_lower = round(res$conf.low, 2),
          CI_upper = round(res$conf.high, 2),
          P_value = round(res$p.value, 3)
        ))
      }
    }
  }
  return(results)
}

# 각 subgroup 분석
subgroup_results <- rbind(
  run_subgroup_cox(period2, "vessel_size", "Vessel Size"),
  run_subgroup_cox(period2 %>% mutate(ACS = ifelse(ACS_n == 1, "ACS", "Stable")), 
                   "ACS", "Clinical Presentation"),
  run_subgroup_cox(period2 %>% mutate(Complex = ifelse(complex_n == 1, "Complex", "Simple")), 
                   "Complex", "Lesion Complexity"),
  run_subgroup_cox(period2 %>% mutate(DM = ifelse(DM_n == 1, "DM", "No DM")), 
                   "DM", "Diabetes")
)

cat("\n=== Subgroup Analysis Results ===\n")
print(subgroup_results)

# Forest plot
if (nrow(subgroup_results) > 0) {
  subgroup_results <- subgroup_results %>%
    mutate(
      label = paste0(Level, " (", N_Orsiro, "/", N_XIENCE, ")"),
      HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, CI_lower, CI_upper)
    )
  
  p_forest_subgroup <- ggplot(subgroup_results, 
                              aes(x = HR, y = reorder(paste(Subgroup, Level, sep = ": "), -row_number()))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
    geom_point(size = 3, color = "#D73027") +
    scale_x_log10(breaks = c(0.5, 1, 2, 4, 8), limits = c(0.3, 10)) +
    labs(
      title = "2014-2016 Orsiro vs XIENCE: Subgroup Analysis",
      subtitle = "HR > 1 favors XIENCE",
      x = "Hazard Ratio (95% CI)",
      y = ""
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  print(p_forest_subgroup)
  ggsave("Figure_2014-2016_Orsiro_Subgroup_Forest.png", p_forest_subgroup,
         width = 10, height = 6, dpi = 300)
}

write.csv(subgroup_results, "Table_2014-2016_Orsiro_Subgroup_Analysis.csv", row.names = FALSE)


################################################################################
# 12. Summary
################################################################################

cat("\n", rep("=", 80), "\n")
cat("SUMMARY\n")
cat(rep("=", 80), "\n")

cat("\n=== KEY FINDINGS ===\n")
cat("\n1. 2014-2016 Orsiro TLF rate:", 
    round(sum(period2$stent_group == "Orsiro" & period2$event_tlf == 1) / 
            sum(period2$stent_group == "Orsiro") * 100, 1), "%")
cat("\n2. 2014-2016 XIENCE TLF rate:", 
    round(sum(period2$stent_group == "XIENCE" & period2$event_tlf == 1) / 
            sum(period2$stent_group == "XIENCE") * 100, 1), "%")
cat("\n3. 2017-2021 Orsiro TLF rate:", 
    round(sum(orsiro_all$era == "2017-2021" & orsiro_all$event_tlf == 1) / 
            sum(orsiro_all$era == "2017-2021") * 100, 1), "%\n")

cat("\n", rep("=", 80), "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n")