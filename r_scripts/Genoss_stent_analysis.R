################################################################################
# Genoss DES Event Rate Analysis
# 2014-2017 period에서 BP-DES > DP-DES 현상이 Genoss DES 때문인지 확인
################################################################################

library(dplyr)
library(survival)
library(ggplot2)
library(tableone)

################################################################################
# 1. Stent type 변수 확인 및 분류
################################################################################

cat("\n", rep("=", 70), "\n")
cat("GENOSS DES EVENT ANALYSIS\n")
cat(rep("=", 70), "\n\n")

# 데이터에서 stent 관련 변수 확인
cat("=== Available stent-related variables ===\n")
stent_vars <- grep("stent|STENT|Stent|DES|des", names(matched_data), value = TRUE)
print(stent_vars)

# Study period 정의 (이미 있다면 skip)
if (!"study_period" %in% names(matched_data)) {
  matched_data <- matched_data %>%
    mutate(
      year = as.numeric(format(CAG_date, "%Y")),
      study_period = case_when(
        year >= 2010 & year <= 2013 ~ "2010-2013",
        year >= 2014 & year <= 2017 ~ "2014-2017",
        year >= 2018 & year <= 2021 ~ "2018-2021",
        TRUE ~ NA_character_
      )
    )
}

################################################################################
# 2. BP-DES 환자 중 Stent type별 분포 확인
################################################################################

# BP-DES 환자만 추출
bp_des_data <- matched_data %>% filter(BP == 1)

cat("\n=== BP-DES Stent Distribution by Period ===\n")

# stent_name 변수가 있다고 가정 (실제 변수명에 맞게 수정 필요)
# 가능한 변수명: stent_name, stent_type, DES_type, stent_1, etc.

# 변수명 찾기 시도
possible_stent_vars <- c("stent_names", "stent_type", "DES_type", "stent_1", 
                         "Stent_name", "STENT_NAME", "stent")
stent_var <- NULL

for (v in possible_stent_vars) {
  if (v %in% names(bp_des_data)) {
    stent_var <- v
    break
  }
}

if (is.null(stent_var)) {
  cat("\n[WARNING] Stent name variable not found. Please specify the correct variable name.\n")
  cat("Available variables:\n")
  print(names(bp_des_data))
} else {
  cat("Using stent variable:", stent_var, "\n\n")
}

################################################################################
# 3. Genoss DES 식별 및 그룹 생성
################################################################################

# Genoss 스텐트 식별 (대소문자 무관)
bp_des_data <- bp_des_data %>%
  mutate(
    is_genoss = grepl("GENOSS|Genoss|genoss", get(stent_var), ignore.case = TRUE),
    stent_group = case_when(
      is_genoss ~ "Genoss DES",
      TRUE ~ "Other BP-DES"
    )
  )

cat("\n=== Genoss DES Distribution ===\n")
genoss_by_period <- bp_des_data %>%
  group_by(study_period, stent_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(study_period) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

print(genoss_by_period)

################################################################################
# 4. 2014-2017 기간 집중 분석
################################################################################

cat("\n", rep("=", 70), "\n")
cat("FOCUS: 2014-2017 PERIOD ANALYSIS\n")
cat(rep("=", 70), "\n")

period2_bp <- bp_des_data %>% filter(study_period == "2014-2017")

cat("\n=== Stent Distribution in 2014-2017 (BP-DES only) ===\n")
stent_dist <- period2_bp %>%
  group_by(!!sym(stent_var)) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

print(stent_dist)

################################################################################
# 5. TLF Event Rate by Stent Type (2014-2017)
################################################################################

cat("\n=== TLF Event Rate by Stent Type (2014-2017) ===\n")

# composite_event 변수 사용 (TLF)
event_by_stent <- period2_bp %>%
  group_by(!!sym(stent_var)) %>%
  summarise(
    n = n(),
    events = sum(composite_event, na.rm = TRUE),
    event_rate = round(events / n * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(event_rate))

print(event_by_stent)

# Genoss vs Others 비교
cat("\n=== Genoss vs Other BP-DES: TLF Comparison (2014-2017) ===\n")

genoss_comparison <- period2_bp %>%
  group_by(stent_group) %>%
  summarise(
    n = n(),
    TLF_events = sum(composite_event, na.rm = TRUE),
    TLF_rate = round(TLF_events / n * 100, 1),
    cardiac_death = sum(cardiac_death, na.rm = TRUE),
    cardiac_death_rate = round(cardiac_death / n * 100, 1),
    TVMI = sum(TVMI, na.rm = TRUE),
    TVMI_rate = round(TVMI / n * 100, 1),
    TLR = sum(TLR, na.rm = TRUE),
    TLR_rate = round(TLR / n * 100, 1),
    any_death = sum(any_death, na.rm = TRUE),
    any_death_rate = round(any_death / n * 100, 1),
    .groups = "drop"
  )

print(genoss_comparison)

################################################################################
# 6. Cox Regression: Genoss vs Other BP-DES (2014-2017)
################################################################################

cat("\n=== Cox Regression: Genoss vs Other BP-DES (2014-2017) ===\n")

# Time to event 계산
period2_bp <- period2_bp %>%
  mutate(
    time_to_tlf = as.numeric(pmin(
      ifelse(composite_event == 1, 
             as.numeric(composite_date - CAG_date), 
             fu_days),
      365 * 5
    )),
    status_tlf = ifelse(time_to_tlf >= 365 * 5 & composite_event == 1, 
                        ifelse(as.numeric(composite_date - CAG_date) <= 365*5, 1, 0),
                        composite_event)
  )

# Genoss를 reference로
period2_bp$stent_group <- factor(period2_bp$stent_group, 
                                 levels = c("Other BP-DES", "Genoss DES"))

cox_genoss <- coxph(Surv(time_to_tlf, composite_event) ~ stent_group, 
                    data = period2_bp)

cat("\nCox Model Summary:\n")
print(summary(cox_genoss))

# HR 추출
cox_summ <- summary(cox_genoss)
hr <- round(cox_summ$conf.int[1, "exp(coef)"], 2)
hr_lower <- round(cox_summ$conf.int[1, "lower .95"], 2)
hr_upper <- round(cox_summ$conf.int[1, "upper .95"], 2)
p_val <- round(cox_summ$coefficients[1, "Pr(>|z|)"], 3)

cat("\n=== Genoss DES vs Other BP-DES (2014-2017) ===\n")
cat(sprintf("HR: %.2f (95%% CI: %.2f - %.2f), P = %.3f\n", hr, hr_lower, hr_upper, p_val))

################################################################################
# 7. 전체 기간 분석: Genoss vs Other BP-DES
################################################################################

cat("\n", rep("=", 70), "\n")
cat("ALL PERIODS: GENOSS VS OTHER BP-DES\n")
cat(rep("=", 70), "\n")

# 전체 기간 Genoss 비교
overall_genoss <- bp_des_data %>%
  group_by(study_period, stent_group) %>%
  summarise(
    n = n(),
    TLF_events = sum(composite_event, na.rm = TRUE),
    TLF_rate = round(TLF_events / n * 100, 1),
    .groups = "drop"
  )

print(overall_genoss)

################################################################################
# 8. Kaplan-Meier Plot: Genoss vs Other BP-DES (2014-2017)
################################################################################

cat("\n=== Generating KM Plot ===\n")

# Survfit
fit_genoss <- survfit(Surv(time_to_tlf, composite_event) ~ stent_group, 
                      data = period2_bp)

# ggsurvplot 사용 (survminer 패키지)
library(survminer)

km_plot <- ggsurvplot(
  fit_genoss,
  data = period2_bp,
  fun = "event",  # Cumulative incidence
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "Stent Type",
  legend.labs = c("Other BP-DES", "Genoss DES"),
  xlab = "Days since PCI",
  ylab = "Cumulative Incidence of TLF",
  title = "TLF: Genoss DES vs Other BP-DES (2014-2017)",
  palette = c("#2E9FDF", "#E7B800"),
  xlim = c(0, 1825),
  break.time.by = 365,
  risk.table.height = 0.25,
  ggtheme = theme_bw()
)

print(km_plot)

ggsave("Figure_Genoss_vs_OtherBPDES_2014-2017.png", km_plot$plot, 
       width = 8, height = 6, dpi = 300)

################################################################################
# 9. 추가 분석: 개별 스텐트별 Event Rate (2014-2017)
################################################################################

cat("\n", rep("=", 70), "\n")
cat("INDIVIDUAL STENT EVENT RATES (2014-2017)\n")
cat(rep("=", 70), "\n")

# 스텐트별 상세 분석
individual_stent_analysis <- period2_bp %>%
  group_by(!!sym(stent_var)) %>%
  summarise(
    N = n(),
    TLF = sum(composite_event, na.rm = TRUE),
    TLF_pct = round(TLF / N * 100, 1),
    Cardiac_Death = sum(cardiac_death, na.rm = TRUE),
    CD_pct = round(Cardiac_Death / N * 100, 1),
    TVMI = sum(TVMI, na.rm = TRUE),
    TVMI_pct = round(TVMI / N * 100, 1),
    TLR = sum(TLR, na.rm = TRUE),
    TLR_pct = round(TLR / N * 100, 1),
    Any_Death = sum(any_death, na.rm = TRUE),
    AD_pct = round(Any_Death / N * 100, 1),
    ST = sum(ST_outcome, na.rm = TRUE),
    ST_pct = round(ST / N * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(TLF_pct))

cat("\nEvent Rates by Individual Stent (2014-2017):\n")
print(individual_stent_analysis)

# CSV 저장
write.csv(individual_stent_analysis, "Stent_Event_Rates_2014-2017.csv", row.names = FALSE)
cat("\nSaved: Stent_Event_Rates_2014-2017.csv\n")

################################################################################
# 10. Summary
################################################################################

cat("\n", rep("=", 70), "\n")
cat("SUMMARY\n")
cat(rep("=", 70), "\n\n")

cat("Key Findings:\n")
cat("1. Genoss DES usage in 2014-2017:", 
    sum(period2_bp$is_genoss), "patients (",
    round(sum(period2_bp$is_genoss) / nrow(period2_bp) * 100, 1), "%)\n")

genoss_tlf <- genoss_comparison %>% filter(stent_group == "Genoss DES")
other_tlf <- genoss_comparison %>% filter(stent_group == "Other BP-DES")

cat("2. TLF rate - Genoss DES:", genoss_tlf$TLF_rate, "%\n")
cat("3. TLF rate - Other BP-DES:", other_tlf$TLF_rate, "%\n")
cat("4. HR (Genoss vs Other):", sprintf("%.2f (%.2f-%.2f), P=%.3f\n", 
                                        hr, hr_lower, hr_upper, p_val))

if (hr > 1 && p_val < 0.1) {
  cat("\n[CONCLUSION] Genoss DES shows HIGHER TLF rate compared to other BP-DES in 2014-2017.\n")
  cat("This may partially explain the worse outcomes of BP-DES during this transitional period.\n")
} else if (hr < 1) {
  cat("\n[CONCLUSION] Genoss DES shows LOWER or SIMILAR TLF rate compared to other BP-DES.\n")
  cat("The worse BP-DES outcomes in 2014-2017 are likely NOT due to Genoss DES.\n")
} else {
  cat("\n[CONCLUSION] No significant difference between Genoss DES and other BP-DES.\n")
}

################################################################################
# 11. 추가: Matched cohort 전체에서 Period별 BP-DES vs DP-DES 비교 (Genoss 제외 시)
################################################################################

cat("\n", rep("=", 70), "\n")
cat("SENSITIVITY ANALYSIS: EXCLUDING GENOSS DES\n")
cat(rep("=", 70), "\n")

# Genoss 제외한 데이터
data_no_genoss <- matched_data %>%
  filter(!(BP == 1 & grepl("GENOSS|Genoss|genoss", get(stent_var), ignore.case = TRUE)))

cat("\nPatients after excluding Genoss DES:\n")
cat("Original:", nrow(matched_data), "\n")
cat("After exclusion:", nrow(data_no_genoss), "\n")
cat("Excluded (Genoss BP-DES):", nrow(matched_data) - nrow(data_no_genoss), "\n")

# 2014-2017 기간만
period2_no_genoss <- data_no_genoss %>% filter(study_period == "2014-2017")

cat("\n=== 2014-2017 Period: BP-DES vs DP-DES (Genoss excluded) ===\n")

# TLF 비교
tlf_comparison <- period2_no_genoss %>%
  group_by(BP) %>%
  summarise(
    n = n(),
    TLF = sum(composite_event, na.rm = TRUE),
    TLF_rate = round(TLF / n * 100, 1),
    .groups = "drop"
  ) %>%
  mutate(Group = ifelse(BP == 1, "BP-DES (excl. Genoss)", "DP-DES"))

print(tlf_comparison)

# Cox model
period2_no_genoss <- period2_no_genoss %>%
  mutate(
    time_to_tlf = pmin(
      ifelse(composite_event == 1, 
             as.numeric(composite_date - CAG_date), 
             fu_days_correct),
      365 * 5
    )
  )

cox_no_genoss <- coxph(Surv(time_to_tlf, composite_event) ~ BP, 
                       data = period2_no_genoss)

cat("\nCox Model (2014-2017, Genoss excluded):\n")
print(summary(cox_no_genoss))

cox_summ2 <- summary(cox_no_genoss)
hr2 <- round(cox_summ2$conf.int[1, "exp(coef)"], 2)
hr2_lower <- round(cox_summ2$conf.int[1, "lower .95"], 2)
hr2_upper <- round(cox_summ2$conf.int[1, "upper .95"], 2)
p_val2 <- round(cox_summ2$coefficients[1, "Pr(>|z|)"], 3)

cat("\n=== BP-DES vs DP-DES (2014-2017, Genoss excluded) ===\n")
cat(sprintf("HR: %.2f (95%% CI: %.2f - %.2f), P = %.3f\n", hr2, hr2_lower, hr2_upper, p_val2))

cat("\n[INTERPRETATION]\n")
cat("Original (with Genoss): HR 1.48 (0.99-2.22), P=0.056\n")
cat(sprintf("Sensitivity (without Genoss): HR %.2f (%.2f-%.2f), P=%.3f\n", 
            hr2, hr2_lower, hr2_upper, p_val2))

if (hr2 < 1.48) {
  cat("\n-> Excluding Genoss ATTENUATES the HR, suggesting Genoss DES contributed to worse BP-DES outcomes.\n")
} else {
  cat("\n-> Excluding Genoss does NOT substantially change the HR.\n")
}