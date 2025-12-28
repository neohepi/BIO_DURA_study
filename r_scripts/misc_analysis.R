# 년도별로 DP / BP 사용수를 비교
yearly_summary <- df.lesion.all %>%
  # 1. 연도 추출
  mutate(procedure_year = year(CAG_date)) %>%
  
  # 2. [핵심 수정] BP(0/1) 값을 명시적인 문자열 라벨로 변환
  # BP == 1 -> "biodegradable", BP == 0 -> "durable"
  mutate(stent_class = if_else(BP == 1, "biodegradable", "durable")) %>%
  
  # 3. 그룹화 및 카운트
  group_by(procedure_year, stent_class) %>%
  summarize(n = n(), .groups = 'drop') %>%
  
  # 4. Pivot (Wide format으로 변환)
  # 이제 stent_class에 있는 "biodegradable", "durable"이 컬럼명이 됩니다.
  pivot_wider(names_from = stent_class, 
              values_from = n, 
              values_fill = 0) %>%
  
  # 5. 합계 및 비율 계산
  mutate(total = biodegradable + durable,
         pct_biodegradable = round(biodegradable / total * 100, 1))

# 결과 확인
print(yearly_summary)


# 0. 전처리: 분석용 데이터셋(df.study) 생성
# 기존 df.lesion.all에서 필요한 파생 변수들을 먼저 만듭니다.
df.study <- df.lesion.all %>%
  mutate(
    # procedure_year = year(CAG_date),
    # BP: 1 -> Biodegradable, 0 -> Durable
    stent_type = if_else(BP == 1, "Biodegradable", "Durable"),
    # 혹시 모를 NA나 빈 문자열 처리를 위해 stent_name 정리 (선택사항)
    stent_name = ifelse(is.na(stent_name) | stent_name == "", "Unknown", stent_name)
  )

# 1. 연도별 구체적인 스텐트 모델 분포 (Table)
stent_by_year <- df.study %>%
  group_by(procedure_year, stent_name) %>% 
  summarize(n = n(), .groups = 'drop') %>%
  arrange(procedure_year, desc(n))

print(stent_by_year, n = 100)

# 2. BP stent별 연도 분포 (Pivot Table)
bp_stents <- df.study %>%
  filter(stent_type == "Biodegradable") %>%
  group_by(stent_name, procedure_year) %>%
  summarize(n = n(), .groups = 'drop') %>%
  pivot_wider(names_from = procedure_year, 
              values_from = n, 
              values_fill = 0) %>%
  # (Tip) 총합(Total) 컬럼을 추가하면 보기 좋습니다
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(Total))

print(bp_stents)

# 3. DP stent별 연도 분포 (Pivot Table)
dp_stents <- df.study %>%
  filter(stent_type == "Durable") %>%
  group_by(stent_name, procedure_year) %>%
  summarize(n = n(), .groups = 'drop') %>%
  pivot_wider(names_from = procedure_year,
              values_from = n,
              values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(Total))

print(dp_stents)

# 4. 시각화 (ggplot)
# 스텐트 종류가 많을 경우 색상이 너무 많아질 수 있으므로 범례 위치 등을 조정했습니다.
ggplot(df.study, aes(x = factor(procedure_year), fill = stent_name)) +
  geom_bar(position = "stack", color = "white", linewidth = 0.2) + # 바 테두리 추가로 구분감 향상
  facet_wrap(~stent_type, scales = "free_y") + # y축 스케일을 풀어주어(DP가 압도적으로 많을 때 유용) 비교 용이
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",       # 범례를 아래로 이동
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90") # Facet 제목 배경색
  ) +
  labs(title = "Distribution of Stent Models by Year",
       subtitle = "Comparison between Biodegradable (BP) and Durable (DP) Polymer Stents",
       x = "Procedure Year", 
       y = "Number of Lesions (Patients)",
       fill = "Stent Model")



################################################################################
library(tidyverse)

# 1. Create study cohort with 3-era split
# (stent_type 변수가 없어도 BP 변수를 이용해 생성하도록 안전장치 추가)
df.study <- df.lesion.all %>%
  filter(procedure_year >= 2010) %>%
  mutate(
    # BP(1/0)를 이용해 명시적인 stent_type 생성 (이미 있다면 덮어씌움)
    stent_type = if_else(BP == 1, "biodegradable", "durable"),
    
    # Era 구분
    era = case_when(
      procedure_year <= 2013 ~ "2010-2013",
      procedure_year <= 2016 ~ "2014-2016",
      # 2017년 이후는 모두 포함 (혹시 2022년 데이터가 섞여도 에러 안 나게 TRUE 사용)
      TRUE ~ "2017-2021" 
    ),
    # Factor Level 순서 지정 (그래프 그릴 때 순서 꼬임 방지)
    era = factor(era, levels = c("2010-2013", "2014-2016", "2017-2021"))
  )

# 2. Verify era distribution (with Margins)
cat("\n=== Era Distribution (Count) ===\n")
# addmargins를 쓰면 행/열 합계를 같이 보여주어 데이터 누락 확인에 좋습니다.
table(df.study$era, df.study$stent_type) %>% 
  addmargins() %>% 
  print()

# 3. Era별 BP percentage Summary
era_summary <- df.study %>%
  group_by(era) %>%
  summarize(
    n_total = n(),
    n_bp = sum(stent_type == "biodegradable"),
    n_dp = sum(stent_type == "durable"),
    # BP 사용 비율 계산
    bp_pct = round(n_bp / n_total * 100, 1),
    .groups = 'drop'
  )

print(era_summary)


################################################################################
library(ggplot2)
library(dplyr)

plot_grouped_trend <- function(data, 
                               era_var,        
                               group_var,      
                               event_var,      
                               title = "Event Rate Trend by Era and Stent Type",
                               y_lab = "Event Rate (%)",
                               group_labels = c("DP-DES", "BP-DES")) { # [New] 범례 이름 지정 옵션
  
  # [Step 0] 변수명 유효성 검사 (Safety Check)
  # 입력받은 변수명이 실제 데이터에 있는지 먼저 체크합니다.
  required_vars <- c(era_var, group_var, event_var)
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(paste0("\n[Error] 다음 변수명을 데이터에서 찾을 수 없습니다:\n -> ", 
                paste(missing_vars, collapse = ", "), 
                "\n 철자를 확인하거나 colnames(data)를 체크해보세요."))
  }
  
  # [Step 1] 데이터 요약 (Summary)
  summary_data <- data %>%
    group_by(Era = .data[[era_var]], Group = .data[[group_var]]) %>%
    summarise(
      N = n(),
      Events = sum(.data[[event_var]], na.rm = TRUE),
      Rate = mean(.data[[event_var]], na.rm = TRUE),
      SE = sqrt((Rate * (1 - Rate)) / N),
      # 95% Confidence Interval
      Lower = (Rate - 1.96 * SE) * 100,
      Upper = (Rate + 1.96 * SE) * 100,
      Rate_Pct = Rate * 100,
      .groups = 'drop'
    )
  
  # 95% CI 하한선이 0보다 작아지지 않게 보정
  summary_data$Lower <- ifelse(summary_data$Lower < 0, 0, summary_data$Lower)
  
  # [Step 2] Group 변수 라벨링 (시각화를 위해 Factor 변환)
  # 만약 0, 1로 되어있다면 group_labels 인자를 이용해 이름 부여
  if(is.numeric(summary_data$Group) || is.character(summary_data$Group)) {
    summary_data$Group <- factor(summary_data$Group, labels = group_labels)
  }
  
  # [Step 3] 그래프 그리기
  p <- ggplot(summary_data, aes(x = Era, y = Rate_Pct, group = Group, color = Group)) +
    
    # Trend Line
    geom_line(linewidth = 1.2) +
    
    # Points
    geom_point(size = 4, shape = 21, fill = "white", stroke = 2) +
    
    # Error Bars (95% CI)
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1, linewidth = 0.8, alpha = 0.7) +
    
    # Value Labels (텍스트 표시)
    geom_text(aes(label = sprintf("%.1f%%", Rate_Pct)), 
              vjust = -1.5, size = 4, show.legend = FALSE, fontface = "bold") +
    
    # Colors & Labels
    scale_color_manual(values = c("#E41A1C", "#377EB8")) + 
    
    # Titles & Axes
    labs(
      title = title,
      subtitle = "Error bars indicate 95% Confidence Intervals",
      y = y_lab,
      x = "Era"
    ) +
    
    # Theme
    theme_classic() +
    theme(
      axis.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 13, face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40")
    )
  
  return(list(plot = p, data = summary_data))
}

# ==============================================================================
# 사용 예시
# ==============================================================================

# 1. Era 순서 지정 (필수)
# 데이터 내 era3 변수가 Factor가 아니거나 순서가 없으면 그래프 x축이 엉킵니다.
if("era" %in% names(result_composite$processed_data)) {
  result_composite$processed_data$era <- factor(
    result_composite$processed_data$era, 
    levels = c("2010-2013", "2014-2016", "2017-2021") # 실제 라벨에 맞춰 수정하세요
  )
}

# 2. 함수 실행
trend_res <- plot_grouped_trend(
  data = result_composite$processed_data, 
  era_var = "era",             
  group_var = "BP",             
  event_var = "analyzed_status",
  title = "Trend of Composite Outcome: BP-DES vs. DP-DES",
  group_labels = c("DP-DES", "BP-DES") # 0=DP, 1=BP 순서대로 이름 지정
)

# 3. 결과 출력
print(trend_res$plot)
print(trend_res$data)


################################################################################
# Publication-Ready Trend Plot Function
################################################################################

plot_publication_trend <- function(data, 
                                   era_var,        
                                   group_var,      
                                   event_var,      
                                   title = "Target Lesion Failure by Era",
                                   y_lab = "Event Rate (%)",
                                   x_lab = "Era",
                                   group_labels = c("Durable Polymer", "Biodegradable Polymer"),
                                   y_limits = c(0, 16),
                                   y_breaks = seq(0, 16, by = 2),
                                   show_n = TRUE,  # Show sample size
                                   p_interaction = NULL) {  # P for interaction
  
  library(ggplot2)
  library(dplyr)
  
  # Validate variables
  required_vars <- c(era_var, group_var, event_var)
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(paste0("Missing variables: ", paste(missing_vars, collapse = ", ")))
  }
  
  # Summarize data
  summary_data <- data %>%
    group_by(Era = .data[[era_var]], Group = .data[[group_var]]) %>%
    summarise(
      N = n(),
      Events = sum(.data[[event_var]], na.rm = TRUE),
      Rate = mean(.data[[event_var]], na.rm = TRUE),
      SE = sqrt((Rate * (1 - Rate)) / N),
      Lower = pmax(0, (Rate - 1.96 * SE) * 100),  # Lower CI, minimum 0
      Upper = (Rate + 1.96 * SE) * 100,           # Upper CI
      Rate_Pct = Rate * 100,
      .groups = 'drop'
    )
  
  # Assign group labels
  summary_data$Group <- factor(summary_data$Group, labels = group_labels)
  
  # Create plot
  p <- ggplot(summary_data, aes(x = Era, y = Rate_Pct, 
                                group = Group, color = Group, shape = Group)) +
    
    # Lines
    geom_line(linewidth = 1.0, alpha = 0.8) +
    
    # Error bars
    geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                  width = 0.15, linewidth = 0.6, alpha = 0.6) +
    
    # Points
    geom_point(size = 3.5, fill = "white", stroke = 1.2) +
    
    # Value labels
    geom_text(aes(label = sprintf("%.1f", Rate_Pct)), 
              vjust = -1.2, size = 3, 
              show.legend = FALSE, 
              color = "black",
              fontface = "plain") +
    
    # Colors (grayscale-friendly + color)
    scale_color_manual(values = c("#D55E00", "#0072B2")) +  # Colorblind-safe
    scale_shape_manual(values = c(16, 17)) +  # Circle and triangle
    
    # Axes
    scale_y_continuous(limits = y_limits, 
                       breaks = y_breaks,
                       expand = c(0.02, 0.02)) +
    
    # Labels
    labs(
      title = title,
      y = y_lab,
      x = x_lab,
      color = NULL,
      shape = NULL
    ) +
    
    # Theme
    theme_classic(base_size = 11) +
    theme(
      # Title
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      
      # Axes
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      
      # Legend
      legend.position = "top",
      legend.text = element_text(size = 10),
      legend.key.width = unit(1.5, "cm"),
      legend.spacing.x = unit(0.3, "cm"),
      legend.margin = margin(t = 0, b = 5),
      
      # Panel
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      
      # Overall
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Add sample sizes if requested
  if (show_n) {
    n_labels <- summary_data %>%
      select(Era, Group, N) %>%
      pivot_wider(names_from = Group, values_from = N) %>%
      mutate(label = sprintf("%s: n=%d\n%s: n=%d", 
                             group_labels[1], .[[2]],
                             group_labels[2], .[[3]]))
    
    # Add N as caption or annotation
    p <- p + 
      annotate("text", 
               x = 1:nrow(n_labels), 
               y = -1.5,  # Below x-axis
               label = paste0("n=", summary_data$N[seq(1, nrow(summary_data), 2)], 
                              "/", summary_data$N[seq(2, nrow(summary_data), 2)]),
               size = 2.5, 
               color = "gray30")
  }
  
  # Add p-value for interaction if provided
  if (!is.null(p_interaction)) {
    p <- p + 
      labs(caption = sprintf("P for interaction = %.3f", p_interaction))
  }
  
  return(list(plot = p, data = summary_data))
}

################################################################################
# Usage Example
################################################################################

# Create the plot
trend_result <- plot_publication_trend(
  data = matched_data,
  era_var = "era",
  group_var = "BP",
  event_var = "composite_event",  # Your TLF variable
  title = "Composite Endpoint by Era",
  y_lab = "Cumulative Incidence (%)",
  x_lab = "Study Period",
  group_labels = c("Durable Polymer", "Biodegradable Polymer"),
  y_limits = c(0, 22),
  y_breaks = seq(0, 22, by = 2),
  show_n = TRUE,
  p_interaction = 0.234  # Replace with actual p-value from interaction test
)

# View plot
print(trend_result$plot)

# View data
print(trend_result$data)

# Save as high-resolution figure
ggsave("Figure_TLF_by_Era.tiff", 
       plot = trend_result$plot,
       width = 6, 
       height = 5, 
       dpi = 300,
       compression = "lzw")

# Also save as PDF (vector format)
ggsave("Figure_TLF_by_Era.pdf", 
       plot = trend_result$plot,
       width = 6, 
       height = 5)

# Save as EPS for some journals
ggsave("Figure_TLF_by_Era.eps", 
       plot = trend_result$plot,
       width = 6, 
       height = 5,
       device = cairo_ps)


################################################################################
# p-for-interaction


################################################################################
# Test for Era × Treatment Interaction
################################################################################

library(survival)
library(dplyr)

# # 1. Cox model WITH interaction (1-year landmark)
# cox_with_interaction <- coxph(
#   Surv(time_to_composite, composite_event) ~ BP * era,
#   data = matched_data %>% filter(time_to_composite <= 365 | composite_event == 0),
#   cluster = subclass  # Account for matching
# )


matched_data <- add_capped_variables(
  data = matched_data,
  time_var = "time_to_composite",
  status_var = "composite_event",
  years = 5
)

matched_data <- add_capped_variables(
  data = matched_data,
  time_var = "time_to_TLR",
  status_var = "TLR",
  years = 5
)

matched_data <- add_capped_variables(
  data = matched_data,
  time_var = "time_to_MI",
  status_var = "next_MI",
  years = 5
)


# 2. Cox 모델 실행 (Interaction Term 포함)
cox_with_interaction <- coxph(
  formula = Surv(time_to_event, status_event) ~ BP * era, 
  data = matched_data, 
  cluster = subclass
)

summary(cox_with_interaction)

# 2. Main Effect 모델 (대조군) - filtering 없이 data_for_cox 사용
cox_without_interaction <- coxph(
  Surv(time_to_event, status_event) ~ BP + era, # 상호작용 없음 (+)
  data = matched_data,
  cluster = subclass
)

summary(cox_without_interaction)

################################################################################
# Solution: Interaction Test with Robust Variance
################################################################################

library(survival)
library(lmtest)

# Method 1: Wald Test for Interaction Terms (Recommended)
################################################################################

# Refit WITH interaction
cox_with_interaction <- coxph(
  Surv(time_to_event, status_event) ~ BP * era, 
  data = matched_data, 
  cluster = subclass
)

# Test if BOTH interaction terms are jointly zero
# H0: BP:era2014-2016 = 0 AND BP:era2017-2021 = 0

#install.packages("aod", "lmtest")
library(aod)  # For wald.test

# Get coefficient positions for interaction terms
# Position 4: BP1:era2014-2016
# Position 5: BP1:era2017-2021

wald_test <- wald.test(
  Sigma = vcov(cox_with_interaction),  # Robust variance-covariance matrix
  b = coef(cox_with_interaction),       # Coefficients
  Terms = c(4, 5)                       # Test terms 4 and 5 jointly
)

print(wald_test)

# Extract p-value
p_interaction <- wald_test$result$chi2["P"]
cat(sprintf("\n=== P for interaction (Wald test): %.3f ===\n", p_interaction))


################################################################################
# Method 2: Manual Wald Test (if aod package not available)
################################################################################

# Extract coefficients and variance-covariance matrix
coefs <- coef(cox_with_interaction)
vcov_robust <- vcov(cox_with_interaction)  # Already robust

# Interaction terms are positions 4 and 5
interaction_coefs <- coefs[4:5]
interaction_vcov <- vcov_robust[4:5, 4:5]

# Wald statistic: β' * Σ^(-1) * β
wald_stat <- t(interaction_coefs) %*% solve(interaction_vcov) %*% interaction_coefs
wald_stat <- as.numeric(wald_stat)

# P-value (chi-square with 2 df)
p_interaction_manual <- pchisq(wald_stat, df = 2, lower.tail = FALSE)

cat(sprintf("\nWald statistic: %.3f\n", wald_stat))
cat(sprintf("P-value (manual): %.3f\n", p_interaction_manual))


################################################################################
# Method 3: Without Robust Variance (for LRT comparison)
################################################################################

# Refit models WITHOUT cluster (for LRT only)
cox_with_int_simple <- coxph(
  Surv(time_to_event, status_event) ~ BP * era, 
  data = matched_data
  # No cluster
)

cox_without_int_simple <- coxph(
  Surv(time_to_event, status_event) ~ BP + era, 
  data = matched_data
  # No cluster
)

# Now anova works
lrt <- anova(cox_without_int_simple, cox_with_int_simple)
print(lrt)

p_interaction_lrt <- lrt$`Pr(>|Chi|)`[2]
cat(sprintf("\nP for interaction (LRT, non-robust): %.3f\n", p_interaction_lrt))

# NOTE: Use robust Wald test (Method 1) for final reporting


################################################################################
# Summary of Results
################################################################################

cat("\n" , rep("=", 70), "\n")
cat("INTERACTION TEST RESULTS\n")
cat(rep("=", 70), "\n\n")

cat("Individual interaction terms:\n")
cat(sprintf("  BP × era2014-2016: HR=%.2f, P=%.3f\n", 
            exp(coefs[4]), 
            summary(cox_with_interaction)$coefficients[4, "Pr(>|z|)"]))
cat(sprintf("  BP × era2017-2021: HR=%.2f, P=%.3f\n", 
            exp(coefs[5]), 
            summary(cox_with_interaction)$coefficients[5, "Pr(>|z|)"]))

cat(sprintf("\nJoint test for interaction: P = %.3f\n", p_interaction))

if (p_interaction > 0.10) {
  cat("\nConclusion: No significant interaction between treatment and era.\n")
  cat("Treatment effect is consistent across study periods.\n")
} else if (p_interaction > 0.05) {
  cat("\nConclusion: Borderline interaction (P=0.05-0.10).\n")
  cat("Consider reporting era-specific effects.\n")
} else {
  cat("\nConclusion: Significant interaction detected.\n")
  cat("Treatment effect varies by era - report stratified results.\n")
}


################################################################################
# Era-Specific HRs (Regardless of Interaction)
################################################################################

cat("\n" , rep("=", 70), "\n")
cat("ERA-SPECIFIC TREATMENT EFFECTS\n")
cat(rep("=", 70), "\n\n")

era_specific_results <- data.frame()

for (era_level in levels(matched_data$era)) {
  
  data_era <- matched_data %>% filter(era == era_level)
  
  cox_era <- coxph(
    Surv(time_to_event, status_event) ~ BP,
    data = data_era,
    cluster = subclass
  )
  
  hr <- exp(coef(cox_era))
  ci <- exp(confint(cox_era))
  p <- summary(cox_era)$coefficients["BP1", "Pr(>|z|)"]
  
  # Get event counts
  events_bp <- sum(data_era$status_event[data_era$BP == 1])
  events_dp <- sum(data_era$status_event[data_era$BP == 0])
  n_bp <- sum(data_era$BP == 1)
  n_dp <- sum(data_era$BP == 0)
  
  era_specific_results <- rbind(
    era_specific_results,
    data.frame(
      Era = era_level,
      N_DP = n_dp,
      Events_DP = events_dp,
      Rate_DP = round(events_dp/n_dp*100, 1),
      N_BP = n_bp,
      Events_BP = events_bp,
      Rate_BP = round(events_bp/n_bp*100, 1),
      HR = round(hr, 2),
      CI_lower = round(ci[1], 2),
      CI_upper = round(ci[2], 2),
      P_value = round(p, 3)
    )
  )
}

print(era_specific_results)


################################################################################
# Create Final Publication Plot
################################################################################

trend_final <- plot_publication_trend(
  data = matched_data,  # 전체 데이터 사용 (filter 제거)
  era_var = "era",
  group_var = "BP",
  event_var = "status_event",  # This will calculate rates properly
  title = "Target Lesion Failure at 1 Year by Era",
  y_lab = "Event Rate (%)",
  x_lab = "Study Period",
  group_labels = c("Durable Polymer", "Biodegradable Polymer"),
  y_limits = c(0, 12),
  y_breaks = seq(0, 12, by = 2),
  show_n = TRUE,
  p_interaction = 0.833
)

print(trend_final$plot)

# Save
ggsave("Figure_TLF_by_Era_final.tiff", 
       plot = trend_final$plot,
       width = 6.5, 
       height = 5, 
       dpi = 300,
       compression = "lzw")


################################################################################
# Table for Manuscript
################################################################################

manuscript_table <- era_specific_results %>%
  mutate(
    `DP-DES` = sprintf("%d/%d (%.1f%%)", Events_DP, N_DP, Rate_DP),
    `BP-DES` = sprintf("%d/%d (%.1f%%)", Events_BP, N_BP, Rate_BP),
    `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR, CI_lower, CI_upper),
    `P-value` = sprintf("%.3f", P_value)
  ) %>%
  select(Era, `DP-DES`, `BP-DES`, `HR (95% CI)`, `P-value`)

print("\n=== Table: Era-Specific Analysis ===\n")
print(manuscript_table)

# Add footnote about interaction
cat(sprintf("\nP for interaction = %.3f\n", p_interaction))

# Export
write.csv(manuscript_table, 
          "Table_Era_Specific_Analysis.csv", 
          row.names = FALSE)


################################################################################
# Alternative: Logistic Regression (for binary outcome at 1 year)
################################################################################

# Create binary outcome at 1 year
matched_data_1y <- matched_data %>%
  mutate(
    event_1y = if_else(time_to_event <= 365 & composite_outcome == 1, 1, 0)
  )

# Logistic with interaction
logistic_with <- glm(
  event_1y ~ BP * era,
  data = matched_data_1y,
  family = binomial
)

# Logistic without interaction
logistic_without <- glm(
  event_1y ~ BP + era,
  data = matched_data_1y,
  family = binomial
)

# LRT
lrt_logistic <- anova(logistic_without, logistic_with, test = "LRT")
print(lrt_logistic)

p_interaction_logistic <- lrt_logistic$`Pr(>Chi)`[2]
cat(sprintf("\nP for interaction (logistic): %.3f\n", p_interaction_logistic))

################################################################################
# Era-specific Treatment Effects
################################################################################

# Calculate HR for each era
era_specific_hr <- data.frame()

for (e in levels(matched_data$era)) {
  
  data_era <- matched_data %>% 
    filter(era == e, time_to_event <= 365 | composite_outcome == 0)
  
  cox_era <- coxph(
    Surv(time_to_event, composite_outcome) ~ BP,
    data = data_era
  )
  
  hr <- exp(coef(cox_era))
  ci <- exp(confint(cox_era))
  p <- summary(cox_era)$coefficients[,"Pr(>|z|)"]
  
  era_specific_hr <- rbind(
    era_specific_hr,
    data.frame(
      Era = e,
      HR = round(hr, 2),
      Lower_CI = round(ci[1], 2),
      Upper_CI = round(ci[2], 2),
      P_value = round(p, 3)
    )
  )
}

print("\n=== Era-Specific Treatment Effects ===")
print(era_specific_hr)

################################################################################
# Summary Table with Interaction
################################################################################

summary_table <- trend_result$data %>%
  select(Era, Group, N, Events, Rate_Pct) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(N, Events, Rate_Pct)
  ) %>%
  left_join(era_specific_hr, by = "Era") %>%
  mutate(
    HR_95CI = sprintf("%.2f (%.2f-%.2f)", HR, Lower_CI, Upper_CI)
  ) %>%
  select(
    Era,
    N_DP = `N_Durable Polymer`,
    Events_DP = `Events_Durable Polymer`,
    Rate_DP = `Rate_Pct_Durable Polymer`,
    N_BP = `N_Biodegradable Polymer`,
    Events_BP = `Events_Biodegradable Polymer`,
    Rate_BP = `Rate_Pct_Biodegradable Polymer`,
    HR_95CI,
    P_value
  )

print("\n=== Summary by Era ===")
print(summary_table)

cat(sprintf("\n=== Overall P for interaction: %.3f ===\n", p_interaction))

################################################################################
# Create Publication Plot with P-value
################################################################################

trend_final <- plot_publication_trend(
  data = matched_data,
  era_var = "era",
  group_var = "BP",
  event_var = "composite_outcome",
  title = "Target Lesion Failure by Era",
  y_lab = "Cumulative Incidence (%)",
  x_lab = "Study Period",
  group_labels = c("Durable Polymer", "Biodegradable Polymer"),
  y_limits = c(0, 25),
  y_breaks = seq(0, 25, by = 5),
  show_n = TRUE,
  p_interaction = p_interaction  # Add calculated p-value
)

print(trend_final$plot)

# Save
ggsave("Figure_TLF_by_Era.tiff", 
       plot = trend_final$plot,
       width = 6, 
       height = 5, 
       dpi = 300,
       compression = "lzw")






