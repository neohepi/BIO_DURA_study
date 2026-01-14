################################################################################
# Subgroup Analysis - TLF at 5 Years Only (Error Fixed)
################################################################################

library(survival)
library(dplyr)
library(aod)

################################################################################
# 1. Prepare Data with Subgroup Variables
################################################################################

matched_data <- matched_data %>%
  mutate(
    # Age groups
    age_group = if_else(age < 65, "<65 years", ">=65 years"),
    age_group = factor(age_group, levels = c("<65 years", ">=65 years")),
    
    # Sex
    sex_group = if_else(male == 1, "Male", "Female"),
    sex_group = factor(sex_group, levels = c("Male", "Female")),
    
    # Diabetes Mellitus
    dm_group = if_else(DM == 1 | DM == "Yes", "Yes", "No"),
    dm_group = factor(dm_group, levels = c("No", "Yes")),
    
    # Chronic Kidney Disease
    ckd_group = if_else(CKD == 1 | CKD == "Yes", "Yes", "No"),
    ckd_group = factor(ckd_group, levels = c("No", "Yes")),
    
    # Clinical Presentation
    presentation_group = if_else(is_ACS == 1, "ACS", "Stable CAD"),
    presentation_group = factor(presentation_group, levels = c("Stable CAD", "ACS")),
    
    # Number of Diseased Vessels - factor를 numeric으로 변환 후 처리
    num_CAOD_numeric = as.numeric(as.character(num_CAOD)),
    vessel_group = case_when(
      num_CAOD_numeric == 1 ~ "1 vessel",
      num_CAOD_numeric == 2 ~ "2 vessels",
      num_CAOD_numeric >= 3 ~ "3 vessels",
      TRUE ~ NA_character_
    ),
    vessel_group = factor(vessel_group, levels = c("1 vessel", "2 vessels", "3 vessels")),
    
    # Bifurcation Lesion
    bifurcation_group = if_else(is_bifurcation == 1 | is_bifurcation == "Yes", "Yes", "No"),
    bifurcation_group = factor(bifurcation_group, levels = c("No", "Yes")),
    
    # Study Period (era)
    era_group = factor(era, levels = c("2010-2013", "2014-2016", "2017-2021"))
  )

# Verify distributions
cat("\n=== Subgroup Distributions ===\n")
cat("\nAge:\n"); print(table(matched_data$age_group, matched_data$BP))
cat("\nSex:\n"); print(table(matched_data$sex_group, matched_data$BP))
cat("\nDiabetes Mellitus:\n"); print(table(matched_data$dm_group, matched_data$BP))
cat("\nChronic Kidney Disease:\n"); print(table(matched_data$ckd_group, matched_data$BP))
cat("\nClinical Presentation:\n"); print(table(matched_data$presentation_group, matched_data$BP))
cat("\nNumber of Diseased Vessels:\n"); print(table(matched_data$vessel_group, matched_data$BP, useNA = "ifany"))
cat("\nBifurcation Lesion:\n"); print(table(matched_data$bifurcation_group, matched_data$BP))
cat("\nStudy Period:\n"); print(table(matched_data$era_group, matched_data$BP))


################################################################################
# 2. Subgroup Analysis Function (Fixed)
################################################################################

perform_subgroup_analysis <- function(data, 
                                      subgroup_var, 
                                      subgroup_name,
                                      time_var = "time_to_composite",
                                      event_var = "composite_event") {
  
  results <- list()
  
  # NA 제거된 데이터
  data_clean <- data %>% filter(!is.na(!!sym(subgroup_var)))
  
  # Get subgroup levels (실제 존재하는 레벨만)
  subgroup_levels <- levels(droplevels(data_clean[[subgroup_var]]))
  
  for (level in subgroup_levels) {
    
    # Subset data
    data_subgroup <- data_clean %>% filter(!!sym(subgroup_var) == level)
    
    # Skip if too few observations or events
    n_events <- sum(data_subgroup[[event_var]], na.rm = TRUE)
    if (nrow(data_subgroup) < 20 || n_events < 5) {
      cat("Warning: Skipping", level, "- insufficient observations/events\n")
      next
    }
    
    # Cox model
    formula_text <- paste0("Surv(", time_var, ", ", event_var, ") ~ BP")
    
    cox_model <- tryCatch({
      coxph(as.formula(formula_text), data = data_subgroup, cluster = subclass)
    }, error = function(e) {
      cat("Warning: Cox model failed for", level, "\n")
      return(NULL)
    })
    
    if (is.null(cox_model)) next
    
    # Extract results
    hr <- exp(coef(cox_model))[1]
    ci <- exp(confint(cox_model))
    p <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]
    
    # Event counts
    events_bp <- sum(data_subgroup[[event_var]][data_subgroup$BP == 1], na.rm = TRUE)
    events_dp <- sum(data_subgroup[[event_var]][data_subgroup$BP == 0], na.rm = TRUE)
    n_bp <- sum(data_subgroup$BP == 1)
    n_dp <- sum(data_subgroup$BP == 0)
    
    results[[level]] <- data.frame(
      Subgroup = subgroup_name,
      Category = level,
      N_DP = n_dp,
      Events_DP = events_dp,
      N_BP = n_bp,
      Events_BP = events_bp,
      HR = hr,
      CI_lower = ci[1],
      CI_upper = ci[2],
      P_value = p,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results
  if (length(results) == 0) {
    return(NULL)
  }
  
  results_df <- do.call(rbind, results)
  rownames(results_df) <- NULL
  
  # Test for interaction (실제 존재하는 레벨이 2개 이상일 때만)
  p_interaction <- NA
  
  if (length(subgroup_levels) >= 2 && nrow(results_df) >= 2) {
    interaction_formula <- paste0("Surv(", time_var, ", ", event_var, ") ~ BP * ", subgroup_var)
    
    cox_interaction <- tryCatch({
      coxph(as.formula(interaction_formula), data = data_clean, cluster = subclass)
    }, error = function(e) NULL)
    
    if (!is.null(cox_interaction)) {
      interaction_term <- grep("^BP1:", names(coef(cox_interaction)), value = TRUE)
      
      if (length(interaction_term) > 0) {
        if (length(interaction_term) > 1) {
          # Multiple levels - Wald test
          interaction_positions <- grep("^BP1:", names(coef(cox_interaction)))
          wald <- tryCatch({
            wald.test(Sigma = vcov(cox_interaction), b = coef(cox_interaction), Terms = interaction_positions)
          }, error = function(e) NULL)
          
          if (!is.null(wald)) {
            p_interaction <- wald$result$chi2["P"]
          }
        } else {
          p_interaction <- tryCatch({
            summary(cox_interaction)$coefficients[interaction_term, "Pr(>|z|)"]
          }, error = function(e) NA)
        }
      }
    }
  }
  
  results_df$P_interaction <- p_interaction
  
  return(results_df)
}


################################################################################
# 3. Run Subgroup Analysis - TLF at 5 Years
################################################################################

cat("\n\n=== Running Subgroup Analysis: TLF at 5 Years ===\n")

subgroup_list <- list()

# Age
subgroup_list$age <- perform_subgroup_analysis(matched_data, "age_group", "Age")

# Sex  
subgroup_list$sex <- perform_subgroup_analysis(matched_data, "sex_group", "Sex")

# Diabetes Mellitus
subgroup_list$dm <- perform_subgroup_analysis(matched_data, "dm_group", "Diabetes Mellitus")

# Chronic Kidney Disease
subgroup_list$ckd <- perform_subgroup_analysis(matched_data, "ckd_group", "Chronic Kidney Disease")

# Clinical Presentation
subgroup_list$presentation <- perform_subgroup_analysis(matched_data, "presentation_group", "Clinical Presentation")

# Number of Diseased Vessels
subgroup_list$vessel <- perform_subgroup_analysis(matched_data, "vessel_group", "Number of Diseased Vessels")

# Bifurcation Lesion
subgroup_list$bifurcation <- perform_subgroup_analysis(matched_data, "bifurcation_group", "Bifurcation Lesion")

# Study Period
subgroup_list$era <- perform_subgroup_analysis(matched_data, "era_group", "Study Period")

# Combine all results
subgroup_results_tlf <- do.call(rbind, subgroup_list)
rownames(subgroup_results_tlf) <- NULL

# Add overall result
overall_cox <- coxph(
  Surv(time_to_composite, composite_event) ~ BP,
  data = matched_data,
  cluster = subclass
)

overall_hr <- exp(coef(overall_cox))[1]
overall_ci <- exp(confint(overall_cox))
overall_p <- summary(overall_cox)$coefficients[1, "Pr(>|z|)"]
overall_events_bp <- sum(matched_data$composite_event[matched_data$BP == 1], na.rm = TRUE)
overall_events_dp <- sum(matched_data$composite_event[matched_data$BP == 0], na.rm = TRUE)

overall_row <- data.frame(
  Subgroup = "Overall",
  Category = "All Patients",
  N_DP = sum(matched_data$BP == 0),
  Events_DP = overall_events_dp,
  N_BP = sum(matched_data$BP == 1),
  Events_BP = overall_events_bp,
  HR = overall_hr,
  CI_lower = overall_ci[1],
  CI_upper = overall_ci[2],
  P_value = overall_p,
  P_interaction = NA,
  stringsAsFactors = FALSE
)

subgroup_results_tlf <- rbind(overall_row, subgroup_results_tlf)
rownames(subgroup_results_tlf) <- NULL

cat("\n=== Subgroup Analysis Results: TLF at 5 Years ===\n")
print(subgroup_results_tlf)


################################################################################
# 4. Format Table
################################################################################

format_subgroup_table <- function(data) {
  data %>%
    mutate(
      `BP-DES` = sprintf("%d/%d (%.1f%%)", Events_BP, N_BP, Events_BP/N_BP*100),
      `DP-DES` = sprintf("%d/%d (%.1f%%)", Events_DP, N_DP, Events_DP/N_DP*100),
      `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR, CI_lower, CI_upper),
      `P Value` = sprintf("%.3f", P_value),
      `P for Interaction` = if_else(!is.na(P_interaction), sprintf("%.3f", P_interaction), "—")
    ) %>%
    select(Subgroup, Category, `BP-DES`, `DP-DES`, `HR (95% CI)`, `P Value`, `P for Interaction`)
}

table_subgroup_tlf <- format_subgroup_table(subgroup_results_tlf)

cat("\n=== Formatted Table: Subgroup Analysis - TLF at 5 Years ===\n")
print(table_subgroup_tlf)

# Export
write.csv(table_subgroup_tlf, "Table_Subgroup_TLF.csv", row.names = FALSE)
write.csv(subgroup_results_tlf, "Subgroup_Results_TLF_Raw.csv", row.names = FALSE)


################################################################################
# 5. Interaction Summary
################################################################################

cat("\n", rep("=", 70), "\n")
cat("INTERACTION TEST SUMMARY - TLF at 5 Years\n")
cat(rep("=", 70), "\n\n")

interaction_summary <- subgroup_results_tlf %>%
  filter(!is.na(P_interaction)) %>%
  group_by(Subgroup) %>%
  slice(1) %>%
  ungroup() %>%
  select(Subgroup, P_interaction) %>%
  mutate(
    Significant = if_else(P_interaction < 0.10, "Yes*", "No"),
    `P for Interaction` = sprintf("%.3f", P_interaction)
  ) %>%
  select(Subgroup, `P for Interaction`, Significant) %>%
  arrange(Subgroup)

print(interaction_summary)

cat("\n* P < 0.10 suggests potential effect modification\n")
cat(rep("=", 70), "\n")


################################################################################
# 6. Forest Plot 데이터 변환 함수
################################################################################

convert_to_forest_format <- function(subgroup_results) {
  
  forest_data <- data.frame()
  
  # Overall 추가
  overall_row <- subgroup_results %>% filter(Category == "All Patients")
  
  if (nrow(overall_row) > 0) {
    forest_data <- rbind(forest_data, data.frame(
      Subgroup = "Overall",
      Category = "All Patients",
      HR = overall_row$HR,
      CI_lower = overall_row$CI_lower,
      CI_upper = overall_row$CI_upper,
      P_value = overall_row$P_value,
      P_interaction = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  # 각 Subgroup 처리
  subgroups <- unique(subgroup_results$Subgroup[subgroup_results$Subgroup != "Overall"])
  
  for (sg in subgroups) {
    sg_data <- subgroup_results %>% filter(Subgroup == sg)
    p_int <- sg_data$P_interaction[1]
    
    # Header 행
    forest_data <- rbind(forest_data, data.frame(
      Subgroup = sg,
      Category = sg,
      HR = NA,
      CI_lower = NA,
      CI_upper = NA,
      P_value = NA,
      P_interaction = p_int,
      stringsAsFactors = FALSE
    ))
    
    # Category 행들
    for (i in 1:nrow(sg_data)) {
      forest_data <- rbind(forest_data, data.frame(
        Subgroup = sg,
        Category = sg_data$Category[i],
        HR = sg_data$HR[i],
        CI_lower = sg_data$CI_lower[i],
        CI_upper = sg_data$CI_upper[i],
        P_value = sg_data$P_value[i],
        P_interaction = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  rownames(forest_data) <- NULL
  return(forest_data)
}

# Forest plot 형식으로 변환
forest_data_tlf <- convert_to_forest_format(subgroup_results_tlf)

cat("\n=== Forest Plot Data ===\n")
print(forest_data_tlf)