################################################################################
# Comprehensive Subgroup Analysis
################################################################################

library(survival)
library(dplyr)
library(forestplot)
library(ggplot2)

################################################################################
# 1. Prepare Data with Subgroup Variables
################################################################################

# Create subgroup variables if not already exist
matched_data <- matched_data %>%
  mutate(
    # Age groups
    age_group = if_else(age < 65, "<65 years", "≥65 years"),
    age_group = factor(age_group, levels = c("<65 years", "≥65 years")),
    
    # Sex (ensure proper coding)
    sex_group = if_else(male == 1, "Male", "Female"),
    sex_group = factor(sex_group, levels = c("Male", "Female")),
    
    # Diabetes
    dm_group = if_else(DM == 1 | DM == "Yes", "Diabetes", "No Diabetes"),
    dm_group = factor(dm_group, levels = c("No Diabetes", "Diabetes")),
    
    # ACS
    acs_group = if_else(is_ACS == 1, "ACS", "Stable CAD"),
    acs_group = factor(acs_group, levels = c("Stable CAD", "ACS"))
  )

# Verify distributions
cat("\n=== Subgroup Distributions ===\n")
table(matched_data$age_group, matched_data$BP)
table(matched_data$sex_group, matched_data$BP)
table(matched_data$dm_group, matched_data$BP)
table(matched_data$acs_group, matched_data$BP)


################################################################################
# 2. Subgroup Analysis Function
################################################################################
################################################################################
# Fixed Subgroup Analysis Function
################################################################################

perform_subgroup_analysis <- function(data, 
                                      subgroup_var, 
                                      subgroup_name,
                                      time_var = "time_to_composite",
                                      event_var = "composite_event",
                                      endpoint_name = "TLF at 5 years") {
  
  results <- data.frame()
  
  # Get subgroup levels
  subgroup_levels <- levels(data[[subgroup_var]])
  
  for (level in subgroup_levels) {
    
    # Subset data
    data_subgroup <- data %>% filter(!!sym(subgroup_var) == level)
    
    # Create formula dynamically (FIXED)
    formula_text <- paste0("Surv(", time_var, ", ", event_var, ") ~ BP")
    
    # Cox model
    cox_model <- coxph(
      as.formula(formula_text),
      data = data_subgroup,
      cluster = subclass
    )
    
    # Extract results
    hr <- exp(coef(cox_model))
    ci <- exp(confint(cox_model))
    p <- summary(cox_model)$coefficients["BP1", "Pr(>|z|)"]
    
    # Event counts
    events_bp <- sum(data_subgroup[[event_var]][data_subgroup$BP == 1], na.rm = TRUE)
    events_dp <- sum(data_subgroup[[event_var]][data_subgroup$BP == 0], na.rm = TRUE)
    n_bp <- sum(data_subgroup$BP == 1)
    n_dp <- sum(data_subgroup$BP == 0)
    
    results <- rbind(results, data.frame(
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
    ))
  }
  
  # Test for interaction (FIXED)
  interaction_formula <- paste0("Surv(", time_var, ", ", event_var, ") ~ BP * ", subgroup_var)
  
  cox_interaction <- coxph(
    as.formula(interaction_formula),
    data = data,
    cluster = subclass
  )
  
  # Extract interaction p-value
  interaction_term <- grep("BP1:", names(coef(cox_interaction)), value = TRUE)
  
  if (length(interaction_term) > 0) {
    # If multiple levels, use Wald test
    if (length(interaction_term) > 1) {
      library(aod)
      interaction_positions <- grep("BP1:", names(coef(cox_interaction)))
      wald <- wald.test(
        Sigma = vcov(cox_interaction),
        b = coef(cox_interaction),
        Terms = interaction_positions
      )
      p_interaction <- wald$result$chi2["P"]
    } else {
      p_interaction <- summary(cox_interaction)$coefficients[interaction_term, "Pr(>|z|)"]
    }
  } else {
    p_interaction <- NA
  }
  
  # Add interaction p-value to results
  results$P_interaction <- p_interaction
  
  return(results)
}


################################################################################
# Now run the subgroup analyses
################################################################################

library(survival)
library(dplyr)
library(aod)

# 1. TLF at 5 years
subgroup_results_tlf <- bind_rows(
  perform_subgroup_analysis(matched_data, "age_group", "Age", 
                            endpoint_name = "TLF at 5 years"),
  perform_subgroup_analysis(matched_data, "sex_group", "Sex",
                            endpoint_name = "TLF at 5 years"),
  perform_subgroup_analysis(matched_data, "dm_group", "Diabetes",
                            endpoint_name = "TLF at 5 years"),
  perform_subgroup_analysis(matched_data, "acs_group", "Clinical Presentation",
                            endpoint_name = "TLF at 5 years")
)

# Add overall result
overall_cox <- coxph(
  Surv(time_to_composite, composite_event) ~ BP,
  data = matched_data,
  cluster = subclass
)

overall_hr <- exp(coef(overall_cox))
overall_ci <- exp(confint(overall_cox))
overall_p <- summary(overall_cox)$coefficients["BP1", "Pr(>|z|)"]
overall_events_bp <- sum(matched_data$composite_event[matched_data$BP == 1], na.rm = TRUE)
overall_events_dp <- sum(matched_data$composite_event[matched_data$BP == 0], na.rm = TRUE)

subgroup_results_tlf <- rbind(
  data.frame(
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
    P_interaction = NA
  ),
  subgroup_results_tlf
)

print("\n=== Subgroup Analysis: TLF at 5 Years ===\n")
print(subgroup_results_tlf)


################################################################################
# 2. MI at 1 year
################################################################################

# Create 1-year capped MI variables if not already done
if (!"time_mi_1yr" %in% names(matched_data)) {
  matched_data <- matched_data %>%
    mutate(
      time_mi_1yr = pmin(time_to_MI, 365),
      # [수정] next_MI를 문자를 거쳐 숫자로 변환해야 안전합니다.
      event_mi_1yr = if_else(
        time_to_MI > 365, 
        0, 
        as.numeric(as.character(next_MI))
      )
    )
}

subgroup_results_mi <- bind_rows(
  perform_subgroup_analysis(matched_data, "age_group", "Age",
                            time_var = "time_mi_1yr",
                            event_var = "event_mi_1yr",
                            endpoint_name = "MI at 1 year"),
  perform_subgroup_analysis(matched_data, "sex_group", "Sex",
                            time_var = "time_mi_1yr",
                            event_var = "event_mi_1yr",
                            endpoint_name = "MI at 1 year"),
  perform_subgroup_analysis(matched_data, "dm_group", "Diabetes",
                            time_var = "time_mi_1yr",
                            event_var = "event_mi_1yr",
                            endpoint_name = "MI at 1 year"),
  perform_subgroup_analysis(matched_data, "acs_group", "Clinical Presentation",
                            time_var = "time_mi_1yr",
                            event_var = "event_mi_1yr",
                            endpoint_name = "MI at 1 year")
)

# Add overall for MI
overall_cox_mi <- coxph(
  Surv(time_mi_1yr, event_mi_1yr) ~ BP,
  data = matched_data,
  cluster = subclass
)

overall_hr_mi <- exp(coef(overall_cox_mi))
overall_ci_mi <- exp(confint(overall_cox_mi))
overall_p_mi <- summary(overall_cox_mi)$coefficients["BP1", "Pr(>|z|)"]
overall_events_bp_mi <- sum(matched_data$event_mi_1yr[matched_data$BP == 1], na.rm = TRUE)
overall_events_dp_mi <- sum(matched_data$event_mi_1yr[matched_data$BP == 0], na.rm = TRUE)

subgroup_results_mi <- rbind(
  data.frame(
    Subgroup = "Overall",
    Category = "All Patients",
    N_DP = sum(matched_data$BP == 0),
    Events_DP = overall_events_dp_mi,
    N_BP = sum(matched_data$BP == 1),
    Events_BP = overall_events_bp_mi,
    HR = overall_hr_mi,
    CI_lower = overall_ci_mi[1],
    CI_upper = overall_ci_mi[2],
    P_value = overall_p_mi,
    P_interaction = NA
  ),
  subgroup_results_mi
)

print("\n=== Subgroup Analysis: MI at 1 Year ===\n")
print(subgroup_results_mi)


################################################################################
# 3. Format and Export Tables
################################################################################

format_subgroup_table <- function(data, endpoint_name) {
  
  data %>%
    mutate(
      `BP-DES` = sprintf("%d/%d (%.1f%%)", 
                         Events_BP, N_BP, Events_BP/N_BP*100),
      `DP-DES` = sprintf("%d/%d (%.1f%%)", 
                         Events_DP, N_DP, Events_DP/N_DP*100),
      `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR, CI_lower, CI_upper),
      `P Value` = sprintf("%.3f", P_value),
      `P for Interaction` = if_else(!is.na(P_interaction),
                                    sprintf("%.3f", P_interaction),
                                    "—")
    ) %>%
    select(Subgroup, Category, `BP-DES`, `DP-DES`, 
           `HR (95% CI)`, `P Value`, `P for Interaction`)
}

table_subgroup_tlf <- format_subgroup_table(subgroup_results_tlf, 
                                            "Target Lesion Failure at 5 Years")
table_subgroup_mi <- format_subgroup_table(subgroup_results_mi,
                                           "Myocardial Infarction at 1 Year")

print("\n=== Formatted Table: Subgroup Analysis - TLF ===")
print(table_subgroup_tlf)

print("\n=== Formatted Table: Subgroup Analysis - MI ===")
print(table_subgroup_mi)

# Export
write.csv(table_subgroup_tlf, 
          "Table_Subgroup_TLF.csv", 
          row.names = FALSE)

write.csv(table_subgroup_mi,
          "Table_Subgroup_MI.csv",
          row.names = FALSE)

# Export raw data too
write.csv(subgroup_results_tlf, 
          "Subgroup_Results_TLF_Raw.csv", 
          row.names = FALSE)

write.csv(subgroup_results_mi,
          "Subgroup_Results_MI_Raw.csv",
          row.names = FALSE)


################################################################################
# 4. Summary of Interactions
################################################################################

cat("\n", rep("=", 70), "\n")
cat("INTERACTION TEST SUMMARY\n")
cat(rep("=", 70), "\n\n")

cat("Target Lesion Failure at 5 Years:\n")
interaction_summary_tlf <- subgroup_results_tlf %>%
  filter(!is.na(P_interaction)) %>%
  group_by(Subgroup) %>%
  slice(1) %>%
  ungroup() %>%
  select(Subgroup, P_interaction) %>%
  mutate(
    Significant = if_else(P_interaction < 0.10, "Yes*", "No"),
    P_formatted = sprintf("%.3f", P_interaction)
  ) %>%
  select(Subgroup, `P for Interaction` = P_formatted, Significant)

print(interaction_summary_tlf)

cat("\n\nMyocardial Infarction at 1 Year:\n")
interaction_summary_mi <- subgroup_results_mi %>%
  filter(!is.na(P_interaction)) %>%
  group_by(Subgroup) %>%
  slice(1) %>%
  ungroup() %>%
  select(Subgroup, P_interaction) %>%
  mutate(
    Significant = if_else(P_interaction < 0.10, "Yes*", "No"),
    P_formatted = sprintf("%.3f", P_interaction)
  ) %>%
  select(Subgroup, `P for Interaction` = P_formatted, Significant)

print(interaction_summary_mi)

cat("\n* P < 0.10 suggests potential effect modification\n")
cat(rep("=", 70), "\n")

