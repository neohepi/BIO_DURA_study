################################################################################
# ERA ANALYSIS FUNCTIONS
# Author: Your Name
# Purpose: Streamlined era-specific treatment effect analysis
################################################################################

library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(aod)

################################################################################
# Function 1: Create Era Variable
################################################################################

create_era_variable <- function(data, 
                                year_var = "procedure_year",
                                era_breaks = c(2010, 2013, 2016, 2021),
                                era_labels = c("2010-2013", "2014-2016", "2017-2021")) {
  #' Create era variable with custom breaks
  #' 
  #' @param data Data frame
  #' @param year_var Name of year variable
  #' @param era_breaks Numeric vector of break points
  #' @param era_labels Character vector of era labels
  #' @return Data frame with era variable added
  
  data <- data %>%
    mutate(
      era = cut(.data[[year_var]], 
                breaks = c(-Inf, era_breaks[-c(1, length(era_breaks))], Inf),
                labels = era_labels,
                include.lowest = TRUE),
      era = factor(era, levels = era_labels)
    )
  
  # Summary
  cat("\n=== Era Distribution ===\n")
  print(table(data$era, useNA = "ifany"))
  
  return(data)
}


################################################################################
# Function 2: Era Distribution Summary
################################################################################

summarize_era_distribution <- function(data, 
                                       era_var = "era",
                                       group_var = "BP") {
  #' Summarize patient distribution by era and treatment group
  #' 
  #' @param data Data frame
  #' @param era_var Name of era variable
  #' @param group_var Name of treatment group variable (0/1)
  #' @return Summary table
  
  summary_table <- data %>%
    group_by(Era = .data[[era_var]], Group = .data[[group_var]]) %>%
    summarise(N = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Group, 
                values_from = N, 
                names_prefix = "Group_",
                values_fill = 0) %>%
    mutate(Total = Group_0 + Group_1,
           Pct_Group1 = round(Group_1 / Total * 100, 1))
  
  cat("\n=== Era × Treatment Distribution ===\n")
  print(summary_table)
  
  return(summary_table)
}


################################################################################
# Function 3: Test Era × Treatment Interaction
################################################################################

test_era_interaction <- function(data,
                                 time_var,
                                 status_var,
                                 group_var = "BP",
                                 era_var = "era",
                                 cluster_var = "subclass") {
  #' Test for interaction between treatment and era
  #' 
  #' @param data Data frame
  #' @param time_var Name of time-to-event variable
  #' @param status_var Name of event status variable
  #' @param group_var Name of treatment group variable
  #' @param era_var Name of era variable
  #' @param cluster_var Name of cluster variable for robust SE (optional)
  #' @return List with test results
  
  # Build formula
  formula_with <- as.formula(
    paste0("Surv(", time_var, ", ", status_var, ") ~ ", 
           group_var, " * ", era_var)
  )
  
  formula_without <- as.formula(
    paste0("Surv(", time_var, ", ", status_var, ") ~ ", 
           group_var, " + ", era_var)
  )
  
  # Fit models
  if (!is.null(cluster_var) && cluster_var %in% names(data)) {
    cox_with <- coxph(formula_with, data = data, cluster = .data[[cluster_var]])
    cox_without <- coxph(formula_without, data = data, cluster = .data[[cluster_var]])
  } else {
    cox_with <- coxph(formula_with, data = data)
    cox_without <- coxph(formula_without, data = data)
  }
  
  # Wald test for interaction terms
  coef_names <- names(coef(cox_with))
  interaction_indices <- grep(":", coef_names)
  
  if (length(interaction_indices) > 0) {
    wald_result <- wald.test(
      Sigma = vcov(cox_with),
      b = coef(cox_with),
      Terms = interaction_indices
    )
    p_interaction <- wald_result$result$chi2["P"]
  } else {
    p_interaction <- NA
  }
  
  # Print results
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("INTERACTION TEST RESULTS\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  cat("Model with interaction:\n")
  print(summary(cox_with)$coefficients)
  
  cat(sprintf("\nJoint Wald test for interaction: P = %.4f\n", p_interaction))
  
  if (p_interaction > 0.10) {
    cat("\n✓ No significant interaction.\n")
    cat("  → Treatment effect is consistent across eras.\n")
  } else if (p_interaction > 0.05) {
    cat("\n⚠ Borderline interaction (0.05 < P < 0.10).\n")
    cat("  → Consider reporting era-specific effects.\n")
  } else {
    cat("\n✓✓ Significant interaction detected!\n")
    cat("  → Treatment effect varies by era.\n")
    cat("  → Report stratified results.\n")
  }
  
  return(list(
    cox_with_interaction = cox_with,
    cox_without_interaction = cox_without,
    p_interaction = p_interaction,
    wald_test = if(exists("wald_result")) wald_result else NULL
  ))
}


################################################################################
# Function 4: Calculate Era-Specific HRs
################################################################################

calculate_era_specific_hr <- function(data,
                                      time_var,
                                      status_var,
                                      group_var = "BP",
                                      era_var = "era",
                                      cluster_var = NULL) {
  #' Calculate treatment effect (HR) within each era
  #' 
  #' @param data Data frame
  #' @param time_var Name of time-to-event variable
  #' @param status_var Name of event status variable
  #' @param group_var Name of treatment group variable
  #' @param era_var Name of era variable
  #' @param cluster_var Name of cluster variable (optional)
  #' @return Data frame with era-specific HRs
  
  era_levels <- levels(data[[era_var]])
  results <- data.frame()
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("ERA-SPECIFIC HAZARD RATIOS\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  for (era_level in era_levels) {
    
    data_era <- data %>% filter(.data[[era_var]] == era_level)
    
    # Cox model for this era
    formula_cox <- as.formula(
      paste0("Surv(", time_var, ", ", status_var, ") ~ ", group_var)
    )
    
    if (!is.null(cluster_var) && cluster_var %in% names(data)) {
      cox_era <- coxph(formula_cox, data = data_era, 
                       cluster = data_era[[cluster_var]])
    } else {
      cox_era <- coxph(formula_cox, data = data_era)
    }
    
    # Extract results
    hr <- exp(coef(cox_era)[1])
    ci <- exp(confint(cox_era)[1, ])
    p <- summary(cox_era)$coefficients[1, "Pr(>|z|)"]
    
    # Event counts
    group0_data <- data_era %>% filter(.data[[group_var]] == 0)
    group1_data <- data_era %>% filter(.data[[group_var]] == 1)
    
    n0 <- nrow(group0_data)
    n1 <- nrow(group1_data)
    events0 <- sum(group0_data[[status_var]], na.rm = TRUE)
    events1 <- sum(group1_data[[status_var]], na.rm = TRUE)
    
    results <- rbind(results, data.frame(
      Era = era_level,
      N_Group0 = n0,
      Events_Group0 = events0,
      Rate_Group0 = round(events0 / n0 * 100, 1),
      N_Group1 = n1,
      Events_Group1 = events1,
      Rate_Group1 = round(events1 / n1 * 100, 1),
      HR = round(hr, 2),
      CI_Lower = round(ci[1], 2),
      CI_Upper = round(ci[2], 2),
      P_value = round(p, 4)
    ))
  }
  
  print(results)
  
  return(results)
}


################################################################################
# Function 5: Create Publication-Ready Trend Plot
################################################################################

plot_era_trends <- function(data,
                            era_var = "era",
                            group_var = "BP",
                            status_var,
                            title = "Event Rate by Era",
                            y_label = "Event Rate (%)",
                            group_labels = c("Control", "Treatment"),
                            y_limits = NULL,
                            p_interaction = NULL) {
  #' Create publication-ready trend plot
  #' 
  #' @param data Data frame
  #' @param era_var Name of era variable
  #' @param group_var Name of group variable
  #' @param status_var Name of event status variable
  #' @param title Plot title
  #' @param y_label Y-axis label
  #' @param group_labels Labels for groups (length 2)
  #' @param y_limits Y-axis limits (optional)
  #' @param p_interaction P-value for interaction (optional)
  #' @return List with plot and data
  
  # Calculate rates
  plot_data <- data %>%
    group_by(Era = .data[[era_var]], Group = .data[[group_var]]) %>%
    summarise(
      N = n(),
      Events = sum(.data[[status_var]], na.rm = TRUE),
      Rate = mean(.data[[status_var]], na.rm = TRUE),
      SE = sqrt((Rate * (1 - Rate)) / N),
      Lower = pmax(0, (Rate - 1.96 * SE) * 100),
      Upper = (Rate + 1.96 * SE) * 100,
      Rate_Pct = Rate * 100,
      .groups = 'drop'
    ) %>%
    mutate(Group = factor(Group, labels = group_labels))
  
  # Auto y-limits if not provided
  if (is.null(y_limits)) {
    max_val <- max(plot_data$Upper, na.rm = TRUE)
    y_limits <- c(0, ceiling(max_val * 1.1))
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = Era, y = Rate_Pct, 
                             group = Group, color = Group, shape = Group)) +
    geom_line(linewidth = 1.0, alpha = 0.8) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                  width = 0.15, linewidth = 0.6, alpha = 0.6) +
    geom_point(size = 3.5, fill = "white", stroke = 1.2) +
    geom_text(aes(label = sprintf("%.1f%%", Rate_Pct)), 
              vjust = -1.2, size = 3.5, 
              show.legend = FALSE, 
              color = "black") +
    scale_color_manual(values = c("#D55E00", "#0072B2")) +
    scale_shape_manual(values = c(16, 17)) +
    scale_y_continuous(limits = y_limits, expand = c(0.02, 0.02)) +
    labs(title = title, y = y_label, x = "Era", 
         color = NULL, shape = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      legend.position = "top",
      legend.text = element_text(size = 10),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
    )
  
  # Add p-value if provided
  if (!is.null(p_interaction)) {
    p <- p + 
      labs(caption = sprintf("P for interaction = %.4f", p_interaction))
  }
  
  return(list(plot = p, data = plot_data))
}


################################################################################
# Function 6: Complete Era Analysis Pipeline
################################################################################

analyze_era_effects <- function(data,
                                time_var,
                                status_var,
                                group_var = "BP",
                                era_var = "era",
                                cluster_var = NULL,
                                outcome_name = "Event",
                                group_labels = c("Control", "Treatment"),
                                save_plots = FALSE,
                                output_prefix = "era_analysis") {
  #' Complete era analysis pipeline
  #' 
  #' @param data Data frame (must have era variable)
  #' @param time_var Name of time-to-event variable
  #' @param status_var Name of event status variable
  #' @param group_var Name of treatment group variable
  #' @param era_var Name of era variable
  #' @param cluster_var Name of cluster variable (optional)
  #' @param outcome_name Name of outcome for titles
  #' @param group_labels Labels for treatment groups
  #' @param save_plots Save plots to files?
  #' @param output_prefix Prefix for output files
  #' @return List with all results
  
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("ERA-SPECIFIC TREATMENT EFFECT ANALYSIS\n")
  cat(rep("=", 80), "\n", sep = "")
  cat(sprintf("Outcome: %s\n", outcome_name))
  cat(sprintf("Time variable: %s\n", time_var))
  cat(sprintf("Status variable: %s\n", status_var))
  cat(rep("=", 80), "\n\n", sep = "")
  
  # Step 1: Distribution summary
  cat("STEP 1: Era Distribution\n")
  cat(rep("-", 80), "\n", sep = "")
  dist_summary <- summarize_era_distribution(data, era_var, group_var)
  
  # Step 2: Interaction test
  cat("\n\nSTEP 2: Test for Era × Treatment Interaction\n")
  cat(rep("-", 80), "\n", sep = "")
  interaction_results <- test_era_interaction(
    data, time_var, status_var, group_var, era_var, cluster_var
  )
  
  # Step 3: Era-specific HRs
  cat("\n\nSTEP 3: Era-Specific Treatment Effects\n")
  cat(rep("-", 80), "\n", sep = "")
  era_hr <- calculate_era_specific_hr(
    data, time_var, status_var, group_var, era_var, cluster_var
  )
  
  # Step 4: Trend plot
  cat("\n\nSTEP 4: Creating Trend Plot\n")
  cat(rep("-", 80), "\n", sep = "")
  trend_plot <- plot_era_trends(
    data, era_var, group_var, status_var,
    title = paste(outcome_name, "by Era"),
    y_label = paste(outcome_name, "Rate (%)"),
    group_labels = group_labels,
    p_interaction = interaction_results$p_interaction
  )
  
  print(trend_plot$plot)
  
  # Save if requested
  if (save_plots) {
    filename_png <- paste0(output_prefix, "_trend.png")
    filename_pdf <- paste0(output_prefix, "_trend.pdf")
    
    ggsave(filename_png, plot = trend_plot$plot, 
           width = 7, height = 5, dpi = 300)
    ggsave(filename_pdf, plot = trend_plot$plot, 
           width = 7, height = 5)
    
    cat(sprintf("\n✓ Plots saved: %s, %s\n", filename_png, filename_pdf))
  }
  
  # Summary
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("SUMMARY\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  cat(sprintf("P for interaction: %.4f\n", interaction_results$p_interaction))
  cat("\nEra-specific results:\n")
  print(era_hr %>% 
          select(Era, HR, CI_Lower, CI_Upper, P_value) %>%
          mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, CI_Lower, CI_Upper)))
  
  # Return all results
  return(list(
    distribution = dist_summary,
    interaction_test = interaction_results,
    era_specific_hr = era_hr,
    trend_plot = trend_plot,
    p_interaction = interaction_results$p_interaction
  ))
}


################################################################################
# USAGE EXAMPLES
################################################################################

# Example 1: Simple usage
# result_composite_era <- analyze_era_effects(
#   data = matched_data,
#   time_var = "time_to_composite",
#   status_var = "composite_event",
#   group_var = "BP",
#   era_var = "era",
#   cluster_var = "subclass",
#   outcome_name = "Composite Endpoint",
#   group_labels = c("DP-DES", "BP-DES"),
#   save_plots = TRUE,
#   output_prefix = "composite"
# )

# Example 2: Multiple endpoints
# endpoints <- list(
#   list(status = "composite_event", time = "time_to_composite", name = "Composite"),
#   list(status = "TLR", time = "time_to_TLR", name = "TLR"),
#   list(status = "next_MI", time = "time_to_MI", name = "MI"),
#   list(status = "cardiac_death", time = "time_to_death", name = "Cardiac Death")
# )
# 
# era_results <- list()
# for (ep in endpoints) {
#   era_results[[ep$name]] <- analyze_era_effects(
#     data = matched_data,
#     time_var = ep$time,
#     status_var = ep$status,
#     outcome_name = ep$name,
#     group_labels = c("DP-DES", "BP-DES"),
#     save_plots = TRUE,
#     output_prefix = tolower(gsub(" ", "_", ep$name))
#   )
# }


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

# Define endpoints
endpoints <- tribble(
  ~name,            ~time_var,           ~status_var,
  "Composite",      "time_to_composite", "composite_event",
  "TLR",            "time_to_TLR",       "TLR",
  "MI",             "time_to_MI",        "next_MI",
  "Death",          "time_to_death",     "any_death",
  "Cardiac Death",  "time_to_death",     "cardiac_death"
)

# Loop through all
all_results <- list()

for (i in 1:nrow(endpoints)) {
  ep <- endpoints[i, ]
  
  cat("\n\n")
  cat(rep("#", 80), "\n", sep = "")
  cat(sprintf("ANALYZING: %s\n", ep$name))
  cat(rep("#", 80), "\n\n", sep = "")
  
  all_results[[ep$name]] <- analyze_era_effects(
    data = matched_data,
    time_var = ep$time_var,
    status_var = ep$status_var,
    outcome_name = ep$name,
    group_labels = c("DP-DES", "BP-DES"),
    save_plots = TRUE,
    output_prefix = tolower(gsub(" ", "_", ep$name))
  )
}

cat("\n✓ Era analysis functions loaded successfully!\n\n")
cat("Main function: analyze_era_effects()\n")
cat("Run ?analyze_era_effects for help\n\n")

# Extract p-values
p_interactions <- sapply(all_results, function(x) x$p_interaction)

# Summary table
summary_table <- data.frame(
  Endpoint = names(all_results),
  P_interaction = round(p_interactions, 4),
  Significant = ifelse(p_interactions < 0.05, "Yes", 
                       ifelse(p_interactions < 0.10, "Borderline", "No"))
)

print(summary_table)

# Export
write.csv(summary_table, "era_interaction_summary.csv", row.names = FALSE)