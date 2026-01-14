################################################################################
# Combined Survival Analysis: Figures and Table 2
# Clinical Outcomes at 1 year and 5 years (PSM Cohort)
################################################################################

library(dplyr)
library(survival)
library(jskm)
library(cowplot)
library(ggplot2)
library(flextable)
library(officer)

################################################################################
# PART 1: UTILITY FUNCTIONS
################################################################################

safe_add_column <- function(data, col_name, values) {
  if (col_name %in% names(data)) {
    data <- data %>% select(-all_of(col_name))
  }
  data[[col_name]] <- values
  return(data)
}

extract_cox_result <- function(model, label) {
  if (is.null(model) || inherits(try(summary(model), silent = TRUE), "try-error")) {
    return(data.frame(Period = label, HR = NA, CI_95 = NA, P_value = NA, stringsAsFactors = FALSE))
  }
  
  summ <- summary(model)
  if (nrow(summ$coefficients) == 0) {
    return(data.frame(Period = label, HR = NA, CI_95 = NA, P_value = NA, stringsAsFactors = FALSE))
  }
  
  coef <- summ$coefficients
  conf <- summ$conf.int
  
  hr      <- round(conf[1, "exp(coef)"], 2)
  ci_low  <- round(conf[1, "lower .95"], 2)
  ci_high <- round(conf[1, "upper .95"], 2)
  pval    <- format.pval(coef[1, "Pr(>|z|)"], eps = 0.001, digits = 3)
  
  return(data.frame(
    Period = label,
    HR = hr,
    CI_95 = paste0(ci_low, " - ", ci_high),
    P_value = pval,
    stringsAsFactors = FALSE
  ))
}

################################################################################
# PART 2: SURVIVAL ANALYSIS FUNCTIONS (FOR FIGURES)
################################################################################

# Landmark Survival Analysis
analyze_landmark_survival <- function(data, 
                                      group_var, 
                                      date_cag_var, 
                                      date_event_var, 
                                      event_status_var, 
                                      fu_days_var, 
                                      cap_years = 5, 
                                      landmark_days = 365, 
                                      yrange = c(0, 0.20),
                                      plot_title = "",
                                      group_labels = c("DP-DES", "BP-DES")) {
  
  required_vars <- c(group_var, date_cag_var, date_event_var, event_status_var, fu_days_var)
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(paste0("[Error] Missing variables: ", paste(missing_vars, collapse = ", ")))
  }
  
  start_date <- data[[date_cag_var]]
  evt_date   <- data[[date_event_var]]
  evt_status <- as.numeric(data[[event_status_var]])
  fu_days    <- as.numeric(data[[fu_days_var]])
  group      <- data[[group_var]]
  
  diff_event <- as.numeric(evt_date - start_date)
  raw_time <- ifelse(evt_status == 1, diff_event, fu_days)
  raw_time <- ifelse(is.na(raw_time) | raw_time <= 0, 0.1, raw_time)
  
  censor_time <- 365 * cap_years
  time_capped <- pmin(raw_time, censor_time)
  status_capped <- ifelse(raw_time > censor_time, 0, evt_status)
  
  jskm_data <- data.frame(
    time = time_capped,
    status = status_capped,
    group = as.factor(group)
  )
  
  df_early <- jskm_data %>%
    mutate(
      time_lm = pmin(time, landmark_days),
      status_lm = ifelse(time > landmark_days, 0, status)
    )
  
  df_late <- jskm_data %>% filter(time > landmark_days)
  
  cox_overall <- tryCatch(coxph(Surv(time, status) ~ group, data = jskm_data), error = function(e) NULL)
  cox_early <- tryCatch(coxph(Surv(time_lm, status_lm) ~ group, data = df_early), error = function(e) NULL)
  cox_late <- tryCatch(coxph(Surv(time, status) ~ group, data = df_late), error = function(e) NULL)
  
  fit_overall <- eval(bquote(survfit(Surv(time, status) ~ group, data = jskm_data)))
  
  if (landmark_days > 0) {
    plot_obj <- jskm(
      fit_overall, data = jskm_data,
      cumhaz = TRUE, table = TRUE, mark = FALSE, 
      ystratalabs = group_labels, ystrataname = "",
      ylab = "Cumulative incidence (%)", surv.scale = "percent", 
      ylims = yrange, timeby = 365, size.label.nrisk = 11,
      cut.landmark = landmark_days, showpercent = TRUE,
      pval = TRUE, pval.size = 4, pval.testname = FALSE,
      pval.coord = c(150, yrange[2] * 0.9), main = plot_title
    )
  } else {
    plot_obj <- jskm(
      fit_overall, data = jskm_data,
      cumhaz = TRUE, table = TRUE, mark = FALSE, 
      ystratalabs = group_labels, ystrataname = "",
      ylab = "Cumulative incidence (%)", surv.scale = "percent", 
      ylims = yrange, timeby = 365, size.label.nrisk = 11,
      showpercent = TRUE, pval = TRUE, pval.size = 4, 
      pval.testname = FALSE, pval.coord = c(150, yrange[2] * 0.9), 
      main = plot_title
    )
  }
  
  result_table <- rbind(
    extract_cox_result(cox_overall, paste0("Overall (0-", cap_years, "y)")),
    extract_cox_result(cox_early, paste0("Early (0-", landmark_days, "d)")),
    extract_cox_result(cox_late, paste0("Late (", landmark_days, "d-", cap_years, "y)"))
  )
  
  return(list(
    plot = plot_obj,
    cox_table = result_table,
    fit = fit_overall,
    analyzed_time = time_capped,
    analyzed_status = status_capped
  ))
}

# Simple Cumulative Incidence Analysis (No Landmark)
analyze_simple_survival <- function(data, 
                                    group_var, 
                                    date_cag_var, 
                                    date_event_var, 
                                    event_status_var, 
                                    fu_days_var, 
                                    cap_years = 5, 
                                    yrange = c(0, 0.05),
                                    plot_title = "",
                                    group_labels = c("DP-DES", "BP-DES")) {
  
  required_vars <- c(group_var, date_cag_var, date_event_var, event_status_var, fu_days_var)
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(paste0("[Error] Missing variables: ", paste(missing_vars, collapse = ", ")))
  }
  
  start_date <- data[[date_cag_var]]
  evt_date   <- data[[date_event_var]]
  evt_status <- as.numeric(data[[event_status_var]])
  fu_days    <- as.numeric(data[[fu_days_var]])
  group      <- data[[group_var]]
  
  diff_event <- as.numeric(evt_date - start_date)
  raw_time <- ifelse(evt_status == 1, diff_event, fu_days)
  raw_time <- ifelse(is.na(raw_time) | raw_time <= 0, 0.1, raw_time)
  
  censor_time <- 365 * cap_years
  time_capped <- pmin(raw_time, censor_time)
  status_capped <- ifelse(raw_time > censor_time, 0, evt_status)
  
  jskm_data <- data.frame(
    time = time_capped,
    status = status_capped,
    group = as.factor(group)
  )
  
  cox_model <- tryCatch(
    coxph(Surv(time, status) ~ group, data = jskm_data), 
    error = function(e) NULL
  )
  
  fit <- eval(bquote(survfit(Surv(time, status) ~ group, data = jskm_data)))
  
  plot_obj <- jskm(
    fit, data = jskm_data,
    cumhaz = TRUE, table = TRUE, mark = FALSE, 
    ystratalabs = group_labels, ystrataname = "",
    ylab = "Cumulative incidence (%)", surv.scale = "percent", 
    ylims = yrange, timeby = 365, size.label.nrisk = 11,
    showpercent = TRUE, pval = TRUE, pval.size = 4, 
    pval.testname = FALSE, pval.coord = c(150, yrange[2] * 0.9), 
    main = plot_title
  )
  
  result_table <- extract_cox_result(cox_model, paste0("Overall (0-", cap_years, "y)"))
  
  return(list(
    plot = plot_obj,
    cox_table = result_table,
    fit = fit,
    analyzed_time = time_capped,
    analyzed_status = status_capped
  ))
}

# Batch Analysis Function
run_batch_analysis <- function(data, group_var = "BP", date_cag_var = "CAG_date",
                               fu_days_var = "fu_days", outcomes,
                               cap_years = 5, landmark_days = 365,
                               group_labels = c("DP-DES", "BP-DES")) {
  
  results <- list()
  output_data <- data
  
  for (outcome_name in names(outcomes)) {
    cfg <- outcomes[[outcome_name]]
    cat("\n=== Analyzing:", outcome_name, "===\n")
    
    result <- analyze_landmark_survival(
      data = data, group_var = group_var, date_cag_var = date_cag_var,
      date_event_var = cfg$date_var, event_status_var = cfg$status_var,
      fu_days_var = fu_days_var, cap_years = cap_years,
      landmark_days = landmark_days, yrange = cfg$yrange,
      plot_title = if (!is.null(cfg$title)) cfg$title else "",
      group_labels = group_labels
    )
    
    results[[outcome_name]] <- result
    time_col_name <- paste0("time_to_", outcome_name)
    output_data <- safe_add_column(output_data, time_col_name, result$analyzed_time)
    print(result$cox_table)
  }
  
  return(list(results = results, data = output_data))
}

# Grid Figure Function
create_grid_figure <- function(plot_list, 
                               labels = c("A", "B", "C", "D"),
                               ncol = 2, nrow = 2,
                               label_size = 14,
                               label_fontface = "bold") {
  
  labeled_plots <- mapply(function(p, lab) {
    p + 
      ggtitle(lab) +
      theme(
        plot.title = element_text(
          size = label_size, 
          face = label_fontface,
          hjust = 0,
          margin = margin(b = 5)
        )
      )
  }, plot_list, labels, SIMPLIFY = FALSE)
  
  grid_plot <- plot_grid(
    plotlist = labeled_plots,
    ncol = ncol,
    nrow = nrow,
    align = "hv",
    axis = "tblr"
  )
  
  return(grid_plot)
}

################################################################################
# PART 3: TABLE 2 FUNCTIONS
################################################################################

extract_outcome_results_table <- function(data, 
                                          group_var = "BP",
                                          date_cag_var = "CAG_date",
                                          date_event_var,
                                          event_status_var,
                                          fu_days_var = "fu_days_correct",
                                          outcome_name = "Outcome") {
  
  start_date <- data[[date_cag_var]]
  evt_date <- data[[date_event_var]]
  evt_status <- as.numeric(data[[event_status_var]])
  fu_days <- as.numeric(data[[fu_days_var]])
  group <- data[[group_var]]
  
  diff_event <- as.numeric(evt_date - start_date)
  raw_time <- ifelse(evt_status == 1, diff_event, fu_days)
  raw_time <- ifelse(is.na(raw_time) | raw_time <= 0, 0.1, raw_time)
  
  results_list <- list()
  
  for (cap_years in c(1, 5)) {
    censor_time <- 365 * cap_years
    time_capped <- pmin(raw_time, censor_time)
    status_capped <- ifelse(raw_time > censor_time, 0, evt_status)
    
    analysis_df <- data.frame(
      time = time_capped,
      status = status_capped,
      group = as.factor(group)
    )
    
    fit <- survfit(Surv(time, status) ~ group, data = analysis_df)
    
    cox_model <- tryCatch(
      coxph(Surv(time, status) ~ group, data = analysis_df),
      error = function(e) NULL
    )
    
    time_point <- censor_time
    summ <- summary(fit, times = time_point, extend = TRUE)
    
    # DP-DES (group=0)
    dp_surv <- summ$surv[1]
    dp_lower <- summ$lower[1]
    dp_upper <- summ$upper[1]
    dp_cum_inc <- (1 - dp_surv) * 100
    dp_ci_lower <- (1 - dp_upper) * 100
    dp_ci_upper <- (1 - dp_lower) * 100
    
    # BP-DES (group=1)
    bp_surv <- summ$surv[2]
    bp_lower <- summ$lower[2]
    bp_upper <- summ$upper[2]
    bp_cum_inc <- (1 - bp_surv) * 100
    bp_ci_lower <- (1 - bp_upper) * 100
    bp_ci_upper <- (1 - bp_lower) * 100
    
    if (!is.null(cox_model)) {
      cox_summ <- summary(cox_model)
      hr <- cox_summ$conf.int[1, "exp(coef)"]
      hr_lower <- cox_summ$conf.int[1, "lower .95"]
      hr_upper <- cox_summ$conf.int[1, "upper .95"]
      p_value <- cox_summ$coefficients[1, "Pr(>|z|)"]
    } else {
      hr <- hr_lower <- hr_upper <- p_value <- NA
    }
    
    results_list[[paste0("year_", cap_years)]] <- list(
      dp_estimate = sprintf("%.1f (%.1f–%.1f)", dp_cum_inc, dp_ci_lower, dp_ci_upper),
      bp_estimate = sprintf("%.1f (%.1f–%.1f)", bp_cum_inc, bp_ci_lower, bp_ci_upper),
      hr = sprintf("%.2f (%.2f–%.2f)", hr, hr_lower, hr_upper),
      p_value = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
    )
  }
  
  return(results_list)
}

generate_table2 <- function(matched_data, fu_days_var = "fu_days_correct") {
  
  outcomes <- list(
    "Target lesion failure" = list(
      date_var = "composite_date",
      status_var = "composite_event"
    ),
    "  Cardiac death" = list(
      date_var = "cardiac_death_date",
      status_var = "cardiac_death"
    ),
    "  Target vessel MI" = list(
      date_var = "TVMI_date",
      status_var = "TVMI"
    ),
    "  Clinically driven TLR" = list(
      date_var = "TLR_date",
      status_var = "TLR"
    ),
    "All-cause death" = list(
      date_var = "death_date",
      status_var = "any_death"
    ),
    "Stent thrombosis" = list(
      date_var = "ST_date",
      status_var = "ST_outcome"
    )
  )
  
  table2_data <- data.frame(
    Outcome = character(),
    DP_1yr = character(),
    BP_1yr = character(),
    HR_1yr = character(),
    P_1yr = character(),
    DP_5yr = character(),
    BP_5yr = character(),
    HR_5yr = character(),
    P_5yr = character(),
    stringsAsFactors = FALSE
  )
  
  for (outcome_name in names(outcomes)) {
    cfg <- outcomes[[outcome_name]]
    
    cat("Processing:", outcome_name, "\n")
    
    results <- extract_outcome_results_table(
      data = matched_data,
      group_var = "BP",
      date_cag_var = "CAG_date",
      date_event_var = cfg$date_var,
      event_status_var = cfg$status_var,
      fu_days_var = fu_days_var,
      outcome_name = outcome_name
    )
    
    row_data <- data.frame(
      Outcome = outcome_name,
      DP_1yr = results$year_1$dp_estimate,
      BP_1yr = results$year_1$bp_estimate,
      HR_1yr = results$year_1$hr,
      P_1yr = results$year_1$p_value,
      DP_5yr = results$year_5$dp_estimate,
      BP_5yr = results$year_5$bp_estimate,
      HR_5yr = results$year_5$hr,
      P_5yr = results$year_5$p_value,
      stringsAsFactors = FALSE
    )
    
    table2_data <- rbind(table2_data, row_data)
  }
  
  return(table2_data)
}

create_table2_flextable <- function(table2_data) {
  
  ft <- flextable(table2_data) %>%
    set_header_labels(
      Outcome = "Outcome",
      DP_1yr = "DP-DES",
      BP_1yr = "BP-DES",
      HR_1yr = "HR (95% CI)",
      P_1yr = "P-value",
      DP_5yr = "DP-DES",
      BP_5yr = "BP-DES",
      HR_5yr = "HR (95% CI)",
      P_5yr = "P-value"
    ) %>%
    add_header_row(
      values = c("", "1-Year", "1-Year", "1-Year", "1-Year", "5-Year", "5-Year", "5-Year", "5-Year"),
      top = TRUE
    ) %>%
    merge_at(i = 1, j = 2:5, part = "header") %>%
    merge_at(i = 1, j = 6:9, part = "header") %>%
    theme_booktabs() %>%
    autofit() %>%
    align(align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "body") %>%
    fontsize(size = 9, part = "all") %>%
    bold(part = "header") %>%
    border_inner_h(border = fp_border(color = "gray80", width = 0.5)) %>%
    border_inner_v(border = fp_border(color = "gray80", width = 0.5))
  
  return(ft)
}

################################################################################
# PART 4: EXECUTION
################################################################################

cat("\n")
cat(rep("=", 80), "\n")
cat("SURVIVAL ANALYSIS: FIGURES AND TABLE 2\n")
cat(rep("=", 80), "\n\n")

# ------------------------------------------------------------------------------
# 4.1 Outcome Configuration for Figures
# ------------------------------------------------------------------------------

outcomes_config <- list(
  cardiac_death = list(
    date_var = "cardiac_death_date",
    status_var = "cardiac_death",
    yrange = c(0, 0.04),
    title = "Cardiac Death"
  ),
  TLR = list(
    date_var = "TLR_date",
    status_var = "TLR",
    yrange = c(0, 0.05),
    title = "TLR"
  ),
  TVMI = list(
    date_var = "TVMI_date",
    status_var = "TVMI",
    yrange = c(0, 0.04),
    title = "TV-MI"
  ),
  any_death = list(
    date_var = "death_date",
    status_var = "any_death",
    yrange = c(0, 0.10),
    title = "All-cause Death"
  )
)

# ------------------------------------------------------------------------------
# 4.2 Run Landmark Analysis (4 outcomes)
# ------------------------------------------------------------------------------

cat("\n--- Running Landmark Survival Analysis ---\n")

batch_results <- run_batch_analysis(
  data = matched_data,
  group_var = "BP",
  date_cag_var = "CAG_date",
  fu_days_var = "fu_days",
  outcomes = outcomes_config,
  cap_years = 5,
  landmark_days = 365
)

matched_data <- batch_results$data

# ------------------------------------------------------------------------------
# 4.3 Run ST Analysis (No Landmark)
# ------------------------------------------------------------------------------

cat("\n=== Analyzing: Stent Thrombosis (No Landmark) ===\n")

result_ST <- analyze_simple_survival(
  data = matched_data,
  group_var = "BP",
  date_cag_var = "CAG_date",
  date_event_var = "ST_date",
  event_status_var = "ST_outcome",
  fu_days_var = "fu_days",
  cap_years = 5,
  yrange = c(0, 0.01),
  plot_title = "Stent Thrombosis",
  group_labels = c("DP-DES", "BP-DES")
)

print(result_ST$cox_table)
matched_data <- safe_add_column(matched_data, "time_to_ST", result_ST$analyzed_time)

# ------------------------------------------------------------------------------
# 4.4 Create and Save Figures
# ------------------------------------------------------------------------------

cat("\n--- Creating Figures ---\n")

# 2x2 Grid (Cardiac Death, TLR, TV-MI, All-cause Death)
plot_list_4 <- list(
  batch_results$results$cardiac_death$plot,
  batch_results$results$TLR$plot,
  batch_results$results$TVMI$plot,
  batch_results$results$any_death$plot
)

grid_fig_4 <- create_grid_figure(
  plot_list = plot_list_4,
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2,
  label_size = 16, label_fontface = "bold"
)

# Save 2x2 grid
ggsave("Figure_secondary_outcomes_2x2.png", grid_fig_4, width = 14, height = 12, dpi = 300)
ggsave("Figure_secondary_outcomes_2x2.tiff", grid_fig_4, width = 14, height = 12, dpi = 300, compression = "lzw")
ggsave("Figure_secondary_outcomes_2x2.pdf", grid_fig_4, width = 14, height = 12)
cat("Saved: Figure_secondary_outcomes_2x2 (png, tiff, pdf)\n")

# Save ST figure separately
ggsave("Figure_ST.png", result_ST$plot, width = 7, height = 6, dpi = 300)
ggsave("Figure_ST.tiff", result_ST$plot, width = 7, height = 6, dpi = 300, compression = "lzw")
ggsave("Figure_ST.pdf", result_ST$plot, width = 7, height = 6)
cat("Saved: Figure_ST (png, tiff, pdf)\n")

# Display figures
print(grid_fig_4)
print(result_ST$plot)

# ------------------------------------------------------------------------------
# 4.5 Generate Table 2
# ------------------------------------------------------------------------------

cat("\n", rep("=", 70), "\n")
cat("GENERATING TABLE 2\n")
cat(rep("=", 70), "\n\n")

table2_data <- generate_table2(matched_data, fu_days_var = "fu_days_correct")

cat("\n=== Table 2. Clinical outcomes at 1 year and 5 years ===\n\n")
print(table2_data, row.names = FALSE)

# Save CSV
write.csv(table2_data, "Table2_Clinical_Outcomes.csv", row.names = FALSE)
cat("\nSaved: Table2_Clinical_Outcomes.csv\n")

# Create and save flextable
ft <- create_table2_flextable(table2_data)
print(ft)

# Save Word document
doc <- read_docx() %>%
  body_add_par("Table 2. Clinical outcomes at 1 year and 5 years (propensity score-matched cohort)", 
               style = "heading 1") %>%
  body_add_flextable(ft) %>%
  body_add_par("") %>%
  body_add_par("Data are presented as Kaplan-Meier estimates (%) with 95% confidence intervals. Hazard ratios (BP-DES vs DP-DES) were calculated using Cox proportional hazards models. Target lesion failure was defined as the composite of cardiac death, target vessel myocardial infarction, or clinically driven target lesion revascularization.", 
               style = "Normal") %>%
  body_add_par("") %>%
  body_add_par("Abbreviations: BP-DES, biodegradable polymer drug-eluting stent; CI, confidence interval; DP-DES, durable polymer drug-eluting stent; HR, hazard ratio; MI, myocardial infarction; TLR, target lesion revascularization.", 
               style = "Normal")

print(doc, target = "Table2_Clinical_Outcomes.docx")
cat("Saved: Table2_Clinical_Outcomes.docx\n")

# Save HTML
save_as_html(ft, path = "Table2_Clinical_Outcomes.html")
cat("Saved: Table2_Clinical_Outcomes.html\n")

# ------------------------------------------------------------------------------
# 4.6 Summary Output
# ------------------------------------------------------------------------------

cat("\n", rep("=", 80), "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n\n")

cat("Generated Files:\n")
cat("  Figures:\n")
cat("    - Figure_secondary_outcomes_2x2.png/tiff/pdf\n")
cat("    - Figure_ST.png/tiff/pdf\n")
cat("  Tables:\n")
cat("    - Table2_Clinical_Outcomes.csv\n")
cat("    - Table2_Clinical_Outcomes.docx\n")
cat("    - Table2_Clinical_Outcomes.html\n")
cat("\n")