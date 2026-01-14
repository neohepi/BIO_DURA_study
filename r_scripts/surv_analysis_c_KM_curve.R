library(dplyr)
library(survival)
library(jskm)
library(cowplot)
library(ggplot2)

# ==============================================================================
# 유틸리티 함수
# ==============================================================================
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

# ==============================================================================
# Landmark Survival 분석 함수
# ==============================================================================
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
    stop(paste0("[Error] 다음 변수를 찾을 수 없습니다: ", paste(missing_vars, collapse = ", ")))
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
      pval.coord = c(150, yrange[2] * 0.9), #main = plot_title
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
      #main = plot_title
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

# ==============================================================================
# Landmark 없는 단순 Cumulative Incidence 분석 함수
# ==============================================================================
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
    stop(paste0("[Error] 다음 변수를 찾을 수 없습니다: ", paste(missing_vars, collapse = ", ")))
  }
  
  # 변수 추출
  start_date <- data[[date_cag_var]]
  evt_date   <- data[[date_event_var]]
  evt_status <- as.numeric(data[[event_status_var]])
  fu_days    <- as.numeric(data[[fu_days_var]])
  group      <- data[[group_var]]
  
  # Time 계산
  diff_event <- as.numeric(evt_date - start_date)
  raw_time <- ifelse(evt_status == 1, diff_event, fu_days)
  raw_time <- ifelse(is.na(raw_time) | raw_time <= 0, 0.1, raw_time)
  
  # Capping
  censor_time <- 365 * cap_years
  time_capped <- pmin(raw_time, censor_time)
  status_capped <- ifelse(raw_time > censor_time, 0, evt_status)
  
  # 분석용 데이터
  jskm_data <- data.frame(
    time = time_capped,
    status = status_capped,
    group = as.factor(group)
  )
  
  # Cox 모델
  cox_model <- tryCatch(
    coxph(Surv(time, status) ~ group, data = jskm_data), 
    error = function(e) NULL
  )
  
  # Survfit
  fit <- eval(bquote(survfit(Surv(time, status) ~ group, data = jskm_data)))
  
  # Plot (landmark 없음)
  plot_obj <- jskm(
    fit, data = jskm_data,
    cumhaz = TRUE, 
    table = TRUE, 
    mark = FALSE, 
    ystratalabs = group_labels, 
    ystrataname = "",
    ylab = "Cumulative incidence (%)", 
    surv.scale = "percent", 
    ylims = yrange, 
    timeby = 365, 
    size.label.nrisk = 11,
    showpercent = TRUE, 
    pval = TRUE, 
    pval.size = 4, 
    pval.testname = FALSE,
    pval.coord = c(150, yrange[2] * 0.9), 
    #main = plot_title
  )
  
  # Cox 결과 테이블
  result_table <- extract_cox_result(cox_model, paste0("Overall (0-", cap_years, "y)"))
  
  return(list(
    plot = plot_obj,
    cox_table = result_table,
    fit = fit,
    analyzed_time = time_capped,
    analyzed_status = status_capped
  ))
}

# ==============================================================================
# 배치 분석 함수
# ==============================================================================
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

# ==============================================================================
# 2x2 Grid Figure 생성 함수
# ==============================================================================
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

# ==============================================================================
# Outcome 설정 (Landmark 분석용)
# ==============================================================================
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
    yrange = c(0, 0.08),
    title = "Any death"
  )
)

# ==============================================================================
# 실행: Landmark 분석 (4개 outcome)
# ==============================================================================
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

# ==============================================================================
# 실행: ST 분석 (Landmark 없음)
# ==============================================================================
cat("\n=== Analyzing: ST (No Landmark) ===\n")

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
print(result_ST$plot)

# ST time 변수 추가
matched_data <- safe_add_column(matched_data, "time_to_ST", result_ST$analyzed_time)

# ==============================================================================
# 2x2 Grid Figure 생성 (Cardiac Death, TLR, TV-MI, Any death)
# ==============================================================================
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

print(grid_fig_4)

# 저장
ggsave("Figure_secondary_outcomes_2x2.png", grid_fig_4, width = 14, height = 12, dpi = 300)
ggsave("Figure_secondary_outcomes_2x2.tiff", grid_fig_4, width = 14, height = 12, dpi = 300, compression = "lzw")
ggsave("Figure_secondary_outcomes_2x2.pdf", grid_fig_4, width = 14, height = 12)

# ==============================================================================
# ST 단독 Figure 저장
# ==============================================================================
ggsave("Figure_ST.png", result_ST$plot, width = 7, height = 6, dpi = 300)
ggsave("Figure_ST.tiff", result_ST$plot, width = 7, height = 6, dpi = 300, compression = "lzw")
ggsave("Figure_ST.pdf", result_ST$plot, width = 7, height = 6)

# ==============================================================================
# (선택) 5개 모두 포함한 Grid (2x3 또는 3x2)
# ==============================================================================
plot_list_5 <- list(
  batch_results$results$cardiac_death$plot,
  batch_results$results$TLR$plot,
  batch_results$results$TVMI$plot,
  batch_results$results$any_death$plot,
  result_ST$plot
)

# 2x3 grid (빈 공간 하나)
grid_fig_5 <- create_grid_figure(
  plot_list = plot_list_5,
  labels = c("A", "B", "C", "D", "E"),
  ncol = 3, nrow = 2,
  label_size = 14, label_fontface = "bold"
)

print(grid_fig_5)

ggsave("Figure_all_outcomes_2x3.png", grid_fig_5, width = 18, height = 12, dpi = 300)