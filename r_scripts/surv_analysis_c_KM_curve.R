library(dplyr)
library(survival)
library(jskm)

library(survival)
library(jskm)

# ==============================================================================
# 함수 정의: analyze_landmark_survival (processed_data 반환 기능 추가)
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
                                      plot_title = "") {
  
  # [Step 0] 변수명 유효성 검사
  required_vars <- c(group_var, date_cag_var, date_event_var, event_status_var, fu_days_var)
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(paste0("\n[Error] 다음 변수명을 데이터에서 찾을 수 없습니다:\n -> ", 
                paste(missing_vars, collapse = ", "), 
                "\n 철자를 확인하거나 colnames(data)를 체크해보세요."))
  }
  
  # [Step 1] 데이터 추출
  start_date <- data[[date_cag_var]]
  evt_date   <- data[[date_event_var]]
  evt_status <- data[[event_status_var]]
  fu_days    <- as.numeric(data[[fu_days_var]])
  group      <- data[[group_var]]
  
  # [Step 2] 날짜 차이 계산 & Time 변수 생성
  diff_event <- as.numeric(evt_date - start_date)
  raw_time <- ifelse(evt_status == 1, diff_event, fu_days)
  
  # NA 처리 및 0 이하 보정
  raw_time[is.na(raw_time)] <- 0.1 
  raw_time[raw_time <= 0] <- 0.1
  
  # [Step 3] Capping (5년)
  CENSOR_TIME <- 365 * cap_years
  
  time_capped <- ifelse(raw_time > CENSOR_TIME, CENSOR_TIME, raw_time)
  status_capped <- ifelse(raw_time > CENSOR_TIME, 0, evt_status)
  
  # [Step 4] 분석용 데이터프레임 생성
  # (A) 내부 플로팅용 (컬럼명 고정)
  working_data <- data.frame(
    time_capped = time_capped,
    status_capped = status_capped,
    group = group
  )
  
  # (B) [New] 외부 반환용 (원본 데이터 + 계산된 결과)
  # 원본 데이터에 분석된 status와 time을 붙여서 반환합니다.
  processed_data <- data %>%
    mutate(
      analyzed_time = time_capped,    # Capping된 시간
      analyzed_status = status_capped # Capping된 이벤트 여부 (0/1)
    )
  
  # [Step 5] Landmark 데이터셋 분리 (내부 분석용)
  # (A) 0 ~ Landmark
  data_0_LM <- working_data
  data_0_LM$time_lm <- ifelse(working_data$time_capped > landmark_days, landmark_days, working_data$time_capped)
  data_0_LM$status_lm <- ifelse(working_data$time_capped > landmark_days, 0, working_data$status_capped)
  
  # (B) Landmark ~ End
  data_LM_End <- working_data[working_data$time_capped > landmark_days, ]
  
  # [Step 6] Survival Fit & Plot
  fit_overall <- survfit(Surv(time_capped, status_capped) ~ group, data = working_data)
  
  # Cox Model
  cox_overall <- coxph(Surv(time_capped, status_capped) ~ group, data = working_data)
  cox_0_LM    <- coxph(Surv(time_lm, status_lm) ~ group, data = data_0_LM)
  cox_LM_End  <- coxph(Surv(time_capped, status_capped) ~ group, data = data_LM_End)
  
  # Plotting with jskm
  if (landmark_days > 0) {
    plot_obj <- jskm(fit_overall, 
                     data = working_data,
                     cumhaz = T, table = T, mark = F, 
                     ystratalabs = c("DP-DES", "BP-DES"), 
                     ylab = "Cumulative incidence (%)", 
                     surv.scale = "percent", 
                     ylims = yrange,
                     timeby = 365,
                     size.label.nrisk = 11,
                     cut.landmark = landmark_days,
                     showpercent = T,
                     pval = T,
                     pval.size = 4,
                     pval.testname = F,
                     pval.coord = c(150, yrange[2]*0.9),
                     main = plot_title
    )
  } else {
    plot_obj <- jskm(fit_overall, 
                     data = working_data,
                     cumhaz = T, table = T, mark = F, 
                     ystratalabs = c("DP-DES", "BP-DES"), 
                     ylab = "Cumulative incidence (%)", 
                     surv.scale = "percent", 
                     ylims = yrange,
                     timeby = 365,
                     size.label.nrisk = 11,
                     showpercent = T,
                     pval = T,
                     pval.size = 4,
                     pval.testname = F,
                     pval.coord = c(150, yrange[2]*0.9),
                     main = plot_title
    )
  }
  
  # [Step 7] 결과 테이블 정리
  extract_cox_res <- function(model, label) {
    if(is.null(model)) return(NULL)
    summ <- summary(model)
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
      stringsAsFactors = F
    ))
  }
  
  res_table <- rbind(
    extract_cox_res(cox_overall, paste0("Overall (0-", cap_years, "y)")),
    extract_cox_res(cox_0_LM,    paste0("Landmark (0-", landmark_days, "d)")),
    extract_cox_res(cox_LM_End,  paste0("Landmark (", landmark_days, "d-", cap_years, "y)"))
  )
  
  return(list(
    plot = plot_obj,
    cox_table = res_table,
    fit = fit_overall,
    processed_data = processed_data # [New] Era 분석을 위해 원본+결과 데이터 반환
  ))
}

################################################################################
# Composite Event
################################################################################
result_composite <- analyze_landmark_survival(
  data = matched_data,
  group_var = "BP",
  date_cag_var = "CAG_date",
  date_event_var = "composite_date",   # <--- 여기만 변경
  event_status_var = "composite_event",      # <--- 여기만 변경
  fu_days_var = "fu_days",
  cap_years = 5,
  landmark_days = 365,
  yrange = c(0, 0.10),
  plot_title = "TLF"
)

print(result_composite$cox_table)
print(result_composite$plot)

# ESC grid figure를 위한 code
outcome <- list()
outcome$composite_plot <- result_composite$plot

matched_data <- add_external_column(
  target_data = matched_data,
  source_data = result_composite$processed_data,
  source_col = "analyzed_time",       # 원래 이름
  new_col_name = "time_to_composite"  # 바꿀 이름 (함수 인자로 지정)
)

################################################################################
# Cardiac death
################################################################################
result_cardiac_death <- analyze_landmark_survival(
  data = matched_data,
  group_var = "BP",
  date_cag_var = "CAG_date",
  date_event_var = "cardiac_death_date",   # <--- 여기만 변경
  event_status_var = "cardiac_death",      # <--- 여기만 변경
  fu_days_var = "fu_days",
  cap_years = 5,
  landmark_days = 365,
  yrange = c(0, 0.04),
  plot_title = "Cardiac Death"
)

print(result_cardiac_death$cox_table)
print(result_cardiac_death$plot)

outcome$cardiac_death_plot <- result_cardiac_death$plot

################################################################################
# TLR
################################################################################
result_TLR <- analyze_landmark_survival(
  data = matched_data,
  group_var = "BP",
  date_cag_var = "CAG_date",
  date_event_var = "TLR_date",   # <--- 여기만 변경
  event_status_var = "TLR",      # <--- 여기만 변경
  fu_days_var = "fu_days",
  cap_years = 5,
  landmark_days = 365,
  yrange = c(0, 0.05),
  plot_title = "TLR"
)

print(result_TLR$cox_table)
print(result_TLR$plot)

outcome$TLR_plot <- result_TLR$plot

################################################################################
# MI
################################################################################
# factor 결과값은 제대로 분석이 안됨
# matched_data$next_MI <- as.numeric(as.character(matched_data$next_MI))
# # 
# result_next_MI <- analyze_landmark_survival(
#   data = matched_data,
#   group_var = "BP",
#   date_cag_var = "CAG_date",
#   date_event_var = "next_MI_date",   # <--- 여기만 변경
#   event_status_var = "next_MI",      # <--- 여기만 변경
#   fu_days_var = "fu_days",
#   cap_years = 5,
#   landmark_days = 365,
#   yrange = c(0, 0.08),
#   plot_title = "MI"
# )
# print(result_next_MI$plot)

result_TVMI <- analyze_landmark_survival(
  data = matched_data,
  group_var = "BP",
  date_cag_var = "CAG_date",
  date_event_var = "TVMI_date",   # <--- 여기만 변경
  event_status_var = "TVMI",      # <--- 여기만 변경
  fu_days_var = "fu_days",
  cap_years = 5,
  landmark_days = 365,
  yrange = c(0, 0.04),
  plot_title = "TVMI"
)


print(result_TVMI$cox_table)
print(result_TVMI$plot)

matched_data <- add_external_column(
  target_data = matched_data,
  source_data = result_composite$processed_data,
  source_col = "analyzed_time",       # 원래 이름
  new_col_name = "time_to_MI"  # 바꿀 이름 (함수 인자로 지정)
)

outcome$MI_plot <-result_TVMI$plot

################################ Early
#install.packages("survminer")
library(survminer)
fit <- survfit(
  Surv(time_to_MI, TVMI) ~ BP,
  data = matched_data
)
#matched_data <- matched_data %>% rename(time_to_composite=time_to_composite...58)


ggsurvplot(
  fit,
  data = matched_data,
  fun = "event",  # Cumulative incidence
  xlim = c(0, 365),  # 첫 1년만
  break.time.by = 30,
  risk.table = TRUE
)

################################################################################
# Any death
################################################################################
result_any_death <- analyze_landmark_survival(
  data = matched_data,
  group_var = "BP",
  date_cag_var = "CAG_date",
  date_event_var = "death_date",   # <--- 여기만 변경
  event_status_var = "any_death",      # <--- 여기만 변경
  fu_days_var = "fu_days",
  cap_years = 5,
  landmark_days = 365,
  yrange = c(0, 0.08),
  plot_title = "Any death"
)

print(result_any_death$cox_table)
print(result_any_death$plot)

matched_data <- add_external_column(
  target_data = matched_data,
  source_data = result_composite$processed_data,
  source_col = "analyzed_time",       # 원래 이름
  new_col_name = "time_to_death"  # 바꿀 이름 (함수 인자로 지정)
)



################################################################################
# TVR
################################################################################
result_TVR <- analyze_landmark_survival(
  data = matched_data,
  group_var = "BP",
  date_cag_var = "CAG_date",
  date_event_var = "TVR_date",   # <--- 여기만 변경
  event_status_var = "TVR",      # <--- 여기만 변경
  fu_days_var = "fu_days",
  cap_years = 5,
  landmark_days = 365,
  yrange = c(0, 0.08),
  plot_title = "TVR"
)

print(result_TVR$cox_table)
print(result_TVR$plot)


# ################################################################################
# # ST
# ################################################################################
# result_ST <- analyze_landmark_survival(
#   data = matched_data,
#   group_var = "BP",
#   date_cag_var = "CAG_date",
#   date_event_var = "ST_date",   # <--- 여기만 변경
#   event_status_var = "ST_outcome",      # <--- 여기만 변경
#   fu_days_var = "fu_days",
#   cap_years = 5,
#   landmark_days = 365,
#   yrange = c(0, 0.01),
#   #  plot_title = "Primary Composite"
# )
# 
# print(result_ST$cox_table)
# print(result_ST$plot)
# 
# matched_data %>% filter(BP == 0 & ST_outcome == 1) %>% count()
# matched_data %>% filter(BP == 1 & ST_outcome == 1) %>% count()
