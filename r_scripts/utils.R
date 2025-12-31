library(dplyr)
library(survival)
library(jskm)

# ==============================================================================
# 함수 정의: analyze_landmark_survival
# ==============================================================================
analyze_landmark_survival <- function(data, 
                                      group_var, 
                                      date_cag_var, 
                                      date_event_var, 
                                      event_status_var, 
                                      fu_days_var, 
                                      cap_years = 5, 
                                      landmark_days = 365, 
                                      yrange = c(0, 0.20)) {
  
  # [Step 0] 변수명 유효성 검사 (안전장치)
  # 입력받은 변수명이 실제 데이터에 있는지 먼저 체크합니다.
  required_vars <- c(group_var, date_cag_var, date_event_var, event_status_var, fu_days_var)
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(paste0("\n[Error] 다음 변수명을 데이터에서 찾을 수 없습니다:\n -> ", 
                paste(missing_vars, collapse = ", "), 
                "\n 철자를 확인하거나 colnames(matched_data)를 체크해보세요."))
  }
  
  # [Step 1] 데이터 추출 (Base R 방식)
  start_date <- data[[date_cag_var]]
  evt_date   <- data[[date_event_var]]
  evt_status <- data[[event_status_var]]
  fu_days    <- as.numeric(data[[fu_days_var]])
  group      <- data[[group_var]]
  
  # [Step 2] 날짜 차이 계산 & Time 변수 생성
  # evt_date가 NA인 경우(이벤트 안 생긴 사람) 결과도 NA가 됨 -> 정상
  diff_event <- as.numeric(evt_date - start_date)
  
  # 이벤트(1)이면 날짜 차이, 아니면(0) fu_days 사용
  # diff_event가 NA여도 evt_status가 0이면 fu_days를 쓰므로 문제 없음
  raw_time <- ifelse(evt_status == 1, diff_event, fu_days)
  
  # NA 처리 (혹시 모를 결측치 방지) 및 0 이하 보정
  raw_time[is.na(raw_time)] <- 0.1 
  raw_time[raw_time <= 0] <- 0.1
  
  # [Step 3] Capping (5년)
  CENSOR_TIME <- 365 * cap_years
  
  time_capped <- ifelse(raw_time > CENSOR_TIME, CENSOR_TIME, raw_time)
  status_capped <- ifelse(raw_time > CENSOR_TIME, 0, evt_status)
  
  # [Step 4] 분석용 데이터프레임 생성
  working_data <- data.frame(
    time_capped = time_capped,
    status_capped = status_capped,
    group = group
  )
  
  # [Step 5] Landmark 데이터셋 분리
  # (A) 0 ~ Landmark
  data_0_LM <- working_data
  data_0_LM$time_lm <- ifelse(working_data$time_capped > landmark_days, landmark_days, working_data$time_capped)
  data_0_LM$status_lm <- ifelse(working_data$time_capped > landmark_days, 0, working_data$status_capped)
  
  # (B) Landmark ~ End (1년 생존자만)
  data_LM_End <- working_data[working_data$time_capped > landmark_days, ]
  
  # [Step 6] Survival Fit & Plot
  fit_overall <- survfit(Surv(time_capped, status_capped) ~ group, data = working_data)
  cox_overall <- coxph(Surv(time_capped, status_capped) ~ group, data = working_data)
  cox_0_LM    <- coxph(Surv(time_lm, status_lm) ~ group, data = data_0_LM)
  cox_LM_End  <- coxph(Surv(time_capped, status_capped) ~ group, data = data_LM_End)
  
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
  )
  
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
    fit = fit_overall
  ))
}

# ==============================================================================
# 함수 정의: count_unique_id
# ==============================================================================
count_unique_id <- function(x) {
  return(length(unique(x)))
}

library(dplyr)
library(rlang) # 동적 변수 생성을 위해 필요
# ==============================================================================
# 함수 정의: add_capped_variables
# ==============================================================================
add_capped_variables <- function(data, 
                                 time_var,       # 원본 시간 변수명
                                 status_var,     # 원본 상태 변수명
                                 years = 5)      # Capping할 연도
{
  # 1. Capping 기준일 계산
  cap_days <- 365 * years
  
  # 2. 변수 생성
  data <- data %>%
    mutate(
      # Time Capping
      time_to_event = pmin(.data[[time_var]], cap_days),
      
      # Status Capping [핵심 수정]
      # .data[[status_var]]가 Factor일 경우를 대비해 
      # as.numeric(as.character(...))를 씌워 숫자로 통일합니다.
      status_event = if_else(
        .data[[time_var]] > cap_days, 
        0, 
        as.numeric(as.character(.data[[status_var]]))
      )
    )
  
  return(data)
}


add_external_column <- function(target_data,   # 원본 데이터 (matched_data)
                                source_data,   # 가져올 데이터 (result_composite$processed_data)
                                source_col,    # 가져올 컬럼명 ("analyzed_time")
                                new_col_name) {# 새로 붙일 이름 ("time_to_composite")
  
  # 1. 행 개수 일치 여부 확인 (안전장치)
  if(nrow(target_data) != nrow(source_data)) {
    stop("Error: 두 데이터프레임의 행(Row) 개수가 달라서 합칠 수 없습니다.")
  }
  
  # 2. 컬럼 추출 및 이름 변경
  # !!sym()을 사용하여 문자열로 된 변수명을 실제 변수처럼 인식시킵니다.
  col_extracted <- source_data %>%
    select(all_of(source_col)) %>%
    rename(!!new_col_name := !!sym(source_col))
  
  # 3. 데이터 합치기 (cbind 대신 안전한 bind_cols 추천)
  # cbind(target_data, col_extracted) 도 가능하지만, dplyr 파이프라인엔 bind_cols가 좋습니다.
  combined_data <- bind_cols(target_data, col_extracted)
  
  return(combined_data)
}
