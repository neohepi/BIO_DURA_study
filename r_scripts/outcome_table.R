library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(stringr)

# --------- helper: KM에서 특정 timepoint의 누적발생(1-S)과 95%CI 문자열 만들기 ---------
km_cuminc_str <- function(df, time_col, status_col, group_col, t_eval, digits = 1) {
  f <- as.formula(paste0("Surv(", time_col, ",", status_col, ") ~ ", group_col))
  fit <- survfit(f, data = df)
  
  s <- summary(fit, times = t_eval)
  
  # strata 이름에서 그룹값 추출 (예: "BP=0", "BP=1")
  grp <- str_replace(s$strata, ".*=", "")
  
  # 누적발생 = 1 - surv
  est <- 1 - s$surv
  lcl <- 1 - s$upper   # survival upper -> incidence lower
  ucl <- 1 - s$lower   # survival lower -> incidence upper
  
  tibble(
    group = grp,
    value = sprintf(
      paste0("%.", digits, "f%% (%.", digits, "f-%.", digits, "f%%)"),
      est * 100, lcl * 100, ucl * 100
    )
  )
}

# --------- helper: 기간별 분석용 time/status 생성 ---------
# 전체(0-5y): 1825에서 행정검열, 1825 이전 event만 event=1
# 0-1y: 365에서 행정검열, 365 이전 event만 event=1
# 1-5y(landmark): 1년까지 event-free & follow-up>=365인 사람만 포함,
#                 365 이후~1825 이내 event만 event=1, 시간은 (time-365)
prep_period_data <- function(data, time_raw, status_raw,
                             admin = 1825, landmark = 365,
                             period = c("overall", "0_1y", "1_5y")) {
  
  period <- match.arg(period)
  
  df <- data %>%
    mutate(
      t = .data[[time_raw]],
      e = .data[[status_raw]]
    )
  
  if (period == "overall") {
    df %>%
      mutate(
        time_p = pmin(t, admin),
        status_p = ifelse(e == 1 & t <= admin, 1, 0)
      )
  } else if (period == "0_1y") {
    df %>%
      mutate(
        time_p = pmin(t, landmark),
        status_p = ifelse(e == 1 & t <= landmark, 1, 0)
      )
  } else { # "1_5y"
    df %>%
      # 1년 시점에서 event-free이며 1년까지 추적이 존재해야 landmark risk set에 들어감
      filter(t > landmark) %>%   # t==365에서 event/censor는 landmark 이후 위험집단이 아님
      mutate(
        t_cap = pmin(t, admin),
        time_p = t_cap - landmark,
        status_p = ifelse(e == 1 & t > landmark & t <= admin, 1, 0)
      )
  }
}

# --------- 메인: 여러 outcome에 대해 wide table 생성 ---------
make_event_rate_wide <- function(data, group_var, outcomes,
                                 admin = 1825, landmark = 365, digits = 1) {
  # outcomes: named list
  # 예) list("TLF"=list(time="time_to_TLF", status="TLF_event"), ...)
  
  periods <- tribble(
    ~period_id, ~period_label, ~t_eval,
    "overall",  "Overall",     admin,
    "0_1y",     "0-1 year",    landmark,
    "1_5y",     "1-5 years",   admin - landmark
  )
  
  res_long <- imap_dfr(outcomes, function(spec, outcome_name) {
    
    bind_rows(lapply(1:nrow(periods), function(i) {
      pid <- periods$period_id[i]
      plabel <- periods$period_label[i]
      teval <- periods$t_eval[i]
      
      dfp <- prep_period_data(
        data = data,
        time_raw = spec$time,
        status_raw = spec$status,
        admin = admin,
        landmark = landmark,
        period = pid
      )
      
      # KM 누적발생률 문자열
      km_cuminc_str(
        df = dfp,
        time_col = "time_p",
        status_col = "status_p",
        group_col = group_var,
        t_eval = teval,
        digits = digits
      ) %>%
        mutate(Outcome = outcome_name, Period = plabel)
    }))
  })
  
  # wide: Period_group 형태로 펼치기 (원하시는 Overall_0, Overall_1 형태)
  res_wide <- res_long %>%
    mutate(col = paste0(Period, "_", group)) %>%
    select(Outcome, col, value) %>%
    pivot_wider(names_from = col, values_from = value)
  
  res_wide
}

# ---------------- 사용 예시 ----------------
# outcome time/status 컬럼명은 실제 데이터에 맞게 바꾸세요.
outcomes <- list(
  "Target lesion failure"           = list(time="time_to_composite", status="composite_event"),
  "Target vessel MI"                = list(time="time_to_MI", status="TVMI"),
  "Target lesion revascularization" = list(time="time_to_TLR", status="TLR"),
  "Target vessel revascularization" = list(time="time_to_TVR", status="TVR"),
  "Cardiac death"                   = list(time="time_to_death", status="cardiac_death"),
  "All-cause death"                 = list(time="time_to_death", status="any_death")
)

event_rate_wide_complete <- make_event_rate_wide(
  data = matched_data,
  group_var = "BP",     # 0/1
  outcomes = outcomes,
  admin = 1825,
  landmark = 365,
  digits = 1
)

event_rate_wide_complete

library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(broom)

# p-value 포맷
fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# HR 포맷
fmt_hr <- function(hr, lo, hi, digits = 2) {
  ifelse(any(is.na(c(hr, lo, hi))), NA_character_,
         sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
                 hr, lo, hi))
}

# 기간별 Cox HR (cluster robust), 0 event이면 NA 반환
cox_hr_period_safe <- function(data, group_var, time_raw, status_raw,
                               admin = 1825, landmark = 365,
                               period = c("overall","0_1y","1_5y"),
                               pair_id = "subclass") {
  period <- match.arg(period)
  
  dfp <- prep_period_data(data, time_raw, status_raw, admin, landmark, period) %>%
    mutate(
      pair = .data[[pair_id]],
      # group을 0이 reference가 되도록 정리 (BP가 factor여도 안전)
      g = factor(as.character(.data[[group_var]]))
    )
  
  # reference level을 "0"으로(있으면)
  if ("0" %in% levels(dfp$g)) dfp$g <- relevel(dfp$g, ref = "0")
  
  # 이벤트 0개면 Cox 불가 -> NA
  if (sum(dfp$status_p, na.rm = TRUE) == 0) {
    return(tibble(HR = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p = NA_real_,
                  n = nrow(dfp), events = 0))
  }
  
  fit <- tryCatch(
    coxph(Surv(time_p, status_p) ~ g + cluster(pair), data = dfp),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(HR = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p = NA_real_,
                  n = nrow(dfp), events = sum(dfp$status_p, na.rm = TRUE)))
  }
  
  td <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  
  # g의 비교항(보통 g1 또는 gBP 등) 첫 줄만 사용
  td1 <- td %>% slice(1)
  
  tibble(
    HR = td1$estimate,
    conf.low = td1$conf.low,
    conf.high = td1$conf.high,
    p = td1$p.value,
    n = nrow(dfp),
    events = sum(dfp$status_p, na.rm = TRUE)
  )
}

# 여러 outcome에 대해 HR wide table 생성
make_hr_wide <- function(data, group_var, outcomes,
                         admin = 1825, landmark = 365,
                         pair_id = "subclass",
                         hr_digits = 2) {
  
  periods <- tribble(
    ~period_id, ~period_label,
    "overall",  "Overall",
    "0_1y",     "0-1 year",
    "1_5y",     "1-5 years"
  )
  
  hr_long <- imap_dfr(outcomes, function(spec, outcome_name) {
    map_dfr(1:nrow(periods), function(i) {
      pid <- periods$period_id[i]
      plb <- periods$period_label[i]
      
      out <- cox_hr_period_safe(
        data = data,
        group_var = group_var,
        time_raw = spec$time,
        status_raw = spec$status,
        admin = admin,
        landmark = landmark,
        period = pid,
        pair_id = pair_id
      )
      
      tibble(
        Outcome = outcome_name,
        Period = plb,
        HR_str = fmt_hr(out$HR, out$conf.low, out$conf.high, digits = hr_digits),
        p_str  = fmt_p(out$p),
        events = out$events,
        n = out$n
      )
    })
  })
  
  # wide로: Overall_HR, Overall_p, ...
  hr_wide <- hr_long %>%
    mutate(
      col_hr = paste0(Period, "_HR"),
      col_p  = paste0(Period, "_p")
    ) %>%
    select(Outcome, col_hr, HR_str, col_p, p_str) %>%
    pivot_longer(cols = c(col_hr, col_p), names_to = "tmp", values_to = "colname") %>%
    mutate(value = ifelse(grepl("_HR$", colname), HR_str, p_str)) %>%
    select(Outcome, colname, value) %>%
    pivot_wider(names_from = colname, values_from = value)
  
  hr_wide
}

hr_wide_complete <- make_hr_wide(
  data = matched_data,
  group_var = "BP",
  outcomes = outcomes,
  admin = 1825,
  landmark = 365,
  pair_id = "subclass",
  hr_digits = 2
)

hr_wide_complete

