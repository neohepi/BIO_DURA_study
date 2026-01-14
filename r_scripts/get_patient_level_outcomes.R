################################################################################
#
#  BP-DES vs DP-DES Study - Data Preprocessing Pipeline (Updated)
#  
#  목적: Lesion 데이터에서 Patient 단위 분석 데이터셋 생성
#  추가: CTO, ISR, LAD/LCX/RCA involvement, Multivessel
#  
################################################################################

library(dplyr)
library(lubridate)
library(data.table)
library(tibble)

################################################################################
# 0. 유틸리티 함수
################################################################################

count_unique_id <- function(x) {
  cat("Unique IDs:", length(unique(x)), "\n")
  return(length(unique(x)))
}


################################################################################
# 1. Index Procedure 추출 (CTO, ISR, Vessel involvement 포함)
################################################################################

get_index_proc_from_lesion_data <- function(df.lesions) {
  
  cat("\n", rep("=", 70), "\n")
  cat("STEP 1: EXTRACTING INDEX PROCEDURES\n")
  cat(rep("=", 70), "\n\n")
  
  # Vessel mapping 정의
  target_LAD <- c("plad", "mlad", "dlad", "dg", "septal_br")
  target_LCX <- c("plcx", "dlcx", "om", "ri")
  target_RCA <- c("prca", "mrca", "drca", "plb", "pda")
  
  index_procedure_df <- df.lesions %>%
    # lesion_anatomy 정규화
    mutate(lesion_norm = tolower(trimws(lesion_anatomy))) %>%
    # 1. 시술(Procedure) 단위로 그룹화
    group_by(pt_id, CAG_date, cath_no) %>%
    summarise(
      # Stent type 확인
      has_DP = any(DP == 1, na.rm = TRUE),
      has_BP = any(BP == 1, na.rm = TRUE),
      
      # Stent names 수집 (unique, comma-separated)
      stent_names = paste(unique(na.omit(stent_name)), collapse = ", "),
      
      # Lesion 정보 수집
      n_lesions = n(),
      lesion_anatomies = paste(unique(na.omit(lesion_anatomy)), collapse = ", "),
      
      # === Procedural Characteristics (from add_procedural_data) ===
      # Stent count, length, diameter
      num_stents = n(),
      total_stent_length = sum(stent_length, na.rm = TRUE),
      avg_stent_diameter = mean(stent_diameter, na.rm = TRUE),
      
      # LM involvement
      has_LM = any(lesion_norm == "lm", na.rm = TRUE),
      
      # Vessel involvement (LAD, LCX, RCA)
      has_LAD = any(lesion_norm %in% target_LAD, na.rm = TRUE),
      has_LCX = any(lesion_norm %in% target_LCX, na.rm = TRUE),
      has_RCA = any(lesion_norm %in% target_RCA, na.rm = TRUE),
      
      # Bifurcation
      has_bifurcation = any(bifurcation == 1 | bifurcation == "Yes", na.rm = TRUE),
      
      # Complex lesion (ACC/AHA type B2 or C)
      has_complex_lesion = any(ACC_type %in% c("B2", "C"), na.rm = TRUE),
      
      # CTO (하나라도 CTO면 1)
      has_CTO = any(CTO == 1 | CTO == "Yes" | CTO == TRUE, na.rm = TRUE),
      
      # ISR (하나라도 ISR이면 1)
      has_ISR = any(ISR_lesion == 1 | ISR_lesion == "Yes" | ISR_lesion == TRUE, na.rm = TRUE),
      
      .groups = 'drop'
    ) %>%
    # 2. Mixed Type 여부 및 대표 Stent Type 결정
    mutate(
      is_mixed = has_DP & has_BP,
      representative_BP = ifelse(has_BP, 1, 0)
    ) %>%
    # 3. 환자별 첫 시술 선택
    group_by(pt_id) %>%
    arrange(CAG_date) %>%
    slice(1) %>%
    ungroup() %>%
    # 4. Mixed Type 제외
    filter(!is_mixed) %>%
    # 5. Multivessel 계산 (LAD, LCX, RCA 중 2개 이상)
    mutate(
      n_vessels = as.integer(has_LAD) + as.integer(has_LCX) + as.integer(has_RCA),
      is_multivessel = n_vessels >= 2
    ) %>%
    # 6. NA 처리 (add_procedural_data 방식)
    mutate(
      num_stents = ifelse(is.na(num_stents), 0, num_stents),
      total_stent_length = ifelse(is.na(total_stent_length), 0, total_stent_length),
      avg_stent_diameter = ifelse(is.na(avg_stent_diameter), 0, avg_stent_diameter)
    ) %>%
    # 7. 컬럼 정리
    select(
      pt_id, 
      cath_no, 
      CAG_date, 
      BP = representative_BP,
      stent_names,
      n_lesions,
      lesion_anatomies,
      # Procedural characteristics
      num_stents,
      total_stent_length,
      avg_stent_diameter,
      # Vessel involvement
      is_LM = has_LM,
      is_LAD = has_LAD,
      is_LCX = has_LCX,
      is_RCA = has_RCA,
      n_vessels,
      is_multivessel,
      # Lesion complexity
      is_bifurcation = has_bifurcation,
      is_complex_lesion = has_complex_lesion,
      is_CTO = has_CTO,
      is_ISR = has_ISR
    )
  
  cat("Total index procedures:", nrow(index_procedure_df), "\n")
  cat("BP-DES:", sum(index_procedure_df$BP == 1), "\n")
  cat("DP-DES:", sum(index_procedure_df$BP == 0), "\n\n")
  
  cat("Procedural characteristics:\n")
  cat("  Mean stents per patient:", round(mean(index_procedure_df$num_stents), 2), "\n")
  cat("  Mean total stent length:", round(mean(index_procedure_df$total_stent_length), 1), "mm\n")
  cat("  Mean stent diameter:", round(mean(index_procedure_df$avg_stent_diameter, na.rm = TRUE), 2), "mm\n\n")
  
  cat("Lesion characteristics:\n")
  cat("  LM involvement:", sum(index_procedure_df$is_LM), "\n")
  cat("  LAD involvement:", sum(index_procedure_df$is_LAD), "\n")
  cat("  LCX involvement:", sum(index_procedure_df$is_LCX), "\n")
  cat("  RCA involvement:", sum(index_procedure_df$is_RCA), "\n")
  cat("  Multivessel:", sum(index_procedure_df$is_multivessel), "\n")
  cat("  Bifurcation:", sum(index_procedure_df$is_bifurcation), "\n")
  cat("  Complex lesion (B2/C):", sum(index_procedure_df$is_complex_lesion), "\n")
  cat("  CTO:", sum(index_procedure_df$is_CTO), "\n")
  cat("  ISR:", sum(index_procedure_df$is_ISR), "\n")
  
  return(index_procedure_df)
}


################################################################################
# 2. Follow-up MI 이벤트 추출
################################################################################

get_first_MI_event_from_index_proc <- function(index_proc, all_lesions, cag_data) {
  
  cat("\n", rep("=", 70), "\n")
  cat("STEP 2: EXTRACTING MI EVENTS\n")
  cat(rep("=", 70), "\n\n")
  
  # Index date 기준
  index_base <- index_proc %>%
    select(pt_id, index_date = CAG_date) 
  
  # Index 이후 첫 ACS 이벤트 찾기
  next_acs_event <- index_base %>%
    left_join(cag_data, by = "pt_id") %>%
    filter(CAG_date > index_date, is_ACS == 1) %>%
    group_by(pt_id) %>%
    arrange(CAG_date) %>%
    slice(1) %>%
    ungroup() %>%
    select(pt_id, next_acs_date = CAG_date, index_date)
  
  # Index procedure에 합치기
  df <- index_proc %>% 
    left_join(next_acs_event, by = "pt_id") %>%
    select(-index_date) %>%
    rename(next_MI_date = next_acs_date) %>%
    mutate(
      next_MI = ifelse(is.na(next_MI_date), 0, 1),
      next_MI_date = as.Date(next_MI_date)
    )
  
  # Last follow-up date 추가
  temp <- all_lesions %>% 
    select(pt_id, last_fu_date) %>%
    mutate(last_fu_date = as.Date(last_fu_date)) %>%
    distinct()
  
  df <- df %>% 
    left_join(temp, by = "pt_id") %>%
    mutate(
      CAG_date = as.Date(CAG_date),
      BP = as.factor(BP),
      next_MI = as.factor(next_MI)
    )
  
  # Time to MI 계산
  df <- df %>%
    mutate(
      days_to_MI = ifelse(
        next_MI == 1,
        as.numeric(next_MI_date - CAG_date),
        as.numeric(last_fu_date - CAG_date)
      ),
      days_to_MI = as.integer(days_to_MI)
    )
  
  cat("MI events - DP-DES:", sum(df$next_MI == 1 & df$BP == 0), "\n")
  cat("MI events - BP-DES:", sum(df$next_MI == 1 & df$BP == 1), "\n")
  
  return(df)
}


################################################################################
# 3. Death 이벤트 추출
################################################################################

get_death_events_from_lesion <- function(df_lesions) {
  
  cat("\n", rep("=", 70), "\n")
  cat("STEP 3: EXTRACTING DEATH EVENTS\n")
  cat(rep("=", 70), "\n\n")
  
  death_df <- df_lesions %>%
    select(pt_id, death_date, cardiac_death) %>%
    mutate(
      death_date = as.Date(death_date),
      any_death = if_else(!is.na(death_date), 1, 0),
      cardiac_death_date = case_when(
        cardiac_death == 1 ~ death_date,
        TRUE ~ as.Date(NA)
      )
    ) %>%
    distinct()
  
  cat("Any death:", sum(death_df$any_death, na.rm = TRUE), "\n")
  cat("Cardiac death:", sum(death_df$cardiac_death, na.rm = TRUE), "\n")
  
  return(death_df)
}


################################################################################
# 4. Stent Thrombosis (ST) 추출
################################################################################

calc_ST_after_index_vessel <- function(
    df_pts, 
    df_lesions,
    pt_col = "pt_id",
    date_col = "CAG_date",
    lesion_col = "lesion_anatomy",
    thrombus_col = "thrombus",
    outcome_col = "ST_outcome",
    st_date_col = "ST_date",
    st_vessel_col = "ST_vessel",
    allow_same_day = FALSE
) {
  
  cat("\n", rep("=", 70), "\n")
  cat("STEP 4: CALCULATING STENT THROMBOSIS\n")
  cat(rep("=", 70), "\n\n")
  
  stopifnot(is.data.frame(df_pts), is.data.frame(df_lesions))
  
  # Vessel mapping
  target_LAD <- c("pLAD", "mLAD", "dLAD", "Dg", "septal_br")
  target_LCX <- c("pLCX", "dLCX", "OM", "RI")
  target_RCA <- c("pRCA", "mRCA", "dRCA", "PLB", "PDA")
  target_LM  <- c("LM")
  
  lad_set <- tolower(target_LAD)
  lcx_set <- tolower(target_LCX)
  rca_set <- tolower(target_RCA)
  lm_set  <- tolower(target_LM)
  
  # Data.table 변환
  pts <- as.data.table(copy(df_pts))
  les <- as.data.table(copy(df_lesions))
  
  pts[, idx_day := as.Date(get(date_col))]
  les[, ev_day  := as.Date(get(date_col))]
  
  les[, lesion_norm := tolower(trimws(as.character(get(lesion_col))))]
  
  les[, vessel := fcase(
    lesion_norm %chin% lad_set, "LAD",
    lesion_norm %chin% lcx_set, "LCX",
    lesion_norm %chin% rca_set, "RCA",
    lesion_norm %chin% lm_set,  "LM",
    default = NA_character_
  )]
  
  les[, thrombus_flag := {
    x <- get(thrombus_col)
    if (is.logical(x)) x
    else if (is.numeric(x)) x == 1
    else if (is.character(x)) toupper(trimws(x)) %in% c("1", "TRUE", "T", "YES", "Y")
    else as.logical(x)
  }]
  
  pts[, .row_id := .I]
  
  # Index vessels 식별
  idx_vessels <- les[
    pts,
    on = c(setNames(pt_col, pt_col), "ev_day" = "idx_day"),
    nomatch = 0L,
    .(.row_id = i..row_id, vessel = vessel)
  ][!is.na(vessel)]
  
  idx_vessels <- unique(idx_vessels, by = c(".row_id", "vessel"))
  
  # Index 이후 thrombus 후보
  cmp <- if (allow_same_day) `>=` else `>`
  
  cand <- les[
    pts,
    on = c(setNames(pt_col, pt_col)),
    allow.cartesian = TRUE,
    nomatch = 0L,
    .(.row_id = i..row_id,
      idx_day = i.idx_day,
      ev_day  = ev_day,
      vessel  = vessel,
      thrombus_flag = thrombus_flag)
  ][!is.na(vessel) & thrombus_flag == TRUE & cmp(ev_day, idx_day)]
  
  cand2 <- cand[idx_vessels, on = .(.row_id, vessel), nomatch = 0L]
  
  if (nrow(cand2) == 0) {
    out <- copy(pts)
    out[, (outcome_col) := 0L]
    out[, (st_date_col) := as.Date(NA)]
    out[, (st_vessel_col) := NA_character_]
    out[, .row_id := NULL]
    out[, idx_day := NULL]
    
    cat("ST events: 0\n")
    return(as.data.frame(out))
  }
  
  setorder(cand2, .row_id, ev_day)
  first_ev <- cand2[, .SD[1], by = .(.row_id)]
  st_ev <- first_ev[, .(.row_id, ST_date = ev_day, ST_vessel = vessel)]
  
  out <- merge(pts, st_ev, by = ".row_id", all.x = TRUE)
  out[, (outcome_col) := as.integer(!is.na(ST_date))]
  setnames(out, "ST_date", st_date_col)
  setnames(out, "ST_vessel", st_vessel_col)
  
  out[, .row_id := NULL]
  out[, idx_day := NULL]
  
  cat("ST events:", sum(out[[outcome_col]]), "\n")
  
  return(as.data.frame(out))
}


################################################################################
# 5. Target Vessel MI (TV-MI) 추출
################################################################################

get_TVMI_events <- function(
    df_pts, 
    df_lesions,
    pt_col = "pt_id",
    date_col = "CAG_date",
    lesion_col = "lesion_anatomy",
    mi_col = "is_ACS", 
    outcome_col = "TVMI",
    tvmi_date_col = "TVMI_date"
) {
  
  cat("\n", rep("=", 70), "\n")
  cat("STEP 5: CALCULATING TARGET VESSEL MI\n")
  cat(rep("=", 70), "\n\n")
  
  # Vessel mapping
  target_LAD <- c("pLAD", "mLAD", "dLAD", "Dg", "septal_br")
  target_LCX <- c("pLCX", "dLCX", "OM", "RI")
  target_RCA <- c("pRCA", "mRCA", "dRCA", "PLB", "PDA")
  target_LM  <- c("LM")
  
  lad_set <- tolower(target_LAD)
  lcx_set <- tolower(target_LCX)
  rca_set <- tolower(target_RCA)
  lm_set  <- tolower(target_LM)
  
  pts <- as.data.table(copy(df_pts))
  les <- as.data.table(copy(df_lesions))
  
  pts[, idx_day := as.Date(get(date_col))]
  les[, ev_day  := as.Date(get(date_col))]
  
  les[, lesion_norm := tolower(trimws(as.character(get(lesion_col))))]
  
  les[, vessel := fcase(
    lesion_norm %chin% lad_set, "LAD",
    lesion_norm %chin% lcx_set, "LCX",
    lesion_norm %chin% rca_set, "RCA",
    lesion_norm %chin% lm_set,  "LM",
    default = NA_character_
  )]
  
  les[, mi_flag := {
    x <- get(mi_col)
    if (is.numeric(x)) x == 1
    else if (is.logical(x)) x
    else !is.na(x) & (x == 1 | x == "1") 
  }]
  
  pts[, .row_id := .I]
  
  # Index vessels 식별
  idx_vessels <- les[
    pts,
    on = c(setNames(pt_col, pt_col), "ev_day" = "idx_day"),
    nomatch = 0L,
    .(.row_id = i..row_id, vessel = vessel)
  ][!is.na(vessel)]
  
  idx_vessels <- unique(idx_vessels, by = c(".row_id", "vessel"))
  
  # Index 이후 MI 후보
  mi_candidates <- les[
    pts,
    on = c(setNames(pt_col, pt_col)),
    allow.cartesian = TRUE,
    nomatch = 0L,
    .(.row_id = i..row_id,
      idx_day = i.idx_day,
      ev_day  = ev_day,
      vessel  = vessel,
      mi_flag = mi_flag)
  ][!is.na(vessel) & mi_flag == TRUE & ev_day > idx_day]
  
  tv_mi_events <- mi_candidates[idx_vessels, on = .(.row_id, vessel), nomatch = 0L]
  
  if (nrow(tv_mi_events) == 0) {
    out <- copy(pts)
    out[, (outcome_col) := 0L]
    out[, (tvmi_date_col) := as.Date(NA)]
    out[, c(".row_id", "idx_day") := NULL]
    
    cat("TV-MI events: 0\n")
    return(as.data.frame(out))
  }
  
  setorder(tv_mi_events, .row_id, ev_day)
  first_ev <- tv_mi_events[, .SD[1], by = .(.row_id)]
  final_ev <- first_ev[, .(.row_id, TVMI_date = ev_day)]
  
  out <- merge(pts, final_ev, by = ".row_id", all.x = TRUE)
  out[, (outcome_col) := as.integer(!is.na(TVMI_date))]
  setnames(out, "TVMI_date", tvmi_date_col)
  
  out[, c(".row_id", "idx_day") := NULL]
  
  cat("TV-MI events:", sum(out[[outcome_col]]), "\n")
  
  return(as.data.frame(out))
}


################################################################################
# 6. Revascularization (TLR, TVR) 추출
################################################################################

get_any_revasc_from_lesion_data <- function(df_pts, df_lesions) {
  
  cat("\n", rep("=", 70), "\n")
  cat("STEP 6: EXTRACTING REVASCULARIZATION EVENTS\n")
  cat(rep("=", 70), "\n\n")
  
  temp <- df_lesions %>%
    select(pt_id, TLR, TLR_date, TVR, TVR_date)
  
  # TLR
  temp_tlr_events <- temp %>%
    filter(TLR == 1) %>% 
    group_by(pt_id) %>%
    summarise(first_tlr_date = min(TLR_date, na.rm = TRUE), .groups = 'drop')
  
  # TVR
  temp_tvr_events <- temp %>%
    filter(TVR == 1) %>% 
    group_by(pt_id) %>%
    summarise(first_tvr_date = min(TVR_date, na.rm = TRUE), .groups = 'drop')
  
  df_updated <- df_pts %>%
    left_join(temp_tlr_events, by = "pt_id") %>%
    left_join(temp_tvr_events, by = "pt_id") %>%
    mutate(
      TLR = if_else(!is.na(first_tlr_date), 1, 0),
      TLR_date = as.Date(first_tlr_date),
      TVR = if_else(!is.na(first_tvr_date), 1, 0),
      TVR_date = as.Date(first_tvr_date)
    ) %>%
    select(-first_tlr_date, -first_tvr_date)
  
  cat("TLR events:", sum(df_updated$TLR), "\n")
  cat("TVR events:", sum(df_updated$TVR), "\n")
  
  return(df_updated)
}


################################################################################
# 7. Primary Endpoint (TLF) 계산
################################################################################

get_primary_endpoints_revised <- function(df_all_events) {
  
  cat("\n", rep("=", 70), "\n")
  cat("STEP 7: CALCULATING PRIMARY ENDPOINT (TLF)\n")
  cat(rep("=", 70), "\n\n")
  
  df <- df_all_events %>% 
    mutate(
      TLR_n = as.numeric(as.character(TLR)),
      TVMI_n = as.numeric(as.character(TVMI)),
      CD_n = as.numeric(as.character(cardiac_death)),
      
      # Composite: TLR + TV-MI + Cardiac Death
      composite_event = pmax(TLR_n, TVMI_n, CD_n, na.rm = TRUE),
      composite_date = pmin(TLR_date, TVMI_date, cardiac_death_date, na.rm = TRUE)
    ) %>%
    mutate(
      composite_event = ifelse(is.na(composite_event), 0, composite_event),
      composite_date = if_else(
        composite_event == 0 | is.infinite(composite_date), 
        as.Date(NA), 
        composite_date
      )
    ) %>%
    select(-ends_with("_n"))
  
  cat("TLF events:", sum(df$composite_event), "\n")
  
  return(df)
}


################################################################################
# 8. Baseline 정보 추출
################################################################################

get_baseline_info_from_lesion_data <- function(df_lesions) {
  
  cat("\n", rep("=", 70), "\n")
  cat("STEP 8: EXTRACTING BASELINE CHARACTERISTICS\n")
  cat(rep("=", 70), "\n\n")
  
  baseline_df <- df_lesions %>% 
    select(pt_id, age, HTN, DM, DL, smoking, prior_MI, prior_PCI, 
           prior_CABG, LVEF, TC, TG, HDL, LDL, BMI, male, AFib, CKD, 
           COPD, CAOD, HF, ESRD, PAOD, prior_stroke, TIA) %>%
    mutate(
      across(c("pt_id", "age"), as.integer),
      across(c("HTN", "DM", "DL", "smoking", "prior_MI", "prior_PCI", "prior_CABG", 
               "male", "AFib", "CKD", "COPD", "CAOD", "HF", "ESRD", "PAOD", 
               "prior_stroke", "TIA"), as.factor)
    ) %>%
    distinct()
  
  cat("Baseline records:", nrow(baseline_df), "\n")
  
  return(baseline_df)
}


################################################################################
# 9. Era 구분
################################################################################

add_era_variable <- function(df) {
  
  cat("\n", rep("=", 70), "\n")
  cat("STEP 9: ADDING ERA VARIABLE\n")
  cat(rep("=", 70), "\n\n")
  
  df <- df %>%
    filter(procedure_year >= 2010 & procedure_year <= 2021) %>% 
    mutate(
      era = cut(
        procedure_year,
        breaks = c(2009, 2013, 2016, 2021), 
        labels = c("2010-2013", "2014-2016", "2017-2021")
      ),
      era = as.factor(era)
    )
  
  cat("Era distribution:\n")
  print(table(df$era))
  
  return(df)
}


################################################################################
#
#  MAIN PIPELINE EXECUTION
#
################################################################################

run_preprocessing_pipeline <- function(df.lesion.all, df.cag.all, all_pcis) {
  
  cat("\n")
  cat(rep("#", 80), "\n")
  cat("#  BP-DES vs DP-DES: DATA PREPROCESSING PIPELINE\n")
  cat(rep("#", 80), "\n")
  
  # Step 1: Index procedure 추출 (CTO, ISR, vessel involvement 포함)
  df.index_proc <- get_index_proc_from_lesion_data(df.lesion.all)
  
  # Step 2: MI 이벤트
  df.pt.events <- get_first_MI_event_from_index_proc(df.index_proc, df.lesion.all, df.cag.all)
  
  # Step 3: Death 이벤트
  temp_death <- get_death_events_from_lesion(df.lesion.all)
  df.pt.events <- df.pt.events %>% left_join(temp_death, by = "pt_id")
  
  # Step 4: ST
  df.pt.events <- as_tibble(calc_ST_after_index_vessel(df.pt.events, all_pcis))
  
  # Step 5: TV-MI
  df.pt.events <- as_tibble(get_TVMI_events(
    df_pts = df.pt.events, 
    df_lesions = all_pcis,
    mi_col = "is_ACS",
    outcome_col = "TVMI",
    tvmi_date_col = "TVMI_date"
  ))
  
  # Step 6: Revascularization
  df.pt.events <- get_any_revasc_from_lesion_data(df.pt.events, df.lesion.all)
  
  # Step 7: Primary endpoint (TLF)
  df.primary <- get_primary_endpoints_revised(df.pt.events)
  
  # Step 8: Baseline 정보
  df.baseline <- get_baseline_info_from_lesion_data(df.lesion.all)
  df.pts <- df.primary %>% left_join(df.baseline, by = "pt_id")
  
  # CAG 정보 추가
  temp_cag <- df.cag.all %>% select(cath_no, is_ACS, ACS_type, num_CAOD)
  df.pts <- df.pts %>% left_join(temp_cag, by = "cath_no")
  
  # Follow-up days 계산 (올바른 방식)
  df.pts <- df.pts %>% 
    mutate(
      procedure_year = year(CAG_date),
      fu_days = as.numeric(last_fu_date - CAG_date)  # 올바른 계산
    )
  
  # Step 9: Era 구분
  df.pts <- df.pts %>% filter(procedure_year >= 2010)
  df.pts <- add_era_variable(df.pts)
  
  # 최종 검증
  cat("\n", rep("=", 70), "\n")
  cat("FINAL DATA VALIDATION\n")
  cat(rep("=", 70), "\n\n")
  
  cat("Total patients:", nrow(df.pts), "\n")
  cat("BP-DES:", sum(df.pts$BP == 1), "\n")
  cat("DP-DES:", sum(df.pts$BP == 0), "\n\n")
  
  cat("Follow-up:\n")
  cat("  Median:", round(median(df.pts$fu_days, na.rm = TRUE)/365.25, 2), "years\n")
  cat("  IQR:", round(quantile(df.pts$fu_days, 0.25, na.rm = TRUE)/365.25, 2), "-",
      round(quantile(df.pts$fu_days, 0.75, na.rm = TRUE)/365.25, 2), "years\n\n")
  
  cat("Procedural characteristics:\n")
  cat("  Mean stents per patient:", round(mean(df.pts$num_stents, na.rm = TRUE), 2), "\n")
  cat("  Mean total stent length:", round(mean(df.pts$total_stent_length, na.rm = TRUE), 1), "mm\n")
  cat("  Mean stent diameter:", round(mean(df.pts$avg_stent_diameter, na.rm = TRUE), 2), "mm\n\n")
  
  cat("Lesion characteristics:\n")
  cat("  LM:", sum(df.pts$is_LM), "\n")
  cat("  LAD:", sum(df.pts$is_LAD), "\n")
  cat("  LCX:", sum(df.pts$is_LCX), "\n")
  cat("  RCA:", sum(df.pts$is_RCA), "\n")
  cat("  Multivessel:", sum(df.pts$is_multivessel), "\n")
  cat("  Bifurcation:", sum(df.pts$is_bifurcation), "\n")
  cat("  Complex lesion (B2/C):", sum(df.pts$is_complex_lesion), "\n")
  cat("  CTO:", sum(df.pts$is_CTO), "\n")
  cat("  ISR:", sum(df.pts$is_ISR), "\n\n")
  
  cat("Events:\n")
  cat("  TLF (Composite):", sum(df.pts$composite_event, na.rm = TRUE), "\n")
  cat("  Cardiac death:", sum(as.numeric(as.character(df.pts$cardiac_death)), na.rm = TRUE), "\n")
  cat("  TV-MI:", sum(df.pts$TVMI, na.rm = TRUE), "\n")
  cat("  TLR:", sum(df.pts$TLR, na.rm = TRUE), "\n")
  cat("  ST:", sum(df.pts$ST_outcome, na.rm = TRUE), "\n")
  cat("  Any death:", sum(as.numeric(as.character(df.pts$any_death)), na.rm = TRUE), "\n\n")
  
  cat("Stent names (top 10):\n")
  print(head(sort(table(df.pts$stent_names), decreasing = TRUE), 10))
  
  cat("\n", rep("#", 80), "\n")
  cat("#  PREPROCESSING COMPLETE\n")
  cat(rep("#", 80), "\n\n")
  
  return(df.pts)
}


################################################################################
# 실행
################################################################################

# 파이프라인 실행
df.pts <- run_preprocessing_pipeline(df.lesion.all, df.cag.all, all_pcis)

# 결과 확인
cat("\n=== Stent Names Distribution ===\n")
print(head(sort(table(df.pts$stent_names), decreasing = TRUE), 20))

cat("\n=== Lesion Characteristics by Stent Type ===\n")
cat("\nLM Disease:\n")
print(table(df.pts$is_LM, df.pts$BP, dnn = c("is_LM", "BP")))

cat("\nLAD Involvement:\n")
print(table(df.pts$is_LAD, df.pts$BP, dnn = c("is_LAD", "BP")))

cat("\nLCX Involvement:\n")
print(table(df.pts$is_LCX, df.pts$BP, dnn = c("is_LCX", "BP")))

cat("\nRCA Involvement:\n")
print(table(df.pts$is_RCA, df.pts$BP, dnn = c("is_RCA", "BP")))

cat("\nMultivessel:\n")
print(table(df.pts$is_multivessel, df.pts$BP, dnn = c("is_multivessel", "BP")))

cat("\nBifurcation:\n")
print(table(df.pts$is_bifurcation, df.pts$BP, dnn = c("is_bifurcation", "BP")))

cat("\nCTO:\n")
print(table(df.pts$is_CTO, df.pts$BP, dnn = c("is_CTO", "BP")))

cat("\nISR:\n")
print(table(df.pts$is_ISR, df.pts$BP, dnn = c("is_ISR", "BP")))

cat("\nComplex Lesion (B2/C):\n")
print(table(df.pts$is_complex_lesion, df.pts$BP, dnn = c("is_complex_lesion", "BP")))

cat("\n=== Procedural Characteristics by Stent Type ===\n")
proc_summary <- df.pts %>%
  group_by(BP) %>%
  summarise(
    N = n(),
    Mean_Stents = round(mean(num_stents, na.rm = TRUE), 2),
    SD_Stents = round(sd(num_stents, na.rm = TRUE), 2),
    Mean_Length = round(mean(total_stent_length, na.rm = TRUE), 1),
    SD_Length = round(sd(total_stent_length, na.rm = TRUE), 1),
    Mean_Diameter = round(mean(avg_stent_diameter, na.rm = TRUE), 2),
    SD_Diameter = round(sd(avg_stent_diameter, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(Group = ifelse(BP == 1, "BP-DES", "DP-DES"))
print(proc_summary)