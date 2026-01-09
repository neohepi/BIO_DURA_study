#df.lesions
df.lesion.all$ISR_outcome

################################################################################
# patient 단위로 재구성
count_unique_id(df.lesion.all$pt_id) # 9,586

################################################################################
# Index procedure 확인
get_index_proc_from_lesion_data <- function(df.lesions) {
  index_procedure_df <- df.lesions %>%
    # 1. 시술(Procedure) 단위로 그룹화 (환자 + 날짜 + 시술번호)
    # 같은 시술 내에 여러 병변(row)이 있을 수 있으므로 이를 요약해야 합니다.
    group_by(pt_id, CAG_date, cath_no) %>%
    summarise(
      # 해당 시술 내에 DP가 하나라도 있는지, BP가 하나라도 있는지 확인
      has_DP = any(DP == 1, na.rm = TRUE),
      has_BP = any(BP == 1, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    # 2. Mixed Type 여부 및 대표 Stent Type 결정
    mutate(
      # DP와 BP가 둘 다 있으면 Mixed
      is_mixed = has_DP & has_BP,
      
      # Mixed가 아닐 때: BP가 있으면(has_BP=TRUE) 1, 없으면(DP만 있으면) 0
      # (연구 목적에 따라 BP=1, DP=0으로 매핑한다고 가정)
      representative_BP = ifelse(has_BP, 1, 0)
    ) %>%
    # 3. 환자별로 정렬하여 가장 빠른 날짜 선택
    group_by(pt_id) %>%
    arrange(CAG_date) %>%
    slice(1) %>% # 환자별 가장 첫 번째 시술만 남김
    ungroup() %>%
    # 4. Mixed Type인 경우(Index procedure가 mixed인 환자) 제외
    filter(!is_mixed) %>%
    # 5. 최종 컬럼 선택 및 이름 정리
    select(pt_id, cath_no, CAG_date, BP = representative_BP)
  
  return(index_procedure_df)
}

df.index_proc.all <- get_index_proc_from_lesion_data(df.lesion.all)

# index procedure 선택
print(df.index_proc.all) # 9,289

table(df.index_proc.all$BP) # BP 3,667; DP 5,622

################################################################################
# MI 확인하기
# cag_data에서 pt_id에 일치하는 row를 확인, CAG_date보다 큰 CAG_date를 갖는 CAG data를 찾고
# 그 중에 is_ACS==1인 가장 빠른 날짜를 확인한다
get_first_MI_event_from_index_proc <- function(index_proc, all_lesions, cag_data) {
  # 1. index_procedure_df에서 필요한 컬럼만 선택하고 날짜 이름 변경 (기준 날짜)
  index_base <- index_proc %>%
    select(pt_id, index_date = CAG_date) 
  
  # 2. 전체 데이터(cag_data)와 조인하여 조건에 맞는 Row 찾기
  next_acs_event <- index_base %>%
    # pt_id를 기준으로 전체 CAG 데이터와 결합 (Left Join)
    left_join(cag_data, by = "pt_id") %>%
    
    # 조건 1: Index Date보다 미래에 발생한 건이어야 함 (Strictly Greater)
    # 조건 2: ACS 이벤트여야 함 (is_ACS == 1)
    filter(CAG_date > index_date, is_ACS == 1) %>%
    
    # 조건 3: 그 중에서 가장 빠른 날짜 찾기
    group_by(pt_id) %>%
    arrange(CAG_date) %>%
    slice(1) %>% # 환자별로 가장 상단의 행(가장 빠른 날짜)만 선택
    ungroup() %>%
    
    # 4. 결과 정리 (필요한 컬럼만 남김)
    select(pt_id, next_acs_date = CAG_date, index_date)

  df <- index_proc %>% 
    left_join(next_acs_event, by="pt_id") %>%
    select(-index_date) %>%
    rename(next_MI_date=next_acs_date) %>%
    mutate(next_MI=ifelse(is.na(next_MI_date), 0, 1),
           next_MI_date=as.Date(next_MI_date))
  
  temp <- all_lesions %>% select(pt_id, last_fu_date) %>%
    mutate(last_fu_date=as.Date(last_fu_date)) %>%
    distinct()
  
  df <- df %>% 
    left_join(temp, by="pt_id") %>%
    mutate(CAG_date=as.Date(CAG_date),
           across(c("BP", "next_MI"), as.factor))
  
  # adding variables to use in the survival analysis
  df <- df %>%
    mutate(
      days_to_MI = ifelse(
        next_MI == 1,
        as.numeric(next_MI_date - CAG_date), 
        # Event 발생 시
        as.numeric(last_fu_date - CAG_date)  # Censored 시
      ),
      days_to_MI = as.integer(days_to_MI)
    )
  
  return(df)
}

df.pt.events <- get_first_MI_event_from_index_proc(df.index_proc.all, df.lesion.all, df.cag.all)

df.pt.events %>% filter(is.na(next_MI_date) == FALSE & BP == 0) %>% count() # DP, MI 629
df.pt.events %>% filter(is.na(next_MI_date) == FALSE & BP == 1) %>% count() # BP, MI 192


################################################################################
# Death/cardiac death

get_death_events_from_lesion <- function(df_lesions) {
  return(df_lesions %>%
           select(pt_id, death_date, cardiac_death) %>%
           mutate(death_date = as.Date(death_date),
                  # 1. any_death: 사망 날짜가 있으면(!is.na) 1, 없으면 0
                  any_death = if_else(!is.na(death_date), 1, 0),
                  
                  # 2. cardiac_death_date: 심인성 사망(1)인 경우에만 날짜를 가져오고, 나머지는 NA 처리
                  # case_when을 쓰면 날짜 형식(Date class)을 유지하면서 NA를 넣기에 안전합니다.
                  cardiac_death_date = case_when(
                    cardiac_death == 1 ~ death_date
                  )) %>%
           distinct()
    
  )
}
temp <- get_death_events_from_lesion(df.lesion.all)

df.pt.events <- df.pt.events %>% # left_join with death events (cardiac and any death)
  left_join(temp, by="pt_id")
df.pt.events


################################################################################
# ST 계산하기
# index procedure 이후에 동일한 lesion anatomy에 시술이 있고 당시 lesion 특성에 thrombus가 1이면 ST로 계산
calc_ST_after_index_vessel <- function(
    df_pts, df_lesions,
    pt_col = "pt_id",
    date_col = "CAG_date",
    lesion_col = "lesion_anatomy",
    thrombus_col = "thrombus",
    outcome_col = "ST_outcome",
    st_date_col = "ST_date",
    st_vessel_col = "ST_vessel",
    allow_same_day = FALSE  # Date만 있으면 TRUE는 index thrombus를 ST로 오인할 수 있어 기본 FALSE 권장
) {
  stopifnot(is.data.frame(df_pts), is.data.frame(df_lesions))
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install with install.packages('data.table').")
  }
  library(data.table)
  
  # --- user-defined anatomy groups ---
  target_LAD <- c("pLAD", "mLAD", "dLAD", "Dg", "septal_br")
  target_LCX <- c("pLCX", "dLCX", "OM", "RI")
  target_RCA <- c("pRCA", "mRCA", "dRCA", "PLB", "PDA")
  target_LM  <- c("LM")
  
  lad_set <- tolower(target_LAD)
  lcx_set <- tolower(target_LCX)
  rca_set <- tolower(target_RCA)
  lm_set  <- tolower(target_LM)
  
  # --- copy to data.table ---
  pts <- as.data.table(copy(df_pts))
  les <- as.data.table(copy(df_lesions))
  
  # --- dates (Date) ---
  pts[, idx_day := as.Date(get(date_col))]
  les[, ev_day  := as.Date(get(date_col))]
  
  # --- normalize lesion_anatomy + map to vessel ---
  les[, lesion_norm := tolower(trimws(as.character(get(lesion_col))))]
  
  les[, vessel := fcase(
    lesion_norm %chin% lad_set, "LAD",
    lesion_norm %chin% lcx_set, "LCX",
    lesion_norm %chin% rca_set, "RCA",
    lesion_norm %chin% lm_set,  "LM",
    default = NA_character_
  )]
  
  # --- thrombus flag ---
  les[, thrombus_flag := {
    x <- get(thrombus_col)
    if (is.logical(x)) x
    else if (is.numeric(x)) x == 1
    else if (is.character(x)) toupper(trimws(x)) %in% c("1", "TRUE", "T", "YES", "Y")
    else as.logical(x)
  }]
  
  # row id to preserve multiple rows in df_pts (if any)
  pts[, .row_id := .I]
  
  # -------------------------------------------
  # 1) index day에 해당 환자의 index_vessels 세트 만들기 (row_id 단위)
  # -------------------------------------------
  idx_vessels <- les[
    pts,
    on = c(setNames(pt_col, pt_col), "ev_day" = "idx_day"),
    nomatch = 0L,
    .(.row_id = i..row_id, vessel = vessel)
  ][!is.na(vessel)]
  
  idx_vessels <- unique(idx_vessels, by = c(".row_id", "vessel"))
  
  # -------------------------------------------
  # 2) index 이후 thrombus==1 후보 (pt_id join 후 필터)
  # -------------------------------------------
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
  
  # index_vessels에 포함되는 vessel만 남김
  cand2 <- cand[idx_vessels, on = .(.row_id, vessel), nomatch = 0L]
  
  # -------------------------------------------
  # 3) row_id별 earliest event
  # -------------------------------------------
  if (nrow(cand2) == 0) {
    out <- copy(pts)
    out[, (outcome_col) := 0L]
    out[, (st_date_col) := as.Date(NA)]
    out[, (st_vessel_col) := NA_character_]
    out[, .row_id := NULL]
    out[, idx_day := NULL]
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
  
  return(as.data.frame(out))
}

# 사용:
# temp <- calc_ST_after_index_vessel(df.pts, df.lesion.all)
# table(temp$ST_outcome)
# head(temp[, c("pt_id", "CAG_date", "ST_outcome", "ST_date", "ST_vessel")])


# 사용 예:
# df_pts2 <- calc_ST_after_index(df.pts, df.lesion.all)
# head(df_pts2[, c("pt_id", "CAG_date", "ST_outcome", "ST_date")])

df.pt.events <- as_tibble(calc_ST_after_index_vessel(df.pt.events, all_pcis))


################################################################################
# TV-MI
library(data.table)

get_TVMI_events <- function(
    df_pts, 
    df_lesions,
    pt_col = "pt_id",
    date_col = "CAG_date",
    lesion_col = "lesion_anatomy",
    # MI 여부를 판단하는 컬럼 (is_ACS 혹은 별도 MI 컬럼)
    # 데이터에 따라 "is_ACS" 혹은 "diagnosys" 등으로 수정 필요
    mi_col = "is_ACS", 
    outcome_col = "TVMI_outcome",
    tvmi_date_col = "TVMI_date"
) {
  
  # --- 1. Data Conversion (data.table) ---
  pts <- as.data.table(copy(df_pts))
  les <- as.data.table(copy(df_lesions))
  
  # 날짜 변환
  pts[, idx_day := as.Date(get(date_col))]
  les[, ev_day  := as.Date(get(date_col))]
  
  # --- 2. Anatomy Mapping (세부 위치 -> 혈관 단위) ---
  # 교수님의 ST 코드와 동일한 매핑 로직 적용
  target_LAD <- c("pLAD", "mLAD", "dLAD", "Dg", "septal_br")
  target_LCX <- c("pLCX", "dLCX", "OM", "RI")
  target_RCA <- c("pRCA", "mRCA", "dRCA", "PLB", "PDA")
  target_LM  <- c("LM")
  
  lad_set <- tolower(target_LAD)
  lcx_set <- tolower(target_LCX)
  rca_set <- tolower(target_RCA)
  lm_set  <- tolower(target_LM)
  
  les[, lesion_norm := tolower(trimws(as.character(get(lesion_col))))]
  
  les[, vessel := fcase(
    lesion_norm %chin% lad_set, "LAD",
    lesion_norm %chin% lcx_set, "LCX",
    lesion_norm %chin% rca_set, "RCA",
    lesion_norm %chin% lm_set,  "LM",
    default = NA_character_
  )]
  
  # --- 3. MI Flag 설정 ---
  # 데이터에 따라 1, TRUE, "Yes" 등 다양할 수 있어 유연하게 처리
  les[, mi_flag := {
    x <- get(mi_col)
    if (is.numeric(x)) x == 1
    else if (is.logical(x)) x
    # 문자열인 경우 "1"이나 "TRUE" 등 확인 (필요시 수정)
    else !is.na(x) & (x == 1 | x == "1") 
  }]
  
  # 환자 식별을 위한 row_id
  pts[, .row_id := .I]
  
  # -------------------------------------------
  # 4) Index Vessel 식별 (환자별 Index 날짜에 시술한 혈관 목록)
  # -------------------------------------------
  idx_vessels <- les[
    pts,
    on = c(setNames(pt_col, pt_col), "ev_day" = "idx_day"),
    nomatch = 0L,
    .(.row_id = i..row_id, vessel = vessel)
  ][!is.na(vessel)]
  
  idx_vessels <- unique(idx_vessels, by = c(".row_id", "vessel"))
  
  # -------------------------------------------
  # 5) Follow-up 중 MI 발생 & Vessel 일치 건 찾기
  # -------------------------------------------
  # 조건: Index 날짜 이후(>) 발생한 MI
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
  
  # Index Vessel 목록과 Join (Target Vessel MI 조건 충족)
  tv_mi_events <- mi_candidates[idx_vessels, on = .(.row_id, vessel), nomatch = 0L]
  
  # -------------------------------------------
  # 6) 환자별 최초 TV-MI 선택
  # -------------------------------------------
  if (nrow(tv_mi_events) == 0) {
    # 이벤트가 하나도 없는 경우 빈 결과 반환
    out <- copy(pts)
    out[, (outcome_col) := 0L]
    out[, (tvmi_date_col) := as.Date(NA)]
    out[, c(".row_id", "idx_day") := NULL]
    return(as.data.frame(out))
  }
  
  setorder(tv_mi_events, .row_id, ev_day)
  first_ev <- tv_mi_events[, .SD[1], by = .(.row_id)]
  
  final_ev <- first_ev[, .(.row_id, TVMI_date = ev_day)]
  
  # 원본 데이터에 병합
  out <- merge(pts, final_ev, by = ".row_id", all.x = TRUE)
  out[, (outcome_col) := as.integer(!is.na(TVMI_date))] # 날짜가 있으면 1, 없으면 0
  setnames(out, "TVMI_date", tvmi_date_col)
  
  # 정리
  out[, c(".row_id", "idx_day") := NULL]
  
  return(as.data.frame(out))
}

# 1. TV-MI 계산 (함수 적용)
# 주의: df.lesion.all에 MI 여부를 나타내는 컬럼이 "is_ACS"라고 가정했습니다. 
# 만약 명확한 'MI' 컬럼이 있다면 mi_col = "MI"로 변경하세요.
df.pt.events <- as_tibble(get_TVMI_events(
  df_pts = df.pt.events, 
  df_lesions = all_pcis,
  mi_col = "is_ACS",    # 혹은 "MI"
  outcome_col = "TVMI",
  tvmi_date_col = "TVMI_date"
))





################################################################################
# any revascularization
# temp (pt_id, TLR, TLR_date)가 있을 때
# df에서 특정 pt_id와 일치하는 temp의 records 들 중에 TLR_date가 가장 빠른 날짜를 TLR_date, TLR 중 하나라도 1이면 1, 다 0이면 0으로 설정

get_any_revasc_from_lesion_data <- function(df_pts, df_lesions) {
  # 필요한 컬럼만 선택
  temp <- df_lesions %>%
    select(pt_id, TLR, TLR_date, TVR, TVR_date)
  
  # ---------------------------------------------------------
  # 1. TLR 요약: TLR=1인 건들 중 가장 빠른 날짜 추출
  # ---------------------------------------------------------
  temp_tlr_events <- temp %>%
    filter(TLR == 1) %>% 
    group_by(pt_id) %>%
    summarise(
      first_tlr_date = min(TLR_date, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # ---------------------------------------------------------
  # 2. TVR 요약: TVR=1인 건들 중 가장 빠른 날짜 추출 (Added)
  # ---------------------------------------------------------
  temp_tvr_events <- temp %>%
    filter(TVR == 1) %>% 
    group_by(pt_id) %>%
    summarise(
      first_tvr_date = min(TVR_date, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # ---------------------------------------------------------
  # 3. 메인 데이터프레임(df)에 병합 및 정리
  # ---------------------------------------------------------
  df_updated <- df_pts %>%
    # TLR 정보 Join
    left_join(temp_tlr_events, by = "pt_id") %>%
    # TVR 정보 Join
    left_join(temp_tvr_events, by = "pt_id") %>%
    mutate(
      # --- TLR 처리 ---
      # Join된 날짜가 있으면 이벤트가 있는 것(1), 없으면(NA) 0
      TLR = if_else(!is.na(first_tlr_date), 1, 0),
      TLR_date = as.Date(first_tlr_date), # 날짜 형식 보장
      
      # --- TVR 처리 (New) ---
      # Join된 날짜가 있으면 이벤트가 있는 것(1), 없으면(NA) 0
      TVR = if_else(!is.na(first_tvr_date), 1, 0),
      TVR_date = as.Date(first_tvr_date)  # 날짜 형식 보장
    ) %>%
    # 임시로 만든 first_... 컬럼들 제거
    select(-first_tlr_date, -first_tvr_date)
  
  return(df_updated)
}

df.pt.events.all <- get_any_revasc_from_lesion_data(df.pt.events, df.lesion.all)

# 결과 확인
df.pt.events.all %>% filter(is.na(TLR_date) == FALSE & BP == 0) %>% count() # DP, TLR 393
df.pt.events.all %>% filter(is.na(TLR_date) == FALSE & BP == 1) %>% count() # BP, TLR 165
df.pt.events.all %>% filter(is.na(TVR_date) == FALSE & BP == 0) %>% count() # DP, TVR 693
df.pt.events.all %>% filter(is.na(TVR_date) == FALSE & BP == 1) %>% count() # BP, TVR 301


################################################################################
# 2. Primary Endpoint (TLF) 재계산 함수 (Any MI -> TV-MI 교체)
get_primary_endpoints_revised <- function(df_all_events) {
  df_all_events %>% 
    mutate(
      # 숫자형 변환 (안전장치)
      TLR_n = as.numeric(as.character(TLR)),
      TVMI_n = as.numeric(as.character(TVMI)), # Any MI 대신 TVMI 사용
      CD_n = as.numeric(as.character(cardiac_death)),
      
      # 1. Composite Event Status (하나라도 1이면 1)
      composite_event = pmax(TLR_n, TVMI_n, CD_n, na.rm = TRUE),
      
      # 2. Composite Event Date (가장 빠른 날짜)
      composite_date = pmin(TLR_date, TVMI_date, cardiac_death_date, na.rm = TRUE)
    ) %>%
    mutate(
      # 이벤트가 없는 경우(0) 날짜를 NA로, 날짜가 Inf면 NA로 처리
      composite_event = ifelse(is.na(composite_event), 0, composite_event),
      composite_date = if_else(composite_event == 0 | is.infinite(composite_date), 
                               as.Date(NA), composite_date)
    ) %>%
    select(-ends_with("_n"))
}

# 3. 최종 데이터 생성
df.primary.endpoints <- get_primary_endpoints_revised(df.pt.events.all)

# 결과 확인 (BP vs DP 비교)
cat("=== Primary Endpoint (TLF) Counts ===\n")
df.primary.endpoints %>% group_by(BP) %>% count(composite_event)


# ################################################################################
# # Composite of cardiac death+MI+TLR
# get_primary_endpoints <- function(df_all_events) {
#   return (
#     df_all_events %>% 
#       mutate(
#         # [핵심 수정] Factor -> Character -> Numeric 변환
#         # pmax 계산을 위해 임시로 숫자형 컬럼을 만듭니다 (혹은 덮어씁니다).
#         TLR_n = as.numeric(as.character(TLR)),
#         next_MI_n = as.numeric(as.character(next_MI)),
#         cardiac_death_n = as.numeric(as.character(cardiac_death)),
#         
#         # 1. Composite Event Status
#         # 변환된 숫자형 변수(_n)들을 사용하여 pmax 계산
#         composite_event = pmax(TLR_n, next_MI_n, cardiac_death_n, na.rm = TRUE),
#         
#         # 2. Composite Event Date (날짜는 그대로 사용)
#         composite_date = pmin(TLR_date, next_MI_date, cardiac_death_date, na.rm = TRUE)
#       ) %>%
#       # 3. 데이터 정합성 유지
#       mutate(
#         composite_date = if_else(composite_event == 0, as.Date(NA), composite_date)
#       ) %>%
#       # (선택사항) 계산에 썼던 임시 변수(_n) 제거
#       select(-ends_with("_n"))
#   )
# }
# 
# df.primary.endpoints <- get_primary_endpoints(df.pt.events.all)
# 
# # 결과 확인
# df.primary.endpoints %>% filter(composite_event == 1 & BP == 0) %>% count() # DP, 1' endpoint 1,004
# df.primary.endpoints %>% filter(composite_event == 1 & BP == 1) %>% count() # BP, 1' endpoint 397

################################################################################
# Background data와 join

# index procedure가 ACS인가 확인
get_baseline_info_from_lesion_data <- function(df_lesions) {
  return(
    df_lesions %>% 
      select(pt_id, age, HTN, DM, DL, smoking, prior_MI, prior_PCI, 
             prior_CABG, LVEF, TC, TG, HDL, LDL, BMI, male, AFib, CKD, 
             COPD, CAOD, HF, ESRD, PAOD, prior_stroke, TIA) %>%
      mutate(across(c("pt_id", "age"), as.integer),
             across(c("HTN", "DM", "DL", "smoking", "prior_MI", "prior_PCI", "prior_CABG", 
                      "male", "AFib", "CKD", "COPD", "CAOD", "HF", "ESRD", "PAOD", 
                      "prior_stroke", "TIA"), as.factor)) %>%
      distinct()
  )
}

df.pt.baseline.all <- get_baseline_info_from_lesion_data(df.lesion.all)

df.pts <- df.primary.endpoints %>% left_join(df.pt.baseline.all, by="pt_id")

temp <- df.cag.all %>% select(cath_no, is_ACS, ACS_type, num_CAOD)
df.pts <- df.pts %>% left_join(temp, by = "cath_no")

library(lubridate)
# procedure year 추가
df.pts <- df.pts %>% 
  mutate(procedure_year = year(CAG_date),
         fu_days = as.numeric(last_fu_date, CAG_date))

################################################################################
# 연구기간 변경
df.pts %>% filter(procedure_year<2010) %>% count() # 897 (yr<2010)

df.pts <- df.pts %>% 
  filter(procedure_year >= 2010) # 8,392 (yr>=2010)

################################################################################
# ERA 구분
# 1) 2010-2013 (early-generation era): characterized by predominant use of BioMatrix and Nobori stents; 
# 2) 2014-2016 (modern-generation era): marked by the introduction of Orsiro and SYNERGY stents; and 
# 3) 2017-2021 (contemporary era): featuring diversified newer-generation biodegradable polymer platforms 
#                                  including SYNERGY, DESyne X2, GENOSS DES, and Firehawk. 
df.pts <- df.pts %>%
  # 1. 분석 대상 연도(2010~2021) 필터링 (선택사항, 데이터 정제를 위해 권장)
  filter(procedure_year >= 2010 & procedure_year <= 2021) %>% 
  
  # 2. Era 구분 (cut 함수 사용)
  mutate(era = cut(procedure_year,
                   # Breaks: 2009초과~2013이하 / 2013초과~2016이하 / 2016초과~2021이하
                   breaks = c(2009, 2013, 2016, 2021), 
                   labels = c("2010-2013", "2014-2016", "2017-2021")),
         # Factor로 변환 (cut이 이미 factor를 반환하지만 명시적으로 작성)
         era = as.factor(era))

# Result
# 2010-2013: 2,804
# 2014-2016: 2,093
# 2017-2021: 3,495
# table(df.pts$procedure_year, df.pts$era)
df.pts %>% count(era)


