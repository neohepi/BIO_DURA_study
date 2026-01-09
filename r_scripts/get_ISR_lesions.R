# lesion data에서 ISR data를 재계산
# ISR_lesion == 1이면 기존 ISR lesion -> 이후에 동일한 lesion_anatomy에 lesion procedure가 있는지 확인
# ISR_lesion == 0이면 이후에 동일한 lesion_anatomy에 있으면 이를 ISR lesion으로 설정

# write_xlsx(df.lesion.all, "lesion_data.xlsx")

library(tidyverse)

df <- df.lesion.all

# 2. ISR 재정의 및 TVR 변수 정리
df_analysis <- df %>%
  # 날짜 형식 변환 (엑셀 로드 시 문자로 인식될 경우를 대비)
  mutate(
    CAG_date = as.Date(CAG_date),
    TVR_date = as.Date(TVR_date),
    last_fu_date = as.Date(last_fu_date)
  ) %>%
  # 환자 및 병변 위치별로 정렬 (시간 순서 중요)
  arrange(pt_id, lesion_anatomy, CAG_date) %>%
  group_by(pt_id, lesion_anatomy) %>%
  mutate(
    # --- [1. ISR Lesion 재정의] ---
    # 조건 1: 원래 ISR_lesion이 1인 경우
    # 조건 2: (0이었더라도) 해당 혈관의 2번째 이상 시술인 경우 (De Novo 후 재발)
    is_ISR_lesion = if_else(ISR_lesion == 1 | row_number() > 1, 1, 0),
    
    # (참고) 몇 번째 시술인지 표시
    proc_seq = row_number()
  ) %>%
  ungroup() %>%
  # --- [2. TVR Outcome 계산 (Survival Analysis용)] ---
  mutate(
    # Event: TVR 컬럼이 1이면 Event(1), 아니면 Censored(0)
    # (데이터에 따라 TVR 컬럼이 'Yes'/'No'라면 ifelse로 수정 필요)
    status_TVR = as.numeric(TVR), 
    
    # Time: TVR이 발생했으면 (TVR날짜 - 시술날짜), 안 했으면 (마지막추적일 - 시술날짜)
    time_TVR = if_else(status_TVR == 1,
                       as.numeric(TVR_date - CAG_date),
                       as.numeric(last_fu_date - CAG_date))
  )

# 3. 결과 확인 (ISR 병변만 따로 보고 싶다면 filter 사용)
# ISR 병변들의 TVR 발생 현황 예시

isr_summary <- df_analysis %>%
  filter(is_ISR_lesion == 1) 
#%>%
  #select(pt_id, lesion_anatomy, CAG_date, is_ISR_lesion, TVR, time_TVR)

isr_summary <- isr_summary %>%
  select(pt_id, lesion_anatomy, CAG_date, age, stent_diameter, stent_length, stent_name, BP, BMI, 
         LM, L_LAD, L_LCX, L_RCA, pre_DS_percent, pre_MLD, post_DS_percent, post_MLD, L_length, ref_diameter,
         thrombus.x, bifurcation, CTO, os_dz, ACS, ACS_type, num_CAOD, HTN, DM, DL, smoking, prior_MI,
         prior_PCI, prior_CABG, LVEF, TC, TG, HDL, LDL, male, AFib, CKD, COPD, CAOD, HF, ESRD, PAOD,
         prior_stroke, TIA, last_fu_date, fu_days, death_date, cardiac_death, ACC_type, procedure_year, is_ACS,
         is_ISR_lesion, ISR_outcome, ST_outcome, TLR, TLR_date, TVR, TVR_date)

isr_summary <- as_tibble(calc_ST_after_index_vessel(isr_summary, all_pcis))
isr_summary <- as_tibble(get_TVMI_events(
  df_pts = isr_summary, 
  df_lesions = all_pcis,
  mi_col = "is_ACS",    # 혹은 "MI"
  outcome_col = "TVMI",
  tvmi_date_col = "TVMI_date"
))

write_xlsx(isr_summary, "ISR_lesion_all.xlsx")

# isr_summary에서 join key를 기준으로 중복을 제거하여 유니크한 '참조표'로 만듭니다.
# (필요한 Outcome 변수들만 select에 포함시키세요)
# isr_summary_unique <- isr_summary %>%
#   select(pt_id, lesion_anatomy, CAG_date, is_ISR_lesion, TVR, time_TVR) %>%
#   distinct(pt_id, lesion_anatomy, CAG_date, .keep_all = TRUE)


# 이제 join을 수행하면 Warning 없이 1:N 매칭이 됩니다. (하나의 Outcome이 해당 병변의 모든 스텐트 행에 붙음)
df.lesion.isr <- isr_summary_unique %>%
  left_join(df.lesion.all, by = c("pt_id", "lesion_anatomy", "CAG_date"))


era_check <- function(df) {
  return (df %>%
    # 1. 분석 대상 연도(2010~2021) 필터링 (선택사항, 데이터 정제를 위해 권장)
    filter(procedure_year >= 2010 & procedure_year <= 2021) %>% 
    
    # 2. Era 구분 (cut 함수 사용)
    mutate(era = cut(procedure_year,
                     # Breaks: 2009초과~2013이하 / 2013초과~2016이하 / 2016초과~2021이하
                     breaks = c(2009, 2013, 2016, 2021), 
                     labels = c("2010-2013", "2014-2016", "2017-2021")),
           # Factor로 변환 (cut이 이미 factor를 반환하지만 명시적으로 작성)
           era = as.factor(era)))
}

df.lesion.isr <- df.lesion.isr %>%
  mutate(is_complex_lesion = case_when(
    # "B2" 또는 "C"인 경우 1
    ACC_type %in% c("B2", "C") ~ 1,
    
    # 그 외 (A, B1 등)는 0
    TRUE ~ 0
  ))

df.lesion.isr <- era_check(df.lesion.isr)

# 1. TLR 데이터 전처리 (NA 제거 및 숫자형 변환)
df.lesion.isr <- df.lesion.isr %>%
  # TLR 값이 있는 행만 남김 (NA 제거)
  filter(!is.na(TLR)) %>%
  # 혹시 모를 문자/Factor 형식을 숫자(0, 1)로 확실히 변환
  mutate(TLR = as.numeric(as.character(TLR)))

df.isr.index <- get_index_proc_from_lesion_data(df.lesion.isr)

df.isr.event <- get_first_MI_event_from_index_proc(df.isr.index, df.lesion.all, df.cag.all)
df.isr.event <- as_tibble(calc_ST_after_index_vessel(df.isr.event, all_pcis))

temp <- get_death_events_from_lesion(df.lesion.all)

df.isr.event <- df.isr.event %>% # left_join with death events (cardiac and any death)
  left_join(temp, by="pt_id")

df.isr.event <- as_tibble(get_TVMI_events(
  df_pts = df.isr.event, 
  df_lesions = all_pcis,
  mi_col = "is_ACS",    # 혹은 "MI"
  outcome_col = "TVMI",
  tvmi_date_col = "TVMI_date"
))

df.isr.event <- get_any_revasc_from_lesion_data(df.isr.event, df.lesion.all)



df.pt.isr <- df.isr.event %>% left_join(df.pt.baseline.all, by="pt_id")

temp <- df.cag.all %>% select(cath_no, is_ACS, ACS_type, num_CAOD)
df.pt.isr <- df.pt.isr %>% left_join(temp, by = "cath_no")

library(lubridate)
# procedure year 추가
df.pt.isr <- df.pt.isr %>% 
  mutate(procedure_year = year(CAG_date),
         fu_days = as.numeric(last_fu_date, CAG_date))

df.pt.isr <- df.pt.isr %>%
  # 1. 분석 대상 연도(2010~2021) 필터링 (선택사항, 데이터 정제를 위해 권장)
  filter(procedure_year >= 2010 & procedure_year <= 2021) %>% 
  
  # 2. Era 구분 (cut 함수 사용)
  mutate(era = cut(procedure_year,
                   # Breaks: 2009초과~2013이하 / 2013초과~2016이하 / 2016초과~2021이하
                   breaks = c(2009, 2013, 2016, 2021), 
                   labels = c("2010-2013", "2014-2016", "2017-2021")),
         # Factor로 변환 (cut이 이미 factor를 반환하지만 명시적으로 작성)
         era = as.factor(era))
df.pt.isr %>% left_join(df.lesion.isr %>% select(pt_id, CAG_date, stent_name, stent_diameter, stent_length),
                        by=c("pt_id", "CAG_date"))
# Result
# 2010-2013: 2,804
# 2014-2016: 2,093
# 2017-2021: 3,495
# table(df.pts$procedure_year, df.pts$era)
df.pt.isr %>% count(era)

write_xlsx(df.pt.isr, "ISR_lesion_data.xlsx")
