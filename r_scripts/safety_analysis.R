# 1. matched_data에 'death_30d' 변수를 먼저 생성하여 업데이트(저장)합니다.
matched_data <- matched_data %>%
  mutate(
    time_to_death = as.numeric(death_date - CAG_date),
    # 사망했고(1), 사망 시점이 30일 이내인 경우만 1, 나머지는 0
    # (NA 처리 안전장치: time_to_death가 NA여도 0이 되도록)
    death_30d = if_else(any_death == 1 & !is.na(time_to_death) & time_to_death <= 30, 1, 0)
  )

# 2. 결과 요약 (생성된 변수 활용)
result_30d_mortality <- matched_data %>%
  group_by(BP) %>%
  summarize(
    n = n(),
    n_deaths = sum(death_30d, na.rm = TRUE),
    rate = n_deaths / n * 100
  )
print(result_30d_mortality)

# 3. Statistical Test (이제 matched_data 안에 death_30d가 있으므로 작동합니다)
table_30d <- table(matched_data$BP, matched_data$death_30d)

print(table_30d)

# Chi-square test (혹은 Fisher's exact test)
chisq_result <- chisq.test(table_30d)
print(chisq_result)

# 만약 빈도가 적어 Fisher test가 필요하다면:
# fisher.test(table_30d)

# Odds ratio
library(epitools)
oddsratio(table_30d)

################################################################################
# Competing risk analysis with death
install.packages("cmprsk" ,"tidycmprsk")
library(cmprsk)
library(tidycmprsk)

library(survival)
library(tidycmprsk)
library(dplyr)

# 1. Competing Risk용 데이터 전처리
matched_data_cr <- matched_data %>%
  mutate(
    # (A) 사망까지의 시간 계산 (사망 안 했으면 fu_days)
    time_to_death_clean = if_else(any_death == 1, as.numeric(death_date - CAG_date), as.numeric(fu_days)),
    time_to_death_clean = if_else(time_to_death_clean <= 0, 0.1, time_to_death_clean), # 0 보정
    
    # (B) 365일 이내 발생 여부 확인 (1년 시점 Landmark)
    # MI가 1년 내 발생했으면 그 날짜, 아니면 NA
    t_mi_1yr = if_else(next_MI == 1 & time_to_MI <= 365, time_to_MI, NA_real_),
    # 사망이 1년 내 발생했으면 그 날짜, 아니면 NA
    t_death_1yr = if_else(any_death == 1 & time_to_death_clean <= 365, time_to_death_clean, NA_real_),
    
    # (C) Event Type 결정 (Time to FIRST Event)
    # MI와 사망이 둘 다 기록된 경우, 날짜가 더 빠른 것을 우선시함
    final_event = case_when(
      # 1. 둘 다 1년 내 발생 -> 더 빨리 온 것이 진짜 이벤트
      !is.na(t_mi_1yr) & !is.na(t_death_1yr) & t_mi_1yr <= t_death_1yr ~ "MI",
      !is.na(t_mi_1yr) & !is.na(t_death_1yr) & t_death_1yr < t_mi_1yr ~ "Death",
      
      # 2. 하나만 발생
      !is.na(t_mi_1yr) ~ "MI",
      !is.na(t_death_1yr) ~ "Death",
      
      # 3. 아무 일도 없음 (Censored)
      TRUE ~ "Censored"
    ),
    
    # (D) Time 결정 (결정된 이벤트의 시간 사용)
    final_time = case_when(
      final_event == "MI" ~ t_mi_1yr,
      final_event == "Death" ~ t_death_1yr,
      # Censored인 경우: 1년(365일) 혹은 그 전 추적종료일 중 빠른 날
      TRUE ~ pmin(fu_days, 365, na.rm = TRUE)
    ),
    
    # [핵심] Factor 변환 (반드시 Censored가 첫 번째 레벨이어야 함)
    event_factor = factor(final_event, levels = c("Censored", "MI", "Death"))
  )

# 2. Cumulative Incidence Function 실행
cif <- cuminc(Surv(final_time, event_factor) ~ BP, 
              data = matched_data_cr)

# 3. 결과 출력 (Table)
print(cif)

# 1. 시각화 전용 패키지 설치 및 로드
if(!require("ggsurvfit")) install.packages("ggsurvfit")
library(ggsurvfit)

# 2. 그래프 그리기 (옵션 조금 더 추가해서 예쁘게)
cif %>%
  ggcuminc(outcome = "MI") +             # 관심 Outcome 지정
  labs(
    title = "Cumulative Incidence of MI (1-year Landmark)",
    subtitle = "Competing Risk Analysis (Death as competing event)",
    x = "Days since PCI",
    y = "Cumulative Incidence Probability"
  ) +
  add_confidence_interval() +            # 95% 신뢰구간 (음영)
  add_risktable(                         # Risk Table 추가
    risktable_stats = "n.risk"           # 남은 환자 수만 표시
  ) +
  add_pvalue(caption = "Gray's Test") +  # P-value 자동 추가 (중요!)
  scale_ggsurvfit()                      # 깔끔한 테마 적용
