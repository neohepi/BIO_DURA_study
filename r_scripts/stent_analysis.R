library(dplyr)
library(tidyr)
library(ggplot2)

# 만약 패키지가 없다면 설치하고 로드하는 안전장치 추가
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

install.packages("scales")

library(patchwork) # 그래프 합치기용
library(RColorBrewer)
library(viridis)
library(scales) # for percent format
library(dplyr)
library(ggplot2)
library(RColorBrewer)

analyze_stent_distribution_pretty <- function(data, 
                                              date_col = "CAG_date", 
                                              group_col = "BP", 
                                              stent_col = "stent_name", 
                                              top_n_cnt = 10,
                                              color_palette = "Set3") {
  
  # 1. 컬럼 유효성 검사
  if (!all(c(date_col, group_col, stent_col) %in% names(data))) {
    stop("Error: 지정하신 컬럼 이름이 데이터프레임에 존재하지 않습니다.")
  }
  
  # 2. 데이터 전처리
  stent_analysis_data <- data %>%
    mutate(Year = as.numeric(format(as.Date(.data[[date_col]]), "%Y"))) %>% 
    filter(Year >= 2010) %>%
    mutate(
      Era = case_when(
        Year <= 2013 ~ "Study Period 1\n(2010-2013)",
        Year <= 2017 ~ "Study Period 2\n(2014-2017)", 
        TRUE         ~ "Study Period 3\n(2018-2021)"
      ),
      Group = ifelse(.data[[group_col]] == 1, "BP-DES", "DP-DES"),
      stent_clean = toupper(trimws(.data[[stent_col]])) 
    ) %>%
    filter(!is.na(Year))
  
  # 3. Top Stents 선정
  top_stents <- stent_analysis_data %>%
    count(stent_clean, sort = TRUE) %>%
    top_n(top_n_cnt, n) %>%
    pull(stent_clean)
  
  # 4. Plot용 데이터 집계
  plot_data <- stent_analysis_data %>%
    mutate(
      Stent_Label = ifelse(stent_clean %in% top_stents, stent_clean, "Others")
    ) %>%
    count(Era, Group, Stent_Label) %>%
    group_by(Era, Group) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ungroup()
  
  # 5. Factor Level 순서 정리
  stent_order_list <- plot_data %>%
    filter(Stent_Label != "Others") %>%
    group_by(Stent_Label) %>%
    summarise(total = sum(n)) %>%
    arrange(desc(total)) %>%
    pull(Stent_Label)
  
  # Others를 맨 마지막 레벨로 지정
  final_levels <- c(stent_order_list, "Others")
  plot_data$Stent_Label <- factor(plot_data$Stent_Label, levels = final_levels)
  
  # =========================================================
  # 6. 색상 설정 (수정됨: Set3 팔레트 내의 회색 제거 로직 추가)
  # =========================================================
  
  # Top Stent 개수 (Others 제외)
  n_top_stents <- length(stent_order_list)
  
  # 기본 Set3 (12색) 가져오기
  raw_palette <- brewer.pal(12, "Set3")
  
  # Set3의 9번째 색상이 회색(#D9D9D9)이므로 이를 목록에서 제거
  # 이렇게 하면 스텐트에는 절대 회색이 배정되지 않음
  clean_palette <- raw_palette[-9] 
  
  # 스텐트 개수에 맞춰 색상 생성 (회색 빠진 팔레트 사용)
  if (n_top_stents <= length(clean_palette)) {
    stent_colors <- clean_palette[1:n_top_stents]
  } else {
    stent_colors <- colorRampPalette(clean_palette)(n_top_stents)
  }
  
  # 이름과 색상 매핑
  my_colors <- setNames(stent_colors, stent_order_list)
  
  # Others에만 별도로 '진한 회색' 지정
  my_colors["Others"] <- "#D3D3D3"
  
  # =========================================================
  
  # 7. ggplot 그리기
  p <- ggplot(plot_data, aes(x = Group, y = Percent, fill = Stent_Label)) +
    
    geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.3) +
    
    geom_text(data = subset(plot_data, Percent > 4),
              aes(label = paste0(round(Percent, 1), "%")),
              position = position_stack(vjust = 0.5),
              size = 3, 
              fontface = "bold", 
              color = "black") + 
    
    facet_wrap(~Era) +
    
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
    
    labs(
      #title = "Distribution of Stent Types by Study Period",
      #subtitle = paste("Top", top_n_cnt, "stents vs. Others"),
      y = "Percentage (%)",
      x = NULL,
      fill = "Stent Name"
    ) +
    
    theme_classic(base_size = 14) + 
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 12),
      
      axis.text.x = element_text(face = "bold", color = "black", size = 11, margin = margin(t = 5)),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      
      strip.background = element_rect(fill = "#F0F0F0", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 11),
      
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.5, "cm")
    )
  
  print(p)
  
  return(invisible(list(plot = p, data = plot_data)))
}

analyze_stent_distribution_pretty(df.lesion.all, top_n_cnt = 10)


# 옵션 변경 예시 (만약 컬럼명이 다르거나 Top 10만 보고 싶다면)
# analyze_stent_distribution(
#   data = my_other_data, 
#   date_col = "procedure_date",   # 날짜 컬럼명 지정
#   group_col = "group_var",       # 그룹 컬럼명 지정
#   stent_col = "device_name",     # 스텐트 컬럼명 지정
#   top_n_cnt = 10                 # 상위 10개만 표시
# )

################################################################################
# Table로 출력
library(dplyr)
#install.packages("knitr")
library(knitr) # 표 출력을 위해 필요 (없으면 install.packages("knitr"))

get_top5_stents_table <- function(data, 
                                  date_col = "CAG_date", 
                                  group_col = "BP", 
                                  stent_col = "stent_name") {
  
  # 1. 데이터 전처리 (시기 구분 및 그룹 정의)
  # 기존 로직 유지: 2010-2013 / 2014-2017 / 2018-2021
  summary_data <- data %>%
    mutate(Year = as.numeric(format(as.Date(.data[[date_col]]), "%Y"))) %>% 
    filter(Year >= 2010) %>%
    mutate(
      Era = case_when(
        Year <= 2013 ~ "Period 1 (2010-2013)",
        Year <= 2017 ~ "Period 2 (2014-2017)", 
        TRUE         ~ "Period 3 (2018-2021)"
      ),
      Group = ifelse(.data[[group_col]] == 1, "BP-DES", "DP-DES"),
      Stent_Name = toupper(trimws(.data[[stent_col]])) # 대문자 통일 및 공백 제거
    ) %>%
    filter(!is.na(Year))
  
  # 2. 집계 및 랭킹 산정
  # 각 Era & Group 내에서 스텐트별 개수 세기 -> 백분율 계산 -> 랭킹 매기기
  rank_data <- summary_data %>%
    group_by(Era, Group, Stent_Name) %>%
    summarise(Count = n(), .groups = "drop_last") %>% # Era, Group별로 그룹화 유지
    mutate(
      Total_in_Subgroup = sum(Count),
      Percent = (Count / Total_in_Subgroup) * 100
    ) %>%
    arrange(Era, Group, desc(Count)) %>% # 개수 많은 순 정렬
    mutate(Rank = row_number()) %>%      # 순위 부여 (1, 2, 3...)
    filter(Rank <= 5) %>%                # Top 5만 남기기
    ungroup()
  
  # 3. 출력용 표 다듬기 (가독성 향상)
  final_table <- rank_data %>%
    select(Era, Group, Rank, Stent_Name, Count, Percent) %>%
    mutate(
      Percent = sprintf("%.1f%%", Percent) # 소수점 1자리 + % 기호 붙이기
    )
  
  # 4. 결과 출력
  cat("\n======================================================\n")
  cat("   Top 5 Stents by Study Period & Group (BP vs. DP)\n")
  cat("======================================================\n\n")
  
  # Period 별로 나누어 출력 (가독성을 위해)
  periods <- unique(final_table$Era)
  
  for (p in periods) {
    cat(paste0(">>> ", p, " <<<\n"))
    
    sub_df <- final_table %>% filter(Era == p)
    
    print(kable(sub_df %>% select(-Era), # Era 컬럼은 제목에 있으니 제외
                format = "simple", 
                align = "cclrc")) # 정렬 (Center, Center, Left, Right, Center)
    cat("\n")
  }
  
  # 5. 데이터 반환 (변수에 저장용)
  return(invisible(final_table))
}

# === 실행 예시 ===
# 결과를 변수에 저장하면 엑셀로 내보내기도 쉽습니다.
top5_result <- get_top5_stents_table(df.lesion.all)

# csv 저장 필요시:
# write.csv(top5_result, "top5_stents_list.csv", row.names = FALSE)











