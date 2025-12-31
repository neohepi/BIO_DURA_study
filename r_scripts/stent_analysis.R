library(dplyr)
library(tidyr)
library(ggplot2)

# 만약 패키지가 없다면 설치하고 로드하는 안전장치 추가
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

library(patchwork) # 그래프 합치기용
library(RColorBrewer)

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)

analyze_stent_distribution <- function(data, 
                                       date_col = "CAG_date", 
                                       group_col = "BP", 
                                       stent_col = "stent_name", 
                                       top_n_cnt = 10,
                                       color_palette = "Set3") {
  
  # 1. 컬럼 이름 유효성 검사
  if (!all(c(date_col, group_col, stent_col) %in% names(data))) {
    stop("지정하신 컬럼 이름이 데이터프레임에 존재하지 않습니다.")
  }
  
  # 2. 데이터 전처리
  stent_analysis_data <- data %>%
    mutate(
      Year = as.numeric(format(as.Date(.data[[date_col]]), "%Y"))
    ) %>% 
    filter(Year >= 2010) %>%
    mutate(
      Era = case_when(
        Year <= 2013 ~ "Study Period 1 (2010-2013)",
        Year <= 2016 ~ "Study Period 2 (2014-2016)",
        TRUE         ~ "Study Period 3 (2017-2021)"
      ),
      Group = ifelse(.data[[group_col]] == 1, "BP-DES", "DP-DES"),
      stent_clean = toupper(trimws(.data[[stent_col]])) 
    ) %>%
    filter(!is.na(Year))
  
  # 3. 전체 리스트 집계
  full_stent_list <- stent_analysis_data %>%
    count(Era, Group, stent_clean, name = "Count") %>%
    group_by(Era, Group) %>%
    mutate(
      Total_In_Group = sum(Count),
      Percent = round(Count / Total_In_Group * 100, 2)
    ) %>%
    arrange(Era, Group, desc(Count))
  
  # 4. Top stents 선택
  top_stents <- stent_analysis_data %>%
    count(stent_clean, sort = TRUE) %>%
    top_n(top_n_cnt, n) %>%
    pull(stent_clean)
  
  plot_data <- stent_analysis_data %>%
    mutate(
      Stent_Label = ifelse(stent_clean %in% top_stents, stent_clean, "Others")
    ) %>%
    count(Era, Group, Stent_Label) %>%
    group_by(Era, Group) %>%
    mutate(Percent = n / sum(n) * 100)
  
  plot_data$Era_Label <- gsub(" \\(", "\n(", plot_data$Era)
  
  # 5. 각 막대별로 정렬 및 누적 계산
  plot_data <- plot_data %>%
    group_by(Era, Group) %>%
    arrange(Era, Group, 
            Stent_Label == "Others",
            desc(Percent)) %>%
    mutate(
      ymax = cumsum(Percent),
      ymin = lag(ymax, default = 0),
      y_position = (ymin + ymax) / 2,
      draw_order = row_number()
    ) %>%
    ungroup()
  
  # 6. Stent 순서 (전체 빈도 기준)
  stent_color_order <- plot_data %>%
    group_by(Stent_Label) %>%
    summarise(Total_Percent = sum(Percent), .groups = 'drop') %>%
    arrange(desc(Total_Percent)) %>%
    pull(Stent_Label)
  
  if ("Others" %in% stent_color_order) {
    stent_color_order <- c(setdiff(stent_color_order, "Others"), "Others")
  }
  
  plot_data$Stent_Label <- factor(plot_data$Stent_Label, levels = stent_color_order)
  
  # ========================================================================
  # 7. 색상 팔레트 선택 (R 표준 라이브러리)
  # ========================================================================
  
  colourCount <- length(unique(plot_data$Stent_Label))
  
  # RColorBrewer 팔레트
  if (color_palette %in% c("Set1", "Set2", "Set3", "Pastel1", "Pastel2", 
                           "Paired", "Dark2", "Accent")) {
    # Qualitative palettes
    max_colors <- brewer.pal.info[color_palette, "maxcolors"]
    if (colourCount <= max_colors) {
      colors <- brewer.pal(max_colors, color_palette)[1:colourCount]
    } else {
      colors <- colorRampPalette(brewer.pal(max_colors, color_palette))(colourCount)
    }
    
  } else if (color_palette %in% c("Spectral", "RdYlBu", "RdYlGn", "RdBu", 
                                  "PiYG", "PRGn", "BrBG", "PuOr")) {
    # Diverging palettes
    max_colors <- brewer.pal.info[color_palette, "maxcolors"]
    colors <- colorRampPalette(brewer.pal(max_colors, color_palette))(colourCount)
    
  } else if (color_palette %in% c("Blues", "Greens", "Reds", "Oranges", 
                                  "Purples", "Greys", "YlOrRd", "YlGnBu")) {
    # Sequential palettes
    max_colors <- brewer.pal.info[color_palette, "maxcolors"]
    colors <- colorRampPalette(brewer.pal(max_colors, color_palette))(colourCount)
    
  } else if (color_palette == "viridis") {
    # Viridis palette
    colors <- viridis(colourCount, option = "D")
    
  } else if (color_palette == "plasma") {
    colors <- viridis(colourCount, option = "C")
    
  } else if (color_palette == "inferno") {
    colors <- viridis(colourCount, option = "B")
    
  } else if (color_palette == "magma") {
    colors <- viridis(colourCount, option = "A")
    
  } else if (color_palette == "cividis") {
    colors <- viridis(colourCount, option = "E")
    
  } else if (color_palette == "turbo") {
    colors <- viridis(colourCount, option = "H")
    
  } else if (color_palette == "rainbow") {
    # Base R rainbow
    colors <- rainbow(colourCount)
    
  } else if (color_palette == "heat") {
    # Base R heat.colors
    colors <- heat.colors(colourCount)
    
  } else if (color_palette == "terrain") {
    # Base R terrain.colors
    colors <- terrain.colors(colourCount)
    
  } else if (color_palette == "topo") {
    # Base R topo.colors
    colors <- topo.colors(colourCount)
    
  } else if (color_palette == "cm") {
    # Base R cm.colors
    colors <- cm.colors(colourCount)
    
  } else {
    # Default: Set3
    colors <- colorRampPalette(brewer.pal(12, "Set3"))(colourCount)
  }
  
  # Others는 회색으로
  if ("Others" %in% levels(plot_data$Stent_Label)) {
    others_idx <- which(levels(plot_data$Stent_Label) == "Others")
    colors[others_idx] <- "#CCCCCC"
  }
  
  # 8. 그래프 생성
  plot_data <- plot_data %>%
    arrange(Era, Group, draw_order)
  
  p <- ggplot(plot_data) +
    
    geom_rect(aes(xmin = as.numeric(factor(Group)) - 0.35,
                  xmax = as.numeric(factor(Group)) + 0.35,
                  ymin = ymin,
                  ymax = ymax,
                  fill = Stent_Label),
              color = "white", linewidth = 0.3) +
    
    geom_text(data = plot_data %>% filter(Percent > 5),
              aes(x = Group, y = y_position, 
                  label = paste0(round(Percent, 1), "%")),
              size = 3.5, fontface = "bold", color = "white") +
    
    facet_wrap(~Era_Label) +
    
    scale_fill_manual(values = colors) +
    scale_x_discrete() +
    
    labs(
      title = "Distribution of Stent Types by Study Period",
      subtitle = paste("Top", top_n_cnt, "stents ordered by percentage (largest at bottom)"),
      y = "Percentage (%)",
      x = "",
      fill = "Stent Name"
    ) +
    
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 9, face = "bold"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.key.size = unit(0.6, "cm"),
      strip.text = element_text(size = 11, face = "bold", color = "#333333"),
      strip.background = element_rect(fill = "#F5F5F5", color = NA),
      axis.text.x = element_text(size = 11, face = "bold", color = "#333333"),
      axis.text.y = element_text(size = 9, color = "#666666"),
      axis.title.y = element_text(size = 10, face = "bold", color = "#333333"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.3),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#666666"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(p)
  
  # 9. Top 5 Summary
  cat("\n=== Top 5 Stents by Era & Group ===\n")
  top5_summary <- plot_data %>%
    group_by(Era, Group) %>%
    filter(Stent_Label != "Others") %>%
    arrange(Era, Group, desc(Percent)) %>%
    slice_head(n = 5) %>%
    select(Era, Group, Stent_Label, Percent) %>%
    mutate(Rank = row_number())
  
  print(as.data.frame(top5_summary))
  
  cat("\n=== Color Palette Used:", color_palette, "===\n")
  
  return(invisible(list(
    plot = p,
    full_data = full_stent_list,
    plot_data = plot_data,
    top5 = top5_summary,
    colors = colors
  )))
}

# 기본값 사용 (컬럼명이 CAG_date, BP, stent_name 인 경우)
analyze_stent_distribution(df.lesion.all, top_n_cnt = 10, color_palette = "Set1")


# 옵션 변경 예시 (만약 컬럼명이 다르거나 Top 10만 보고 싶다면)
# analyze_stent_distribution(
#   data = my_other_data, 
#   date_col = "procedure_date",   # 날짜 컬럼명 지정
#   group_col = "group_var",       # 그룹 컬럼명 지정
#   stent_col = "device_name",     # 스텐트 컬럼명 지정
#   top_n_cnt = 10                 # 상위 10개만 표시
# )