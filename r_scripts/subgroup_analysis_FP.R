################################################################################
# Subgroup Forest Plot - 기존 분석 결과 사용
################################################################################

library(ggplot2)
library(dplyr)

################################################################################
# Function: 기존 subgroup_results를 forest plot 형식으로 변환
################################################################################

convert_to_forest_format <- function(subgroup_results) {
  
  # 결과 저장용
  forest_data <- data.frame()
  
  # 1. Overall 추가
  overall_row <- subgroup_results %>% 
    filter(Category == "All Patients")
  
  if (nrow(overall_row) > 0) {
    forest_data <- rbind(forest_data, data.frame(
      Subgroup = "Overall",
      Category = "All Patients",
      HR = overall_row$HR,
      CI_lower = overall_row$CI_lower,
      CI_upper = overall_row$CI_upper,
      P_value = overall_row$P_value,
      P_interaction = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  # 2. 각 Subgroup 처리
  subgroups <- unique(subgroup_results$Subgroup[subgroup_results$Subgroup != "Overall"])
  
  for (sg in subgroups) {
    sg_data <- subgroup_results %>% filter(Subgroup == sg)
    
    # P_interaction 값 추출 (첫 행에서)
    p_int <- sg_data$P_interaction[1]
    
    # Header 행 추가
    forest_data <- rbind(forest_data, data.frame(
      Subgroup = sg,
      Category = sg,  # Header: Subgroup == Category
      HR = NA,
      CI_lower = NA,
      CI_upper = NA,
      P_value = NA,
      P_interaction = p_int,
      stringsAsFactors = FALSE
    ))
    
    # Category 행 추가
    for (i in 1:nrow(sg_data)) {
      forest_data <- rbind(forest_data, data.frame(
        Subgroup = sg,
        Category = sg_data$Category[i],
        HR = sg_data$HR[i],
        CI_lower = sg_data$CI_lower[i],
        CI_upper = sg_data$CI_upper[i],
        P_value = sg_data$P_value[i],
        P_interaction = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(forest_data)
}

################################################################################
# Function: Create Publication-Style Forest Plot
################################################################################

create_forest_plot <- function(data,
                               x_breaks = c(0.5, 0.75, 1.0, 1.5, 2.0),
                               x_limits = c(0.35, 2.5),
                               filename = NULL) {
  
  # 데이터 준비
  n_rows <- nrow(data)
  
  plot_data <- data %>%
    mutate(
      row_id = row_number(),
      y_pos = n_rows - row_id + 1,
      
      is_header = (Category == Subgroup | Category == "All Patients"),
      is_overall = (Category == "All Patients"),
      
      display_label = case_when(
        Category == "All Patients" ~ "Overall",
        Category == Subgroup ~ Category,
        TRUE ~ paste0("  ", Category)
      ),
      
      hr_text = if_else(
        !is.na(HR),
        sprintf("%.2f (%.2f-%.2f)", HR, CI_lower, CI_upper),
        ""
      ),
      
      p_int_text = case_when(
        !is.na(P_interaction) & is_header & !is_overall ~ 
          sprintf("P int = %.3f", P_interaction),
        TRUE ~ ""
      ),
      
      label_color = if_else(is_overall, "red", "black")
    )
  
  # 플롯 생성
  p <- ggplot(plot_data, aes(y = y_pos)) +
    
    # 배경 줄무늬
    geom_rect(
      data = plot_data %>% filter(row_id %% 2 == 0),
      aes(ymin = y_pos - 0.5, ymax = y_pos + 0.5),
      xmin = -Inf, xmax = Inf,
      fill = "gray96", inherit.aes = FALSE
    ) +
    
    # Reference line
    geom_vline(xintercept = 1, linetype = "dashed", color = "darkgreen", linewidth = 0.7) +
    
    # CI lines (NA 제외)
    geom_segment(
      data = plot_data %>% filter(!is.na(HR)),
      aes(x = CI_lower, xend = CI_upper, yend = y_pos),
      linewidth = 0.7, color = "darkblue"
    ) +
    
    # Points - Overall (diamond)
    geom_point(
      data = plot_data %>% filter(is_overall & !is.na(HR)),
      aes(x = HR),
      shape = 23, size = 4, fill = "red", color = "darkred"
    ) +
    
    # Points - Subgroups (square)
    geom_point(
      data = plot_data %>% filter(!is_overall & !is.na(HR)),
      aes(x = HR),
      shape = 22, size = 2.5, fill = "steelblue", color = "darkblue"
    ) +
    
    # 왼쪽 라벨
    geom_text(
      aes(x = x_limits[1] * 0.55, label = display_label, 
          fontface = if_else(is_header, "bold", "plain"),
          color = label_color),
      hjust = 0, size = 3.3
    ) +
    
    # 오른쪽 HR 텍스트
    geom_text(
      aes(x = x_limits[2] * 1.05, label = hr_text,
          fontface = if_else(is_overall, "bold", "plain"),
          color = label_color),
      hjust = 0, size = 3
    ) +
    
    # P interaction
    geom_text(
      aes(x = x_limits[2] * 1.55, label = p_int_text),
      hjust = 0, size = 2.8, color = "gray40", fontface = "italic"
    ) +
    
    # Scales
    scale_x_continuous(
      trans = "log10",
      breaks = x_breaks,
      labels = x_breaks,
      limits = c(x_limits[1] * 0.45, x_limits[2] * 2),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0.03, 0.03)) +
    scale_color_identity() +
    
    # Favors 라벨
    annotate("text", x = 0.6, y = 0, label = "Favors\nBP-DES", 
             size = 2.8, fontface = "italic", color = "gray50", lineheight = 0.9) +
    annotate("text", x = 1.7, y = 0, label = "Favors\nDP-DES", 
             size = 2.8, fontface = "italic", color = "gray50", lineheight = 0.9) +
    
    # X축 라벨
    labs(x = "Hazard Ratio (BP-DES vs DP-DES)") +
    
    # Theme
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 8)),
      axis.text.x = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3),
      plot.margin = margin(10, 120, 10, 10),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    
    coord_cartesian(clip = "off")
  
  # 저장
  if (!is.null(filename)) {
    height_calc <- n_rows * 0.32 + 1.5
    
    ggsave(paste0(filename, ".png"), p, width = 10, height = height_calc, dpi = 300)
    ggsave(paste0(filename, ".tiff"), p, width = 10, height = height_calc, dpi = 300, compression = "lzw")
    ggsave(paste0(filename, ".pdf"), p, width = 10, height = height_calc)
    
    cat("Saved:", filename, "(png/tiff/pdf)\n")
  }
  
  return(p)
}

################################################################################
# 실행: 기존 subgroup_results_tlf 사용
################################################################################

# 1. 기존 결과를 forest plot 형식으로 변환
forest_data_tlf <- convert_to_forest_format(subgroup_results_tlf)

# 2. 결과 확인
print(forest_data_tlf)

# 3. Forest plot 생성
forest_tlf <- create_forest_plot(
  data = forest_data_tlf,
  x_breaks = c(0.5, 0.75, 1.0, 1.5, 2.0),
  x_limits = c(0.4, 2.2),
  filename = "Figure_Subgroup_TLF_5years"
)

print(forest_tlf)

################################################################################
# MI at 1 Year도 동일하게 처리
################################################################################

# forest_data_mi <- convert_to_forest_format(subgroup_results_mi)
# 
# forest_mi <- create_forest_plot(
#   data = forest_data_mi,
#   x_breaks = c(0.5, 0.75, 1.0, 1.5, 2.0),
#   x_limits = c(0.3, 2.5),
#   filename = "Figure_Subgroup_MI"
# )
# 
# print(forest_mi)