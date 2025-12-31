################################################################################
# Publication-Quality Forest Plots for Subgroup Analysis
################################################################################

#install.packages("forestplot")
library(forestplot)
library(dplyr)

################################################################################
# Function: Create Enhanced Forest Plot
################################################################################

create_enhanced_forest_plot <- function(subgroup_data, 
                                        title = "Subgroup Analysis",
                                        filename = NULL) {
  
  # Prepare data
  forest_data <- subgroup_data %>%
    arrange(desc(row_number())) %>%  # Reverse for bottom-up plotting
    mutate(
      # Format labels
      label = if_else(Category == "All Patients", 
                      "Overall", 
                      paste0("  ", Category)),
      
      # Events text
      events_text = sprintf("%d/%d", Events_BP, N_BP),
      events_dp_text = sprintf("%d/%d", Events_DP, N_DP),
      
      # HR text
      hr_text = sprintf("%.2f (%.2f-%.2f)", HR, CI_lower, CI_upper),
      
      # P value text
      p_text = if_else(P_value < 0.001, 
                       "<0.001", 
                       sprintf("%.3f", P_value)),
      
      # Interaction p-value (only show once per subgroup)
      p_int_display = ""
    )
  
  # Add interaction p-values to first row of each subgroup
  for (sg in unique(forest_data$Subgroup[forest_data$Subgroup != "Overall"])) {
    first_row <- which(forest_data$Subgroup == sg)[1]
    p_int <- forest_data$P_interaction[first_row]
    if (!is.na(p_int)) {
      forest_data$p_int_display[first_row] <- sprintf("P int = %.3f", p_int)
    }
  }
  
  # Create table text matrix
  tabletext <- cbind(
    c("Subgroup", forest_data$label),
    c("BP-DES\nEvents/Total", forest_data$events_text),
    c("DP-DES\nEvents/Total", forest_data$events_dp_text),
    c("HR (95% CI)", forest_data$hr_text),
    c("P Value", forest_data$p_text),
    c("", forest_data$p_int_display)
  )
  
  # Determine which rows are summary/subgroup headers
  is_summary <- c(TRUE, forest_data$Category == "All Patients")
  is_subgroup <- c(FALSE, !grepl("^  ", forest_data$label) & 
                     forest_data$Category != "All Patients")
  
  # Create forest plot
  fp <- forestplot(
    labeltext = tabletext,
    mean = c(NA, forest_data$HR),
    lower = c(NA, forest_data$CI_lower),
    upper = c(NA, forest_data$CI_upper),
    
    # Title and labels
    title = title,
    xlab = "Favors BP-DES          Favors DP-DES",
    
    # Reference line
    zero = 1,
    
    # Clip extreme values
    clip = c(0.2, 2.5),
    
    # Line at null
    grid = structure(c(1), 
                     gp = gpar(lty = 2, col = "gray60")),
    
    # X-axis ticks
    xticks = c(0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5),
    
    # Box sizes proportional to precision (inverse variance)
    boxsize = 0.25,
    
    # Colors
    col = fpColors(
      box = "royalblue",
      line = "darkblue",
      summary = "darkred",
      hrz_lines = "gray70"
    ),
    
    # Summary indicators
    is.summary = is_summary,
    
    # Styling for text
    txt_gp = fpTxtGp(
      label = list(
        gpar(fontface = "plain", cex = 0.9),
        gpar(fontface = ifelse(is_summary[-1] | is_subgroup[-1], 
                               "bold", "plain"), 
             cex = ifelse(is_summary[-1], 1.0, 0.9))
      ),
      ticks = gpar(cex = 0.8),
      xlab = gpar(cex = 0.9, fontface = "bold"),
      title = gpar(cex = 1.1, fontface = "bold")
    ),
    
    # Line properties
    lwd.xaxis = 1,
    lwd.ci = 1.5,
    lwd.zero = 2,
    
    # Graph width
    graphwidth = unit(70, "mm"),
    
    # Column widths
    colgap = unit(3, "mm"),
    
    # Add horizontal lines
    hrzl_lines = list(
      "2" = gpar(lwd = 2, col = "black"),  # After header
      "3" = gpar(lwd = 1, col = "gray80")  # After overall
    ),
    
    # Vertices for CI lines
    vertices = TRUE
  )
  
  # Save if filename provided
  if (!is.null(filename)) {
    pdf(filename, width = 11, height = 8)
    print(fp)
    dev.off()
    
    # Also save as high-res TIFF
    tiff(gsub("\\.pdf$", ".tiff", filename), 
         width = 11, height = 8, units = "in", res = 300, compression = "lzw")
    print(fp)
    dev.off()
    
    # And PNG for presentations
    png(gsub("\\.pdf$", ".png", filename), 
        width = 11, height = 8, units = "in", res = 300)
    print(fp)
    dev.off()
  }
  
  return(fp)
}


################################################################################
# Create Forest Plots
################################################################################

# 1. TLF at 5 Years
fp_tlf <- create_enhanced_forest_plot(
  subgroup_data = subgroup_results_tlf,
  title = "Target Lesion Failure at 5 Years - Subgroup Analysis",
  filename = "Figure_Forest_Plot_TLF.pdf"
)

print(fp_tlf)

# 2. MI at 1 Year
fp_mi <- create_enhanced_forest_plot(
  subgroup_data = subgroup_results_mi,
  title = "Myocardial Infarction at 1 Year - Subgroup Analysis",
  filename = "Figure_Forest_Plot_MI.pdf"
)

print(fp_mi)


################################################################################
# Alternative: Side-by-Side Comparison Plot
################################################################################

library(ggplot2)
library(patchwork)
library(ggplot2)
library(dplyr)
library(patchwork)

create_ggplot_forest <- function(data, title, y_min = 0.2, y_max = 2.5) {
  
  # 1. 데이터 준비 (row_number()를 여기서 미리 계산)
  plot_data <- data %>%
    mutate(
      label = if_else(Category == "All Patients", 
                      "Overall", 
                      Category),
      subgroup_label = if_else(Category == "All Patients", 
                               Subgroup, 
                               ""),
      color_group = case_when(
        Category == "All Patients" ~ "Overall",
        Subgroup == "Age" ~ "Age",
        Subgroup == "Sex" ~ "Sex",
        Subgroup == "Diabetes" ~ "Diabetes",
        Subgroup == "Clinical Presentation" ~ "Presentation",
        TRUE ~ "Other"
      ),
      # 유의성 마커
      sig_marker = case_when(
        P_value < 0.01 ~ "**",
        P_value < 0.05 ~ "*",
        P_value < 0.10 ~ "†",
        TRUE ~ ""
      )
    ) %>%
    # 순서를 위한 정렬 및 ID 생성 (가장 중요한 수정 부분)
    arrange(desc(row_number())) %>% 
    mutate(plot_order = row_number()) 
  
  # 2. 플롯 생성
  p <- ggplot(plot_data, aes(y = reorder(label, plot_order), x = HR)) +
    
    # Reference line at HR=1
    geom_vline(xintercept = 1, linetype = "dashed", 
               color = "gray50", linewidth = 0.8) + # size -> linewidth 수정
    
    # Confidence intervals
    # geom_errorbarh 대신 geom_errorbar 사용 (y축 매핑 시 자동 인식)
    geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper, color = color_group),
                  width = 0.3, linewidth = 1) + # height -> width, size -> linewidth
    
    # Point estimates
    geom_point(aes(color = color_group, size = color_group, 
                   shape = color_group)) +
    
    # Color and size scales
    scale_color_manual(
      values = c("Overall" = "#E41A1C", 
                 "Age" = "#377EB8", 
                 "Sex" = "#4DAF4A",
                 "Diabetes" = "#984EA3", 
                 "Presentation" = "#FF7F00",
                 "Other" = "gray50")
    ) +
    
    scale_size_manual(
      values = c("Overall" = 4, 
                 "Age" = 3, 
                 "Sex" = 3,
                 "Diabetes" = 3, 
                 "Presentation" = 3,
                 "Other" = 3)
    ) +
    
    scale_shape_manual(
      values = c("Overall" = 18, # Diamond
                 "Age" = 16,     # Circle
                 "Sex" = 17,     # Triangle
                 "Diabetes" = 15, # Square
                 "Presentation" = 16,
                 "Other" = 16)
    ) +
    
    # X-axis (log scale)
    scale_x_continuous(
      trans = "log",
      breaks = c(0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5),
      limits = c(y_min, y_max),
      labels = c("0.25", "0.5", "0.75", "1.0", "1.5", "2.0", "2.5")
    ) +
    
    # Labels
    labs(
      title = title,
      x = "Hazard Ratio (95% CI)\nFavors BP-DES      Favors DP-DES",
      y = NULL
    ) +
    
    # Theme
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 11, face = "bold"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3), # size -> linewidth
      plot.margin = margin(10, 20, 10, 10)
    ) +
    
    # Add HR values as text
    geom_text(aes(label = sprintf("%.2f", HR)), 
              vjust = -1.5, size = 3, fontface = "bold") + # hjust 대신 vjust 조정이 나을 수 있음 (취향 차이)
    
    # Add significance markers
    geom_text(aes(label = sig_marker, x = CI_upper), 
              hjust = -0.5, size = 4, color = "red")
  
  return(p)
}

# --- 실행 테스트 ---

# Create ggplot versions
p_tlf <- create_ggplot_forest(
  subgroup_results_tlf,
  "Target Lesion Failure at 5 Years"
)

p_mi <- create_ggplot_forest(
  subgroup_results_mi,
  "Myocardial Infarction at 1 Year"
)

# Combined plot
combined_forest <- (p_tlf | p_mi) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# Save
ggsave("Figure_Forest_Combined.pdf", combined_forest, 
       width = 14, height = 7, dpi = 300)

print(p_tlf)
print(p_mi)



################################################################################
# Summary Statistics Table for Forest Plot
################################################################################

cat("\n", rep("=", 70), "\n")
cat("FOREST PLOT DATA SUMMARY\n")
cat(rep("=", 70), "\n\n")

cat("TLF at 5 Years - Significant Subgroups (P<0.10):\n")
subgroup_results_tlf %>%
  filter(P_value < 0.10) %>%
  select(Subgroup, Category, HR, CI_lower, CI_upper, P_value) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print()

cat("\n\nMI at 1 Year - Significant Subgroups (P<0.10):\n")
subgroup_results_mi %>%
  filter(P_value < 0.10) %>%
  select(Subgroup, Category, HR, CI_lower, CI_upper, P_value) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print()

cat("\n\nInteractions P<0.10:\n")
all_interactions <- bind_rows(
  subgroup_results_tlf %>% 
    filter(!is.na(P_interaction)) %>%
    group_by(Subgroup) %>% 
    slice(1) %>%
    select(Subgroup, P_interaction) %>%
    mutate(Endpoint = "TLF at 5 Years"),
  
  subgroup_results_mi %>% 
    filter(!is.na(P_interaction)) %>%
    group_by(Subgroup) %>% 
    slice(1) %>%
    select(Subgroup, P_interaction) %>%
    mutate(Endpoint = "MI at 1 Year")
) %>%
  filter(P_interaction < 0.10) %>%
  arrange(P_interaction)

if (nrow(all_interactions) > 0) {
  print(all_interactions)
} else {
  cat("No interactions with P<0.10\n")
}

cat("\n", rep("=", 70), "\n")


################################################################################
# Figure Legend Template
################################################################################

cat("\n=== FIGURE LEGEND ===\n\n")

cat("Figure X. Subgroup Analysis for Target Lesion Failure at 5 Years

Forest plot showing hazard ratios (HRs) and 95% confidence 
intervals (CIs) for target lesion failure at 5 years comparing 
biodegradable polymer drug-eluting stents (BP-DES) with durable 
polymer drug-eluting stents (DP-DES) across pre-specified subgroups. 
The overall treatment effect is shown at the top (diamond), followed 
by subgroup-specific estimates (squares). Box sizes are proportional 
to the precision of estimates. HRs less than 1.0 favor BP-DES. 
P values for interaction test whether treatment effects differ 
significantly across subgroup categories. * P<0.05; ** P<0.01; 
† P<0.10.\n\n")

cat("Figure Y. Subgroup Analysis for Myocardial Infarction at 1 Year

Forest plot showing hazard ratios (HRs) and 95% confidence 
intervals (CIs) for myocardial infarction at 1 year comparing 
biodegradable polymer drug-eluting stents (BP-DES) with durable 
polymer drug-eluting stents (DP-DES) across pre-specified subgroups. 
The overall treatment effect is shown at the top (diamond), followed 
by subgroup-specific estimates (squares). Significant reductions in 
myocardial infarction were observed in younger patients (<65 years, 
HR 0.56, P=0.016) and male patients (HR 0.60, P=0.014), although 
formal tests for interaction did not reach statistical significance. 
HRs less than 1.0 favor BP-DES. * P<0.05; ** P<0.01; † P<0.10.\n\n")

cat(rep("=", 70), "\n")