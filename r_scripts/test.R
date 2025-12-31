#install.packages("gridExtra", "grid")
library(gridExtra)
library(grid)

plots_list = list(outcome$composite_plot, outcome$cardiac_death_plot, outcome$MI_plot, outcome$TLR_plot)

# Arrange plots (top half)
plots <- arrangeGrob(
  grobs = plots_list,
  ncol = 2,
  nrow = 2
)

# Display
grid.newpage()
grid.draw(plots)
