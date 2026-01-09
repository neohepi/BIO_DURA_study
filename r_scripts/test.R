library(readxl)

temp.death <- read_excel("excel_files/cause_of_death.xlsx")

library(dplyr)

temp <- matched_data %>%
  rows_update(
    temp.death %>% select(pt_id, cardiac_death), 
    by = "pt_id",
    unmatched = "ignore"  # 이 옵션을 추가하면 에러 없이 일치하는 행만 업데이트됩니다.
  )
matched_data <- temp
