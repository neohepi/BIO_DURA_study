install.packages(c("usethis", "gitcreds"))
library(usethis)
library(gitcreds)

usethis::use_git_config(
  user.name = "neohepi",
  user.email = "neohepi@gmail.com"
)
usethis::create_github_token()
gitcreds::gitcreds_set()

#usethis::use_git()
usethis::use_github()
