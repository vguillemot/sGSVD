library(testthat)
detach("package:sGSVD", unload = TRUE)
devtools::install_git(
  "https://github.com/vguillemot/sGSVD.git",
  credentials = git2r::cred_user_pass("juchiyu", getPass::getPass())
)
library(sGSVD)

test_check("sGSVD")

