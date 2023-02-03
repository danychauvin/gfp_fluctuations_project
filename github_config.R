library(usethis)
usethis::use_git_config(user.name = "danychauvin", user.email = "dany.chauvin@unibas.ch")
## create a personal access token for authentication:
usethis::create_github_token()
## in case usethis version < 2.0.0: usethis::browse_github_token() (or even better: update usethis!)
## set personal access token:
credentials::set_github_pat("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")
## or store it manually in '.Renviron':
usethis::edit_r_environ()
## store your personal access token with: GITHUB_PAT=xxxyyyzzz
## and make sure '.Renviron' ends with a newline

