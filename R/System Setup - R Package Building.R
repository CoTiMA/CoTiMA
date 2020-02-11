#start new Project
#copy file to project folder R


install.packages(c("devtools",
                   "roxygen2",
                   "testthat",
                   "knitr",
                   "usethis",
                   "rmarkdown",
                   "pacman"))

library(pacman)
p_load(devtools,
       knitr,
       roxygen2,
       testthat,
       usethis,
       rmarkdown)

has_devel()


usethis::use_build_ignore("System Setup - R Package Building.R")

load_all()
