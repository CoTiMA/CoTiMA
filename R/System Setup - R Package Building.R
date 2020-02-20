#start new Project
#copy file to project folder R


install.packages(c("devtools",
                   "roxygen2",
                   "testthat",
                   "knitr",
                   "usethis",
                   "rmarkdown",
                   "pacman",
                   "tidyverse"
                   "fs"))

library(pacman)
p_load(devtools,
       knitr,
       roxygen2,
       testthat,
       usethis,
       rmarkdown,
       tidyverse,
       fs)

has_devel()


usethis::use_build_ignore("System Setup - R Package Building.R")
usethis::use_build_ignore("Testfile.R")

load_all()
check()
#use_mit_license("Markus Homberg")
document()
install()
use_testthat()
use_test("functionname")
test()
use_package("packagename")
use_readme_rmd()
rmarkdown::render("README.Rmd")
