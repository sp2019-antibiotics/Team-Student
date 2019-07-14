##Creating data-raw/ by entering the command "devtools::use_data_raw()".
##We moved the raw Dataset "ZD.csv" into the now existent "data-raw" file
##and into the general t.mix file.
## code to prepare `DATASET` dataset goes here

library(dplyr)
library(magrittr)

experiment <-
  read.csv('ZD.csv') %>%
  mutate(experiment = 1)
devtools::use_data(experiment)

##Even though it constantly says that "devtools" is outdated, it works.

##usethis::use_data("DATASET")
