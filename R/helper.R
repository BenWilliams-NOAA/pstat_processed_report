library(tidyverse)
library(lubridate)
library(mgcViz)
library(funcr) # devtools::install_github("ben-williams/FNGr")
library(mgcv)
library(readxl)
library(ggridges)
library(EnvStats)
theme_set(theme_report())
library(PBSmapping)
data("nepacLL")
nepacLL %>% 
  dplyr::select(group=PID, POS=POS, long=X, lat=Y) -> ak_map