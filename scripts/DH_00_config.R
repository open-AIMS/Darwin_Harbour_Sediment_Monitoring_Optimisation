## ---- packages
library(knitr)
library(tidyverse)
library(readxl)
library(sf)
library(ggmap)
library(ggspatial)
library(INLA)
library(inlabru)
library(fields)
library(patchwork)
library(gstat)
library(mgcv)
library(cli)
library(spsann)
source('../scripts/DH_functions.R')
## ----end

## ---- Paths
DATA_PATH <<- '../data/'
OUTPUT_PATH <<- '../output/'
fig.width <<- 6
if (!dir.exists(paste0(DATA_PATH, 'processed'))) dir.create(paste0(DATA_PATH, 'processed'))
## ----end

