## ---- packages
library(foreach)
library(parallel)
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
PACKAGES <- c('tidyverse', 'sf', 'gstat', 'spsann', 'cli')
## ----end

## ---- Paths
DATA_PATH <<- '../data/'
OUTPUT_PATH <<- '../output/'
fig.width <<- 6
if (!dir.exists(DATA_PATH)) dir.create(DATA_PATH)
if (!dir.exists(OUTPUT_PATH)) dir.create(OUTPUT_PATH)
if (!dir.exists(paste0(DATA_PATH, 'processed'))) dir.create(paste0(DATA_PATH, 'processed'))
if (!dir.exists(paste0(OUTPUT_PATH, 'figures'))) dir.create(paste0(OUTPUT_PATH, 'figures'))
if (!dir.exists(paste0(OUTPUT_PATH, 'tables'))) dir.create(paste0(OUTPUT_PATH, 'tables'))
## ----end
