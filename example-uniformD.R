### Simple example generating an optimal SCR design
### (uniform animal density)

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(stringr)
library(secr)
library(secrdesign)
library(oSCR)
library(raster)
library(kofnGA)

source("oSCR/SCRdesignGAenrm.R")
source("oSCR/SCRdesignOFenrm.R")

# input file contains data frames with mask locations, and all potential camera locations
load("data/TostExample.RData")

# design parameters
ndesigns <- 1
nT <- 30

# assumptions about animal density and movement
lambda0 <- 1  # beta0 = log(lambda0) = log(K * 'p0')
dens_per_100km2 <- 2 # mean animal density per 100km2, SLs are ~1
D <- dens_per_100km2 / 10000
sigma <- 3000

dt <- "count"

mask <- read.mask(data = bigmask_df)
plot(mask)
plot(read.mask(data = alltraps_df), col = 'blue', add = TRUE, axes =T)              

# traps = subset of mask points that are potential camera locations
alltraps <- alltraps_df[,1:2] %>% as.matrix()

# generate optimal design
mnr <- scrdesignGAenrm(statespace = mask,
                       alltraps = alltraps,
                       ntraps = nT,
                       beta0 = log(lambda0),
                       sigma = sigma,
                       D_per_mask_cell = D,
                       occasions = 1,
                       detector = "count",
                       ngen = 5,
                       popsize = 500,
                       crit = 3,
                       seed = 700)

# extract camera locations
mnr_traps <- as.data.frame(mnr$optimaltraps)

# plot
mnr_traps %>%
  ggplot(aes(x = x, y = y)) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  scale_fill_viridis_c() + coord_equal()
