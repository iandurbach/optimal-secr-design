### Simple example generating an optimal SCR design
### (non-uniform baseline encounter rate)

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

source("oSCR/scrdesignGAenrm.R")
source("oSCR/scrdesignOFenrm.R")

# requires updating secrdesign's En and Er functions, these are C fns so need compiling, see repo readme
# # https://www.r-bloggers.com/three-ways-to-call-cc-from-r/
source("oSCR/LambdaL.R")
source("oSCR/utils_for_enrmL.R")
dyn.load("oSCR/mysecrdesign.so")

# input file contains data frames with mask locations, and all potential camera locations
load("data/TostExample.RData")

# design parameters
ndesigns <- 1
nT <- 30

# assumptions about animal density and movement
lambda0 <- 1  
dens_per_100km2 <- 2 # mean animal density per 100km2, SLs are ~1
D <- dens_per_100km2 / 10000
sigma <- 3000

dt <- "count"

b_det <- 1 # D = exp(b_0 + b_ac * covariate), see below

mask <- read.mask(data = bigmask_df)
plot(mask)
plot(read.mask(data = alltraps_df), col = 'blue', add = TRUE, axes =T)              

# detectability assumed to be some function of a covariate
# 1) user needs to specify strength of this relationship (b_det above)
# 2) these need to be defined for all possible detectors
splam0 <- exp(as.numeric(b_det * scale(-alltraps_df$x)))
# now standardize the Dac so that its average detectability over traps is lambda0
splam0 <- splam0/sum(splam0) * (lambda0 * nrow(alltraps_df))
mean(splam0) # should be D

# traps = subset of mask points that are potential camera locations
alltraps <- alltraps_df[,1:2] %>% as.matrix()

# generate optimal design
mnr <- scrdesignGAenrm(statespace = mask,
                       alltraps = alltraps,
                       ntraps = nT,
                       lambda0 = splam0,
                       sigma = sigma,
                       D = D,
                       occasions = 1,
                       detector = "count",
                       ngen = 30,
                       popsize = 1000,
                       crit = 3,
                       seed = 700)

# extract camera locations
mnr_traps <- as.data.frame(mnr$optimaltraps)

# plot
mnr_traps %>%
  ggplot(aes(x = x, y = y)) + 
  geom_tile(data = mask_df, aes(fill = scale(-x)), colour = "black") +
  geom_point(size = 1, colour = "red") + 
  scale_fill_viridis_c() + coord_equal()
