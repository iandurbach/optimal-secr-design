### Simple example generating an optimal SCR design
### (non-uniform animal density)

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(stringr)
library(secr)
library(secrdesign)
#library(oSCR)
library(raster)
library(kofnGA)

source("oSCR/e2dist.R")
source("oSCR/scrdesignGAenrm.R")
source("oSCR/scrdesignOFenrm.R")

# # https://www.r-bloggers.com/three-ways-to-call-cc-from-r/
# source("oSCR/LambdaL.R")
# source("oSCR/utils_for_enrmL.R")
# dyn.load("oSCR/mysecrdesign.so")

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

b_ac <- 1 # D = exp(b_0 + b_ac * covariate), see below

mask <- read.mask(data = bigmask_df)
plot(mask)
plot(read.mask(data = alltraps_df), col = 'blue', add = TRUE, axes =T)              

# activity center density assumed to be some function of a covariate
# 1) user needs to specify strength of this relationship (b_ac above)
# 2) these need to be defined for all mask points
Dac <- exp(b_ac * bigmask_df$stdGC)
# now standardize the Dac so that its average density is D
Dac <- Dac/sum(Dac) * (D * nrow(bigmask_df))
mean(Dac) # should be D

# plot
bigmask_df$Dac <- Dac
bigmask_df %>%
  ggplot(aes(x = x, y = y)) + 
  geom_tile(aes(fill = Dac), colour = "black") +
  scale_fill_viridis_c() + coord_equal()

# traps = subset of mask points that are potential camera locations
alltraps <- alltraps_df[,1:2] %>% as.matrix()

# generate optimal design
mnr <- scrdesignGAenrm(statespace = mask,
                       alltraps = alltraps,
                       ntraps = nT,
                       lambda0 = lambda0,
                       sigma = sigma,
                       D = Dac,
                       occasions = 1,
                       detector = "count",
                       ngen = 20,
                       popsize = 500,
                       seed = 700)

# extract camera locations
mnr_traps <- as.data.frame(mnr$optimaltraps)

# plot
mnr_traps %>%
  ggplot(aes(x = x, y = y)) + 
  geom_tile(data = mask_df, aes(fill = stdGC), colour = "black") +
  geom_point(size = 1, colour = "red") + 
  scale_fill_viridis_c() + coord_equal()
