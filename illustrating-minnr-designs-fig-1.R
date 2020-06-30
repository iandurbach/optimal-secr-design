### Toy example illustrating optimizing SCR designs by exchange algorithm or genetic algorithm, and 
### comparing simulated CV(D) of these designs with those of regular 2sigma grid design. Generates 
### Fig. 1 of paper. Main point is to show that despite smaller sample sizes (n,r), regular grids 
### can have similar or better CV(D) than min(n,r) designs, because of more precise estimates of sigma.

library(dplyr)
library(tidyr)
library(ggplot2)
library(secr)
library(secrdesign)
library(oSCR)
library(kofnGA)
library(sf)

# design functions
source("oSCR/SCRdesignGAenrm.R")
source("oSCR/SCRdesignOFenrm.R")
source("oSCR/SCRdesignEAenrm.R")
source("fn-simulating-cvD.R")

# detection function parameters
lambda0 <- 0.2
sigma <- 20
D <- 20
my_origin <- c(500,500) 

# make grid of possible camera locations
sp.ratio <- 0.5 # spacing as proportion of sigmas
detector <- "count"
grids_2s <- read.traps(data = expand.grid(x = seq(from = 500, by = 2 * sigma, length.out = 12),
                                          y = seq(from = 500, by = 2 * sigma, length.out = 12)), detector = detector)
masks <- make.mask(grids_2s, spacing = sigma / 2, type = 'traprect', buffer = 4*sigma)
grids <- read.traps(data = expand.grid(x = seq(from = min(grids_2s$x), to = max(grids_2s$x), by = sp.ratio * sigma),
                                       y = seq(from = min(grids_2s$y), to = max(grids_2s$y), by = sp.ratio * sigma)), detector = detector)
plot(masks)
plot(grids, add = T)

# convert mask and trap locations to matrices for design functions
statespace <- masks[,1:2] %>% as.matrix()
traps <- grids[,1:2] %>% as.matrix()

#################################
#### generating designs
#################################

# generate regular 5x5 grid with 2sigma spacing
eqs_grids <- make.grid(nx = 5, ny = 5, spacex = 2*sigma, detector = detector, originxy = c(640,640))
plot(masks)
plot(eqs_grids, add = T)    

# generate min(n,r) design with exchange algorithm, no penalty

set.seed(700)
optSCR_ex0 <- SCRdesignEAenrm(statespace = statespace,
                        all.traps = grids,
                        ntraps = 25, # number of cameras available
                        ndesigns = 1, # number of random starting points
                        beta0 = log(lambda0),
                        sigma = sigma, # SCR pars (log(lambda0) and usual sigma)
                        D = D,
                        occasions = 5,
                        detector = detector,
                        crit = 6, # criterion = 6 for min(n,r)
                        use.secr = TRUE, 
                        pen_wt = 0)

# extract camera locations
opt_traps_ex0 <- as.data.frame(optSCR_ex0$Xlst[[1]])

#plot
opt_grid_ex0 <- read.traps(data = opt_traps_ex0, detector = "count")
plot(masks, axes=T)
plot(opt_grid_ex0, add = T)                

# generate min(n,r) design with exchange algorithm, with penalty
set.seed(700)
optSCR_ex <- SCRdesignEAenrm(statespace = statespace,
                       all.traps = grids,
                       ntraps = 25, # number of cameras available
                       ndesigns = 1, # number of random starting points
                       beta0 = log(lambda0),
                       sigma = sigma, # SCR pars (log(lambda0) and usual sigma)
                       D = D,
                       occasions = 5,
                       detector = detector,
                       crit = 6,
                       use.secr = TRUE, 
                       pen_wt = 50) # added penalty

# extract camera locations
opt_traps_ex <- as.data.frame(optSCR_ex$Xlst[[1]])

#plot
opt_grid_ex <- read.traps(data = opt_traps_ex, detector = "count")
plot(masks, axes=T)
plot(opt_grid_ex, add = T)                

# generate min(n,r) design with genetic algorithm and no penalty
optSCR0 <- scrdesignGAenrm(statespace = masks,
                           alltraps = grids,
                           ntraps = 25, # number of cameras available
                           beta0 = log(lambda0),
                           sigma = sigma, 
                           D = D, # per mask cell not per ha!
                           occasions = 5,
                           detector = detector,
                           ngen = 20,
                           popsize = 1000,
                           crit = 3, # criterion 3 for min(n,r)
                           pen_wt = 0, sum_wt = 0, seed = 700)

# extract camera locations
opt_traps0 <- as.data.frame(optSCR0$optimaltraps)

#plot
opt_grid0 <- read.traps(data = opt_traps0, detector = "count")
plot(masks, axes=T)
plot(opt_grid0, add = T)                

# generate min(n,r) design with genetic algorithm and with penalty
optSCR <- scrdesignGAenrm(statespace = masks,
                          alltraps = grids,
                          ntraps = 25, # number of cameras available
                          beta0 = log(lambda0),
                          sigma = sigma, 
                          D = D, # per mask cell not per ha!
                          occasions = 5,
                          detector = "count",
                          ngen = 20,
                          popsize = 1000,
                          crit = 3,
                          pen_wt = 50, sum_wt = 0, seed = 700)

# extract camera locations
opt_traps <- as.data.frame(optSCR$optimaltraps)

#plot
opt_grid <- read.traps(data = opt_traps, detector = "count")
plot(masks, axes=T)
plot(opt_grid, add = T)                

#################################
#### simulating CV(D) of designs
#################################

set.seed(111)
nsim <- 100
eq_sim <- simulating_cvD(grids = eqs_grids, masks = masks, D = D, lambda0 = lambda0, sigma = sigma, noccasions = 5, my_pop = NULL, nrepl = nsim)
set.seed(111)
opt_sim_ex0<- simulating_cvD(grids = opt_grid_ex0, masks = masks, D = D, lambda0 = lambda0, sigma = sigma, noccasions = 5, my_pop = NULL, nrepl = nsim)
set.seed(111)
opt_sim_ex<- simulating_cvD(grids = opt_grid_ex, masks = masks, D = D, lambda0 = lambda0, sigma = sigma, noccasions = 5, my_pop = NULL, nrepl = nsim)
set.seed(111)
opt_sim0 <- simulating_cvD(grids = opt_grid0, masks = masks, D = D, lambda0 = lambda0, sigma = sigma, noccasions = 5, my_pop = NULL, nrepl = nsim)
set.seed(111)
opt_sim <- simulating_cvD(grids = opt_grid, masks = masks, D = D, lambda0 = lambda0, sigma = sigma, noccasions = 5, my_pop = NULL, nrepl = nsim)

# theoretical enrm 
enrm_eq <- Enrm(D = D, eqs_grids, masks, list(lambda0 = lambda0, sigma = sigma), noccasions = 5)
enrm_opt_ex0 <- Enrm(D = D, opt_grid_ex0, masks, list(lambda0 = lambda0, sigma = sigma), noccasions = 5)
enrm_opt_ex <- Enrm(D = D, opt_grid_ex, masks, list(lambda0 = lambda0, sigma = sigma), noccasions = 5)
enrm_opt0 <- Enrm(D = D, opt_grid0, masks, list(lambda0 = lambda0, sigma = sigma), noccasions = 5)
enrm_opt <- Enrm(D = D, opt_grid, masks, list(lambda0 = lambda0, sigma = sigma), noccasions = 5)

# extract
eq_sim_sum <- eq_sim$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "grid2s",  En = enrm_eq[1], Er = enrm_eq[2], Em = enrm_eq[3])
opt_sim_ex0_sum <- opt_sim_ex0$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "opt_sim_ex0", En = enrm_opt_ex0[1], Er = enrm_opt_ex0[2], Em = enrm_opt_ex0[3])
opt_sim_ex_sum <- opt_sim_ex$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "opt_sim_ex", En = enrm_opt_ex[1], Er = enrm_opt_ex[2], Em = enrm_opt_ex[3])
opt_sim0_sum <- opt_sim0$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "opt_sim0", En = enrm_opt0[1], Er = enrm_opt0[2], Em = enrm_opt0[3])
opt_sim_sum <- opt_sim$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "opt_sim", En = enrm_opt[1], Er = enrm_opt[2], Em = enrm_opt[3])

# combine
sim_sum <- rbind(eq_sim_sum, opt_sim_ex0_sum, opt_sim_ex_sum, opt_sim0_sum, opt_sim_sum)

# compare
all_sum <- sim_sum %>% group_by(design) %>% summarize(approx_cv = 1 / sqrt(min(mean(En), mean(Er))),
                                                      sim_cv = sd(E.N) / mean(E.N),
                                                      n = mean(n), En = mean(En),
                                                      r = mean(r), Er = mean(Er),
                                                      m = mean(moves), Em = mean(Em),
                                                      esa = mean(esa),
                                                      cvN = mean(cvN), cva = mean(cva), cvD = mean(cvD),
                                                      nsim = n())

#################################
#### plot for paper
#################################

# save(grids, eqs_grids, opt_grid_ex0, opt_grid_ex, opt_grid0, opt_grid, eq_sim, opt_sim_ex0, 
# opt_sim_ex, opt_sim0, opt_sim, sim_sum, all_sum, file = "output/toyex-variable-opt-algos")
# 
# load("output/toyex-variable-opt-algos")

mask_df <- data.frame(x = masks$x, y = masks$y) %>% mutate(ptt = 1)
allgrid_df <- data.frame(x = grids$x, y = grids$y) %>% mutate(method = "Possible locations", ptt = 1)
grid_df <- data.frame(x = eqs_grids$x, y = eqs_grids$y) %>% mutate(method = "Grid", ptt = 2)
opt_grid_df <- data.frame(x = opt_grid$x, y = opt_grid$y) %>% mutate(method = "GA", ptt = 2)
opt_grid_ex_df <- data.frame(x = opt_grid_ex$x, y = opt_grid_ex$y) %>% mutate(method = "Exchange + penalty", ptt = 2)
opt_grid_ex0_df <- data.frame(x = opt_grid_ex0$x, y = opt_grid_ex0$y) %>% mutate(method = "Exchange no penalty", ptt = 2)
opt_grid0_df <- data.frame(x = opt_grid0$x, y = opt_grid0$y) %>% mutate(method = "GA no penalty", ptt = 2)
traps_all <- rbind(allgrid_df, grid_df, opt_grid_ex_df, opt_grid_ex0_df, opt_grid_df)

p1 <- traps_all %>% 
  mutate(method = factor(method, levels = c("Possible locations", "Exchange no penalty", "Exchange + penalty", "GA",  "Grid"))) %>%
  ggplot(aes(x = x, y = y)) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black", alpha = 0.3) +
  geom_point(aes(size = factor(ptt)), colour = "red") + 
  scale_size_manual(values = c(.2,.7)) +
  facet_grid(. ~ method) + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
p1

# ggsave("output/setup.png", p1, width=9, height=3, dpi = 300)

