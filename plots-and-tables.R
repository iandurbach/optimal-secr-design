## 

library(dplyr)
library(ggplot2)
library(gridExtra)
library(kableExtra)
library(oSCR)
library(tidyr)

######################
## Maps of example designs for uniform D, euclidean
######################

# load results
load("output/Tost_examples_all_D0.Rdata")  

# split out 
opt_traps_all$buffer <- opt_traps_all$buffer / 1000

# New facet label names for nT
nT.labs <- c("20 detectors", "40 detectors", "60 detectors")
names(nT.labs) <- c("20", "40", "60")

m_abc <- data.frame(nT = c(20,40,60), x = min(mask_df$x), y = min(mask_df$y), labs = c("(a)", "(b)", "(c)"))
m0 <- opt_traps_all %>% filter(sigma == 6000, beta0 == -0.8, nT != 50, buffer == 18, trap_id == 1) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ nT, labeller = labeller(nT = nT.labs)) + 
  labs(title = bquote("Effect of varying number of detectors (baseline"~ sigma~"and"~lambda[0]~","~3*sigma~"buffer)")) +
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

m_abc <- data.frame(nT = c(20,40,60), x = min(mask_df$x), y = min(mask_df$y), labs = c("(d)", "(e)", "(f)"))
m1 <- opt_traps_all %>% filter(sigma == 6000, beta0 == -0.8, nT != 50, buffer == 18, trap_id == 2) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ nT, labeller = labeller(nT = nT.labs)) + 
  labs(title = "Effect of random starts (conditions as above)") +
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

m_abc <- data.frame(nT = c(20,40,60), x = min(mask_df$x), y = min(mask_df$y), labs = c("(g)", "(h)", "(i)"))
m2 <- opt_traps_all %>% filter(sigma == 6000, beta0 == -0.8, nT != 50, buffer == 0, trap_id == 1) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ nT, labeller = labeller(nT = nT.labs)) + 
  labs(title = bquote("Effect of zero buffer (baseline"~ sigma~"and"~lambda[0]~")")) +
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

opt_traps_all$beta02 <- factor(opt_traps_all$beta0, levels = c("-1.5","-0.8","-0.4", "-0.113"),
                               ordered = TRUE, labels = c(expression(paste("50% of baseline ", lambda[0])),
                                                          expression(paste("baseline ", lambda[0])),
                                                          expression(paste("150% of baseline ", lambda[0])),
                                                          expression(paste("200% of baseline ", lambda[0]))))

m_abc <- data.frame(beta02 = levels(opt_traps_all$beta02)[-2], x = min(mask_df$x), y = min(mask_df$y), labs = c("(j)", "(k)", "(l)"))
m3 <- opt_traps_all %>% filter(sigma == 6000, buffer == 18, nT == 40, trap_id == 1, beta0 != -0.8) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ beta02, labeller = label_parsed) + 
  labs(title = bquote("Effect of varying encounter rate"~lambda[0]~"(40 detectors, baseline"~ sigma~","~3*sigma~"buffer)")) +
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

opt_traps_all$sigma2 <- factor(opt_traps_all$sigma, levels = c("3000","6000", "9000", "12000"),
                               ordered = TRUE, labels = c(expression(paste("50% of baseline ", sigma)),
                                                          expression(paste("baseline ", sigma)),
                                                          expression(paste("150% of baseline ", sigma)),
                                                          expression(paste("200% of baseline ", sigma))))

m_abc <- data.frame(sigma2 = levels(opt_traps_all$sigma2)[-2], x = min(mask_df$x), y = min(mask_df$y), labs = c("(m)", "(n)", "(o)"))
m4 <- opt_traps_all %>% filter(nT == 40, trap_id == 1, beta0 == -0.8, sigma != 6000) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ sigma2, labeller = label_parsed) + 
  labs(title = bquote("Effect of varying animal movement parameter"~sigma~"(40 detectors, baseline"~ lambda[0]~","~3*sigma~"buffer)")) +
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


# New facet label names for dt
opt_traps_all$dt <- factor(opt_traps_all$dt, levels = c("count", "proximity", "multi"),
                           labels = c("Count proximity", "Binary proximity", "Multi-catch"))

m_abc <- data.frame(dt = c("Count proximity", "Binary proximity", "Multi-catch"), x = min(mask_df$x), y = min(mask_df$y), labs = c("(p)", "(q)", "(r)"))
m4b <- opt_traps_all %>% filter(nT == 30, beta0 == -0.1, sigma == 2000) %>%
  ggplot(aes(x = x, y = y)) + 
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ dt) + 
  labs(title = bquote("Effect of detector type (30 detectors, 5 occasions, "~ lambda[0] == e^-0.1~","~sigma==2~"km, "~3*sigma~"buffer)")) +
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("output/tost-uni-designs.png", 
       arrangeGrob(m0,m1,m2,m3,m4,m4b, ncol = 1), width=10, height=12, dpi = 300)

######################
## Maps of example designs for non-Uniform D, euclidean
######################

load("output/Tost_examples_nonuniD.Rdata")  

opt_traps_all$nT2 <- factor(opt_traps_all$nT, levels = c("20", "40", "60"), 
                            labels = c(expression(paste("20 detectors")), 
                                       expression(paste("40 detectors")), 
                                       expression(paste("60 detectors"))))
opt_traps_all$b_ac2 <- factor(opt_traps_all$b_ac, levels = c("-1","1","3"),
                              ordered = TRUE, labels = c(expression(paste(alpha[1]," = -1")),
                                                         expression(paste(alpha[1]," = 1")),
                                                         expression(paste(alpha[1]," = 3"))))

m_abc <- expand.grid(b_ac2 = levels(opt_traps_all$b_ac2), 
                     nT2 = levels(opt_traps_all$nT2)) %>%
  mutate(x = min(mask_df$x), y = min(mask_df$y), 
         labs = c("(a)", "(d)", "(g)", "(b)", "(e)", "(h)", "(c)", "(f)", "(i)"))
m5 <- opt_traps_all %>% filter(trap_id == 1) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_tile(data = mask_df, aes(fill = stdGC), colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(b_ac2 ~ nT2, labeller = label_parsed) +
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


ggsave("output/tost-designs-nonU.png", 
       m5, width=10, height=6, dpi = 300)

######################
## Maps of example designs for non-euclidean, uniform D
######################

load("output/Tost_examples_noneuc.Rdata")  

opt_traps_all$nT2 <- factor(opt_traps_all$nT, levels = c("20", "40", "60"), 
                            labels = c(expression(paste("20 detectors")), 
                                       expression(paste("40 detectors")), 
                                       expression(paste("60 detectors"))))
opt_traps_all$alpha22 <- factor(opt_traps_all$alpha2, levels = c("-1","1","3"),
                                ordered = TRUE, labels = c(expression(paste(alpha[2]," = -1")),
                                                           expression(paste(alpha[2]," = 1")),
                                                           expression(paste(alpha[2]," = 3"))))

m_abc <- expand.grid(alpha22 = levels(opt_traps_all$alpha22), 
                     nT2 = levels(opt_traps_all$nT2)) %>%
  mutate(x = min(mask_df$x), y = min(mask_df$y), 
         labs = c("(a)", "(d)", "(g)", "(b)", "(e)", "(h)", "(c)", "(f)", "(i)"))
m6 <- opt_traps_all %>% filter(trap_id == 1) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_tile(data = mask_df, aes(fill = stdGC), colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(alpha22 ~ nT2, labeller = label_parsed) +
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("output/tost-designs-nonEuc.png", 
       m6, width=10, height=6, dpi = 300)

######################
## Line graph comparing CV(D) of optimal design vs non-optimal designs 
######################

load("output/Tost_Enrm_opt.Rdata")
optEnrm <- optEnrm %>% mutate(cvD = sqrt(1/pmin(En, Er)), method = "opt") %>% rename(ntraps = V1, b_ac = V2, beta0 = V3)

load("output/Tost_Enrm_nonopt.Rdata")
surveyEnrm <- surveyEnrm %>% mutate(cvD = sqrt(1/pmin(En, Er)), method = "survey") %>% rename(ntraps = V1, b_ac = V2, beta0 = V3)
gridEnrm <- gridEnrm %>% mutate(cvD = sqrt(1/pmin(En, Er)), method = "grid") %>% rename(ntraps = V1, b_ac = V2, beta0 = V3)
optgridEnrm <- optgridEnrm %>% mutate(cvD = sqrt(1/pmin(En, Er)), method = "opt_grid") %>% rename(ntraps = V1, b_ac = V2, beta0 = V3) %>% dplyr::select(-V4)

all_res <- rbind(optEnrm, surveyEnrm, gridEnrm, optgridEnrm)

summ_res <- all_res %>% group_by(beta0, b_ac, ntraps, method) %>% 
  summarize(mincvD = min(cvD), mediancvD = median(cvD)) %>%
  filter(!(method == "survey" & ntraps > 40))
summ_res <- summ_res %>% filter(b_ac == 1, (method == "grid")|(method == "survey")) %>% mutate(cvD = mediancvD) %>%
  rbind(summ_res %>% filter(b_ac != 1 | method %in% c("opt_grid","opt")) %>% mutate(cvD = mincvD))

summ_res$beta02 <- factor(summ_res$beta0, levels = c("-1.5","-0.8","-0.4"),
                          ordered = TRUE, labels = c(expression(paste("50% of baseline ", lambda[0])),
                                                     expression(paste("baseline ", lambda[0])),
                                                     expression(paste("150% of baseline ", lambda[0]))))
summ_res$b_ac2 <- factor(summ_res$b_ac, levels = c("0", "1"), 
                         labels = c(expression(paste("Uniform D")), 
                                    expression(paste("Non-uniform D (", alpha[1]," = 1)"))))

summ_res$method <- factor(summ_res$method, levels = c("survey", "grid", "opt_grid", "opt"),
                          labels = c("Actual survey", "Regular grid", "Regular grid + optimal.Spacing", "Optimized"))
m7 <- summ_res %>% 
  ggplot(aes(x = ntraps, y = cvD, colour = method, shape = method)) +
  facet_grid(b_ac2 ~ beta02, labeller = label_parsed) +
  geom_line() + geom_point() +
  scale_color_brewer(palette = "Dark2") +
  xlab("Number of detectors") + ylab(expression(paste("Approximate CV(", hat(D),")"))) + 
  theme_bw(base_size = 14) +
  theme(legend.position="bottom",legend.title = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("output/tost-opt-vs-nonopt.png", 
       m7, width=10, height=6, dpi = 300)

######################
## Distributions of inter-detector distances in optimal designs
######################

load("output/Tost_many_D0_for_dbns.Rdata")

# for now opt traps only available for beta0 > -1.5, for D~0, run these later if needed
all_traps <- all_traps %>% filter(beta0 != -1.5) %>% filter(b_ac == 0)

# extract distance to closest detector for each detector 
myseq <- seq(2000,16000,1000)
closest_detectors <- all_traps %>% 
  group_by(method, nT, beta0, trap_id) %>% 
  mutate(nearest_detdist = apply(e2dist(cbind(x,y), cbind(x,y)), 1, function(x){sort(x, FALSE)[2]}),
         nearest_cat = cut(nearest_detdist, breaks = 8, right = FALSE),
         mean_detdist = apply(e2dist(cbind(x,y), cbind(x,y)), 1, mean),
         mean_cat = cut(mean_detdist, breaks = 8, right = FALSE)) %>%
  arrange(method, nT, beta0, nearest_detdist) %>%
  ungroup() %>% group_by(method, nT, beta0) %>% 
  mutate(Fx = row_number() / n())

closest_detectors$nT2 <- factor(closest_detectors$nT, levels = c("20", "40", "60"), 
                                labels = c(expression(paste("20 detectors")), 
                                           expression(paste("40 detectors")), 
                                           expression(paste("60 detectors"))))
closest_detectors$beta02 <- factor(closest_detectors$beta0, levels = c("-1.5","-0.8","-0.4"),
                                   ordered = TRUE, labels = c(expression(paste("50% of baseline ", lambda[0])),
                                                              expression(paste("baseline ", lambda[0])),
                                                              expression(paste("150% of baseline ", lambda[0]))))
closest_detectors$method2 <- factor(closest_detectors$method, levels = c("survey", "grid", "opt_grid", "opt"),
                                    labels = c("Actual survey", "Regular grid", "Regular grid + optimal.Spacing", "Optimized"))


m8 <- closest_detectors %>% filter(nT %in% c(20,40,60)) %>%
  ggplot(aes(x = nearest_detdist/1000, y = Fx, colour = method2)) + 
  geom_path() + facet_grid(beta02 ~ nT2, labeller = label_parsed) +
  xlab("Nearest-detector distance (km)") + ylab("Cumulative probability") + 
  theme_bw(base_size = 14) +
  coord_cartesian(xlim = c(0, 12)) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position="bottom",legend.title = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

m8a <- closest_detectors %>% filter(nT %in% c(20,40,60)) %>%
  ggplot(aes(x = nearest_detdist/1000, y = Fx, colour = method2)) + 
  geom_path() + facet_grid(beta02 ~ nT2, labeller = label_parsed) +
  xlab("Nearest-detector distance (km)") + ylab("Cumulative probability") + 
  theme_bw(base_size = 14) +
  coord_cartesian(xlim = c(0, 12)) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position="none",legend.title = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

m_abc <- data.frame(nT2 = levels(closest_detectors$nT2), x = 12, y = 0, labs = c("(a)", "(b)", "(c)"))
m8b <- closest_detectors %>% filter(nT %in% c(20,40,60), beta0 == -0.8) %>%
  ggplot(aes(x = nearest_detdist/1000, y = Fx)) + 
  geom_path(aes(colour = method2)) + facet_grid(. ~ nT2, labeller = label_parsed) +
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 1, vjust = 0) + 
  xlab("Nearest-detector distance (km)") + ylab("Cumulative probability") + 
  theme_bw(base_size = 14) +
  coord_cartesian(xlim = c(0, 12)) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks=seq(0,12,2)) +
  theme(legend.position="none",legend.title = element_blank(), 
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


cdf_data <- data.frame(dists = as.numeric(), Fx = as.numeric(), nT = as.integer(), beta0 = as.numeric(), method = as.character())
for(m in unique(all_traps$method)){
  for(i in c(20,40,60)){
    if(!((m == "survey")&(i == 60))){
      for(b in c(-0.4, -0.8)){
        allxt <- c()
        for(j in 1:10){
          xx <- all_traps %>% filter(nT == i, trap_id == j, beta0 == b, method == m)
          ds <- e2dist(cbind(xx$x, xx$y),cbind(xx$x, xx$y))
          diag(ds) <- NA
          allxt <- c(allxt, unlist(ds[!is.na(ds)]))
        }
        
        cdf_data <- rbind(cdf_data, 
                          data.frame(dists = sort(allxt), 
                                     Fx = seq(0, 1, length.out = length(allxt)),
                                     nT = i,
                                     beta0 = b,
                                     method = m)) 
      }
    }
  }
}

cdf_data$nT2 <- factor(cdf_data$nT, levels = c("20", "40", "60"), 
                       labels = c(expression(paste("20 detectors")), 
                                  expression(paste("40 detectors")), 
                                  expression(paste("60 detectors"))))
cdf_data$beta02 <- factor(cdf_data$beta0, levels = c("-1.5","-0.8","-0.4"),
                          ordered = TRUE, labels = c(expression(paste("50% of baseline ", lambda[0])),
                                                     expression(paste("baseline ", lambda[0])),
                                                     expression(paste("150% of baseline ", lambda[0]))))
cdf_data$method2 <- factor(cdf_data$method, levels = c("survey", "grid", "opt_grid", "opt"),
                           labels = c("Actual survey", "Regular grid", "Regular grid + optimal.Spacing", "Optimized"))

m9 <- cdf_data %>% ggplot(aes(x = dists/1000, y = Fx, colour = method2)) + 
  geom_path() + facet_grid(beta02 ~ nT2, labeller = label_parsed) +
  xlab("Between-detector distance (km)") + ylab("Cumulative probability") + 
  theme_bw(base_size = 14) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position="bottom",legend.title = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

m_abc <- data.frame(nT2 = levels(cdf_data$nT2), x = 100, y = 0, labs = c("(d)", "(e)", "(f)"))
m9b <- cdf_data %>% filter(beta0 == -0.8) %>% 
  ggplot(aes(x = dists/1000, y = Fx)) + 
  geom_path(aes(colour = method2)) + facet_grid(. ~ nT2, labeller = label_parsed) +
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 1, vjust = 0) + 
  xlab("Between-detector distance (km)") + ylab("Cumulative probability") + 
  theme_bw(base_size = 14) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  theme(legend.position="bottom",legend.title = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("output/tost-closest-trap-dists.png", 
       m8, width=10, height=6, dpi = 300)

ggsave("output/tost-all-trap-dists.png", 
       m9, width=10, height=6, dpi = 300)

ggsave("output/tost-trap-dists.png", 
       arrangeGrob(m8b, m9b, ncol = 1, heights = c(5,6)), width=10, height=6, dpi = 300)

######################
## Covariate coverage plots
######################

load("output/Tost_examples_nonuniD.Rdata")  
load("data/Tost.RData")
mask <- TostMask 
b_ac <- 1

opt_traps_all <- opt_traps_all %>% left_join(mask_df, by = c("x", "y"))

linedata <- expand.grid(stdGC = seq(from = min(mask_df$stdGC), to = max(mask_df$stdGC), length.out = 100),
                        b_ac = c(-1, 1, 3)) %>%
  group_by(b_ac) %>%
  mutate(y = exp(b_ac * stdGC)) %>%
  mutate(y = y / max(y)) %>%
  ungroup()

opt_traps_all <- opt_traps_all %>% mutate(y = (nT - 20) / 200)

opt_traps_all$nT2 <- factor(opt_traps_all$nT, levels = c("20", "40", "60"), 
                            labels = c("20 detectors", 
                                       "40 detectors", 
                                       "60 detectors"))
opt_traps_all$b_ac2 <- factor(opt_traps_all$b_ac, levels = c("-1","1","3"),
                              ordered = TRUE, labels = c(expression(paste(alpha[1]," = -1")),
                                                         expression(paste(alpha[1]," = 1")),
                                                         expression(paste(alpha[1]," = 3"))))

linedata$b_ac2 <- factor(linedata$b_ac, levels = c("-1","1","3"),
                         ordered = TRUE, labels = c(expression(paste(alpha[1]," = -1")),
                                                    expression(paste(alpha[1]," = 1")),
                                                    expression(paste(alpha[1]," = 3"))))

m10 <- linedata %>% ggplot(aes(x = stdGC, y = y)) + 
  geom_line() + 
  geom_point(data = opt_traps_all, aes(pch = factor(nT2), colour = factor(nT2))) +
  facet_grid(. ~ b_ac2, labeller = label_parsed) +
  xlab("Terrain ruggedness (standardized)") + ylab("Density (standardized)") + 
  theme_bw(base_size = 14) +
  scale_color_brewer(palette = "Dark2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme(legend.position="bottom",legend.title = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("output/tost-coverage.png", 
       m10, width=10, height=3.5, dpi = 300)

######################
## tables of Enrm
######################

load("output/Tost_Enrm_opt.Rdata")
optEnrm <- optEnrm %>% mutate(cvD = sqrt(1/pmin(En, Er)), method = "opt") %>% rename(ntraps = V1, b_ac = V2, beta0 = V3)

load("output/Tost_Enrm_nonopt.Rdata")
surveyEnrm <- surveyEnrm %>% mutate(cvD = sqrt(1/pmin(En, Er)), method = "survey") %>% rename(ntraps = V1, b_ac = V2, beta0 = V3)
gridEnrm <- gridEnrm %>% mutate(cvD = sqrt(1/pmin(En, Er)), method = "grid") %>% rename(ntraps = V1, b_ac = V2, beta0 = V3)
optgridEnrm <- optgridEnrm %>% mutate(cvD = sqrt(1/pmin(En, Er)), method = "opt_grid") %>% rename(ntraps = V1, b_ac = V2, beta0 = V3) %>% dplyr::select(-V4)

all_res <- rbind(optEnrm, surveyEnrm, gridEnrm, optgridEnrm)

summ_res <- all_res %>% group_by(beta0, b_ac, ntraps, method) %>% 
  mutate(rankcvD = rank(cvD, ties.method = "random")) %>%
  filter(!(method == "survey" & ntraps > 40)) %>%
  mutate(EnEr = paste0(round(En, 1), " (", round(Er, 1), ")")) %>% ungroup()

summ_res <- summ_res %>% 
  filter(b_ac == 1, (method == "grid")|(method == "survey"), rankcvD == 50) %>% 
  rbind(summ_res %>% filter(b_ac != 1 | method %in% c("opt_grid","opt"), rankcvD == 1))

summ_res %>% 
  dplyr::select(beta0, b_ac, ntraps, method, EnEr) %>%
  spread(method, EnEr) %>%
  arrange(b_ac, beta0, ntraps) %>%
  filter(ntraps %in% c(20, 40, 60)) %>%
  dplyr::select(beta0, ntraps, opt, grid, opt_grid, survey) %>%
  kable("latex", caption = "Group Rows", booktabs = T) %>% 
  kable_styling() %>%
  pack_rows("Uniform activity centre density", 1, 3) %>%
  pack_rows("", 4, 6) %>%
  pack_rows("", 7, 9) %>%
  pack_rows("Non-uniform activity centre density", 10, 12, latex_gap_space = "1.3em") %>%
  pack_rows("", 13, 15) %>%
  pack_rows("", 16, 18)

######################
## supplementary material -- plot of sensitivity to misspecification 
######################

load("output/posthoc-check-robustness-to-misspec.RData")

Enrm_results_df <- data.frame(Enrm_results)
Enrm_results_df$method <- factor(Enrm_results_df$method, levels = 1:4, labels = c("opt", "grid", "opt_grid", "survey"))
Enrm_results_df$varying <- factor(Enrm_results_df$varying, levels = 1:3, labels = c("vary_lambda0", "vary_sigma", "vary_b_ac"))


Enrm_results_df$nT2 <- factor(Enrm_results_df$nT, levels = c("20", "40", "60"), 
                              labels = c(expression(paste("20 detectors")), 
                                         expression(paste("40 detectors")), 
                                         expression(paste("60 detectors"))))
Enrm_results_df$b_ac2 <- factor(Enrm_results_df$b_ac, levels = c("-1","1","3"),
                                ordered = TRUE, labels = c(expression(paste(alpha[1]," = -1")),
                                                           expression(paste(alpha[1]," = 1")),
                                                           expression(paste(alpha[1]," = 3"))))

Enrm_results_df$method <- factor(Enrm_results_df$method, levels = c("survey", "grid", "opt_grid", "opt"),
                                 labels = c("Actual survey", "Regular grid", "Regular grid + optimal.Spacing", "Optimized"))

sm1 <- Enrm_results_df %>% filter(varying == "vary_lambda0") %>% 
  mutate(cvD = 1/sqrt(pmin(En, Er))) %>% group_by(method, error, nT2, beta0, sigma, b_ac2) %>% 
  summarize(median_cvD = median(cvD), min_cvD = min(cvD)) %>%
  ggplot(aes(x = error, y = min_cvD, colour = method, shape = method)) + geom_line() + geom_point() +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  xlab(expression(paste("True log ", lambda[0], " (as multiple of assumed log ", lambda[0],")"))) + ylab(expression(paste("Approximate CV(", hat(D),")"))) + 
  facet_grid(b_ac2 ~ nT2, labeller = label_parsed, scales = "free") + theme_bw(base_size = 14) +
  theme(legend.position="bottom",legend.title = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

sm2 <- Enrm_results_df %>% filter(varying == "vary_sigma") %>% 
  mutate(cvD = 1/sqrt(pmin(En, Er))) %>% group_by(method, error, nT2, beta0, sigma, b_ac2) %>% 
  summarize(median_cvD = median(cvD), min_cvD = min(cvD)) %>%
  ggplot(aes(x = error, y = min_cvD, colour = method, shape = method)) + geom_line() + geom_point() +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  xlab(expression(paste("True ", sigma, " (as multiple of assumed ", sigma,")"))) + ylab(expression(paste("Approximate CV(", hat(D),")"))) + 
  facet_grid(b_ac2 ~ nT2, labeller = label_parsed, scales = "free") + theme_bw(base_size = 14) +
  theme(legend.position="bottom",legend.title = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

sm3 <- Enrm_results_df %>% filter(varying == "vary_b_ac") %>% 
  mutate(cvD = 1/sqrt(pmin(En, Er))) %>% group_by(method, error, nT2, beta0, sigma, b_ac2) %>% 
  summarize(median_cvD = median(cvD), min_cvD = min(cvD)) %>%
  ggplot(aes(x = error, y = min_cvD, colour = method, shape = method)) + geom_line() + geom_point() +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  xlab(expression(paste("True ", alpha[1], " (as multiple of assumed ", alpha[1],")"))) + ylab(expression(paste("Approximate CV(", hat(D),")"))) + 
  facet_grid(b_ac2 ~ nT2, labeller = label_parsed, scales = "free") + theme_bw(base_size = 14) +
  theme(legend.position="bottom",legend.title = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("output/tost-misspec-lambda.png", sm1, width=10, height=6, dpi = 300)
ggsave("output/tost-misspec-sigma.png", sm2, width=10, height=6, dpi = 300)
ggsave("output/tost-misspec-alpha.png", sm3, width=10, height=6, dpi = 300)

######################
## supplementary material -- table of bias 
######################

load("output/posthoc-check-bias-nonuniD-noneuc.RData")
results_noneucBADSTART <- results_noneuc

load("output/posthoc-check-bias-noneuc-3.RData")
results1 <- results_noneuc
load("output/posthoc-check-bias-noneuc-4.RData")
results2 <- results_noneuc
results_noneuc <- rbind(results1, results2)

results <- results %>% 
  mutate(nT = recode(i, `1` = 20L, `2` = 20L, `3` = 20L, `4` = 40L, `5` = 40L, `6` = 40L,
                     `7` = 60L, `8` = 60L, `9` = 60L)) %>% 
  mutate(b_ac = recode(i, `1` = -1L, `2` = 1L, `3` = 3L, `4` = -1L, `5` = 1L, `6` = 3L,
                     `7` = -1L, `8` = 1L, `9` = 3L))

results_noneuc <- results_noneuc %>% 
  mutate(nT = recode(i, `1` = 20L, `2` = 20L, `3` = 20L, `4` = 40L, `5` = 40L, `6` = 40L,
                     `7` = 60L, `8` = 60L, `9` = 60L)) %>% 
  mutate(b_ac = recode(i, `1` = 1L, `2` = -1L, `3` = 3L, `4` = 1L, `5` = -1L, `6` = 3L,
                       `7` = 1L, `8` = -1L, `9` = 3L)) %>% 
  mutate(RB = 100 * (E.N - true.N) / true.N)

bias_nonu <- results %>% filter(E.N < 100) %>% 
  mutate(RB = 100 * (E.N - true.N) / true.N) %>% group_by(nT, b_ac) %>% 
  summarize(medianRB = round(median(RB),1), meanRB = round(mean(RB),1), 
            minRB = round(min(RB),1), maxRB = round(max(RB),1),
            seRB = round(sd(RB) / sqrt(n()),1)) %>% ungroup() %>%
  arrange(b_ac, nT) 

bias_noneuc <- results_noneuc %>% filter(E.N < 100) %>% 
  mutate(RB = 100 * (E.N - true.N) / true.N) %>% group_by(nT, b_ac) %>% 
  summarize(medianRB = round(median(RB),1), meanRB = round(mean(RB),1), 
            minRB = round(min(RB),1), maxRB = round(max(RB),1),
            seRB = round(sd(RB) / sqrt(n()), 1),
            mean_n = mean(n), mean_dets = mean(detections),
            mean_detsused = mean(dets_visited)) %>% ungroup() %>%
  arrange(b_ac, nT)  

bias_table <- bias_nonu %>% dplyr::select(b_ac, nT, medianRB, meanRB, seRB) %>%
  cbind(bias_noneuc %>% dplyr::select(medianRB, meanRB, seRB))

bias_table %>%
  kable("latex", caption = "Group Rows", booktabs = T) %>% 
  kable_styling() %>%
  pack_rows("", 1, 3) %>%
  pack_rows("", 4, 6) %>%
  pack_rows("", 7, 9) 
