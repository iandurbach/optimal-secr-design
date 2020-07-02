## 

library(dplyr)
library(ggplot2)
library(gridExtra)
library(kableExtra)
library(oSCR)
library(tidyr)

######################
## Figure 1: approximation accuracy for spatially-varying designs
######################

load("output/approx-checks-D123.RData")

plot1.cols <- c( "#F5940F", "#0AAA2E", "#B061C2") 
plot2.cols <- c("#EC142B", "#1435EC", "#B061C2") 
sims_sum <- sims_sum %>% mutate(gridsize = factor(gridsize, levels = c("6x6", "8x8", "10x10")))

sims_sum$cov2 <- factor(sims_sum$cov, levels = c("D1","D3", "D2"),
                        ordered = TRUE, labels = c(expression(paste(D %~% 1)),
                                                   expression(paste(D %~% -sqrt((x-bar(x))^2+(y-bar(y))^2))),
                                                   expression(paste(D %~% sqrt(xy)))))

sims_sum$lambda02 <- factor(sims_sum$lambda0, levels = c("0.05","0.1", "0.2"),
                            ordered = TRUE, labels = c(expression(paste(lambda[0],"=0.05")),
                                                       expression(paste(lambda[0],"=0.1")),
                                                       expression(paste(lambda[0],"=0.2"))))

m_abc <- data.frame(cov2 = levels(sims_sum$cov2), x = 0.1, y = 0, labs = c("(a)", "(b)", "(c)"))
p1 <- sims_sum %>% filter(lambda0 == 0.2, Dmod == "D23") %>%
  ggplot(aes(x = sp.ratio, y = 100*emp_cvD, colour = gridsize)) + 
  geom_point() + geom_line() +
  facet_grid(.~cov2, labeller = label_parsed) +
  geom_text(data = m_abc, inherit.aes = FALSE, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_point(aes(y = 100*th_cvD), shape = 2) +
  geom_line(aes(y = 100*th_cvD), linetype = 2) +
  scale_color_manual(values = c("6x6" = plot1.cols[1], "8x8" = plot1.cols[2], "10x10" = plot1.cols[3])) +
  xlab(bquote("Detector spacing ("*sigma~"units)")) + ylab(bquote("CV("*hat(D)*") %")) + theme_bw() + xlim(c(0,3.5)) +
  coord_cartesian(xlim = c(0,3.5), ylim = c(0, 70)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.margin=margin(0,0,5,0),
        legend.box.margin=margin(-10,-10,-10,-10))
p1

m_abc <- data.frame(cov2 = levels(sims_sum$cov2), x = 0.1, y = 0, labs = c("(d)", "(e)", "(f)"))
p2 <- sims_sum %>% filter(gridsize == "10x10", Dmod == "D23") %>%
  mutate(lambda0 = factor(lambda0)) %>%
  ggplot(aes(x = sp.ratio, y = 100*emp_cvD, colour = lambda0)) + 
  geom_point() + geom_line() +
  facet_grid(.~cov2, labeller = label_parsed) +
  geom_text(data = m_abc, inherit.aes = FALSE, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_point(aes(y = 100*th_cvD), shape = 2) +
  geom_line(aes(y = 100*th_cvD), linetype = 2) +
  scale_color_manual(values = plot2.cols, labels = c(expression(paste(lambda[0],"=0.05")),
                                                     expression(paste(lambda[0],"=0.1")),
                                                     expression(paste(lambda[0],"=0.2")))) +
  xlab(bquote("Detector spacing ("*sigma~"units)")) + ylab(bquote("CV("*hat(D)*") %")) + theme_bw() + xlim(c(0,3.5)) +
  coord_cartesian(xlim = c(0,3.5), ylim = c(0, 70)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.margin=margin(0,0,5,0),
        legend.box.margin=margin(-10,-10,-10,-10))
p2

ggsave("output/approx-nonuniform.png", p1, width=9, height=3, dpi = 300)
ggsave("output/approx-nonuniform2.png", p2, width=9, height=3, dpi = 300)

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
m0 <- opt_traps_all %>% filter(sigma == 3000, row_number() < 720, lambda0 == 1, buffer != 0, trap_id == 1) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ nT, labeller = labeller(nT = nT.labs)) + 
  labs(title = bquote("Effect of varying number of detectors (baseline"~ sigma~"and"~lambda[0]*","~3*sigma~"buffer)")) +
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

m_abc <- data.frame(nT = c(20,40,60), x = min(mask_df$x), y = min(mask_df$y), labs = c("(d)", "(e)", "(f)"))
m1 <-  opt_traps_all %>% filter(sigma == 3000, lambda0 == 1, buffer != 0, trap_id == 2) %>% 
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
m2 <- opt_traps_all %>% filter(sigma == 3000, lambda0 == 1, buffer == 0, trap_id == 1) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ nT, labeller = labeller(nT = nT.labs)) + 
  labs(title = bquote("Effect of zero buffer (baseline"~ sigma~"and"~lambda[0]*")")) +
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

opt_traps_all$lambda02 <- factor(opt_traps_all$lambda0, levels = c("0.5","1","1.5", "2"),
                               ordered = TRUE, labels = c(expression(paste("50% of baseline ", lambda[0])),
                                                          expression(paste("baseline ", lambda[0])),
                                                          expression(paste("150% of baseline ", lambda[0])),
                                                          expression(paste("200% of baseline ", lambda[0]))))

m_abc <- data.frame(lambda02 = levels(opt_traps_all$lambda02)[-2], x = min(mask_df$x), y = min(mask_df$y), labs = c("(j)", "(k)", "(l)"))
m3 <- opt_traps_all %>% filter(sigma == 3000, nT == 60, trap_id == 1, lambda0 != 1) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ lambda02, labeller = label_parsed) + 
  labs(title = bquote("Effect of varying encounter rate"~lambda[0]~"(60 detectors, baseline"~ sigma*","~3*sigma~"buffer)")) +
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

opt_traps_all$sigma2 <- factor(opt_traps_all$sigma, levels = c("1500","3000", "4500", "6000"),
                               ordered = TRUE, labels = c(expression(paste("50% of baseline ", sigma)),
                                                          expression(paste("baseline ", sigma)),
                                                          expression(paste("150% of baseline ", sigma)),
                                                          expression(paste("200% of baseline ", sigma))))

m_abc <- data.frame(sigma2 = levels(opt_traps_all$sigma2)[-2], x = min(mask_df$x), y = min(mask_df$y), labs = c("(m)", "(n)", "(o)"))
m4 <- opt_traps_all %>% filter(nT == 60, trap_id == 1, lambda0 == 1, sigma != 3000) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ sigma2, labeller = label_parsed) + 
  labs(title = bquote("Effect of varying animal movement parameter"~sigma~"(60 detectors, baseline"~ lambda[0]*","~3*sigma~"buffer)")) +
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
m4b <- opt_traps_all %>% filter(row_number() > 720, nT == 60, lambda0 == 1, sigma == 3000) %>%
  ggplot(aes(x = x, y = y)) + 
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_tile(data = mask_df, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ dt) + 
  labs(title = bquote("Effect of detector type (60 detectors, 5 occasions, "~ lambda[0]==1/5*", baseline"~sigma*","~3*sigma~"buffer)")) +
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
m4b

ggsave("output/tost-uni-designs.png", 
       arrangeGrob(m0,m1,m2,m3,m4,m4b, ncol = 1), width=10, height=12, dpi = 300)

######################
## Maps of example designs for non-Uniform D, euclidean
######################

load("output/Tost_examples_nonuniD.Rdata")  

opt_traps_D1_all$nT2 <- factor(opt_traps_D1_all$nT, levels = c("40", "60"), 
                            labels = c(expression(paste("40 detectors")), 
                                       expression(paste("60 detectors"))))
opt_traps_D1_all$b_ac2 <- factor(opt_traps_D1_all$b_ac, levels = c("-1", "0", "1","3"),
                              ordered = TRUE, labels = c(expression(paste(alpha[D]," = -1")),
                                                         expression(paste(alpha[D]," = 0")),
                                                         expression(paste(alpha[D]," = 1")),
                                                         expression(paste(alpha[D]," = 3"))))
opt_traps_D1_all$b_det2 <- factor(opt_traps_D1_all$b_det, levels = c("-0.75", "0", "0.75", "1.5"),
                                 ordered = TRUE, labels = c(expression(paste(alpha[lambda[0]]," = -0.75")),
                                                            expression(paste(alpha[lambda[0]]," = 0")),
                                                            expression(paste(alpha[lambda[0]]," = 0.75")),
                                                            expression(paste(alpha[lambda[0]]," = 1.5"))))
                                                            
x_m5a <- opt_traps_D1_all
x_m5a$xx <- "xx"
x_m5a$xx <- factor(x_m5a$xx, levels = c("xx"), labels = expression(paste("uniform ", lambda[0])))
m_abc <- data.frame(b_ac2 = levels(x_m5a$b_ac2)[-2], x = min(mask_df$x), y = min(mask_df$y), labs = c("(a)", "(b)", "(c)"))
m5a <- x_m5a %>% filter(nT == 40, trap_id == 1, b_ac != 0, b_det == 0) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_tile(data = mask_df, aes(fill = stdGC), colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ b_ac2, labeller = label_parsed) +
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  labs(title = bquote("Designs with non-uniform activity center density (40 detectors, baseline conditions)")) +
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
m5a

x_m5b <- opt_traps_D1_all
x_m5b$xx <- "xx"
x_m5b$xx <- factor(x_m5b$xx, levels = c("xx"), labels = expression(paste("uniform ", D)))
m_abc <- data.frame(b_det2 = levels(x_m5b$b_det2)[-2], x = min(mask_df$x), y = min(mask_df$y), labs = c("(d)", "(e)", "(f)"))
m5b <- x_m5b %>% filter(nT == 40, trap_id == 1, b_ac == 0, b_det != 0) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_tile(data = alltraps_df, aes(fill = -x), colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ b_det2, labeller = label_parsed) +
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  labs(title = bquote("Designs with spatial detector covariate and uniform density ("*alpha[D]==0*")")) +
  scale_fill_distiller(palette = "BrBG") + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
m5b

x_m5c <- opt_traps_D1_all
x_m5c$xx <- "xx"
x_m5c$xx <- factor(x_m5c$xx, levels = c("xx"), labels = expression(paste(alpha[1],"=1")))
m_abc <- data.frame(b_det2 = levels(x_m5c$b_det2)[-2], x = min(mask_df$x), y = min(mask_df$y), labs = c("(g)", "(h)", "(i)"))
m5c <- x_m5c %>% filter(nT == 40, trap_id == 1, b_ac != 0, b_det != 0) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_tile(data = mask_df, aes(fill = stdGC), colour = "black") +
  geom_point(size = 1, colour = "red") + 
  facet_grid(. ~ b_det2, labeller = label_parsed) +
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  labs(title = bquote("Designs with spatial detector covariate and non-uniform density ("*alpha[D]==1*")")) +
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
m5c

ggsave("output/tost-designs-nonU.png", 
       arrangeGrob(m5a,m5b,m5c, ncol = 1), width=10, height=6, dpi = 300)

######################
## Distributions of inter-detector distances in optimal designs
######################

load("output/Tost_examples_all_D0.Rdata")  
load("output/Tost_examples_nonopt_D0.Rdata")  

all_traps <- rbind(opt_traps_all %>% filter(row_number() <= 120) %>% mutate(method = "mnr"),
                   grid_traps_all %>% filter(row_number() <= 120) %>% mutate(method = "grid") %>% dplyr::select(-dist2pt), 
                   opt_grid_traps_all%>% filter(row_number() <= 120) %>% mutate(method = "optgrid") %>% dplyr::select(-dist2pt))

# add a camera id to each array
all_traps <- all_traps %>% group_by(method, nT) %>% mutate(cam_id = row_number()) %>% ungroup()

# extract distance to closest detector for each detector 
myseq <- seq(2000,16000,1000)
closest_detectors <- all_traps %>% 
  group_by(method, nT) %>% 
  mutate(nearest_detdist = apply(e2dist(cbind(x,y), cbind(x,y)), 1, function(x){sort(x, FALSE)[2]}),
         nearest_cat = cut(nearest_detdist, breaks = 8, right = FALSE),
         mean_detdist = apply(e2dist(cbind(x,y), cbind(x,y)), 1, mean),
         mean_cat = cut(mean_detdist, breaks = 8, right = FALSE)) %>%
  arrange(method, nT, nearest_detdist) %>%
  ungroup() %>% group_by(method, nT) %>% 
  mutate(Fx = row_number() / n()) %>% ungroup()

closest_detectors$nT2 <- factor(closest_detectors$nT, levels = c("20", "40", "60"), 
                                labels = c(expression(paste("20 detectors")), 
                                           expression(paste("40 detectors")), 
                                           expression(paste("60 detectors"))))
closest_detectors$method2 <- factor(closest_detectors$method, levels = c("grid", "optgrid", "mnr"),
                                    labels = c("Regular grid", "Regular grid + optimal.Spacing", "Min(n,r)"))


m_abc <- data.frame(nT2 = levels(closest_detectors$nT2), x = 12, y = 0, labs = c("(a)", "(b)", "(c)"))
m8b <- closest_detectors %>% filter(nT %in% c(20,40,60)) %>%
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
m8b

cdf_data <- data.frame(dists = as.numeric(), Fx = as.numeric(), nT = as.integer(), method = as.character())
for(m in unique(all_traps$method)){
  for(i in c(20,40,60)){
    if(!((m == "survey")&(i == 60))){
        allxt <- c()
        for(j in 1:1){
          xx <- all_traps %>% filter(nT == i, trap_id == j, method == m)
          ds <- e2dist(cbind(xx$x, xx$y),cbind(xx$x, xx$y))
          diag(ds) <- NA
          allxt <- c(allxt, unlist(ds[!is.na(ds)]))
        }
        
        cdf_data <- rbind(cdf_data, 
                          data.frame(dists = sort(allxt), 
                                     Fx = seq(0, 1, length.out = length(allxt)),
                                     nT = i,
                                     method = m)) 
    }
  }
}

cdf_data$nT2 <- factor(cdf_data$nT, levels = c("20", "40", "60"), 
                       labels = c(expression(paste("20 detectors")), 
                                  expression(paste("40 detectors")), 
                                  expression(paste("60 detectors"))))
cdf_data$method2 <- factor(cdf_data$method, levels = c("grid", "optgrid", "mnr"),
                                    labels = c("Regular grid", "Regular grid + optimal.Spacing", "Min(n,r)"))

m_abc <- data.frame(nT2 = levels(cdf_data$nT2), x = 100, y = 0, labs = c("(d)", "(e)", "(f)"))
m9b <- cdf_data %>% 
  ggplot(aes(x = dists/1000, y = Fx)) + 
  geom_path(aes(colour = method2)) + facet_grid(. ~ nT2, labeller = label_parsed) +
  geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 1, vjust = 0) + 
  xlab("Between-detector distance (km)") + ylab("Cumulative probability") + 
  theme_bw(base_size = 14) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  theme(legend.position="bottom",legend.title = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        legend.margin=margin(0,0,5,0), legend.box.margin=margin(-10,-10,-10,-10))
m9b

ggsave("output/tost-trap-dists.png", 
       arrangeGrob(m8b, m9b, ncol = 1, heights = c(5,6)), width=10, height=6, dpi = 300)


######################
## Covariate coverage plots
######################

load("output/Tost_examples_all_D0.Rdata")  
load("output/Tost_examples_nonuniD.Rdata")  

opt_traps_D1_all <- opt_traps_D1_all %>% left_join(mask_df, by = c("x", "y"))

opt_traps_D1_all <- opt_traps_D1_all %>% mutate(y = (nT - 20) / 200)

opt_traps_D1_all$nT2 <- factor(opt_traps_D1_all$nT, levels = c("20", "40", "60"), 
                               labels = c("20 detectors", 
                                          "40 detectors", 
                                          "60 detectors"))
opt_traps_D1_all$b_ac2 <- factor(opt_traps_D1_all$b_ac, levels = c("-1", "0", "1","3"),
                                 ordered = TRUE, labels = c(expression(paste(alpha[D]," = -1")),
                                                            expression(paste(alpha[D]," = 0")),
                                                            expression(paste(alpha[D]," = 1")),
                                                            expression(paste(alpha[D]," = 3"))))

opt_traps_D1_all$b_det2 <- factor(opt_traps_D1_all$b_det, levels = c("-0.75", "0", "0.75","1.5"),
                                 ordered = TRUE, labels = c(expression(paste(alpha[lambda[0]]," = -0.75")),
                                                            expression(paste(alpha[lambda[0]]," = 0")),
                                                            expression(paste(alpha[lambda[0]]," = 0.75")),
                                                            expression(paste(alpha[lambda[0]]," = 1.5"))))
# D covariate coverage

linedata <- expand.grid(stdGC = seq(from = min(mask_df$stdGC), to = max(mask_df$stdGC), length.out = 100),
                        b_ac = c(-1, 1, 3)) %>%
  group_by(b_ac) %>%
  mutate(y = exp(b_ac * stdGC)) %>%
  mutate(y = y / max(y)) %>%
  ungroup()

linedata$b_ac2 <- factor(linedata$b_ac, levels = c("-1", "1","3"),
                         ordered = TRUE, labels = c(expression(paste(alpha[D]," = -1")),
                                                    expression(paste(alpha[D]," = 1")),
                                                    expression(paste(alpha[D]," = 3"))))


m10a <- linedata %>% ggplot(aes(x = stdGC, y = y)) + 
  geom_line() + 
  geom_point(data = opt_traps_D1_all %>% filter(b_det == 0), inherit.aes = FALSE, aes(x = stdGC, y = y, pch = factor(nT2), colour = factor(nT2))) +
  facet_grid(. ~ b_ac2, labeller = label_parsed) +
  xlab("Terrain ruggedness (standardized)") + ylab("Density (standardized)") + 
  theme_bw(base_size = 14) +
  scale_color_brewer(palette = "Dark2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme(legend.position="none",legend.title = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
m10a

# lambda0 covariate coverage

linedata <- expand.grid(x = seq(from = min((mask_df$x)), to = max((mask_df$x)), length.out = 100),
                        b_det = c(-0.75, 0.75, 1.5)) %>%
  group_by(b_det) %>%
  mutate(y = exp(b_det * -scale(x))) %>%
  mutate(y = y / max(y)) %>%
  ungroup()

linedata$b_det2 <- factor(linedata$b_det, levels = c("-0.75", "0.75","1.5"),
                          ordered = TRUE, labels = c(expression(paste(alpha[lambda[0]]," = -0.75")),
                                                     expression(paste(alpha[lambda[0]]," = 0.75")),
                                                     expression(paste(alpha[lambda[0]]," = 1.5"))))

m10b <- linedata %>% ggplot(aes(x = (x - min(mask_df$x))/(max(mask_df$x) - min(mask_df$x)), y = y)) + 
  geom_line() + 
  geom_point(data = opt_traps_D1_all %>% filter(b_ac == 0, b_det != 0), 
             inherit.aes = FALSE, aes(x = (x - min(mask_df$x))/(max(mask_df$x) - min(mask_df$x)), 
                                      y = y, pch = factor(nT2), colour = factor(nT2))) +
  facet_grid(. ~ b_det2, labeller = label_parsed) +
  xlab("Longitude (standardized)") + ylab(bquote(lambda[0]~"(standardized)")) + 
  theme_bw(base_size = 14) +
  scale_color_brewer(palette = "Dark2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme(legend.position="bottom",legend.title = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        legend.margin=margin(0,0,5,0), legend.box.margin=margin(-10,-10,-10,-10))
m10b

ggsave("output/tost-coverage.png", arrangeGrob(m10a, m10b, ncol = 1), width=10, height=6, dpi = 300)

######################
## supplementary material -- tost
######################
library(secr)
library(secrdesign)

load("data/Tost.RData")
traps <- data.frame(x = TostCams$x, y = TostCams$y)
pt <- traps %>% 
  ggplot(aes(x = x, y = y)) + 
  #geom_text(data = m_abc, aes(x = x, y = y, label = labs), hjust = 0.5, vjust = 0) + 
  geom_tile(data = mask_df, aes(fill = stdGC), colour = "black") +
  geom_point(colour = "red") + 
  scale_fill_viridis_c() + coord_equal() + theme_bw(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
pt

ggsave("output/tost.png", pt, width=8, height=3, dpi = 300)

# Enrm for Tost
load("data/Tost.RData")
mask <- all_masks[[1]]
traps <- read.traps(data = data.frame(x = TostCams$x, y = TostCams$y), detector = "count")
Enrm(D = 2/10000, traps, mask, list(lambda0 = 1, sigma = 3000), noccasions = 1)

######################
## supplementary material -- plot of sensitivity to misspecification 
######################

load("output/posthoc-check-robustness-to-misspec.RData")

Enrm_results_df <- data.frame(Enrm_results)
Enrm_results_df$method <- factor(Enrm_results_df$method, levels = 1:3, labels = c("opt", "grid", "opt_grid"))
Enrm_results_df$varying <- factor(Enrm_results_df$varying, levels = 1:3, labels = c("vary_lambda0", "vary_sigma", "vary_b_ac"))


Enrm_results_df$nT2 <- factor(Enrm_results_df$nT, levels = c("40", "60"), 
                              labels = c(expression(paste("40 detectors")), 
                                         expression(paste("60 detectors"))))
Enrm_results_df$b_ac2 <- factor(Enrm_results_df$b_ac, levels = c("-1","1","3"),
                                ordered = TRUE, labels = c(expression(paste(alpha[D]," = -1")),
                                                           expression(paste(alpha[D]," = 1")),
                                                           expression(paste(alpha[D]," = 3"))))

Enrm_results_df$method <- factor(Enrm_results_df$method, levels = c("grid", "opt_grid", "opt"),
                                 labels = c("Regular grid", "Regular grid + optimal.Spacing", "Min(n,r)"))

sm1 <- Enrm_results_df %>% filter(varying == "vary_lambda0") %>% 
  mutate(cvD = 1/sqrt(pmin(En, Er))) %>% group_by(method, error, nT2, beta0, sigma, b_ac2) %>% 
  summarize(median_cvD = median(cvD), min_cvD = min(cvD)) %>%
  ggplot(aes(x = error, y = min_cvD, colour = method, shape = method)) + geom_line() + geom_point() +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  xlab(expression(paste("True log ", lambda[0], " (as multiple of assumed log ", lambda[0],")"))) + ylab(expression(paste("Approximate CV(", hat(D),")"))) + 
  facet_grid(nT2 ~ b_ac2, labeller = label_parsed, scales = "free") + theme_bw(base_size = 14) +
  theme(legend.position="bottom",legend.title = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        legend.margin=margin(0,0,5,0), legend.box.margin=margin(-10,-10,-10,-10))

sm2 <- Enrm_results_df %>% filter(varying == "vary_sigma") %>% 
  mutate(cvD = 1/sqrt(pmin(En, Er))) %>% group_by(method, error, nT2, beta0, sigma, b_ac2) %>% 
  summarize(median_cvD = median(cvD), min_cvD = min(cvD)) %>%
  ggplot(aes(x = error, y = min_cvD, colour = method, shape = method)) + geom_line() + geom_point() +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  xlab(expression(paste("True ", sigma, " (as multiple of assumed ", sigma,")"))) + ylab(expression(paste("Approximate CV(", hat(D),")"))) + 
  facet_grid(nT2 ~ b_ac2, labeller = label_parsed, scales = "free") + theme_bw(base_size = 14) +
  theme(legend.position="bottom",legend.title = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        legend.margin=margin(0,0,5,0), legend.box.margin=margin(-10,-10,-10,-10))

sm3 <- Enrm_results_df %>% filter(varying == "vary_b_ac") %>% 
  mutate(cvD = 1/sqrt(pmin(En, Er))) %>% group_by(method, error, nT2, beta0, sigma, b_ac2) %>% 
  summarize(median_cvD = median(cvD), min_cvD = min(cvD)) %>%
  ggplot(aes(x = error, y = min_cvD, colour = method, shape = method)) + geom_line() + geom_point() +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  xlab(expression(paste("True ", alpha[D], " (as multiple of assumed ", alpha[D],")"))) + ylab(expression(paste("Approximate CV(", hat(D),")"))) + 
  facet_grid(nT2 ~ b_ac2, labeller = label_parsed, scales = "free") + theme_bw(base_size = 14) +
  theme(legend.position="bottom",legend.title = element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        legend.margin=margin(0,0,5,0), legend.box.margin=margin(-10,-10,-10,-10))

ggsave("output/tost-misspec-lambda.png", sm1, width=10, height=5, dpi = 300)
ggsave("output/tost-misspec-sigma.png", sm2, width=10, height=5, dpi = 300)
ggsave("output/tost-misspec-alpha.png", sm3, width=10, height=5, dpi = 300)

######################
## supplementary material -- table of bias 
######################

load("output/mnr-res-D1.RData")

b_ac <- c(-1, -1, 1, 1, 3, 3, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
b_d <- c(0, 0, 0, 0, 0, 0, -0.75, -0.75, 0.75, 0.75, 1.5, 1.5, -0.75, -0.75, 0.75, 0.75, 1.5, 1.5)
nT <- rep(c(40,60), times = 9)
dens_per_100km2 <- 2

mnr_sim_sum <- mnr_sim[[1]]$mod_summary
mnr_sim_sum$id <- 1
mnr_sim_sum$b_ac <- b_ac[1]
mnr_sim_sum$b_d <- b_d[1]
mnr_sim_sum$nT <- nT[1]
mnr_sim_sum$Dt <- dens_per_100km2
for(i in 2:18){
  mnr_sim_sum_t <- mnr_sim[[i]]$mod_summary
  mnr_sim_sum_t$id <- i
  mnr_sim_sum_t$b_ac <- b_ac[i]
  mnr_sim_sum_t$b_d <- b_d[i]
  mnr_sim_sum_t$nT <- nT[i]
  mnr_sim_sum_t$Dt <- dens_per_100km2
  mnr_sim_sum <- rbind(mnr_sim_sum, mnr_sim_sum_t)
}

bias_nonu <- mnr_sim_sum %>%  filter(E.N < 4 * true.N) %>%
  mutate(RB = 100 * (E.N - true.N) / true.N) %>% group_by(id) %>%
  summarize(nT = mean(nT), b_ac = mean(b_ac), b_d = mean(b_d), 
            medianRB = round(median(RB),1), meanRB = round(mean(RB),1), 
            minRB = round(min(RB),1), maxRB = round(max(RB),1),
            seRB = round(sd(RB) / sqrt(n()), 1),
            mean_n = mean(n), mean_r = mean(r), mean_dets = mean(detections),
            mean_detsused = mean(dets_visited), 
            n = n()) %>% ungroup() %>%
  mutate(obs_nr = paste0(round(mean_n, 1), " (", round(mean_r, 1), ")")) 

bias_table <- bias_nonu %>% dplyr::select(b_ac, nT, obs_nr, medianRB, meanRB, seRB) %>% arrange(nT)

bias_table %>%
  kable("latex", caption = "Group Rows", booktabs = T) %>% 
  kable_styling() %>%
  pack_rows("40 detector array", 1, 3) %>%
  pack_rows("", 4, 6) %>%
  pack_rows("", 7, 9) %>%
  pack_rows("60 detector array", 10, 12, latex_gap_space = "1.3em") %>%
  pack_rows("", 13, 15) %>%
  pack_rows("", 16, 18) 

# region.N suspect for varying lambda0, compute RB on D instead
mnr_sim_sum %>% filter(id %in% 7:12) %>% mutate(true.D = log(true.N/maskarea)) %>%
  mutate(RB = 100 * (D - true.D) / true.D) %>% group_by(id) %>%
  summarize(nT = mean(nT), b_ac = mean(b_ac), b_d = mean(b_d), 
            medianRB = round(median(RB),1), meanRB = round(mean(RB),1), 
            minRB = round(min(RB),1), maxRB = round(max(RB),1),
            seRB = round(sd(RB) / sqrt(n()), 1),
            mean_n = mean(n), mean_r = mean(r), mean_dets = mean(detections),
            mean_detsused = mean(dets_visited), 
            n = n()) %>% ungroup() %>%
  mutate(obs_nr = paste0(round(mean_n, 1), " (", round(mean_r, 1), ")")) 


