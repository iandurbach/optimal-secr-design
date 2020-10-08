# statespace: state space points
# alltraps:  set of all possible trap locations
# ntraps:     number of trap locations we WANT
# beta0:      estimated bld
# sigma:      estimated sigma
# D:          density per ha
# occasions:  the number of survey occasions
# detector:   one of "count", "proximity", or "multi" (see secr for details).
# crit: which criterion to use (1 = En, 2 = Em, 3 = min(En,Er))
# ngen:       number of generations to run the GA for
# popsize:    number of designs in each population
# seed:       set a random seed for reproducibility
# USE THE FOLLOWING WITH CAUTION!
# pen_gridsigma:   penalties for bad trap spacing are based on comparing min(n,r) grid
#                  to what you would get from a X-sigma regular grid. X = pen_gridsigma
# pen_wt:          penalty for spacings that are too different from regular grid

scrdesignGAenrm <- function(statespace = NULL,
                            alltraps = NULL, 
                            ntraps = 9, 
                            lambda0 = exp(-0.6), 
                            sigma = 2, 
                            D = 1,
                            occasions = 1,
                            detector = "proximity",
                            crit = 3, 
                            pen_gridsigma = 2,
                            pen_wt = 0, 
                            #N = 100, 
                            verbose=1,
                            ngen = 50,
                            popsize = 50,
                            seed = NULL,
                            ...){
  
  if(!is.null(seed)) set.seed(seed)
  if(is.null(statespace)) stop("Must supply a 'statespace' object (coords of the study area)")
  if(is.null(alltraps))   stop("Must supply a 'alltraps' object (coords of all possible trap locations)")
  
  # find distribution of trap spacings on a close to regular grid, to ensure later optimized grid has spaced
  # enough detectors sufficiently far apart to get low var(sigma)
  
  # hacky way to make a polygon just bigger than boundary of possible detector points
  sp <- diff(sort(unique(alltraps[,1])))[1] 
  pg <- data.frame(alltraps) %>% st_as_sf(coords = c("x", "y")) %>% st_buffer(dist = sp, endCapStyle = "SQUARE") %>% st_union() 
  pg <- pg %>% st_buffer(dist = -sp * 0.4, endCapStyle = "SQUARE")
  
  # place a grid over the area, with cells X sigma apart
  my_grid <- st_make_grid(pg, cellsize = c(pen_gridsigma * sigma, pen_gridsigma * sigma), what = "centers") %>% st_intersection(pg)
  grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
  
  # choose a subset of nT of these, starting from a random trap and then choosing nT nearest neighbours
  all_grid_traps <- list()
  
  grid_traps <- grid_traps_full
  xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
  grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
  grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= ntraps)
  grid_traps <- read.traps(data = grid_traps, detector = "count")
  
  # find out how many detector pairs are between 2-3 and 3-4 sigma apart; these are important for sigma estimation         
  td <- e2dist(grid_traps, grid_traps) 
  n2ish <- sum((td <= 2.499 * sigma) & (td > 1.499 * sigma))
  n3ish <- sum((td <= 3.499 * sigma) & (td > 2.499 * sigma))
  n4ish <- sum((td <= 4.449 * sigma) & (td > 3.499 * sigma))
  
  # set reasonable targets for optimized grids
  g_n2ish <- 0 # n2ish * 0.8
  g_n3ish <- n3ish * 1
  g_n4ish <- n4ish * 1
  
  des <- kofnGA(n = nrow(alltraps), 
                k = ntraps, 
                OF = scrdesignOFenrm,
                verbose = verbose,
                ...,
                alltraps = alltraps,
                statespace = statespace,
                lambda0 = lambda0,
                sigma = sigma,
                D = D,
                occasions = occasions,
                detector = detector,
                ngen = ngen,
                popsize = popsize,
                g_n234 = c(g_n2ish, g_n3ish, g_n4ish),
                pen_wt = pen_wt, 
                crit=crit)
  optimaltraps <- alltraps[des$bestsol,]
  scrdesign <- list(des=des, statespace=statespace, alltraps=alltraps, optimaltraps=optimaltraps, 
                    sigma=sigma, beta0=exp(lambda0*occasions), lambda0=lambda0)
  class(scrdesign) <- "scrdesign"
  return(scrdesign)
}