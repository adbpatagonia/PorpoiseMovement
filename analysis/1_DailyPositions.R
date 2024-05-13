# step 1: choose initial Position:  day t
# step 2: propose new positions: day t + 1 - must stay on the range.
# propose multiple positions
# Step 3: Compute acceptance probability for this candidate
# calculate weight for each member of the forecast ensemble
# Step 4: Choose whether to accept candidate
# select a single member of the ensemble
s1 <- Sys.time()
# libraries ------
library(tidyverse)
library(data.table)
library(sf)
library(raster)
library(lubridate)
# library(CircStats) # use for random draws of von Mises dist
library(circular) # use for random draws of von Mises dist
library(akima) # use for 2D linear interpolation
library(units)
library(mapview)
library(dismo)
library(osmdata)
library(stars)
library(gganimate)
library(proj4)
library(geosphere)
library (rgdal)

mapviewOptions(basemaps = c("CartoDB.Positron", 'Esri.WorldImagery', "OpenStreetMap"))

# functions -----
source("R/split_poly.r")
source("R/compile_rmd.R")

# parameters -----
## nburnin ----
# number of days to discard from each month so that we avoid sampling from the transient positions
# while whales move away from the random initial position
nburnin <- 0
## ncandidates ----
# number of candidate new positions
ncandidates <- 10

## n agents ----
# simulate n_agents whales
n_agents <- 10

## n areas ----
# number of areas in which to divide the North Sea for step 1, i.e. choose initial position of agents
nareas <- 10

## turning angle parameter ----
# kappa parameter of the von Mises distribution. Defines degree to which angles clump around zero
# very low values allow agents to turn anywhere from 0 to 360
kap <- 0.01

## displacements -----
# source the script that extracts mean and sd of displacements (in km) from dataset
# and calculates the
# the script also plots the distributions
source("analysis/0_StepLength.r")
scale_gamma <- pull(dis.sum[Location == "combined",.(scale_gamma)])
shape_gamma <- pull(dis.sum[Location == "combined",.(shape_gamma)])


## days in month ----
# From Gilles et al (2016)
# Survey effort was very low and
# patchy during winter; therefore, data collected
# during December–February were excluded from
# further analysis. The three remaining seasons
# were defined according to the meteorological
# start in north temperate zones (i.e., Mar–May
#                                 = spring; Jun.–Aug. = summer; Sep.–Nov. = fall)
# Start the year in March 1 (day 60)
month_days <-
  data.table(month =3:11,
             ndays = c( 31, 30, 31, 30, 31, 31, 30, 31, 30))
month_days[, firstday := 0]
month_days[1, firstday := 60]
for(i in 2:9){
  month_days$firstday[i] <- month_days$firstday[i-1] + month_days$ndays[i - 1]
}
month_days[1:8, lastday := month_days[2:9, .(firstday)] -1]
month_days[9, lastday := 365 - 31]

# repeat the table including all years
# this wil be used to switch the underlying density rasters at the start of each month
year_month_days <- month_days
year_month_days[,year := 2050]

year_month_days <- year_month_days[month %in% c(6, 9)]

## ndays ------
# ndays is 365 to cover the entire average year
nyears <- 1
# ndays <- 365 * nyears
sum(month_days$ndays)
ndays <- sum(month_days$ndays) * nyears

# read data -----

# step 1: choose initial Position:  day t ----
## divide the area in big blocks ----
# read density raster as sf object
# Given the description from the Gilles et al paper, I downloaded the data that Anita sent to
# data-raw\AnitaGilles\commondata\ramboll_shapes\
# and copied the shapefiles to
# data\shp\
# and renamed the seasons to the month when the seasons begin (so that my loop code below works)
this.density.sf <- st_as_sf(
  st_read("data/shp/hp_prediction_Mar.shp"))
# obtain crs
NS_shp <- readOGR("data/shp/hp_prediction_Mar.shp") #read in file
ns_crs <- proj4string(NS_shp)

# set crs equal to North Sea CRS: ns_crs
this.density.sf <-
  st_transform(this.density.sf, crs = ns_crs)


this.density.sf <- this.density.sf %>%
  rename(ID_shpfile = Id)

# divide the polygon containing the densities into n_areas of equal area
this.density.blocks <- split_poly(this.density.sf, n_areas = nareas)

# mapview(this.density.blocks)

## Calculate proportional abundance in each of the blocks ----
# density is in animals/km2. Area is in m2. Transform to km2
this.density.blocks$area <- units::set_units(st_area(st_transform(this.density.blocks)), km^2)

# the column that contains density is called the same as the raster file
# change it to "density"
colnames(this.density.blocks)[which(grepl(pattern = "Avg", x = colnames(this.density.blocks)))] <- 'density'

this.density.blocks <- this.density.blocks %>%
  mutate(nw = density * area)
# mapview(this.density.blocks)

this.density.summary <-
  this.density.blocks %>%
  group_by(id) %>%
  summarise(nw = sum(nw)) %>%
  mutate(prop = nw / sum(nw))
# mapview(this.density.summary)

## Start nwhales. The number of whales in each of those blocks will be proportional to the proportional abundance in the block  ----
# calculate how many agents to start in each polygon. Use floor to obtain integers
this.density.summary <- this.density.summary %>%
  mutate(agents = floor(prop * n_agents))

# given that I used floor, the numbers allocated above result in less than n_agents
# need to assign the difference
dif_agents <- n_agents - as.numeric(sum(this.density.summary$agents))

# create a vector with these agents that have not been assigned and shuffle them
extra_agents <- sample(
  c(rep(1, dif_agents), rep(0, nareas - dif_agents)))
# assign extra agents to this.density.summary$agents
this.density.summary$agents <- as.numeric(this.density.summary$agents) + extra_agents

this.density.summary <- this.density.summary %>% as.data.table()

## Within the block, choose a starting position based on the density surface, like in Ruth's paper ----
# within each of the areas defined above
# generate random start point based on density surface
# initial condition (initial position)

# create the object from which the positions will be sampled
# we will sample positions from the this.density.blocks object
this.z <- this.density.blocks

# obtain lat and lon
this.z$point <- st_centroid(this.z$geometry)
this.z$x <-  unlist(lapply(this.z$point,"[",1))
this.z$y <-  unlist(lapply(this.z$point,"[",2))

this.z <- as.data.table(this.z)

init.positions <- data.table()
for (i in 1:nareas){
  # number of whales in this block
  this.nwales <- pull(this.density.summary[id == i,.( agents)])
  # subset the z dataset
  this.iter.z <- this.z[id == i,]
  # add nrow column
  this.iter.z$nrow <- 1:nrow(this.iter.z)
  # randomply sample this.nwales row numbers, proportional to column density
  inds <- sample(x = this.iter.z$nrow,
                 size = this.nwales,
                 replace = T,
                 prob = this.iter.z$density)
  # assign the initial positions using the row numbers, and subset columns
  init.pos <- this.iter.z[inds, .(id_block = id, x, y, z = density)]
  # accumulate init positions
  init.positions <- rbindlist(l = list(init.positions, init.pos))
  rm(init.pos)
}
# end step 1 -----


# step 2: propose new positions: day t + 1  -----
# must stay on the range.
# propose multiple positions

## set up positions object -----
# this object will be the final output
positions <- data.table(day = as.numeric(rep(NA, (nburnin  + ndays) * n_agents)),
                        agent = as.numeric(rep(NA, (nburnin  + ndays) * n_agents)),
                        x = as.numeric(rep(NA, (nburnin  + ndays) * n_agents)),
                        y = as.numeric(rep(NA, (nburnin  + ndays) * n_agents)),
                        z = as.numeric(rep(NA, (nburnin  + ndays) * n_agents))
)

# set up positions object with initial positions for all agents
# setting
positions[1:n_agents, day := -nburnin + month_days$firstday[1]]
positions[day == -nburnin + month_days$firstday[1], agent := 1:n_agents]
positions[day == -nburnin + month_days$firstday[1], x := init.positions$x]
positions[day == -nburnin + month_days$firstday[1], y := init.positions$y]
positions[day == -nburnin + month_days$firstday[1], z := init.positions$z]

## read whale density first month -----
# these data will be used to interpolate the density in the candidate locations
# first month
# this.density.raster <- raster::raster("data/Sperm Whales/MonthlyDensityRasters/gomx_sperm_whale_Jan_2003_density_v2.img")
# mapview(this.density.raster)

# transform sf object to data.frame
this.density <- as.data.frame(this.density.sf )
this.density$point <- st_centroid(this.density$geometry)
this.density$x <-  unlist(lapply(this.density$point,"[",1))
this.density$y <-  unlist(lapply(this.density$point,"[",2))
this.density <- this.density %>%
  rename(Layer_1 = colnames(this.density)[which(grepl(pattern = "Avg", x = colnames(this.density)))]) %>%
  dplyr::select(x,y,Layer_1) %>%
  arrange(x,y,Layer_1)

# create all random step lengths and turning angles that will be used in the simulation
# create them outside looping over days to reduce run-time
ran_steplength_all <- array(data = rgamma(n =  ncandidates * n_agents * ndays, shape = shape_gamma, scale = scale_gamma),
                            dim = c(ncandidates, n_agents, ndays))
ran_angles_all <- suppressWarnings(
  array(data = rvonmises(n = ncandidates * n_agents * ndays, mu = 0, kappa = kap, control.circular=list(units="radians")),
        dim = c(ncandidates, n_agents, ndays)))

## start loop over days ----
for (dy in (-nburnin + month_days$firstday[1] + 1):(month_days$firstday[1] + ndays - 1)){
  ###  switch  whale density  ----
  # switch the underlying whale density in the first day of each month
  (print(dy))
  # apply the switch only the first day of each month
  if(any(dy == year_month_days$firstday)){
    mh <- unique(year_month_days$month[which(dy == year_month_days$firstday)])
    this.mh <- switch(mh, "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

    this.yr <- unique(year_month_days$year[which(dy == year_month_days$firstday)])

    # read the new shape file
    this.density.sf <- st_as_sf(
      st_read(paste0("data/shp/hp_prediction_", this.mh, ".shp")))

    # set crs
    this.density.sf <-
      st_transform(this.density.sf, crs = ns_crs)
    this.density.sf <- this.density.sf %>%
      rename(ID_shpfile = Id)

    # transform sf object to data.frame
    this.density <- as.data.frame(this.density.sf )
    this.density$point <- st_centroid(this.density$geometry)
    this.density$x <-  unlist(lapply(this.density$point,"[",1))
    this.density$y <-  unlist(lapply(this.density$point,"[",2))
    this.density <- this.density %>%
      rename(Layer_1 = colnames(this.density)[which(grepl(pattern = "Avg", x = colnames(this.density)))]) %>%
      dplyr::select(x,y,Layer_1) %>%
      arrange(x,y,Layer_1)


    print(cat("day:",dy, "month:",this.mh, "year:",this.yr, "file:", paste0("data/shp/hp_prediction_", this.mh, ".shp")))

  }

  # include a while() loop to ensure that there is at least two valid candidate per agent
  # create a dummy candidates object to initiate the while() loop
  candidates <- data.table(agent = 0, candidate = 0)
  while(!min(candidates[,.N, keyby = agent][[2]]) > 1){
    ### draw steplengths and turning angles ----

    # get the matrix n_agents*ncandidates of random step lengths, and shuffle them randomly
    # the shuffling is done to allow new values in each iteration of the while loop in case there are not enough random candidate positions
    ran_steplength <- matrix(data = sample(ran_steplength_all[,,dy]),
                             nrow = ncandidates, ncol = n_agents)
    # get the matrix n_agents*ncandidates of random turning angles, and shuffle them randomly
    ran_angles <- matrix(data = sample(ran_angles_all[,,dy]),
                         nrow = ncandidates, ncol = n_agents)

    ### create ncandidates positions ----
    # https://math.stackexchange.com/questions/3877337/how-can-i-get-a-new-position-using-a-turning-angle-the-previous-point-and-an-a
    # x dimension
    old_x <- t(as.matrix(positions[day == dy - 1, .(x)]))
    candidate_x <- rep(old_x, each = nrow(ran_steplength)) + ran_steplength * cos(ran_angles)
    # y dimension
    old_y <- t(as.matrix(positions[day == dy - 1, .(y)]))
    candidate_y <- rep(old_y, each = nrow(ran_steplength)) + ran_steplength * sin(ran_angles)

    # put candidate values into data.table format
    ## x dimension
    candidates_x <- as.data.table(candidate_x)
    candidates_x[, candidate := 1:ncandidates]
    candidates_x <- candidates_x %>%
      pivot_longer(cols = !candidate, names_to = 'agent', values_to = 'x1') %>%
      mutate(agent = as.numeric(gsub('[V]', '', agent))) %>%
      as.data.table()
    ## y dimension
    candidates_y <- as.data.table(candidate_y)
    candidates_y[, candidate := 1:ncandidates]
    candidates_y <- candidates_y %>%
      pivot_longer(cols = !candidate, names_to = 'agent', values_to = 'y1') %>%
      mutate(agent = as.numeric(gsub('[V]', '', agent))) %>%
      as.data.table()
    ## merge
    candidates <- merge(candidates_x, candidates_y, by = c("agent", "candidate"))

    ### ensure candidate positions are in the GoM ----
    # set random candidates as sf elements
    point.sf <- st_as_sf(candidates[,.(x1,y1)], coords = c("x1","y1"))
    # set the coordinate reference system of the random candidates identical to the crs of the ns polygon
    st_crs(point.sf) <- ns_crs

    # filter (keep) the random candidates that are inside the North Sea polygon, and extract their coordinates
    valid_candidates <- st_filter(point.sf, this.density.sf) %>%   # this line filters
      st_coordinates() %>%                               # this line extracts the coordinates
      as.data.table() %>%
      rename(x1 = X,
             y1 = Y)
    # inner join to discard the invalid random candidates
    candidates <- inner_join(valid_candidates, candidates, by = c("x1", "y1"))
    # end step 2 -----

    # Step 3: Compute acceptance probability for these candidates ----
    # calculate weight for each member of the forecast ensemble

    ## whale density at candidate/trail position ----
    # interpolation using function from akima package
    candidates[, z1 := akima::interpp(
      x = as.vector(this.density[,'x']),
      y = as.vector(this.density[,'y']),
      z = as.vector(this.density[,'Layer_1']),
      xo = pull(candidates[,.(x1)]),
      yo = pull(candidates[,.(y1)])
    )[['z']]]

    # discard non-valid values
    candidates <- candidates[z1 > 0]

  } # end while loop

  # whale density at previous step
  whale_denom <- positions[day == dy - 1,.(agent, x, y, z)]

  # calculate the ratio of densities in the candidate locations
  # relative to the locations in the previous step
  candidates <- merge(candidates, whale_denom, by = c("agent"))
  candidates[, likeratio := z1/z]

  # get probabilities for each candidate position for each agent

  candidates <-  candidates %>%
    group_by(agent) %>%
    summarise(sumlik = sum(likeratio), .groups = 'drop') %>%
    right_join(candidates, by = "agent") %>%
    mutate(prob = likeratio/sumlik) %>%
    as.data.table()

  # sample the candidate position relative to the probability
  accepted_candidates <- candidates[, .(acc_candidate = sample(x = candidate, size = 1, replace = TRUE, prob = prob)),
                                    keyby = .(agent)]

  # subset the candidate data.table by keeping only the accepted candidates
  candidates <- candidates[accepted_candidates, on = .(agent = agent, candidate = acc_candidate)]
  candidates <- candidates %>%
    rename(acceptedx = x1,
           acceptedy = y1,
           acceptedz = z1)


  # end step 3 ----



  ## accumulate the results ----

  if(dy < 1){
    positions[((nburnin - abs(dy) - 1) * n_agents + 1):((nburnin - abs(dy)) * n_agents), day := dy]
  }
  if(dy > 0){
    positions[((dy - 1) * n_agents + 1 + n_agents * nburnin):(dy * n_agents +  n_agents * nburnin), day := dy]
  }

  # positions[((dy-1) * n_agents + 1):(dy*n_agents), day := dy]
  positions[day == dy, agent := 1:n_agents]
  positions[day == dy, x := candidates$acceptedx]
  positions[day == dy, y := candidates$acceptedy]
  positions[day == dy, z := candidates$acceptedz]

}
s3 <- Sys.time()
# end loops over days, months, -------

# arrange output ----
positions$agent <- as.factor(positions$agent)

# add month and year to the positions df ------
lkp_daymonth <-
  data.table(
    year = 2050,
    month = as.factor(c(
      rep(3, 31),
      rep(4, 30),
      rep(5, 31),
      rep(6, 30),
      rep(7, 31),
      rep(8, 31),
      rep(9, 30),
      rep(10, 31),
      rep(11, 30)
    )),
    day = c(
      60:90,
      91:120,
      121:151,
      152:181,
      182:212,
      213:243,
      244:273,
      274:304,
      305:334
    ))

positions <- merge(positions, lkp_daymonth, by = "day")

positions$size <- ifelse(positions$day %in% range(positions$day), 6, 3)

# convert positions to sf object
positions_sf <-
  st_as_sf(positions, coords = c("x","y"))
# set crs of positions
st_crs(positions_sf) <- ns_crs

# convert positions to lines
positions_trimmed <- positions

lines_sf <- sfheaders::sf_linestring( positions %>%  arrange(agent, day, x,y), x = "x", y = "y",  linestring_id = "agent", keep = TRUE )
lines_trimmed_sf <- sfheaders::sf_linestring( positions_trimmed %>%  arrange(agent, day, x,y), x = "x", y = "y",  linestring_id = "agent", keep = TRUE )
st_crs(lines_sf) <- ns_crs
st_crs(lines_trimmed_sf) <- ns_crs

lines_trimmed_sf$agent <- as.factor(lines_trimmed_sf$agent)

# plot ----
## subset agents, year round -----
ags <- c(1, 64,      52, 90,31, 22, 72)
ags <- c(27, 31, 54, 80, 38, 89, 9, 76, 10)
ags <- 1:5
p.subset.yearround <-
  mapview(lines_trimmed_sf[lines_trimmed_sf$agent %in% ags, ],
          zcol = "agent",
          legend = FALSE,
          alpha = 0.2,
          color = hcl.colors(n = length(ags), palette = "inferno")
          # legend.opacity = 0,
  ) +
  mapview( positions_sf[positions_sf$agent %in% ags,],
           cex = "size",
           alpha.regions = 0.3,
           alpha = 0.3,
           legend = FALSE,
           zcol = "month",
           burst = F)



s2 <- Sys.time()


# source plot animation -----
# source('analysis/2_animate_plot.r')

# calculate ditances travelled ------
# Transformed data
pj <- proj4::project(xy = positions[,.(x, y)], proj =  ns_crs, inverse = TRUE)
latlon <- data.table( lon=pj$x, lat=pj$y, positions)
distances <-
  latlon %>%
  arrange(agent, day) %>%
  group_by(agent) %>%
  mutate(distance = distGeo(
    p1 = cbind(lon, lat),
    p2 = cbind(lag(lon), lag(lat))) /1000
  ) %>%
  as.data.table()

p.dists.travelled <- ggplot(data = distances %>%
                              na.omit()) +
  # geom_density(data = femdists %>%
  #                filter(type == "all"), aes(x = disps/1000), color = 'black', linewidth = 2) +
  geom_density(aes(x = distance, color = as.factor(agent)), linetype = 1, alpha = 0.2) +
  xlab("Daily displacement (km)") +
  theme_classic() +
  theme(legend.position = "none")

# output ------
# positions_trimmed[, size := NULL]
# write.csv(x = positions_trimmed,
#           file = "output/MonthlyDensities_simulation_100agents_100candidates.csv",
#           row.names = FALSE)

save.image(file = 'rdata/MonthlyDensities_simulation_800agents_100candidates.rdata')
# compile report
# compile_rmd("SpermWhale_movement")

