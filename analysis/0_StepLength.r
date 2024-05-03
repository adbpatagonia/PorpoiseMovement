# ADB
# 2024-05-03


# libraries -----
library(circular) # use for random draws of von Mises dist
library(CircStats) # von Mises Density function
library(tidyverse)
library(data.table)
library(here)

# read data -----
displacements <- fread(here('data', 'PP_displacement.csv')) %>%
  na.omit()

# parameters ----
## turning angle parameter ----
# kappa parameter of the von Mises distribution. Defines degree to which angles clump around zero
# very low values allow agents to turn anywhere from 0 to 360
kap <- 0.01

# wrangle data ----
# displacements[, sex := tolower(sex)]
# displacements[, aggsex := ifelse(sex == 'male', 'male', 'other')]

# drop one record where a male swam over 250 km in one day
# displacements <- displacements[daily.displacement < 200]

# convert displacement to meters
displacements[, daily.displacement := daily.displacement * 1000]



## obtain mean, sd, and parameters for gamma distribution ------
## combined female+juv, and combined males

# scale_par = Var(x)/E[x]
# shape_par = E[x]^2/Var(x)

dis.sum.location <-
  displacements %>%
  group_by(Location) %>%
  reframe(mu_displacement = mean(daily.displacement),
          sd_displacement = sd(daily.displacement)) %>%
  mutate(scale_gamma = sd_displacement^2 / mu_displacement,
         shape_gamma =  mu_displacement^2 / sd_displacement^2) %>%
  as.data.table()
dis.sum <-
  displacements %>%
  reframe(mu_displacement = mean(daily.displacement),
          sd_displacement = sd(daily.displacement)) %>%
  mutate(scale_gamma = sd_displacement^2 / mu_displacement,
         shape_gamma =  mu_displacement^2 / sd_displacement^2) %>%
  mutate(Location = "combined") %>%
  as.data.table()

dis.sum <- rbindlist(l = list(dis.sum.location, dis.sum), use.names = TRUE)


# create gamma distributions ------
n_random_candidates <- 100000

dists <-
  rbindlist(l = list(
    data.table(disps = rgamma(n = n_random_candidates,
                              shape = dis.sum$shape_gamma[dis.sum$Location == "combined"],
                              scale = dis.sum$scale_gamma[dis.sum$Location == "combined"]),
               Location = "combined"),

    data.table(disps = rgamma(n = n_random_candidates,
                              shape = dis.sum$shape_gamma[dis.sum$Location == "West Greenland"],
                              scale = dis.sum$scale_gamma[dis.sum$Location == "West Greenland"]),
               Location = "West Greenland"),

    data.table(disps = rgamma(n = n_random_candidates,
                              shape = dis.sum$shape_gamma[dis.sum$Location == "NE coast of USA and Canada"],
                              scale = dis.sum$scale_gamma[dis.sum$Location == "NE coast of USA and Canada"]),
               Location = "NE coast of USA and Canada")
  ))

# turning angles ------
angles <- data.table(angle = rvonmises(n = 10000, mu = 0, k = kap, control.circular=list(units="radians")))
angles[, density := dvm(angles$angle, mu = 0, kappa = kap)]


# plot -----
p.dists <- ggplot(data = dists, aes(x = disps/1000, color = Location)) +
  geom_density() +
  xlab("Daily displacement (km)") +
  theme_classic() +
  theme(legend.position = "bottom")

p.angles <- ggplot(data = angles, aes(x = angle, y = density)) +
  geom_line() +
  xlab("Turning angle (radians)") +
  # xlab(expression("Turning angle ("*pi*")")) +
  ylim(0, max(angles$density) * 1.1) +
  theme_classic()

# output -----
ggsave('output/StepLengthDistributions.png', plot = p.dists)
ggsave('output/TurningAngleDistributions.png', plot = p.angles)

