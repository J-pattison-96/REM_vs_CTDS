## Code for all species_area combinations

#library
library(dplyr)
library(lubridate)
library(sbd)
library(activity)
library(Distance)
library(knitr)
library(tidyverse)


### load data ###
# Observation data
setwd("/Users/jessepattison/Desktop/Ecology-Masters/Masters_Thesis/Raw_data/File_results/combined_tables/combined/Final_tables")
## Load observation table
obs.table <- read.csv("combined_obs_fixed_angle.csv")
# subset to species-area that you are analyzing
obs.table <- subset(obs.table, species_area == "fallow_deer_skane")

# Load deployment table
dep.table <- read.csv("fixed_times_combined_dep.csv")
# subset to species-area that you are analyzing
dep.table <- subset(dep.table, species_area == "fallow_deer_skane")





## Parameters to estimate 
## y = number of observations
## t = time cameras were active
## v = day range
## r = radius (distance from camera) of detection zone
## tau = angle of detection zone

# estimate y

# make copy of original observation table
obs.table2<-obs.table

## create new file with only unique ids
obs.table3 <- obs.table2[!duplicated(obs.table2$uniqueID),] ## selects only the unique ids

## create right number of individuals
obs.table2$ids <- as.character(obs.table2$uniqueID)
an.list <- split(obs.table2, obs.table2$uniqueID, drop=T)
seq.count.df <- data.frame(uniqueID = NA, seq_count = NA)
for(i in 1:length(an.list)){
  mediaID <- an.list[[i]][1,"mediaID"]
  count <- sum(an.list[[i]][an.list[[i]]$mediaID==mediaID,"count"])
  count.new <- sum(an.list[[i]][2:length(an.list),"countNew"],na.rm=T)
  seq.count.df[i,1] <- an.list[[i]][1,"ids"]
  seq.count.df[i,2] <- count + count.new
  seq.count.df[i,"locationID"]<-an.list[[i]][1,"locationID"]
}

seq.count.df$uniqueID <- as.factor(seq.count.df$uniqueID)
obs.table3$seq_count <- seq.count.df[match(obs.table3$uniqueID,seq.count.df$uniqueID),"seq_count"] ## adds the correct number to the animals2 object

sum(obs.table3$seq_count)

# sum the total number of individuals
y<-sum(obs.table3$seq_count)
y


#### Estimate t ####
# Calculate total sampling effort in days

dep.table$effort_days<- as.numeric(difftime(
  dep.table$deploymentEnd, 
  dep.table$deploymentStart, 
  units = "days"
))

# remove number of days cameras were covered by snow
dep.table$effort_days <- dep.table$effort_days - dep.table$snow.cover.days

t<-sum(as.numeric(dep.table$effort_days))
t

## encounter rate
encounter.rate<-y/t
encounter.rate


#check - the total number of sequences
sum(dep.table$seq_count)


# use Trigger image to estimate effective detection distance & effective detection angle

## Get only trigger image (images are already in chronological order)
obs.table <- obs.table[order(obs.table$seq, obs.table$eventStart), ]

trigger.image <- obs.table[!duplicated(obs.table$seq), ]


## Estimate r & tau for trigger image

## Estimate r & tau for trigger image
## create needed columns from classifications
classinfo <- do.call(rbind,strsplit(trigger.image$observationTags,"\\|"))
colnames(classinfo) <- c("eng","Distance","HorProp")
classinfo <- as.data.frame(classinfo)
classinfo$Distance <- as.numeric(gsub("\\D","",classinfo$Distance))
classinfo$HorProp2 <- as.numeric(do.call(rbind,strsplit(classinfo$HorProp,":"))[,2])
classinfo$Angle <- classinfo$HorProp2 * 29

trigger.image$individualPositionRadius <- classinfo$Distance
trigger.image$individualPositionAngle <- classinfo$Angle
trigger.image$individualPositionRadians <- (trigger.image$individualPositionAngle/180)*pi

breakpoints<-c(seq(0,18,2))

effective_detect_radius <- subset(trigger.image,individualPositionRadius >=0 & individualPositionRadius<=18) # selecting data to estimate EDD
hist(effective_detect_radius$individualPositionRadius,breaks = breakpoints,ylim = c(0,200),col = "blue")

trunc.list<-list(left=0,right=18)
mybreaks<-c(seq(0,18,2))

hist(effective_detect_radius$individualPositionRadius,breaks = mybreaks, main = "trigger image distance",
     xlab = "Radial distance (m)", ylim = c(0,100),xlim = c(0,18), col = "blue")

table(effective_detect_radius$individualPositionRadius)

source("REM_functions.R") # importing some key functions to run the analysis, from vignette


# to avoid overfitting i will just try hn0, hn1, hn2, hr0, hr1
## half-normal functions
hn0<-ds(effective_detect_radius$individualPositionRadius, transect = "point", key="hn", adjustment = NULL,cutpoints = mybreaks,truncation = trunc.list)


hn1<-ds(effective_detect_radius$individualPositionRadius, transect = "point", key="hn", adjustment = "cos",nadj = 1,cutpoints = mybreaks,truncation = trunc.list)


hn2<-ds(effective_detect_radius$individualPositionRadius, transect = "point", key="hn", adjustment = "cos",nadj = 2,cutpoints = mybreaks,truncation = trunc.list)


# hazard-rate

hr0<-ds(effective_detect_radius$individualPositionRadius, transect = "point", key="hr", adjustment = NULL, cutpoints = mybreaks,truncation = trunc.list)


hr1 <-ds(effective_detect_radius$individualPositionRadius, transect = "point", key="hr", adjustment = "poly",nadj = 1,cutpoints = mybreaks,truncation = trunc.list)

## model comparison
library(kableExtra)
AIC_Table <- AIC(hn0, hn1, hn2, hr0, hr1)
AIC_Table[order(AIC_Table$AIC), ]


chats <- chi2_select (hn1, hr1)$criteria
modnames <- unlist(lapply(list(hn1, hr1), function(x) x$ddf$name.message))
results <- data.frame(modnames, chats)
results.sort <- results[order(results$chats),]
knitr::kable(results.sort, digits=2, row.names = FALSE,
             caption="Compare with Table S5 of Howe et al. (2018)") %>%
  kable_paper(full_width = FALSE) %>%
  row_spec(1, bold=TRUE,  background = "#4da6ff")

best_modRad<-hn1 ## select best most parsimonious model ### CHANGE EVERYTIME



r_key <- "hn"
r_adjustment <- "cos"
r_nadj <- 1

# For hn
# r_key <- "hn"
# r_adjustment <- NULL
# r_nadj <- 0

# For hn1
# r_key <- "hn"
# r_adjustment <- "cos"
# r_nadj <- 1

# For hn2
#r_key <- "hn"
#r_adjustment <- "cos"
#r_nadj <- 2

# For hr0
#r_key <- "hr"
#r_adjustment <- NULL
#r_nadj <- 0

# For hr1
#r_key <- "hr"
#r_adjustment <- "poly"
#r_nadj <- 1

# Estimating effective detection radius and (SE)
EfecRad <- EDRtransform(best_modRad)
EfecRad$EDR

# effective detectio radius standard error
EfecRad$se.EDR<-EfecRad$se.EDR # SE (m)
EfecRad$se.EDR[1]


par(mfrow=c(1,2))

### plot detection curve # CHANGE KEY EVERYTIME
library(ggplot2)
plot(hn1,main="detection probability \nby distance of trigger image", xlab= "Distance (m)",ylim=c(0,3.0),showpoints=FALSE,lwd=3,breaks=mybreaks)

plot(hn1, main="Probability density \nof Daytime activity", xlab="Distance (m)", pdf=TRUE,
     showpoints=FALSE, lwd=3, breaks=mybreaks, ylim= c(0,0.25))

r <- EfecRad$EDR
r
# put into correct units 
r<-r/1000
r

# standard error and unit conversion
EfecRad$se.EDR[1]
r.se<-EfecRad$se.EDR[1]/1000
r.se


# Effective Detection Angle

effective_detect_radians <- abs(effective_detect_radius$individualPositionRadians)  # selecting data to estimate EDA
FOV<-58 #field of view of the cameara (degrees)
w_ang<-(FOV/2)*0.0174533 # truncation angle (in radians)

## fitting models 

# half-normal
hn_Ang <- ds(effective_detect_radius$individualPositionRadians, transect = "line", key="hn", adjustment = NULL, truncation=
               w_ang)
## fitting half-normal key functions
## AIC = -347.998


uni_Ang<-ds(effective_detect_radius$individualPositionRadians, transect = "line", key="unif", adjustment = NULL, truncation=
              w_ang)
# Fitting uniform key function


AIC(hn_Ang,uni_Ang)


# model selection 
best_modAng <- uni_Ang  # select best most parsmonious model


tau_key <- "unif"
tau_adjustment <- NULL

#tau_key <- "unif"
#tau_adjustment <- NULL

#tau_key <- "hn"
#tau_adjustment <- NULL


# Estimating effective detection radius and (SE)
summary_ang<- summary(best_modAng$ddf)

EfecAng_mean <- summary_ang$average.p*w_ang # mean (radians)
EfecAng_mean
### 0.5061457
EfecAng_SE <- summary_ang$average.p.se*w_ang # SE (radians)
EfecAng_SE<-EfecAng_SE[1]*2

EfecAng_SE
tau.SE<-EfecAng_SE
tau.SE

plot(uni_Ang,main="detection probability \nby distance of trigger image", xlab= "Distance (m)",ylim=c(0,3),showpoints=FALSE,lwd=3)

tau <- EfecAng_mean
tau<-tau*2 # multiple by 2 because the angle taken from 2 sides.
tau


#### Estimate v  (day range)####
## can be done for all images as properly classified from the start

## formula for calculating distance moved (c) based on radius at t-1 (a), radius at t (b) 
## and angle between them (gamma) calculated as angle at t-1 minus angle at t
## c^2 = a^2 + b^2 - 2ab*cos(gamma)
a = obs.table2[1,"individualPositionRadius"] # Distance from camera in first image
b = obs.table2[2,"individualPositionRadius"] # Distance from camera in second image
gamma = obs.table2[2,"individualPositionRadians"]-obs.table2[1,"individualPositionRadians"] # difference in angle
c = sqrt(a^2 + b^2 - (2*a*b*cos(gamma))) # calculation for law of cosines
# c will give the hypotenuse, which is the distance travelled between images.

## get distance walked and time per sequence
obs.table2$uniqueID2 <- paste0(obs.table2$seq,obs.table2$scientificName,obs.table2$lifeStage,obs.table2$sex,obs.table2$individualID)

## generate distance column
seq.list <- split(obs.table2,obs.table2$uniqueID2,drop=T)
n <- length(seq.list)


for(i in 1:n){
  seq.list[[i]][1,"distanceMoved"] <- 0
  if(nrow(seq.list[[i]]) > 1)
    for(j in 2:nrow(seq.list[[i]])){
      seq.list[[i]][j,"distanceMoved"] <- sqrt(seq.list[[i]][j,"individualPositionRadius"]^2 
                                               + seq.list[[i]][j-1,"individualPositionRadius"]^2 
                                               - 2*seq.list[[i]][j,"individualPositionRadius"]*
                                                 seq.list[[i]][j-1,"individualPositionRadius"]*
                                                 cos(seq.list[[i]][j,"individualPositionRadians"]-
                                                       seq.list[[i]][j-1,"individualPositionRadians"]))
    }
  
}

obs.table6 <- do.call(rbind,seq.list)
rownames(obs.table6) <- c()

#View(obs.table6)
head(obs.table6)

## generate sum distance walked and time per sequence
## create new file with only unique ids
obs.table7 <- obs.table6[!duplicated(obs.table6$uniqueID2),] ## selects only the unique ids


## create right number of individuals
obs.table6$ids2 <- as.character(obs.table6$uniqueID2)
speed.list <- split(obs.table6, obs.table6$uniqueID2, drop=T)
speed.df <- data.frame(uniqueID2 = NA, distanceMoved = NA, timeInFrame = NA)
for(i in 1:length(speed.list)){
  distance.moved <- sum(speed.list[[i]]$distanceMoved,na.rm=T)
  time.in.frame <- as.numeric(difftime(speed.list[[i]][nrow(speed.list[[i]]),"eventStart"],
                                       speed.list[[i]][1,"eventStart"],
                                       units="secs"))
  speed.df[i,1] <- speed.list[[i]][1,"ids2"]
  speed.df[i,2] <- distance.moved
  speed.df[i,3] <- time.in.frame
}

speed.df$uniqueID2 <- as.factor(speed.df$uniqueID2)
speed.df$speed <- speed.df$distanceMoved/speed.df$timeInFrame

#obs.nor7 is a data.frame which has
#a single row per unique inidvidual with its speed in meters per second
obs.table7$speed <- speed.df[match(obs.table7$uniqueID2,speed.df$uniqueID2),"speed"]

summary(obs.table7) ## checking the results


## analysis based on mean speed
# subset to speeds greater than 0.01 and use harmonic mean as there is a high variance
library(sbd)
m.s.speed<- hmean(obs.table7[obs.table7$scientificName=="Dama dama" & obs.table7$speed>0.01&!is.na(obs.table7$speed),"speed"])
mean.m.s.speed<-m.s.speed$mean # need for dayrange function

mean.m.s.speed


se.m.s.speed<-m.s.speed$se # need for dayrange function
se.m.s.speed


# number of speed obs, need for dayrange function
speed.length<- (obs.table7[obs.table7$scientificName=="Dama dama"& obs.table7$speed>0.01 & !is.na(obs.table7$speed),"speed"])
num_observations<-length(speed.length)
num_observations


speed<- hmean(obs.table7[obs.table7$scientificName=="Dama dama" & obs.table7$speed>0.01&!is.na(obs.table7$speed),"speed"])$mean*3.6
speed

m.s.speed


speed.data.frame <- data.frame(
  mean_m_s_speed_ = mean.m.s.speed,
  se_m_s_speed_ = se.m.s.speed,
  num_observations = num_observations)



# time in frame speed
library(dplyr)

time.in.frame <- speed.df %>%
  filter(timeInFrame != 0, speed >= 0.01) %>%
  select(timeInFrame)

 
## Estimating Activity Patterns

## estimate activity level ## looks like Tim uses all the images where the distance moved is greater than 0
library(activity)
obs.table7$radianTime <- gettime(obs.table7$eventStart)


# need to update scientific name
t.obs<-obs.table7[obs.table7$scientificName=="Dama dama","radianTime"] 
t.obs
## anchor using solartime
act.pat<-obs.table7[obs.table7$scientificName=="Dama dama",] 
act.pat
# solar time
library(lubridate)
# first need to convert time to proper format
# Your date-time data
eventStart<- act.pat$eventStart

eventStart.converted <- format(ymd_hms(eventStart, tz = "Europe/Stockholm"), "%Y/%m/%d %H:%M")

act.pat$eventStart.converted<-eventStart.converted
act.pat
times<-solartime(act.pat$eventStart.converted,act.pat$latitude,act.pat$longitude,0) 
times ## it looks like the output is in radians
rawAct <- fitact(times$clock)
rawAct # fits activity model using unanchored times
avgAct <- fitact(times$solar)
avgAct # fits activity model using anchored times



par(mfrow=c(1,1))
plot(rawAct, ylim=c(0,100) , main="anchored activity patterns")
plot(avgAct, add=TRUE, data="n", tline=list(col="cyan"))


### activity patterns
actmod<-fitact(times$solar,sample = "data",reps = 100)
actmod@act

act<-actmod@act[1]
act

act_se<-actmod@act[2]
act_se

# manually calculate day range
v <- speed*as.numeric(act)*24 # put velocity into day range
v

# check against day range function
library(trappingmotion)
dayrange(act,act_se,speed.data.frame)

#### REM estimate ####
#D = (y/t)*(pi/(v*r*(2+tau)))

D = (y/t)*(pi/(v*r*(2+tau)))
D


# Bootstrapping
## bootstrap --- Jesse Code -----
# June 14, 2026 ### WORKING CODE FOR R AND TAU

### Reproducible bootstrap code
### table moose

library(Distance)
library(parallel)
library(sbd)

RNGkind("L'Ecuyer-CMRG")
set.seed(1996)

n_cores <- 1
n_iter <- 500

FOV <- 58
w_ang <- (FOV / 2) * 0.0174533

trunc.list <- list(left = 0, right = 18)
mybreaks <- c(seq(0, 18, 2))




sample_size_yt <- nrow(dep.table)
sample_size_speed <- num_observations

## -----------------------------
## Pre-generate bootstrap indices
## -----------------------------

idx_r_tau <- replicate(
  n_iter,
  sample.int(nrow(effective_detect_radius),
             size = nrow(effective_detect_radius),
             replace = TRUE),
  simplify = FALSE
)

idx_yt <- replicate(
  n_iter,
  sample.int(nrow(dep.table),
             size = sample_size_yt,
             replace = TRUE),
  simplify = FALSE
)

valid_speeds <- obs.table7$speed[
  obs.table7$scientificName == "Dama dama" &
    !is.na(obs.table7$speed) &
    is.finite(obs.table7$speed) &
    obs.table7$speed > 0.01
]

idx_speed <- replicate(
  n_iter,
  sample.int(length(valid_speeds),
             size = sample_size_speed,
             replace = TRUE),
  simplify = FALSE
)

## -----------------------------
## Bootstrap r and tau
## -----------------------------

# add percentage bar
library(pbmcapply)

boot_list <- pbmclapply(seq_len(n_iter), function(i) {
  
  sampled_data <- effective_detect_radius[idx_r_tau[[i]], ]
  
  r_sample <- sampled_data$individualPositionRadius
  tau_sample <- abs(sampled_data$individualPositionRadians)
  
  
  r.model <- try(
    ds(
      r_sample,
      transect = "point",
      key = r_key,
      adjustment = r_adjustment,
      nadj = r_nadj,
      cutpoints = mybreaks,
      truncation = trunc.list
    ),
    silent = TRUE
  )
  

  tau.model <- try(
    ds(
      tau_sample,
      transect = "line",
      key = tau_key,
      adjustment = tau_adjustment,
      truncation = w_ang
    ),
    silent = TRUE
  )
  
  edr <- if (!inherits(r.model, "try-error") && !is.null(r.model)) {
    EDRtransform(r.model)$EDR
  } else {
    NA_real_
  }
  
  eda <- if (!inherits(tau.model, "try-error") && !is.null(tau.model)) {
    
    summary_ang_boot <- summary(tau.model)
    EfecAng_mean_boot <- summary_ang_boot$ds$average.p * w_ang
    EfecAng_mean_boot
    
  } else {
    NA_real_
  }
  
  data.frame(r = edr, tau = eda)
  
}, mc.cores = n_cores, mc.set.seed = TRUE, mc.preschedule = TRUE)

boot_results <- do.call(rbind, boot_list)

## Convert r from metres to kilometres
boot_results$r <- boot_results$r / 1000

## -----------------------------
## Bootstrap y and t
## -----------------------------

t_sums <- numeric(n_iter)
y_sums <- numeric(n_iter)

for (i in seq_len(n_iter)) {
  
  sampled_rows <- dep.table[idx_yt[[i]], ]
  
  t_sums[i] <- sum(sampled_rows$effort_days, na.rm = TRUE)
  y_sums[i] <- sum(sampled_rows$ind_count, na.rm = TRUE)
}

boot_results$t <- t_sums
boot_results$y <- y_sums

## -----------------------------
## Bootstrap speed
## -----------------------------

boot_hmean <- numeric(n_iter)

for (i in seq_len(n_iter)) {
  
  sample_speeds <- valid_speeds[idx_speed[[i]]]
  
  boot_hmean[i] <- hmean(sample_speeds)$mean
}

boot_results$speed <- boot_hmean

## Convert speed to day range using activity
boot_results$speed <- boot_results$speed * 3.6 * 24 * act

## -----------------------------
## Density calculation
## -----------------------------

boot_results$d <- (
  (boot_results$y / boot_results$t) *
    (pi / (
      boot_results$r *
        boot_results$speed *
        (2 + (boot_results$tau*2))
    ))
)

## -----------------------------
## Summary statistics
## -----------------------------

d.median <- median(boot_results$d, na.rm = TRUE)
d.median
d.lower.C.I. <- as.numeric(quantile(boot_results$d, probs = 0.025, na.rm = TRUE))
d.lower.C.I.
d.upper.C.I. <- as.numeric(quantile(boot_results$d, probs = 0.975, na.rm = TRUE))
d.upper.C.I.

y
y.median <- median(boot_results$y, na.rm = TRUE)
y.median
y.low.ci <- as.numeric(quantile(boot_results$y, probs = 0.025, na.rm = TRUE))
y.low.ci
y.high.ci <- as.numeric(quantile(boot_results$y, probs = 0.975, na.rm = TRUE))
y.high.ci

t
t.median <- median(boot_results$t, na.rm = TRUE)
t.median
t.low.ci <- as.numeric(quantile(boot_results$t, probs = 0.025, na.rm = TRUE))
t.low.ci
t.high.ci <- as.numeric(quantile(boot_results$t, probs = 0.975, na.rm = TRUE))
t.high.ci


# fix errors in estimated r


r
r.median <- median(boot_results$r, na.rm = TRUE)
r.median
r.low.ci <- as.numeric(quantile(boot_results$r, probs = 0.025, na.rm = TRUE))
r.low.ci
r.high.ci <- as.numeric(quantile(boot_results$r, probs = 0.975, na.rm = TRUE))
r.high.ci

tau
tau.median <- median(boot_results$tau, na.rm = TRUE)
tau.median*2
tau.low.ci <- as.numeric(quantile(boot_results$tau, probs = 0.025, na.rm = TRUE))
tau.low.ci*2
tau.high.ci <- as.numeric(quantile(boot_results$tau, probs = 0.975, na.rm = TRUE))
tau.high.ci*2

v
v.median <- median(boot_results$speed, na.rm = TRUE)
v.median
v.low.ci <- as.numeric(quantile(boot_results$speed, probs = 0.025, na.rm = TRUE))
v.low.ci
v.high.ci <- as.numeric(quantile(boot_results$speed, probs = 0.975, na.rm = TRUE))
v.high.ci

summary_table <- data.frame(
  parameter = c("density", "y", "t", "r", "tau", "speed"),
  median = c(d.median, y.median, t.median, r.median, tau.median, v.median),
  lower_95_CI = c(d.lower.C.I., y.low.ci, t.low.ci, r.low.ci, tau.low.ci, v.low.ci),
  upper_95_CI = c(d.upper.C.I., y.high.ci, t.high.ci, r.high.ci, tau.high.ci, v.high.ci)
)

print(summary_table)



SE_D <- sd(boot_results$d, na.rm = TRUE)
SE_D

## -----------------------------
## Plots
## -----------------------------

hist(boot_results$r,
     main = "Bootstrap EDR Estimates",
     xlab = "Radial distance (km)",
     col = "red",
     ylim = c(0, 150))

hist(boot_results$tau,
     main = "Bootstrap EDA Estimates",
     xlab = "Tau",
     col = "lightblue")

hist(boot_results$speed,
     breaks = 30,
     col = "lightblue",
     main = "Bootstrapped Harmonic Mean Speeds",
     xlab = "Speed / day range")

hist(boot_results$d,
     breaks = 10,
     col = "lightblue",
     main = "table Bootstrapped Density Estimates",
     xlab = "Density",
     border = "white",
     ylim = c(0, 120))

abline(v = d.upper.C.I., col = "red", lwd = 2, lty = 2)
abline(v = d.lower.C.I., col = "red", lwd = 2, lty = 2)
abline(v = d.median, col = "blue", lwd = 2, lty = 2)
abline(v = D, col = "black", lwd = 2, lty = 2)


sessionInfo()

summary(boot_results$d)
quantile(
  boot_results$d,
  probs = c(0.99, 0.995, 0.999, 1),
  na.rm = TRUE
)


