## CTDS July 6 2026


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
setwd("### Set working directory")
## Load observation table
obs.table <- read.csv("combined_obs_fixed_angle.csv")
# subset to species-area that you are analyzing
obs.table <- subset(obs.table, species_area == "fallow_deer_skane") # select species area combination
# Add region label
obs.table$Region.Label <- obs.table$species_area
# set Area to 0
obs.table$Area <- 0




# Load deployment table
dep.table <- read.csv("fixed_times_combined_dep.csv")
# subset to species-area that you are analyzing
dep.table <- subset(dep.table, species_area == "fallow_deer_skane") # select species area combination


# remove number of days cameras were covered by snow
dep.table$effort_days <- dep.table$effort_days - dep.table$snow.cover.days


# make copy of observation table

obs.table2<-obs.table

### change Vertical_Distance into distance ## required for CT_DS package
colnames(obs.table2)[colnames(obs.table2)=="individualPositionRadius"]<-"distance"
# check to see that column name has been changed
obs.table2

### only take every second image from each sequence from obs.table2
# Ensure eventStart is in POSIXct format
obs.table2$eventStart <- as.POSIXct(obs.table2$eventStart)

# Split data by sequence
split_obs <- split(obs.table2, obs.table2$seq)



obs.table2 <- obs.table2 %>%
  # Step 1: Only keep even seconds (2-second photoperiods)
  filter(second(eventStart) %% 2 == 0) %>%
  arrange(seq, eventStart) %>%
  group_by(seq, eventStart) %>%
  # Step 2: Keep only rows that share the first mediaID in each timestamp
  filter(mediaID == first(mediaID)) %>%
  ungroup()


# show more decimal places
options(digits = 15)  # Show more decimal places

# Calculate sampling effort by photoperiod

dep.table$Effort <- dep.table$effort_days*43200 # 43200 is the number of 2 sec time steps in a 24 hr period

t.seconds<-sum(as.numeric(dep.table$Effort))
t.seconds # in seconds


# remove any sequences that have over 1000 images (animal is sleeping in front of camera)
library(dplyr)

obs.table2 <- obs.table2 %>%
  group_by(seq) %>%
  filter(n() <= 1000) %>%
  ungroup()

# Adjustment for temporal availability
# using activity package
## We use the first image of each unique individual
# start by creating a unique ID

obs.table3<-obs.table2

obs.table3$uniqueID <- paste0(obs.table3$seq,obs.table3$scientificName,obs.table3$lifeStage,obs.table3$sex,obs.table3$individualID)

obs.table3 <- obs.table3[!duplicated(obs.table3$uniqueID),] ## selects only the unique ids

library(activity)
obs.table3$radianTime <- gettime(obs.table3$eventStart)



# the data frame below selects via scientific name from 

t.working<-obs.table3[obs.table3$scientificName=="Dama dama","radianTime"] 
t.working

## anchor using solartime ### MAKE SURE TO UPDATE SCIENTIFIC NAME ###
act.pat<-obs.table3[obs.table3$scientificName=="Dama dama",] 
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


# plot activity patterns
par(mfrow=c(1,1))
plot(rawAct,main= "anchored activity patterns")
plot(avgAct, add=TRUE, data="n", tline=list(col="cyan"))


### activity patterns
actmod<-fitact(times$solar,sample = "data",reps = 100)
actmod@act

act<-actmod@act[1]
act

act_se<-actmod@act[2]
act_se



avail <- list(creation=data.frame(rate = actmod@act[1],
                                  SE   = actmod@act[2]))

avail # Estimated proportion of a 24-hour period during which the species is active (activity level).


## add object field which uniquely identifies each observation
obs.table2$object[!is.na(obs.table2$distance)] <- 1:sum(!is.na(obs.table2$distance))



## add Sample.Label field = it is locationID
obs.table2$Sample.Label<-obs.table2$locationID

# Remove old Effort column if it already exists
obs.table2$Effort <- NULL

obs.table2 <- merge(
  obs.table2,
  dep.table[, c("locationID", "Effort")],
  by = "locationID",
  all.x = TRUE
)



obs.table4 <- obs.table2

#remove columns
obs.table4 <- obs.table4[, !(names(obs.table4) %in% c(
  "observationID",
  "deploymentID",
  "mediaID",
  "eventID",
  "eventStart",
  "eventEnd",
  "observationLevel",
  "cameraSetupType",
  "count",
  "countNew",
  "lifeStage",
  "sex",
  "behavior",
  "individualID",
  "bboxes",
  "classificationMethod",
  "classifiedBy",
  "classificationTimestamp",
  "classificationProbability",
  "observationTags",
  "observationComments",
  "X_id",
  "ids",
  "uniqueID",
  "individualPositionAngle",
  "individualPositionRadians",
  "uniqueID2",
  "distanceMoved",
  "ids2"
))]

obs.table4 <- obs.table4[obs.table4$distance <= 18, ]# remove distances over 19 meters


table(obs.table4$distance)


## Distance recording 
breakpoints<- c(seq(0,18,2)) ## we have not yet truncated the distance
hist(obs.table4$distance, breaks = breakpoints, main = "distance from camera (m)", xlab = "Radial Distance (m)")

# if left truncation is needed, remove images from 1 meter
obs.table4 <- obs.table4[obs.table4$distance > 1, ]
table(obs.table4$distance)



## truncation and function fits
library(Distance)
mybreaks<-c(seq(0,18,2)) ### Adjust if left truncation is needed.
trunc.list<-list(left=0, right=18) ### Adjust if left truncation is needed.
conversion <- convert_units("meter", NULL, "square kilometer")


# once you determine if left truncation is needed, then
# then you can define obs.table4 again and add object column


# set n
sum(
  obs.table4$distance > trunc.list$left &
    obs.table4$distance < trunc.list$right,
  na.rm = TRUE
) # check number of observations

# set n
n<-sum(
  obs.table4$distance > trunc.list$left &
    obs.table4$distance < trunc.list$right,
  na.rm = TRUE
)

# fitting half normal function

hn0 <- ds(obs.table4$distance, transect = "point", key = "hn", adjustment = NULL,
          convert_units = conversion, cutpoints = mybreaks, truncation = trunc.list)

hn1 <- ds(obs.table4$distance, transect = "point", key = "hn", adjustment = "cos",
          nadj = 1, convert_units = conversion,
          cutpoints = mybreaks, truncation = trunc.list)

hn2 <- ds(obs.table4$distance, transect = "point", key = "hn", adjustment = "cos",
          nadj = 2, convert_units = conversion,
          cutpoints = mybreaks, truncation = trunc.list)

# hazard-rate

hr0<-ds(obs.table4$distance, transect = "point", key="hr", adjustment = NULL,
        cutpoints = mybreaks,truncation = trunc.list)


hr1 <-ds(obs.table4$distance, transect = "point", key="hr", adjustment = "poly",
         nadj = 1,cutpoints = mybreaks,truncation = trunc.list)


## Model selection adjustments from overdispersion
AIC(hn0, hn1,hn2,hr0,hr1) 

## QAIC # use for model selection

QAIC(hn0,hn1,hn2)
QAIC(hr0,hr1)

library(kableExtra)

chats <- chi2_select(hn0,hr0)$criteria # Compare candidate models and select the most parsimonious model from each key function.
modnames <- unlist(lapply(list(hn0,hr0), function(x) x$ddf$name.message))
results <- data.frame(modnames, chats)
results.sort <- results[order(results$chats),]
knitr::kable(results.sort, digits=2, row.names = FALSE,
             caption="Compare with Table S5 of Howe et al. (2018)") %>%
  kable_paper(full_width = FALSE) %>%
  row_spec(1, bold=TRUE,  background = "#4da6ff")


p_a <- hr0$ddf$fitted[1]
p_a
w <- range(mybreaks)[2] - range(mybreaks)[1]
rho <- sqrt(p_a * w^2)
rho


## Selected detection function
par(mfrow=c(1,2))

plot(hr0, main="Daytime activity", xlab="Distance (m)",
     showpoints=FALSE, lwd=3, xlim=c(0, 18))
plot(hr0, main="Daytime activity", xlab="Distance (m)", pdf=TRUE,
     showpoints=FALSE, lwd=3, xlim=c(0, 18))

viewangle <- 29 # degrees is it the total angle or half?
samfrac <- viewangle / 360



##adding missing camera-trap locations to obs.table4
# remove any existing camera-only rows first
obs.table4 <- obs.table4[!is.na(obs.table4$observationType), ]

missing_locations <- setdiff(dep.table$locationID, obs.table4$locationID)

new_rows <- dep.table[dep.table$locationID %in% missing_locations, ]

missing_cols <- setdiff(names(obs.table4), names(new_rows))
for (col in missing_cols) {
  new_rows[[col]] <- NA
}

extra_cols <- setdiff(names(new_rows), names(obs.table4))
new_rows <- new_rows[, !(names(new_rows) %in% extra_cols)]

new_rows <- new_rows[, names(obs.table4)]

obs.table4 <- rbind(obs.table4, new_rows)

obs.table4$object[is.na(obs.table4$observationType)] <- NA
obs.table4$Area <- 0
obs.table4$Region.Label <- obs.table4$species_area

# add numbers to object column
# object IDs only for real detections
obs.table4$object <- NA_integer_
obs.table4$object[!is.na(obs.table4$distance)] <- seq_len(sum(!is.na(obs.table4$distance)))


# fill NAs in Sample.Label column
obs.table4$Sample.Label[is.na(obs.table4$Sample.Label)] <-
  obs.table4$locationID[is.na(obs.table4$Sample.Label)]


# calculate density

peak.uni.dens <- dht2(hr0, flatfile= obs.table4, strat_formula = ~Region.Label,
                      sample_fraction = samfrac, er_est = "P2", multipliers = avail,
                      convert_units = conversion)
print(peak.uni.dens, report="density")



 
## Bootstrap
set.seed(1996)

mysummary <- function(ests, fit){
  return(data.frame(Label = ests$individuals$D$Label,
                    Dhat = ests$individuals$D$Estimate))
}


mult <- list(availability= make_activity_fn(act.pat$radianTime, sample="data",
                                            detector_daily_duration=24))

n.cores <- 1

# MAKE SURE TO UPDATE MODEL IF NEEDED
daytime.boot.uni <- bootdht(model=hr0, flatfile=obs.table4,
                            resample_transects = TRUE, nboot=500, 
                            cores = n.cores,
                            summary_fun=mysummary, sample_fraction = samfrac,
                            convert_units = conversion, multipliers=avail)

print(summary(daytime.boot.uni))

par(mfrow=c(1,1))

hist(daytime.boot.uni$Dhat, breaks = 20, 
     xlab="Estimated density", main="D-hat estimates bootstraps")
abline(v=quantile(daytime.boot.uni$Dhat, probs = c(0.025,0.975), na.rm=TRUE), lwd=2, lty=3)



t<-2
n <-671
tau<-1.01229
w<-0.018
t.seconds
act
T.total<-(t.seconds*act)
T.total
P<- p_a

D<-(2*t*n)/(tau*(w^2)*T.total*P)
D


