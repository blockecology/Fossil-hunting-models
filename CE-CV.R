setwd("~/Dropbox/Block et al 2015 Fossil Hunting") # Change this path to match the location of the directory in your computer

# Load dplyr package to process the data 
library(dplyr)

# Loading and filtering data
## Load fossil database
fossils <- read.csv("Data/fossils.csv")
grid_id <- read.table("Data/Sahul_gridcellID.txt", quote="\"")
fossils$cellid <- paste(grid_id$V4, grid_id$V5, sep=",")
rm(grid_id)

## Correct wrong coordinates for Millicent, Ningaloo and TEC fossils
fossils[fossils$Cave.Site == "Millicent", 5] <- 140.4

fossils[fossils$Cave.Site == "Ningaloo", 4] <- -22.7 
fossils[fossils$Cave.Site == "Ningaloo", 5] <- 113.674

fossils[fossils$Cave.Site == "Tight Entrance Cave", 4] <- -34.067 
fossils[fossils$Cave.Site == "Tight Entrance Cave", 5] <- 115.167

## Select variables and rows that we will use
fossils <- select(fossils, Latitude:Class, Status, Megafauna, Age.calib., Sd.calib., LR, Context, cellid)

## Assign correct classes to variables
fossils$Latitude <- as.numeric(as.character(fossils$Latitude))
fossils$Longitude <- as.numeric(as.character(fossils$Longitude))
fossils$Temp <- as.numeric(as.character(fossils$Temp))
fossils$Rain <- as.numeric(as.character(fossils$Rain))
fossils$cellid <- as.factor(fossils$cellid)
fossils$Sd.calib. <- as.numeric(as.character(fossils$Sd.calib.))

## Assign 0 to all fossils without values for dating uncertainity (Sd.calib)
fossils[is.na(fossils$Sd.calib.), 13] <- 0

## Divide Age.calib. and Sd.calib. between 1000
fossils$Age.calib. <- fossils$Age.calib. / 1000
fossils$Sd.calib. <- fossils$Sd.calib. / 1000

## Filter missing and erroneous data
fossils <- filter(fossils, !is.na(Latitude))
fossils <- filter(fossils, Longitude != 0)

## Filter only Australian fossils
fossils <- filter(fossils, Latitude < -5)

## Drop Homo rows
fossils <- filter(fossils, Species != "Homo")

## Filter bad quality fossils
fossB <- filter(fossils, LR != "A" & LR != "A (*)" & LR != "A(*)" & LR != "A*/A" & LR != "A*")

## Drop fossils with bad-quality dates, keep just A and A* (this are the ones we can use for the SDM)
fossA <- filter(fossils, LR == "A" | LR == "A (*)" | LR == "A(*)" | LR == "A*/A" | LR == "A*")

## Add a column with the age for which we have climate data that is closest to the fossil age
## Finding the right time
climate_times <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,
                   26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,
                   70,72,74,76,78,80,84,88,92,96,100,104,108,112,116,120)


library(FastImputation)
fossA$Age2 <- LimitToSet(fossA$Age.calib., climate_times)

# To account for date uncertainity, we also create columns with age plus/minus 1 and 2 standard deviations
fossA$Age2plsd <- LimitToSet(fossA$Age.calib. + fossA$Sd.calib., climate_times)
fossA$Age2mnsd <- LimitToSet(fossA$Age.calib. - fossA$Sd.calib., climate_times)
fossA$Age2pl2sd <- LimitToSet(fossA$Age.calib. + 2 * fossA$Sd.calib., climate_times)
fossA$Age2mn2sd <- LimitToSet(fossA$Age.calib. - 2 * fossA$Sd.calib., climate_times)


## Context - Keep only fossils in a context that makes their dates reliable (exclude above and below)
fossA <- filter(fossA, Context != "A" & Context != "B")

## Remove very old fossils for which we don't have paleoclimate data
fossA <- filter(fossA, !is.na(Temp))

## Change order of Latitude and Longitude
fossA <- fossA[,c(2,1,3:21)]



# Load packages 
require(raster)
require(dismo)



genus <- "Thylacoleo"  # For example, Thylacoleo
paleosp <- subset(fossA, gen.indetus == genus) 

# Create a data frame to store the coordinates, climate data and age of each fossil
presences <- data.frame(Longitude = numeric(), Latitude = numeric(),
                        Temp = numeric(), Rain = numeric(), Age2 = factor())

# This loop extracts climate data for each fossil
for(i in 1:nrow(paleosp)){
        
        if(paleosp$Sd.calib.[i] != 0) {
                prespts <- paleosp[i,1:2]
                
                # Mean age
                age <- levels(as.factor(paleosp$Age2[i]))
                if(as.numeric(age) < 10) age <- paste("00", age, sep="")
                if(as.numeric(age) > 9 & as.numeric(age) < 100) age <- paste("0", age, sep="")
                pathclim <- paste("Data/Paleoclimate/Climate/", age, ".grd", sep="")
                paleoclim <- brick(pathclim)
                mean.presvals <- extract(paleoclim, prespts)
                
                # Age plus 1 sd
                age <- levels(as.factor(paleosp$Age2plsd[i]))
                if(as.numeric(age) < 10) age <- paste("00", age, sep="")
                if(as.numeric(age) > 9 & as.numeric(age) < 100) age <- paste("0", age, sep="")
                pathclim <- paste("Data/Paleoclimate/Climate/", age, ".grd", sep="")
                paleoclim <- brick(pathclim)
                pl1sd.presvals <- extract(paleoclim, prespts)  
                
                # Age minus 1 sd
                age <- levels(as.factor(paleosp$Age2mnsd[i]))
                if(as.numeric(age) < 10) age <- paste("00", age, sep="")
                if(as.numeric(age) > 9 & as.numeric(age) < 100) age <- paste("0", age, sep="")
                pathclim <- paste("Data/Paleoclimate/Climate/", age, ".grd", sep="")
                paleoclim <- brick(pathclim)
                mn1sd.presvals <- extract(paleoclim, prespts)
                
                # Age plus 2 sd
                age <- levels(as.factor(paleosp$Age2pl2sd[i]))
                if(as.numeric(age) < 10) age <- paste("00", age, sep="")
                if(as.numeric(age) > 9 & as.numeric(age) < 100) age <- paste("0", age, sep="")
                pathclim <- paste("Data/Paleoclimate/Climate/", age, ".grd", sep="")
                paleoclim <- brick(pathclim)
                pl2sd.presvals <- extract(paleoclim, prespts)
                
                # Age minus 2 sd
                age <- levels(as.factor(paleosp$Age2mn2sd[i]))
                if(as.numeric(age) < 10) age <- paste("00", age, sep="")
                if(as.numeric(age) > 9 & as.numeric(age) < 100) age <- paste("0", age, sep="")
                pathclim <- paste("Data/Paleoclimate/Climate/", age, ".grd", sep="")
                paleoclim <- brick(pathclim)
                mn2sd.presvals <- extract(paleoclim, prespts)
                
                # Weighted average
                w <- dnorm(c(paleosp$Age2[i], 
                             paleosp$Age2plsd[i],
                             paleosp$Age2mnsd[i],
                             paleosp$Age2pl2sd[i],
                             paleosp$Age2mn2sd[i]), paleosp$Age.calib.[i], paleosp$Sd.calib.[i])
                
                t <- c(mean.presvals[,1], pl1sd.presvals[,1], mn1sd.presvals[,1], 
                       pl2sd.presvals[,1], mn2sd.presvals[,1])
                
                r <- c(mean.presvals[,2], pl1sd.presvals[,2], mn1sd.presvals[,2], 
                       pl2sd.presvals[,2], mn2sd.presvals[,2])
                
                if (sum(w) > 0) {
                        # Final values
                        presx <- data.frame(Longitude = prespts[,1], Latitude = prespts[,2], 
                                            Temp = weighted.mean(t, w), Rain = weighted.mean(r, w), 
                                            Age2 = paleosp$Age2[i])
                        presences <- rbind(presences, presx)
                        
                        rm(pathclim, paleoclim, prespts, presx, w, t, r, mean.presvals, pl1sd.presvals, 
                           mn1sd.presvals, pl2sd.presvals, mn2sd.presvals, age)
                } else {
                        # Final values
                        presx <- data.frame(Longitude = prespts[,1], Latitude = prespts[,2], 
                                            Temp = mean.presvals[,1], Rain = mean.presvals[,2], 
                                            Age2 = paleosp$Age2[i])
                        presences <- rbind(presences, presx)
                        
                        rm(pathclim, paleoclim, prespts, presx, mean.presvals, age)
                }
                
                
        }
        else {
                prespts <- paleosp[i,1:2]
                
                # Mean age
                age <- levels(as.factor(paleosp$Age2[i]))
                if(as.numeric(age) < 10) age <- paste("00", age, sep="")
                if(as.numeric(age) > 9 & as.numeric(age) < 100) age <- paste("0", age, sep="")
                pathclim <- paste("Data/Paleoclimate/Climate/", age, ".grd", sep="")
                paleoclim <- brick(pathclim)
                mean.presvals <- extract(paleoclim, prespts)
                
                # Final values
                presx <- data.frame(Longitude = prespts[,1], Latitude = prespts[,2], 
                                    Temp = mean.presvals[,1], Rain = mean.presvals[,2], 
                                    Age2 = paleosp$Age2[i])
                presences <- rbind(presences, presx)
                
                rm(pathclim, paleoclim, prespts, presx, mean.presvals, age)
        }
}

# Remove any duplicate presences
presences <- unique(presences)

# Paleogenus' climate envelope (pseudo-absence will be taken from outside this envelope)
ce.temp <- quantile(presences$Temp, c(0.05, 0.95))
ce.rain <- quantile(presences$Rain, c(0.05, 0.95))

# Fossil pseudoabs
cellids <- as.character(unique(fossils$cellid)) # Pixels with fossils, from which to take pseudo-absences

# Let's extract the geographic coordinates of the pseudo-absences by taken the Lat and Long of first fossil in each pixel
psabs.xy <- data.frame(x = numeric(length(cellids)), y = numeric(length(cellids))) 
for (i in 1:length(cellids)) {
        cid <- subset(fossils, cellid == cellids[i])
        psabs.xy[i,1:2] <- cid[1,c(2,1)]
        rm(cid)
}

# Pseudoabsence selection
paleo_ages <- levels(as.factor(paleosp$Age2))

ages <- paleo_ages

for (i in 1:length(ages)){
        if(as.numeric(ages[i]) < 10) ages[i] <- paste("00", ages[i], sep="")
        if(as.numeric(ages[i]) > 9 & as.numeric(ages[i]) < 100) ages[i] <- paste("0", ages[i], sep="")
}

pseudoabs <- data.frame(Longitude = numeric(), Latitude = numeric(),
                        Temp = numeric(), Rain = numeric(), Age2 = factor())


# Pseudoabsences
for(i in 1:length(ages)){
        pathclim <- paste("Data/Paleoclimate/Climate/", ages[i], ".grd", sep="")
        paleoclim <- brick(pathclim)
        pspagex <- subset(presences, Age2 == paleo_ages[i])
        n <- nrow(pspagex) * 10
        absvals <- as.data.frame(extract(paleoclim, psabs.xy))
        absvals <- subset(absvals, Temp < ce.temp[1] | Temp > ce.temp[2] | Rain < ce.rain[1] | Rain > ce.rain[2])
        if (n <= nrow(absvals)) {absvals <- absvals[sample(nrow(absvals), n), ]}
        pts <- psabs.xy[as.numeric(row.names(absvals)), ]
        pseudoabsx <- data.frame(Longitude = pts[,1], Latitude = pts[,2], Temp = absvals[,1], 
                                 Rain = absvals[,2], Age2 = rep(paleo_ages[i], nrow(absvals)))
        pseudoabs <- rbind(pseudoabs, pseudoabsx)
        rm(pathclim, paleoclim, pspagex, n, absvals, pseudoabsx)
}





## Making training and test groups with k-folding
sbtss <- numeric()
sxtss <- numeric()
sgtss <- numeric()

sbauc <- numeric()
sxauc <- numeric()
sgauc <- numeric()

for (r in 1:100) {
        k <- 5
        pr.group <- kfold(presences, k)
        pa.group <- kfold(pseudoabs, k)
        eb <- list()
        ex <- list()
        eg <- list()
        
        for(i in 1:k){
                pres_train <- presences[pr.group != i, 3:4]
                pres_test <- presences[pr.group == i, 3:4]
                
                backg_train <- pseudoabs[pa.group != i, 3:4]
                backg_test <- pseudoabs[pa.group == i, 3:4]
                
                train <- rbind(pres_train, backg_train)
                pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
                envtrain <- data.frame( cbind(pa=pb_train, train) )
                
                ## Bioclim
                bc <- bioclim(pres_train)
                eb[[i]] <- evaluate(pres_test, backg_test, bc)
                
                ## Maxent
                xm <- maxent(train, envtrain$pa)
                ex[[i]] <- evaluate(pres_test, backg_test, xm)
                
                ## GLM binomial logit
                gm1 <- glm(pa ~ Temp * Rain, family = binomial(link = "logit"), data = envtrain)
                eg[[i]] <- evaluate(pres_test, backg_test, gm1)
        }
        
        # Calculate the True Skill Statistics for each model
        sbtss.i <- sapply( eb, function(x){ max(x@TPR + x@TNR) - 1} )
        sxtss.i <- sapply( ex, function(x){ max(x@TPR + x@TNR) - 1} )
        sgtss.i <- sapply( eg, function(x){ max(x@TPR + x@TNR) - 1} )
        
        sbtss <- c(sbtss, sbtss.i)
        sxtss <- c(sxtss, sxtss.i)
        sgtss <- c(sgtss, sgtss.i)
        
        # Calculate the True Skill Statistics for each model
        sbauc.i <- sapply( eb, function(x){ x@auc } )
        sxauc.i <- sapply( ex, function(x){ x@auc } )
        sgauc.i <- sapply( eg, function(x){ x@auc } )
        
        sbauc <- c(sbauc, sbauc.i)
        sxauc <- c(sxauc, sxauc.i)
        sgauc <- c(sgauc, sgauc.i)
        
        print(r)
}

summary(sbtss)
summary(sxtss)
summary(sgtss)

summary(sbauc)
summary(sxauc)
summary(sgauc)

# Temporal validation
teb <- list()
tex <- list()
teg <- list()

tbtss <- numeric()
txtss <- numeric()
tgtss <- numeric()

tbauc <- numeric()
txauc <- numeric()
tgauc <- numeric()

for (r in 1:100) {
        for(i in 1:length(ages)){
                pres_train <- presences[presences$Age2 != paleo_ages[i], 3:4]
                pres_test <- presences[presences$Age2 == paleo_ages[i], 3:4]
                
                backg_train <- pseudoabs[pseudoabs$Age2 != paleo_ages[i], 3:4]
                backg_test <- pseudoabs[pseudoabs$Age2 == paleo_ages[i], 3:4]
                
                train <- rbind(pres_train, backg_train)
                pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
                envtrain <- data.frame( cbind(pa=pb_train, train) )
                
                ## Bioclim
                bc <- bioclim(pres_train)
                teb[[i]] <- evaluate(pres_test, backg_test, bc)
                
                ## Maxent
                xm <- maxent(train, envtrain$pa)
                tex[[i]] <- evaluate(pres_test, backg_test, xm)
                
                ## GLM binomial logit
                gm1 <- glm(pa ~ Temp * Rain, family = binomial(link = "logit"), data = envtrain)
                teg[[i]] <- evaluate(pres_test, backg_test, gm1)
        }
        
        # Calculate True Skill Statistics for each model
        tbtss.i <- sapply( teb, function(x){ max(x@TPR + x@TNR) - 1} ) # Bioclim
        txtss.i <- sapply( tex, function(x){ max(x@TPR + x@TNR) - 1} ) # MaxEnt
        tgtss.i <- sapply( teg, function(x){ max(x@TPR + x@TNR) - 1} ) # GLM
        
        tbtss <- c(tbtss, tbtss.i)
        txtss <- c(txtss, txtss.i)
        tgtss <- c(tgtss, tgtss.i)
        
        # Calculate the True Skill Statistics for each model
        tbauc.i <- sapply( teb, function(x){ x@auc } )
        txauc.i <- sapply( tex, function(x){ x@auc } )
        tgauc.i <- sapply( teg, function(x){ x@auc } )
        
        tbauc <- c(tbauc, tbauc.i)
        txauc <- c(txauc, txauc.i)
        tgauc <- c(tgauc, tgauc.i)
        
        print(r)
}

summary(tbtss)
summary(txtss)
summary(tgtss)

summary(tbauc)
summary(txauc)
summary(tgauc)


genus.ce.tcv <- data.frame(BC_TSS = tbtss, MX_TSS = txtss, GLM_TSS = tgtss,
                         BC_AUC = tbauc, MX_AUC = txauc, GLM_AUC = tgauc)

genus.ce.st.cv <- data.frame(BC_TSS = sbtss, MX_TSS = sxtss, GLM_TSS = sgtss,
                          BC_AUC = sbauc, MX_AUC = sxauc, GLM_AUC = sgauc)

path.tcv <- paste("Data/CE_CV/", genus,"_temporal_crossval.csv", sep = "")
path.stcv <- paste("Data/CE_CV/", genus,"_spatiotemporal_crossval.csv", sep = "")

write.csv(genus.ce.tcv, path.tcv)
write.csv(genus.ce.st.cv, path.stcv)
