
# Instructions for Problem Set 1 ------------------------------------------

# Using the ‘twindata’ dataset packaged with OpenMx, evaluate the role of genetic and
# environmental factors in explaining variation in body mass index in young males (MZ: zyg=2;
#                                                                                  DZ: zyg=4).
# Provide a summary of results in tables (goodness-of-fit statistics and estimates) and write a
# summary paragraph(s) with key results, addressing model assumptions.
# Maximum 1 page including tables, due Tuesday 9/1/20 at 9am to hmaes@vcu.edu.
#
#go to https://hermine-maes.squarespace.com/#/one/ for sample code to run the various twin models



# Part 1: Make Saturated model (perfect model fit) ----------------------------------------------------------------------------------------------------------------------
# Program: oneSATc.R  
#  Author: Hermine Maes
#    Date: 10 22 2018 
#
# Twin Univariate Saturated model to estimate means and (co)variances across multiple groups
# Matrix style model - Raw data - Continuous data
# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|

# Load Libraries & Options
rm(list=ls())
library(OpenMx)
library(psych); library(polycor)
source("miFunctions.R")
#source('https://openmx.ssri.psu.edu/software/getOpenMx.R')

# Create Output 
filename    <- "oneSATc"
sink(paste(filename,".Ro",sep=""), append=FALSE, split=TRUE)

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA

# Load Data
data(twinData)
dim(twinData)
describe(twinData[,1:12], skew=F)

# Select Variables for Analysis
vars      <- 'bmi'                     # list of variables names
nv        <- 1                         # number of variables
ntv       <- nv*2                      # number of total variables
selVars   <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")

# Select Data for Analysis
mzData    <- subset(twinData, zyg==2, selVars)  #selecting males
dzData    <- subset(twinData, zyg==4, selVars)  #selecting males

# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")

# Set Starting Values
svMe      <- 20                        # start value for means
svVa      <- .8                        # start value for variance
lbVa      <- .0001                     # lower bound for variance

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL

# Create Algebra for expected Mean Matrices
meanMZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mMZ1","mMZ2"), name="meanMZ" )
meanDZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZ1","mDZ2"), name="meanDZ" )

# Create Algebra for expected Variance/Covariance Matrices
covMZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=c("vMZ1","cMZ21","vMZ2"), name="covMZ" )
covDZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=c("vDZ1","cDZ21","vDZ2"), name="covDZ" )

# Create Data Objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="covMZ", means="meanMZ", dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="covDZ", means="meanDZ", dimnames=selVars )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
modelMZ   <- mxModel( meanMZ, covMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( meanDZ, covDZ, dataDZ, expDZ, funML, name="DZ" )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Confidence Interval Objects
ciCov     <- mxCI( c('MZ.covMZ','DZ.covDZ') )
ciMean    <- mxCI( c('MZ.meanMZ','DZ.meanDZ') )

# Build Saturated Model with Confidence Intervals
modelSAT  <- mxModel( "oneSATc", modelMZ, modelDZ, multi, ciCov, ciMean )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run Saturated Model
fitSAT    <- mxRun( modelSAT, intervals=F )
sumSAT    <- summary( fitSAT )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitSAT)
fitEsts(fitSAT)
mxGetExpected( fitSAT, c("means","covariance") )summary(fitSAT)

# ----------------------------------------------------------------------------------------------------------------------
# RUN SUBMODELS

# Constrain expected Means to be equal across Twin Order
modelEMO  <- mxModel( fitSAT, name="oneEMOc" )
modelEMO  <- omxSetParameters( modelEMO, label=c("mMZ1","mMZ2"), free=TRUE, values=svMe, newlabels='mMZ' )
modelEMO  <- omxSetParameters( modelEMO, label=c("mDZ1","mDZ2"), free=TRUE, values=svMe, newlabels='mDZ' )
fitEMO    <- mxRun( modelEMO, intervals=F )
fitGofs(fitEMO); fitEsts(fitEMO)

summary(fitEMO)

# Constrain expected Means and Variances to be equal across Twin Order
modelEMVO <- mxModel( fitEMO, name="oneEMVOc" )
modelEMVO <- omxSetParameters( modelEMVO, label=c("vMZ1","vMZ2"), free=TRUE, values=svVa, newlabels='vMZ' )
modelEMVO <- omxSetParameters( modelEMVO, label=c("vDZ1","vDZ2"), free=TRUE, values=svVa, newlabels='vDZ' )
fitEMVO   <- mxRun( modelEMVO, intervals=F )
fitGofs(fitEMVO); fitEsts(fitEMVO)

# Constrain expected Means and Variances to be equal across Twin Order and Zygosity
modelEMVZ <- mxModel( fitEMVO, name="oneEMVZc" )
modelEMVZ <- omxSetParameters( modelEMVZ, label=c("mMZ","mDZ"), free=TRUE, values=svMe, newlabels='mZ' )
modelEMVZ <- omxSetParameters( modelEMVZ, label=c("vMZ","vDZ"), free=TRUE, values=svVa, newlabels='vZ' )
fitEMVZ   <- mxRun( modelEMVZ, intervals=F )
fitGofs(fitEMVZ); fitEsts(fitEMVZ)

summary(fitEMVZ)

# Print Comparative Fit Statistics
mxCompare( fitSAT, subs <- list(fitEMO, fitEMVO, fitEMVZ) )

# ----------------------------------------------------------------------------------------------------------------------
sink()
save.image(paste(filename,".Ri",sep=""))


# Part 2 ------------------------------------------------------------------

cor(mzData, use = "complete.obs", method = c("pearson")) #rMZ = 0.77, 1/2 rMZ = 0.385
cor(dzData, use = "complete.obs", method = c("pearson")) #rDZ = 0.32
#since rDZ < 1/2 rMZ, we will use ADE model  


corMZ <- mxEval(cov2cor(covMZ), fitSAT$MZ)
corMZ
corDZ <- mxEval(cov2cor(covDZ), fitSAT$DZ)
corDZ 

summary(corMZ)

# Part 3: Use ADE Model. Estimate Contributions of genetic/environmental effects on total variance of phenotype ------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Program: oneADEvc.R  
#  Author: Hermine Maes
#    Date: 10 22 2018 
#
# Twin Univariate ADE model to estimate causes of variation across multiple groups
# Matrix style model - Raw data - Continuous data
# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|

# Load Libraries & Options
rm(list=ls())
library(OpenMx)
library(psych); library(polycor)
source("miFunctions.R")

# Create Output 
filename    <- "oneADEvc"
sink(paste(filename,".Ro",sep=""), append=FALSE, split=TRUE)

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA

# Load Data
data(twinData)
dim(twinData)
describe(twinData[,1:12], skew=F)

# Select Variables for Analysis
vars      <- 'bmi'                     # list of variables names
nv        <- 1                         # number of variables
ntv       <- nv*2                      # number of total variables
selVars   <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")

# Select Data for Analysis
mzData    <- subset(twinData, zyg==2, selVars)
dzData    <- subset(twinData, zyg==4, selVars)

# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")

# Set Starting Values
svMe      <- 20                        # start value for means
svPa      <- .2                        # start value for path coefficient
svPe      <- .5                        # start value for path coefficient for e
# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL

# Create Algebra for expected Mean Matrices
meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labVars("mean",vars), name="meanG" )

# Create Matrices for Variance Components
covA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VA11", name="VA" ) 
covD      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VD11", name="VD" )
covE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="VE11", name="VE" )

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP      <- mxAlgebra( expression= VA+VD+VE, name="V" )
covMZ     <- mxAlgebra( expression= VA+VD, name="cMZ" )
covDZ     <- mxAlgebra( expression= 0.5%x%VA+ 0.25%x%VD, name="cDZ" )
expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )

# Create Data Objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="meanG", dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="meanG", dimnames=selVars )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars      <- list( meanG, covA, covD, covE, covP )
modelMZ   <- mxModel( pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Algebra for Unstandardized and Standardized Variance Components
rowUS     <- rep('US',nv)
colUS     <- rep(c('VA','VD','VE','SA','SD','SE'),each=nv)
estUS     <- mxAlgebra( expression=cbind(VA,VD,VE,VA/V,VD/V,VE/V), name="US", dimnames=list(rowUS,colUS) )

# Create Confidence Interval Objects
ciADE     <- mxCI( "US[1,1:3]" )

# Build Model with Confidence Intervals
modelADE  <- mxModel( "oneADEvc", pars, modelMZ, modelDZ, multi, estUS, ciADE )
# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ADE Model
fitADE    <- mxRun( modelADE, intervals=T )
sumADE    <- summary( fitADE )

# Compare with Saturated Model
#if saturated model fitted in same session
mxCompare( fitSAT, fitADE )
#if saturated model prior to genetic model
#lrtSAT(fitADE,4055.9346,1767)

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitADE)
fitEstCis(fitADE)
summary(fitADE)
# ----------------------------------------------------------------------------------------------------------------------
# RUN SUBMODELS

# Run AE model
modelAE   <- mxModel( fitADE, name="oneAEvc" )
modelAE   <- omxSetParameters( modelAE, labels="VD11", free=FALSE, values=0 )
fitAE     <- mxRun( modelAE, intervals=T )
fitGofs(fitAE); fitEstCis(fitAE)

# Run E model
modelE    <- mxModel( fitAE, name="oneEvc" )
modelE    <- omxSetParameters( modelE, labels="VA11", free=FALSE, values=0 )
fitE      <- mxRun( modelE, intervals=T )
fitGofs(fitE); fitEstCis(fitE)

# Print Comparative Fit Statistics
mxCompare( fitADE, nested <- list(fitAE, fitE) )
mxCompare(fitAE, fitE)
round(rbind(fitADE$US$result,fitAE$US$result,fitE$US$result ),4)

# ----------------------------------------------------------------------------------------------------------------------
sink()
save.image(paste(filename,".Ri",sep=""))
