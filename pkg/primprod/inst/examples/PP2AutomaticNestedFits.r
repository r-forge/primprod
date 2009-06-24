path = "C:/r-forge/primprod/pkg/R"
setwd(path)
source("PP2Models.R",verbose=F)
source("PP2Functions.R",verbose=F)

require(FME) 
#require(nlme);require(drc) both packages contain routines on (nested) model selection 

#Load data and models
#source("./R/PP2ReadData.R",verbose=F)
path = "C:/r-forge/primprod/pkg/data"
setwd(path)
load("data34b.Rdata")

#######
#A wrapper for the modFit function, since a fixed naming will be used throughout the script
fitWrap <- function(){
modFit(f=Modelcost,p=init,mod=modname,fixpars=fixpars,data=data,observed=observed,datavarnames=datavarnames,modvarnames=modvarnames,logtransform=logtransform,lower=1e-7)
}

########
#0. Estimating alfa and Pm for each individual PE-curve 
########

PEdata <- PP34b[c("E","PPC","time")]
names(PEdata) <- c("E","P","index")

ap <- PEfit(PEdata, "Platt", init = c(a = 0.003,Pm = 5),index="index",
            logtransform=TRUE, replace=TRUE, upper=c(0.15,100))

#Indentical syntax for other PE-models, e.g. Eilers-Peeters
#abc <- PEfit(PEdata,"EP",init = c(a = 0.0005,b = 0.1, c = 50),logtransform=F)

PP34b[c("a","Pm","lin.slope","accept")] <- PEdata[c("a","Pm","lin.slope","accept")]
ZPP34b <- subset(PP34b,distance>90)
apdata <- PP34b[match(unique(PP34b$time),PP34b$time),]
Zapdata <- subset(apdata,distance>90)

########
#A. Seperate two step fits: seperate fit of alpha and Pm
#######

modname <- "alfamod"
logtransform <- F
observed <- "a" #"alfaC"

#Freshwater alpha data and forcings
data <- Zapdata #adata #dfb
datavarnames <- c("Si","Dia","temperature","salinity0")
modvarnames <- c("SiO2","dia","temperature","Cl")
vars <- as.list(data[datavarnames])
names(vars) <- modvarnames

################
#Freshwater alpha nest-scheme

NestScheme <- list(mod0=list(algalfa = 0, kSi = 0, Q10dia = 1, Q10alg = 1,bCl=0),
                    mod1=list(algalfa = 0, kSi = 0, Q10alg = 1, bCl = 0),
                    mod2=list(algalfa = 0, Q10alg = 1, bCl = 0),
                    mod3=list(bCl = 0))

#Freshwater initial values
InitialValues <- c(dialfa = 0.011,algalfa = 0.011, Q10alg = 2.0, Q10dia = 2.0,kSi = 0.05)

AlfaFNested <- nestedFits(model=modname,NestScheme=NestScheme,InitialValues=InitialValues,data=data,observed=observed,datavarnames=datavarnames,modvarnames = modvarnames,logtransform=F)

#################
#Freshwater Pm nested scheme. 
modname <- "Pmmod"
data <- Zapdata[Zapdata$accept,] #dfb[!is.na(dfb$PmC),]
observed <- "Pm" #"PmC"

NestScheme <- list(mod0=list(algP = 0, kSi = 0, Q10dia = 1, Q10alg = 1,bCl=0),
                    mod1=list(algP = 0, kSi = 0, Q10alg = 1, bCl = 0),
                    mod2=list(algP = 0, Q10alg = 1, bCl = 0),
                    mod3=list(bCl = 0))
InitialValues <- c(diaP = 4,algP = 4, Q10alg = 2.0, Q10dia = 2.0,kSi = 0.05)
PmFNested <- nestedFits(model=modname,NestScheme=NestScheme,InitialValues=InitialValues,data=data,observed=observed,datavarnames=datavarnames,modvarnames = modvarnames,logtransform=F)

#################
#FW and brackish alpha nested scheme.
modname <- "alfamod"
data <- apdata #dfClb
observed <- "a" #"alfaC"
NestScheme <- list(mod0=list(algalfa = 0, kSi = 0, Q10dia = 1, Q10alg = 1),
                    mod1=list(algalfa = 0, kSi = 0, Q10alg = 1),
                    mod2=list(algalfa = 0, Q10alg = 1),
                    mod3=NULL)
InitialValues <- c(dialfa = 0.011,algalfa = 0.011, Q10alg = 2.0, Q10dia = 2.0,kSi = 0.05,bCl = 0.05)

AlfaFBNested <- nestedFits(model=modname,NestScheme=NestScheme,InitialValues=InitialValues,data=data,observed=observed,datavarnames=datavarnames,modvarnames = modvarnames,logtransform=F)

#################
#FW and brackish Pm nested scheme.
modname <- "Pmmod"
data <- apdata[apdata$accept,] #dfClb[!is.na(dfClb$PmC),]
observed <- "Pm"        #"PmC"

NestScheme <- list(mod0=list(algP = 0, kSi = 0, Q10dia = 1, Q10alg = 1),
                    mod1=list(algP = 0, kSi = 0, Q10alg = 1),
                    mod2=list(algP = 0, Q10alg = 1),
                    mod3=NULL)
InitialValues <- c(diaP = 4,algP = 4, Q10alg = 2.0, Q10dia = 2.0,kSi = 0.05, bCl = 0.05)
PmFBNested <- nestedFits(model=modname,NestScheme=NestScheme,InitialValues=InitialValues,data=data,observed=observed,datavarnames=datavarnames,modvarnames = modvarnames,logtransform=F)

########
#B. Freshwater joint 2 step nested Scheme
#######
modname1 <- "alfamod"
modname2 <- "Pmmod"
data <- apdata
data$Pm <- replace(data$Pm,!(data$accept),NA)
observed1 <- "a"
observed2 <- "Pm"        #"PmC"

NestScheme <- list(mod0=list(algP = 0, algalfa = 0, kSi = 0, Q10dia = 1, Q10alg = 1),
                    mod1=list(algP = 0, algalfa = 0, kSi = 0, Q10alg = 1),
                    mod2=list(algP = 0, algalfa = 0, Q10alg = 1),
                    mod3=NULL)
InitialValues <- c(diaP = 4,algP = 4, dialfa = 0.011, algalfa = 0.011, Q10alg = 2.0, Q10dia = 2.0,kSi = 0.05, bCl = 0.05)
aPmFBNested <- nestedFits(f=Modelcost,model=list(mod1=modname1,mod2 = modname2),NestScheme=NestScheme,InitialValues=InitialValues,data=data,observed=list(obs1=observed1,obs2=observed2),datavarnames=datavarnames,modvarnames = modvarnames,logtransform=F,weights="sd")

########
#C. Direct fits of P-E data
#######
modname <- "PEmod"
logtransform <- T
observed <- "PPC"

#Freshwater PE data and forcings
data <- ZPP34b  #Zoet34b
datavarnames <- c("E","Si","Dia","temperature","salinity0")
modvarnames <- c("E","Si","Dia","T","Cl")
vars <- as.list(data[datavarnames])
names(vars) <- modvarnames

#Freshwater PE nest scheme
NestScheme <- list(mod0=list(aa = 0, Pma = 1, kSi = 0, Q10d = 1, Q10a = 1, bCl=0),
                    mod1=list(aa = 0, Pma = 1, kSi = 0, Q10a = 1, bCl=0),
                    mod2=list(aa = 0, Pma = 1, Q10a = 1, bCl=0),
                    mod3=list(bCl=0))
                    
InitialValues <- c(aa = 0.011,ad=0.011,Pma=2,Pmd=2,Q10a=2,Q10d=2,kSi=0.01)
PEFNested <- nestedFits(model=modname,NestScheme=NestScheme,InitialValues=InitialValues,data=data,observed=observed,datavarnames=datavarnames,modvarnames = modvarnames,logtransform=logtransform)

#Freshwater and brackish PE data and forcings
data <- PP34b
datavarnames <- c("E","Si","Dia","temperature","salinity0")
modvarnames <- c("E","Si","Dia","T","Cl")
vars <- as.list(data[datavarnames])
names(vars) <- modvarnames

#Freshwater and brackish PE nest scheme
NestScheme <- list(mod0=list(aa = 0, Pma = 1, kSi = 0, Q10d = 1, Q10a = 1),
                    mod1=list(aa = 0, Pma = 1, kSi = 0, Q10a = 1),
                    mod2=list(aa = 0, Pma = 1, Q10a = 1),
                    mod3=NULL)

InitialValues <- c(aa = 0.011,ad=0.011,Pma=2,Pmd=2,Q10a=2,Q10d=2,kSi=0.01,bCl=0.05)
PEFBNested <- nestedFits(model=modname,NestScheme=NestScheme,InitialValues=InitialValues,data=data,observed=observed,datavarnames=datavarnames,modvarnames = modvarnames,logtransform=logtransform)

