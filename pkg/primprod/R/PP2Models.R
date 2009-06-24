##########################
#The models
#########################

#model for fitting alfa on fitted P-E curves
#this is the most complex model in the two step fit, all reduced models are derived by setting a subset of parameters to a fixed value

alfamod <- function(temp, Si, dia, S,
  ad = 0.01, aa = 0.01, Q10d = 2, Q10a = 2, kSi = 0.005, bS = 0.09) {
  LSi <- Si/(Si+kSi)
  fS = 1+bS*S
  alfaDia <- ad*LSi*Q10d**((temp-10)/10)
  alfaAlg <-  aa*Q10a**((temp-10)/10)
  alfa <- alfaDia*pmax(alfaAlg==0,dia/100) + alfaAlg*(100-dia)/100
  return(fS*alfa)
}

#Exactly the same models are derived for Pm, only the parameters are named differently
Pmmod <- function(temp, Si, dia, S,
  Pmd = 5, Pma = 5, Q10d = 2, Q10a = 2, kSi = 0.005, bS = 0.09) {
  LSi <- Si/(Si+kSi)
  fS = 1+bS*S
  PmDia <- Pmd*LSi*Q10d**((temp-10)/10)
  PmAlg <-  Pma*Q10a**((temp-10)/10)
  Pm <- PmDia*pmax(PmAlg==0,dia/100) + PmAlg*(100-dia)/100
  return(fS*Pm)
}

#And for combined alfa, Pm fitting
PEmod <- function(E, temp, Si, dia, S,
  Pmd = 5, Pma = 5, ad = 0.01, aa = 0.01, Q10d = 2,  Q10a = 2, kSi = 0.005, bS = 0.09) {
  LSi <- Si/(kSi+Si)
  fS <- 1+bS*S
  PPd <- Pmd * LSi * Q10d**(temp/10-1) * (1 - exp(-ad*E/Pmd))
  PPa <- Pma * Q10a**(temp/10-1) * (1 - exp(-aa*E/Pma))
  PP <- PPd *pmax(PPa==0,dia/100) + PPa * (1-dia/100)
  PPS <- PP * fS
  return(PPS)
}

PEmod2 <- function(E, temp, Si, dia, S,
  Pmd = 5, Pma = 5, ad = 0.01, aa = 0.01, Q10da = 2,  Q10dP = 2,
  Q10aa = 2, Q10aP = 2, kSia = 0.005, kSiP = 0.005, bS = 0.09) {
  LSia <- Si/(kSia+Si)
  LSiP <- Si/(kSiP+Si)
  fS <- 1+bS*S
  fTda <- Q10da**(temp/10-1)
  fTdP <- Q10dP**(temp/10-1)
  fTaa <- Q10aa**(temp/10-1)
  fTaP <- Q10aP**(temp/10-1)
  PPd <- Pmd * LSiP * fTdP * (1 - exp(-ad*fTda*LSia*E/(Pmd*fTdP*LSiP)))
  PPa <- Pma * fTaP * (1 - exp(-aa*fTaa*E/(Pma*fTaP)))
  PP <- PPd *pmax(PPa==0,dia/100) + PPa * (1-dia/100)
  PPS <- PP * fS
  return(PPS)
}

PEmodFlynn <- function(E, temp, Si, dia, S,
  Pmd = 5, Pma = 5, ad = 0.01, aa = 0.01,
  Q10d = 2, Q10a = 2, kSi = 0.005, bS = 0.09) {
  LSi <- Si/(kSi+Si)
  fTd <-  Q10d**(temp/10-1)
  fTa <-  Q10a**(temp/10-1)
  fS <- 1+bS*S
  PPd <- Pmd * LSi * fTd * (1 - exp(-ad*E/(Pmd*LSi*fTd)))
  PPa <- Pma * fTa * (1 - exp(-aa*E/(Pma*fTa)))
  PP <- PPd *pmax(PPa==0,dia/100) + PPa * (1-dia/100)
  PPS <- PP * fS
  return(PPS)
}

# P(E) (Simple Platt) 
Platt <- function(E, Pm = 5, alpha = 0.005) {
  PP <- Pm*(1-exp(-alpha*E/Pm))
  return(PP)
}

#Eilers Peeters
EP <- function(E,a = 0.0005, b = 0.1, c = 50)  {
  PP <- E*(a*E^2 + b*E + c)^-1
  return(PP)
}

