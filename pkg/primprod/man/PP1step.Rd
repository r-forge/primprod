\name{Onestep}
\alias{Onestep}
\alias{EP}
\alias{Platt}
\alias{PEmod}
\alias{PEmod2}
\alias{PEmodFlynn}

\title{Photosynthesis-irradiance models}
\description{
  Calculates the photosynthesis as a function of temperature, salinity,
  Si, and the relative contribution of diatoms.

  \itemize{
    \item EP estimates the Eilers-Peeters model
    \item Platt implements the Platt formulation
    \item PEmod includes temperature dependency only on Pmax
    \item Pemod2 includes temperature dependency on Pmax and alpha
    \item PemodFlynn includes the model sensu Flynn
  }

  Used in the one-step fit approach
}
\usage{
EP (E, a = 0.0005, b = 0.1, c = 50)

Platt (E, Pm = 5, alpha = 0.005)

PEmod (E, temp, Si, dia, S,
  Pmd = 5, Pma = 5, ad = 0.01, aa = 0.01,
  Q10d = 2,  Q10a = 2, kSi = 0.005, bS = 0.09)

PEmod2 <- function(E, temp, Si, dia, S,
  Pmd = 5, Pma = 5, ad = 0.01, aa = 0.01,
  Q10da = 2,  Q10dP = 2, Q10aa = 2, Q10aP = 2,
  kSia = 0.005, kSiP = 0.005, bS = 0.09)
  
PEmodFlynn (E, temp, Si, dia, S,
  Pmd = 5, Pma = 5, ad = 0.01, aa = 0.01,
  Q10d = 2, Q10a = 2, kSi = 0.005, bS = 0.09)

}
\arguments{
  \item{E }{the incident light intensity, (micromol photons/m2/s),
  }
  \item{temp }{the temperature, dgC,
  }
  \item{Si }{dissolved silicate concentration, mmol/m3,
  }
  \item{dia }{percentage diatoms,
  }
  \item{S }{salinity, -
  }
  \item{a, b, c }{fitting parameters,
  }
  \item{Pm }{maximal photosynthetic rate, Pmax, at reference
    temperature of 10 dgC, (gC/gChl/hr),
  }
  \item{Pmd }{diatom maximal photosynthetic rate, Pmax, at reference
    temperature of 10 dgC, (gC/gChl/hr),
  }
  \item{Pma }{non-diatom maximal photosynthetic rate, Pmax, at reference
    temperature of 10 dgC, (gC/gChl/hr),
  }
  \item{alpha }{photosynthetic efficiency, at reference temperature of
    10 dgC, (gC/gChl/hr)/(micromol photons / m2 / s),
  }
  \item{ad }{diatom photosynthetic efficiency, alpha, at reference temperature of
    10 dgC, (gC/gChl/hr)/(micromol photons / m2 / s),
  }
  \item{aa }{non-diatom photosynthetic efficiency, alpha, at reference temperature of
    10 dgC, (gC/gChl/hr)/(micromol photons / m2 / s),
  }
  \item{Q10d }{diatom Q10, parameter determining the temperature dependence,
    the multiplication factor for a 10 dgC temperature increase, (-),
  }
  \item{Q10a }{non-diatom Q10, parameter determining the temperature dependence,
    the multiplication factor for a 10 dgC temperature increase, (-),
  }
  \item{kSi, kSiP, kSia }{half-saturation silicate concentration, used in
    Monod-limitation formula, (mmol Si/m3),
  }
  \item{bS }{linear factor expressing salinity dependence, f(S) = 1+bS*S, (-).
  }
}
\value{
  
  The estimated photosynthesis rate, gC/gChl/hour.

}
\note{
  This implements the models required for the one-step fitting approach.
}

\author{Tom Cox <t.cox@nioo.knaw.nl>

  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{
pr <- par(mfrow=c(2,2))

E <- seq(0,1500, by=10)

## PE-curves
pemod <- PEmod(E, temp=10, Si =1, dia = 50, S = 0)
plot(E,pemod, xlab="Incident light", main="photosynthesis", type="l")

platt <- Platt(E)
lines(E,platt, col="red")

eilers <- EP(E)
lines(E,eilers, col="blue")

flynn <- PEmodFlynn(E, temp=10, Si =1, dia = 50, S = 0)
lines(E,flynn, col="green")

## Function of temperature, Si, Salinity
P <- PEmod(E=100,temp=1:20, Si =1, dia = 50, S = 0)
plot(1:20, P, xlab="Temperature", main="photosynthesis", type="l")

Si <- seq(0,0.5,0.01)
P <- PEmod(E=100,temp=10, Si =Si, dia = 50, S = 0)
plot(Si, P, xlab="Silicate",  main="photosynthesis", type="l")

Sal <- seq(0,35,1)
P <- PEmod(E=100,temp=10, Si =1, dia = 50, S = Sal)
plot(Sal, P, xlab="Salinity",  main="photosynthesis", type="l")

par(mfrow=pr)
}
\details{

  The Eilers-Peeters model is given by:
  \deqn{P=\frac{E}{a E^2 + b E + c}}{P=E/(aE^2+bE+c)}

  The Webb model is:
  \deqn{P=F(E)}{P=F(E)}

  The PEmod is:
  \deqn{P=(Pmd * F_d(E) * f_d(T) * M(Si)*dia/100 + Pma * F_a(E) * (1-dia/100) *
    f_a(T))*g(S)}{(Pmd*Fd(E)*fd(T)*M(si)*dia/100+Pma*Fa(E)*fa(T)*(1-dia/100))*g(S)}

  The PEmod2 is:
  
  ???
  
  
  and where
  
  \deqn{F(E)=1-exp(\frac{-\alpha E}{Pm})}{F(E)=(1-exp((-alpha E)/Pm)}
  \deqn{f(t)=Q10^{(T-10)/10}}{f(T)=Q10^((T-10)/10)}
  \deqn{g(S)=(1+bS*S)}{f(S)=1+bS*S}
  \deqn{M(Si)=Si/(Si+kSi)}{M(Si)=Si/(Si+kSi)}
}
\seealso{
 \code{\link{alfamod}}, \code{\link{Pmmod}} for two-step models
 }
\keyword{math}

