\name{Twostep}
\alias{Twostep}
\alias{alfamod}
\alias{Pmmod}
\title{Two step fitting models for alfa, Pm}
\description{
  Estimates the values of alfa and Pmax, as a function of temperature,
  salinity, Si and the relative contribution of diatoms.

  \itemize{
    \item alfamod estimates the photosynthetic efficiency
    \item Pmmod estimates the maximal photosynthetic rate
  }

  Used in the two-step fit approach
}
\usage{

alfamod (temp, Si, dia, S,
  ad = 0.01, aa = 0.01, Q10d = 2, Q10a = 2,
  kSi = 0.005, bS = 0.09)

Pmmod (temp, Si, dia, S,
  Pmd = 5, Pma = 5, Q10d = 2, Q10a = 2,
  kSi = 0.005, bS = 0.09)

}
\arguments{
  \item{temp }{the temperature, dgC,
  }
  \item{Si }{dissolved silicate concentration, mmol/m3,
  }
  \item{dia }{percentage diatoms,
  }
  \item{S }{salinity, -
  }
  \item{Pmd }{diatom maximal photosynthetic rate, Pmax, at reference
    temperature of 10 dgC, (gC/gChl/hr),
  }
  \item{Pma }{non-diatom maximal photosynthetic rate, Pmax, at reference
    temperature of 10 dgC, (gC/gChl/hr),
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
  \item{kSi }{half-saturation silicate concentration, used in Monod-limitation
    formula, (mmol Si/m3),
  }
  \item{bS }{linear factor expressing salinity dependence, f(S) = 1+bS*S, (-).
  }


}
\value{
  
  The estimated value of \code{alpha}, community photosynthetic efficiency
  (\code{alfamod}) or \code{Pmax}, community maximal photosynthetic rate
  (\code{Pmmod}).
  
}
\note{
  This implements the models required for the two-step fitting approach.
}

\author{Tom Cox <t.cox@nioo.knaw.nl>

  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{
pr <- par(mfrow=c(2,2))

## Alfa
alfaT <- alfamod(temp=1:20, Si =1, dia = 50, S = 0)
plot(1:20,alfaT, xlab="Temperature", main="alpha")

Si <- seq(0,0.5,0.01)
alfaSi <- alfamod(temp=10, Si =Si, dia = 50, S = 0)
plot(Si, alfaSi, xlab="Silicate", main="alpha")

Dia <- seq(0,100,1)
alfaDia <- alfamod(temp=10, Si =1, dia = Dia, S = 0)
plot(Dia, alfaDia, xlab="\% diatoms", main="alpha")

Sal <- seq(0,35,1)
alfaSal <- alfamod(temp=10, Si =1, dia = 50, S = Sal)
plot(Sal, alfaSal, xlab="Salinity", main="alpha")


## Pmax
PmaxT <- Pmmod(temp=1:20, Si =1, dia = 50, S = 0)
plot(1:20, PmaxT, xlab="Temperature", main="Pmax")

Si <- seq(0,0.5,0.01)
PmaxSi <- Pmmod(temp=10, Si =Si, dia = 50, S = 0)
plot(Si, PmaxSi, xlab="Silicate", main="Pmax")

Dia <- seq(0,100,1)
PmaxDia <- Pmmod(temp=10, Si =1, dia = Dia, S = 0)
plot(Dia, PmaxDia, xlab="\% diatoms", main="Pmax")

Sal <- seq(0,35,1)
PmaxSal <- Pmmod(temp=10, Si =1, dia = 50, S = Sal)
plot(Sal, PmaxSal, xlab="Salinity", main="Pmax")

par(mfrow=pr)
}
\details{

  Alfa is estimated as:

  \deqn{(ad * f_d(T) * M(Si)*dia/100 + aa*(1-dia/100) *
    f_a(T))*g(S)}{(ad*fd(T)*M(si)*dia/100+aa*fa(T)*(1-dia/100))*g(S)}

  Pmax is estimated as:

  \deqn{(Pmd * f_d(T) * M(Si)*dia/100 + Pma*(1-dia/100) *
    f_a(T))*g(S)}{(Pmd*fd(T)*M(si)*dia/100+Pma*fa(T)*(1-dia/100))*g(S)}

  where:
  \deqn{f(t)=Q10^{(T-10)/10}}{f(T)=Q10^((T-10)/10)}
  \deqn{g(S)=(1+bS*S)}{f(S)=1+bS*S}
  \deqn{M(Si)=Si/(Si+kSi)}{M(Si)=Si/(Si+kSi)}

}
\references{%%to be filled in
}
\seealso{
  \code{\link{Platt}} for the Platt model,
  etc
}
\keyword{math}

