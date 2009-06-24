#path = "C:/r-forge/primprod/pkg/data"
#setwd(path)
#load("data34b.Rdata")   # contains PP34b
require(primprod)


Times <- unique(PP34b$time)

PI <- PP34b[months(PP34b$time)=="March" & PP34b$station =="Kruibeke",]

ap <- PEfit(PI, "Platt", init = c(a = 0.003,Pm = 5),index="time",
            logtransform=TRUE, upper=c(0.15,100))



plot(PI$I,PI$P, pch=c(16,17)[factor(PI$mY)])
legend("bottomright",pch=c(16,17),legend=levels(factor(PI$mY)))
