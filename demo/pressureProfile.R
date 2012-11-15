#R
library(lattice)
data(pressureProfile)
pressureProfilePlot(pressureProfile)


op<-par(mfrow=c(1,1), mar=c(5,5,5,1))
pp<-pressureProfileSummary(pressureProfile, time=seq(0,110,by=10))
boxplot(Pc~time, data=pp, main='pressure profiles', xlab='time [in minutes]', ylab='Pc(psi)')
par(op)

pp<-pressureProfileSummary(pressureProfile, time=seq(25,40,by=5))
print(xyplot(Pc ~ as.factor(file) | paste("time =", 
    as.character(time), "minutes"),
    panel = function(x, y){
        m<-sum(y)/length(y)
        m5<-(max(y)-min(y))*0.05
        panel.abline(h=c(m-m5,m,m+m5), 
            col=rep("#ffcccc",3),lwd=c(1,2,1))
        panel.grid(h=-1, v=0)
        panel.xyplot(x, y)
    },
    ylab='Pc [psi]',
    layout=c(1,4),
    sub='The three read lines indicate avg plus min 5%.',
    scales = list(x = list(rot = 45)),
    data=pp))

pp<-pressureProfileSummary(pressureProfile, time=seq(0,140,length=128))
print(levelplot(Pc ~ time * as.factor(file),
    main='pressure profiles',
    sub='color represents Pc(psi)',
    data=pp,
    col.regions=rainbow(100)[1:80]))

pp<-pressureProfileSummary(pressureProfile, time=seq(0,100,length=128))
print(levelplot(Pc ~ time * as.factor(file),
    data=pp,
    main='pressure profiles',
    sub='color represents Pc(psi)',
    col.regions=rainbow(100)[1:80]))
