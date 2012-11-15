#R

# <cp@fgcz.ethz.ch>, Bernd Roschitzki <bernd.roschitzki@fgcz.uzh.ch>, Christian Trachsel <christian.trachsel@fgcz.uzh.ch>

# 2012-08-28 CP,BR
# 2012-09-06 CP,BR
# 2012-09-18 CP,CT
# 2012-11-08 CP,CT


pressureProfilePlot<-function(data){
    # iterate over all files group files having similar rt range
    dd<-round(tapply(data[,2], data[,1], max))

    for (i in unique(sort(dd))){

        op<-par(mfrow=c(1,2))
        s.filter<-as.data.frame(data[data[,1] %in% names(dd[dd==i]),])
        s.filter.names<-as.character(unique(s.filter[,1]))

        cm<-rainbow(length(s.filter.names))

       # split.screen(c(1,2))
       # split.screen(c(2,1),screen=2)

       # screen(4)
       # plot(Qb~time,type="n",cex=0.1,
       #         data=s.filter,
       #         xlab='time', ylab='Qb(nL/min)' , main="flow profile")
       #
       # for (j in c(1:length(s.filter.names))){
       #     lines(Qb ~ time, data=s.filter[s.filter[,1]==s.filter.names[j],], col=cm[j])
       # }
       #
       # idx<-c(1:length(s.filter.names))
       # #legend("bottomleft", s.filter.names[idx], col=cm[idx],pch=22,cex=0.5)

       # screen(3)
        plot(Pc ~ time,type="n",cex=0.1,
                data=s.filter,
                xlab='time', ylab='Pc(psi)' , main="pressure profile")

        for (j in c(1:length(s.filter.names))){
            lines(Pc ~ time,data=s.filter[s.filter[,1]==s.filter.names[j],], col=cm[j])
        }

        idx<-c(1:length(s.filter.names))
       # screen(1)
        plot(0,0,type='n',axes='F',xlab='',ylab='')
        legend("bottomleft", s.filter.names[idx], col=cm[idx],pch=22,cex=0.50)
        #close.screen(all = TRUE)
        par(op)

    }
}

pressureProfileSummary<-function(data, time=seq(min(data$time), max(data$time), by=1)){
    rPc<-numeric()
    rFile<-numeric()
    rTime<-numeric()
    rTimeDiff<-numeric()

    for (i in unique(data$file)){
        s<-data[data$file == i, ]

        idx<-findNN(time,s$time)


        rPc<-c(rPc,s$Pc[idx])

        rTimeDiff<-c(rTimeDiff, s$time[idx]-time)
        rFile<-c(rFile, rep(i, length(idx)))
        rTime<-c(rTime,time)
    }
    return(list(file=rFile, time=rTime, error=rTimeDiff, Pc=rPc))
}

#r<-pressureProfileSummary(pressureProfile, 5)
#par(mfrow=c(2,1))
#boxplot(r$Pc~r$time, xlab='time [in minutes]', ylab='Pc(psi)')
#boxplot(r$Qb~r$time, xlab='time [in minutes]', ylab='Qb(nL/min)')
