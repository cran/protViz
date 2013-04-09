#R

# <cp@fgcz.ethz.ch>, Bernd Roschitzki <bernd.roschitzki@fgcz.uzh.ch>, Christian Trachsel <christian.trachsel@fgcz.uzh.ch>

# 2012-08-28 CP,BR
# 2012-09-06 CP,BR
# 2012-09-18 CP,CT
# 2012-11-08 CP,CT
# 2013-03-14 CP
# 2013-03-19 CP


ppp<-function(data){
    # iterate over all files group files having similar rt range
    dd<-round(tapply(data[,2], data[,1], max))

    for (i in unique(sort(dd))){

        op<-par(mfrow=c(1,2))

        s.filter<-as.data.frame(data[data[,1] %in% names(dd[dd==i]),])
        s.filter.names<-as.character(unique(s.filter[,1]))

         
        cm<-rainbow(length(s.filter.names))
                            
        s.filter.Pc.max <- max(s.filter$Pc, rm.na=TRUE)
        plot(Pc ~ time, type='n', 
                ylim=c(0, s.filter.Pc.max),
                data=s.filter,
                xlab='time', 
                ylab='Pc(psi)', 
                main="pressure profile")

        abline(h=seq(0,s.filter.Pc.max, length=11), col='lightgrey', lty=5)
        axis(4, at=seq(0,s.filter.Pc.max, length=11), paste(seq(0,100, length=11),"%"))

        Pc.max <- max(s.filter$Pc, rm.na=TRUE)
        Qb.max <- max(s.filter$Qb, rm.na=TRUE)
        diff<-numeric()

        for (j in c(1:length(s.filter.names))){

            lines(Pc ~ time, 
                data=s.filter[s.filter[,1]==s.filter.names[j],], 
                col=cm[j])

            lines(Pc.max*(Qb/Qb.max) ~ time, 
                data=s.filter[s.filter[,1]==s.filter.names[j],], 
                col='grey', lty=6)

            data.local<-s.filter[s.filter[,1] == s.filter.names[j],]

            t0 <- data.local$time[data.local$Qb == max(data.local$Qb)] 
            t1.idx <- (data.local$time > t0[1])
            
            t1 <- data.local$time[data.local$Pc == min(data.local$Pc[t1.idx]) & t1.idx]
            
            # abline(v=t0[1],  col='lightgrey') 
            # abline(v=t1[1], col='lightgrey') 
            # diff<-c(diff, round(t1[1]-t0[1]))
        }

        idx<-c(1:length(s.filter.names))
        plot(0,0,type='n',axes='F',xlab='',ylab='')

        legend("bottomleft", 
            paste(s.filter.names[idx], " / t0=", round(t0,2), " / t1=", round(t1, 2), sep=''), 
            pt.bg=cm[idx],
            pt.cex=1.5,
            col=cm[idx],
            pch=22,
            cex=1.25)
        par(op)

    }
}

pps<-function(data, time=seq(min(data$time), max(data$time), by=1)){
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
