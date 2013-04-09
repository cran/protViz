#R
# Mass Deviations Plot
# Christian Panse / Jonas Grossmann <{cp,jg}@fgcz.ethz.ch>
# 2013-02-07

lcmsmap<-function(data){
    s.rtinseconds<-lapply(data, function(x){return (x$rtinseconds)})
    s.pepmass<-lapply(data, function(x){return (x$pepmass)})
    s.intensity<-lapply(data, function(x){return (sum(x$intensity))})
    s.score<-lapply(data, function(x){return (x$mascotScore)})
    s.charge<-as.double(lapply(data, function(x){return (x$charge)}))
    

    plot(s.rtinseconds, s.pepmass,
        type='n',
        main='LC-MS overview',
        sub=data[[1]]$title,
        xlab='rt [seconds]',
        ylab='pepmass')

    points(s.rtinseconds[s.score < 30 & ! is.na(s.score)], s.pepmass[s.score < 30 & ! is.na(s.score)], 
        pch=22, col='grey', bg='grey', cex=0.75)
    points(s.rtinseconds[s.score >= 30 & ! is.na(s.score) ], s.pepmass[s.score >= 30 & ! is.na(s.score)], 
        pch=22, col=rgb(0.1,0.1,0.8,alpha=0.4), bg=rgb(0.1,0.1,0.8,alpha=0.4), cex=1.25)

    legend("topleft", c('cutOffScore30', 'rest'), lwd=c(4), col=c(rgb(0.1,0.1,0.8,alpha=0.9), 'grey')) 

    for (c in 2:max(s.charge)){
        plot(s.rtinseconds, s.pepmass,
            type='n',
            main=paste('LC-MS overview [', c, ' +]', sep=''),
            sub=data[[1]]$title,
            xlab='rt [seconds]',
            ylab='pepmass')

        points(s.rtinseconds[s.score < 30 & ! is.na(s.score) & s.charge == c], s.pepmass[s.score < 30 & ! is.na(s.score) & s.charge == c], 
            pch=22, col='grey', bg='grey', cex=0.75)
        points(s.rtinseconds[s.score >= 30 & ! is.na(s.score) & s.charge == c ], s.pepmass[s.score >= 30 & ! is.na(s.score) & s.charge == c], 
            pch=22, col=rgb(0.1,0.1,0.8,alpha=0.4), bg=rgb(0.1,0.1,0.8,alpha=0.4), cex=1.25)

        legend("topleft", c('cutOffScore30', 'rest'), lwd=c(4), col=c(rgb(0.1,0.1,0.8,alpha=0.9), 'grey')) 
    }
}
