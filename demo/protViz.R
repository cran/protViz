#R

fetuin<-c('MK', 'SFVLLFCLAQLWGCHSIPLDPVAGYK',
'EPACDDPDTEQAALAAVDYINK',
'HLPR', 'GYK', 'HTLNQIDSVK', 'VWPR',
'RPTGEVYDIEIDTLETTCHVLDPTPLANCSVR',
'QQTQHAVEGDCDIHVLK', 'QDGQFSVLFTK',
'CDSSPDSAEDVR', 'K', 'LCPDCPLLAPLNDSR',
'VVHAVEVALATFNAESNGSYLQLVEISR',
'AQFVPLPVSVSVEFAVAATDCIAK',
'EVVDPTK', 'CNLLAEK', 'QYGFCK',
'GSVIQK', 'ALGGEDVR',
'VTCTLFQTQPVIPQPQPDGAEAEAPSAVPDAAGPTPSAAGPPVASVVVGPSVVAVPLPLHR',
'AHYDLR', 'HTFSGVASVESSSGEAFHVGK',
'TPIVGQPSIPGGPVR', 'LCPGR', 'IR', 'YFK', 'I')



plot(pim<-parentIonMass(fetuin), ylab="parent ion mass", main="FETUA_BOVIN")
hist(pim, main="Histogram of parent ion mass of FETUA_BOVIN", xlab="m/Z")

fi<-fragmentIons(fetuin)

peptides<-c('HTLNQIDSVK', 'ALGGEDVR', 'TPIVGQPSIPGGPVR')

pim<-parentIonMass(peptides)
fi<-fragmentIons(peptides)

par(mfrow=c(3,1)); 
for (i in 1:length(peptides)){
    plot(0,0,
        xlab='m/Z',
        ylab='',
        xlim=range(c(fi[i][[1]]$b,fi[i][[1]]$y)),
        ylim=c(0,1),
        type='n',
        axes=FALSE,
        sub=paste("fragment ions of", peptides[i], "/ parent ion mass =", pim[i], "Da"));
    box()
    axis(1,fi[i][[1]]$b,round(fi[i][[1]]$b,2))
    pepSeq<-strsplit(peptides[i],"")
    axis(3,fi[i][[1]]$b,pepSeq[[1]])

    abline(v=fi[i][[1]]$b, col='red',lwd=2) 
    abline(v=fi[i][[1]]$c, col='orange') 
    abline(v=fi[i][[1]]$y, col='blue',lwd=2)
    abline(v=fi[i][[1]]$z, col='cyan')
}

spec<-list(scans=1138,
    title="178: (rt=22.3807) [20080816_23_fetuin_160.RAW]",
    rtinseconds=1342.8402,
    charge=2,
    mZ=c(195.139940, 221.211970, 239.251780, 290.221750, 
    316.300770, 333.300050, 352.258420, 448.384360, 466.348830, 
    496.207570, 509.565910, 538.458310, 547.253380, 556.173940, 
    560.358050, 569.122080, 594.435500, 689.536940, 707.624790, 
    803.509240, 804.528220, 822.528020, 891.631250, 909.544400, 
    916.631600, 973.702160, 990.594520, 999.430580, 1008.583600, 
    1017.692500, 1027.605900),
    intensity=c(931.8, 322.5, 5045, 733.9, 588.8, 9186, 604.6,
    1593, 531.8, 520.4, 976.4, 410.5, 2756, 2279, 5819, 2.679e+05,
    1267, 1542, 979.2, 9577, 3283, 9441, 1520, 1310, 1.8e+04,
    587.5, 2685, 671.7, 3734, 8266, 3309)
)


par(mfrow=c(2,3));
m1<-psm('HTLNQIDSVK', spec,plot=TRUE)
hist(m1$mZ.Da.error,20)
abline(v=c(-5,5),col='grey')
hist(m1$mZ.ppm.error,20)
abline(v=c(-1000,1000),col='grey')

m2<-psm('ENNTHLLVSK', spec,plot=TRUE)
hist(m2$mZ.Da.error,20)
abline(v=c(-5,5),col='grey')
hist(m2$mZ.ppm.error,20)
abline(v=c(-1000,1000),col='grey')


library(lattice)
data(fetuinLFQ)

cv<-1-1:7/10
t<-trellis.par.get("strip.background")
t$col<-(rgb(cv,cv,cv))
trellis.par.set("strip.background",t)

my.xlab="Fetuin concentration spiked into experiment [fmol]"
my.ylab<-"Abundance"

xyplot(abundance~conc|prot*method, data=fetuinLFQ$apex, groups=prot,
    aspect=1,
    main="package:protViz - data(fetuinLFQ)",
    panel = function(x, y, subscripts, groups) {
        if (groups[subscripts][1] == "Fetuin")  {
            panel.fill(col="#ffcccc")
        }
                panel.grid(h=-1,v=-1)
                panel.xyplot(x, y)
                panel.loess(x,y, span=1)
            },
    xlab=my.xlab,
    ylab=my.ylab
)


xyplot(abundance~conc|prot*method,data=fetuinLFQ$empai,groups=prot,
    aspect=1,
    main="package:protViz - data(fetuinLFQ)" ,
    panel = function(x, y, subscripts, groups) {
        if (groups[subscripts][1] == "Fetuin")  {
            panel.fill(col="#ffcccc")
        }
                panel.grid(h=-1,v=-1)
                panel.xyplot(x, y)
                panel.loess(x,y, span=1)
            },
    xlab=my.xlab,
    ylab=my.ylab
)


xyplot(abundance~conc|prot*method,data=fetuinLFQ$t3pq,groups=prot,
    aspect=1,
    main="package:protViz - data(fetuinLFQ)" ,
    panel = function(x, y, subscripts, groups) {
        if (groups[subscripts][1] == "Fetuin")  {
            panel.fill(col="#ffcccc")
        }
                panel.grid(h=-1,v=-1)
                panel.xyplot(x, y)
                panel.loess(x,y, span=1)
        if (groups[subscripts][1] == "Fetuin")  {
            panel.text(min(fetuin.t3pq$conc),
                max(fetuin.t3pq$abundance),
                paste("R-squared:", 
                round(summary(lm(x~y))$r.squared,2)),
                cex=0.75,
                pos=4)
        }
            },
    xlab=my.xlab,
    ylab=my.ylab
)

# iTRAQ
data(iTRAQ)
barchart(pl<-unique(iTRAQ[,1:2]),
    main="package:protViz - data(iTRAQ) / Overview" ,
    col='grey',
    xlab='number of unique peptides',
)

par(mfrow=c(4,4),mar=c(5,5,5,5));
for (i in 3:10){
    hist(asinh(iTRAQ[,i]),
    main=names(iTRAQ)[i])
            qqnorm(asinh(iTRAQ[,i]), main=names(iTRAQ)[i])
            qqline(asinh(iTRAQ[,i]), main=names(iTRAQ)[i])
}

par(mfrow=c(8,8),mar=c(1,1,1,1))
for (i in 3:10){
    for (j in 3:10){
        if (i==j){
            hist(asinh(iTRAQ[,i]), main=names(iTRAQ)[i])
        }
        else if (i<j){
            qqplot(asinh(iTRAQ[,i]) , asinh(iTRAQ[,j]),axes=F);box()
        } else{plot(0,0,type='n',axes=FALSE, xlab="",ylab="")}
        
    }
}

par(mfrow=c(1,1),mar=c(5,5,5,5));
b<-boxplot(asinh(iTRAQ[,c(3:10)]), axes=FALSE,
    main="package:protViz - data(iTRAQ) / QC" )
box()
axis(1,3:10-2,names(iTRAQ)[3:10])
axis(2)


par(mfrow=c(1,5),mar=c(6,5,5,1))
qPeptide<-iTRAQ2GroupAnalysis(data=iTRAQ, group1=c(3,4,5,6), group2=7:10, INDEX=iTRAQ$prot, plot=TRUE, FUN=asinh)

op<-par(mfrow=c(3,5), mar=c(5,4,5,1))
qPeptide<-iTRAQ2GroupAnalysis(data=iTRAQ, group1=c(3,4,5,6), group2=7:10, INDEX=paste(iTRAQ$prot,iTRAQ$peptide), plot=TRUE,FUN=asinh)
par(op)

######################
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
