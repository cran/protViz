#R

pgLFQaov<-function(data, groups, names, idx=1:nrow(data), 
    plot=FALSE, 
    FUN=function(x){return(x)}){

    a<-rep(NA, nrow(data))

    for (i in idx){
        x.aov<-aov(FUN(as.numeric(data[i,])) ~ groups)
        a[i]<-round(summary(x.aov)[[1]][["Pr(>F)"]][1],4)

        boxplot.color='#ffcccc'
        if (a[i] < 0.05)
            boxplot.color='lightgreen'

        if (plot == TRUE){
            boxplot(FUN(as.numeric(data[i,]))~groups, 
                col=boxplot.color,
                sub=paste("ANOVA Pr(>F):",a[i]),
                main=names[i])
         }
    }
    return(a)
}

