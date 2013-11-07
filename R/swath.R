#R

# <christian.trachsel,jg,cp@fgcz.ethz.ch>
# Wed Nov 28 15:12:28 CET 2012
# Wed Jul 17 10:10:45 CEST 2013
# Wed Aug 21 14:57:11 CEST 2013
# Sat Sep  7 16:52:55 CEST 2013 /  refactored 
# Mon Sep  9 11:02:08 CEST 2013
# Wed Sep 18 16:20:29 CEST 2013 / bugfix mapply / topN peptides
# Thu Sep 19 12:38:50 CEST 2013 bugfix
# Tue Sep 24 15:26:41 CEST 2013 added .normalize_iRT_peptides

# TODO Mon Sep 23 15:08:54 CEST 2013
# iso type label  (Jonas)
# group by fileID for iRT normalization; so there is a normalization for each run?!

# ProteinId;    StrippedSequence ;    iRT;     FragmentGroupId; PrecursorCharge;     PrecursorMz;    FragmentCharge;     FragmentType; FragmentNumber;     FragmentMz;     RelativeFragmentIntensity; ModifiedSequence


# this function does separate regressions iRT normalization GROUP by 'fileName'
.normalize_iRT_peptides_sep <- function(data, iRT=iRTpeptides, plot=FALSE){

    rt <- unlist(lapply(data, function(x){x$rt}))
    peptide <- as.character(unlist(lapply(data, function(x){x$peptideSequence})))
    fileName <-  as.factor(unlist(lapply(data, function(x){x$fileName})))
    df <- data.frame(rt, peptide, fileName)

    data<-aggregate(df$rt, by=list(df$peptide, df$fileName), FUN=mean)
    names(data)<-c('peptide', 'fileName', 'aggregateInputRT')

    m <- merge(iRT, data, by.x='peptide', by.y='peptide')

    # build the model
    # We can do this only of we have found iRT peptides!
    fit <- lm(formula = rt ~ aggregateInputRT * fileName, data=m)

    # apply the model to my data 
    fileName[!fileName %in% m$fileName] <- NA
    rt.predicted <- predict(fit, data.frame(aggregateInputRT=rt, fileName=fileName))

    return(rt.predicted)
}

.normalize_iRT_peptides <- function(data, iRT=iRTpeptides, plot=FALSE){
# TODO
# ignore peptides for fitting iff rt range is higher than 4 minutes

    rt <- unlist(lapply(data, function(x){x$rt}))
    peptide <- unlist(lapply(data, function(x){x$peptideSequence}))
    df <- data.frame(rt, peptide)

    # aggreagte rt if we have more than one psm
    data<-aggregate(df$rt, by=list(df$peptide), FUN=mean)
    names(data)<-c('peptide', 'aggregateRT')

    # check if we have  a sufficient number of iRT peptides and if yes do a merge
    m<-merge(iRT, data, by.x='peptide', by.y='peptide')
    if (nrow(m) < 5){
        warning("no iRT peptides found!")
        return(rt)
    }

    measuredRT <- m$aggregateRT 
    outputRT <- m$rt

    # build the model
    fit <- lm(outputRT ~ measuredRT)

    # apply the model to my data 
    rt.predicted <- predict(fit, data.frame(measuredRT=rt))


    if (plot){

    #pdf("/tmp/swath_debug.pdf", 12,12)
     # TODO: finding the parameter is to heuristic
     slp <- round(fit$coefficients[2], 2)
     yaxis <- round(fit$coefficients[1], 2)
     rSqrd <- round(summary(fit)$r.squared, 4)


    plot(m$rt ~ m$aggregateRT,
        main="fitting retention times, using iRTs", 
        sub=paste("Correlation is based on: ", nrow(m), " out of ", nrow(iRT), " datapoints.", sep=""),
        xlab="measured retention time", 
        ylab="independent rt (iRT)")

    # should be a line of points
    points(rt, rt.predicted, col='lightgrey', pch='.')

    text(m$aggregateRT, m$rt, m$peptide, cex=0.5, pos=4, srt=90)

     legend("bottomright", 
        c(paste("y = ", yaxis ," + x * ",slp,sep=""), 
        paste("R-squared: ", rSqrd, sep="")))
     #dev.off()
    }

    return(rt.predicted)
}

genSwathIonLib <- function(data, mascotIonScoreCutOFF=20, proteinIDPattern='', 
    file="genSwathIonLib.txt", 
    max.mZ.Da.error=0.1, 
    ignoreMascotIonScore=TRUE, 
    topN=10,
    fragmentIonMzRange = c(200, 2000),
    fragmentIonRange = c(2,100), 
    fragmentIonTyp = c('b','y'), 
    iRT=iRTpeptides){

    if (fragmentIonRange[1] < 2){
        fragmentIonRange = c(2,100)
        warning("min fragmentIonRange  should be at least set to 2. reset fragmentIonRange = c(2,100).")
    }

    # data(AA)

    # rt.max<-max(unlist(lapply(GasPhaseFract_reconvoluted_all_plates.blib, function(x){x$rtinseconds})))

    x.peptideSeq<-unlist(lapply(data, function(x){x$peptideSequence}))

    x.varModMass<-(lapply(data, function(x){x$varModification}))

    x.AAmass<-.Call("aa2mass_main", x.peptideSeq , AA$Monoisotopic,  AA$letter1, PACKAGE="protViz")$output

    x.AAmodifiedMass <- mapply(function(x,y){x+y}, x.varModMass, x.AAmass, SIMPLIFY = FALSE)

    x.rt <- unlist(lapply(data, function(x){x$rt}))

    if ( length(iRT) > 1 ){
        if (sum(unlist(lapply(data, function(x){exists("x$fileName")}))) > 0){
            x.rt <- .normalize_iRT_peptides_sep(data, iRT, plot=FALSE)
        }else{
            x.rt <- .normalize_iRT_peptides(data, iRT, plot=FALSE)
        }

    }


    # determine b and y fragment ions while considering the var mods
    fi<-lapply(x.AAmodifiedMass, function(x){fragmentIons(x)[[1]]})
    
    # find NN peak
    findNN.idx<-mapply(function(x.fi,y.data){findNN(c(x.fi$b,x.fi$y), y.data$mZ)}, fi, data, SIMPLIFY = FALSE)

    # determine mZ error
    mZ.error<-mapply(function(x, y.findNN.idx, z){abs(x$mZ[y.findNN.idx] - c(z$b,z$y))}, data, findNN.idx, fi, SIMPLIFY = FALSE)

    # prepare table for output
    output<-mapply (function(x, fi, findNN.idx, mZ.error_, rt){
        m <-length(2 * nchar(x$peptideSequence))
        q1 <-rep(x$pepmass, m)
        q3 <-x$mZ[findNN.idx]
        irt <- round(rep(rt, m), 2)
        frg_z <- rep(1, m)
        decoy <- rep(0, m)
        prec_z <- rep(x$charge, m)
        frg_nr  <- rep(1:length(fi$b),2)

        frg_type <- c(rep('b',length(fi$b)), rep('y', length(fi$y)))

        group_id <- rep(paste(x$peptideSequence,".", x$charge,";", x$pepmass, sep=''), m)

        intensity <- 100 * round(x$intensity[findNN.idx] / max(x$intensity[findNN.idx]), 2)

        peptide_sequence <- rep(x$peptideSequence, m)
        peptideModSeq <- rep(x$peptideModSeq, m)

        res<-cbind(group_id, peptide_sequence, q1, q3, decoy, prec_z, frg_type, frg_nr, frg_z, relativeFragmentIntensity=intensity, irt, peptideModSeq, mZ.error_)

         
        massErrorFilter <- ( (mZ.error_ < max.mZ.Da.error) & (frg_type %in% fragmentIonTyp) & (fragmentIonMzRange[1] < q3 & q3 < fragmentIonMzRange[2]) )
        rr <- res[massErrorFilter, ]


        tryCatch(
        if (fragmentIonRange[1] <= sum(massErrorFilter) & sum(massErrorFilter) <= fragmentIonRange[2]){

            intensity.idx <- rev( order(as.double(rr[,10]) ) )

            if (length(intensity.idx) >= topN){
               return( rr[intensity.idx[1:topN], ] )
            }else{
                return( rr[intensity.idx, ] )
            }
        #}else if (sum(massErrorFilter) == 1){
        #    return (rr)

        }, error=function(e){
            print(e);
            print(rr); 
            print(sum(massErrorFilter)); 
            print("-----")})

    }, data, fi, findNN.idx, mZ.error, x.rt, SIMPLIFY = FALSE)

    xx<-cbind(group_id='', peptide_sequence='', q1='', q3='', decoy='', prec_z='', frg_type='', frg_nr='', frg_z='', relativeFragmentIntensity='', irt='', peptideModSeq='', mZ.error='')

    write.table(xx, file=file, col.names = TRUE, 
        row.names = FALSE, quote = FALSE, sep = "\t", append=FALSE)

    res<-lapply(output, function(x){
            write.table(x, file=file, row.names = FALSE, 
                col.names = FALSE, quote = FALSE, sep = "\t", append=TRUE)
            })

    return(output)
}
