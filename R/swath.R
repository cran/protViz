#R

# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/swath.R $
# $Id: swath.R 6452 2014-05-16 12:14:59Z cpanse $
# $Date: 2014-05-16 14:14:59 +0200 (Fri, 16 May 2014) $


# this function is for normalizing the rt on data
# for building the model data.fit  is used
.normalize_rt <- function(data, data.fit, iRT=iRTpeptides, plot=FALSE){

    message("normalizing RT ...")

    rt <- unlist(lapply(data, function(x){x$rt}))
    rt.fit <- unlist(lapply(data.fit, function(x){x$rt}))

    peptide <- as.character(unlist(lapply(data, function(x){x$peptideSequence})))
    peptide.fit <- as.character(unlist(lapply(data.fit, function(x){x$peptideSequence})))

    fileName <-  as.factor(unlist(lapply(data, function(x){x$fileName})))
    fileName.fit <-  as.factor(unlist(lapply(data.fit, function(x){x$fileName})))

    df <- data.frame(rt, peptide, fileName)
    df.fit <- data.frame(rt=rt.fit, peptide=peptide.fit, fileName=fileName.fit)

    data<-aggregate(df$rt, by=list(df$peptide, df$fileName), FUN=mean)
    data.fit<-aggregate(df.fit$rt, by=list(df.fit$peptide, df.fit$fileName), FUN=mean)

    names(data)<-c('peptide', 'fileName', 'aggregateInputRT')
    names(data.fit)<-c('peptide', 'fileName', 'aggregateInputRT')

    m <- merge(iRT, data.fit, by.x='peptide', by.y='peptide')

    message(paste("found ", length(unique(m$peptide)), " and ", length(unique(m$fileName)), 
	'filename(s).', sep=''))

    # build the model
    # We can do this only of we have found iRT peptides!

    nFileName <- length(unique (as.factor(fileName)))

    if ( nFileName == 1){

	message('model with only one file.')
    	fit <- lm(formula = rt ~ aggregateInputRT, data=m)

    } else if (nFileName > 1){

    	fit <- lm(formula = rt ~ aggregateInputRT * fileName, data=m)

    } else {
	stop("problem in .normalize_rt.")
	}

    # apply the model to my data 
    fileName[!fileName %in% m$fileName] <- NA
    rt.predicted <- predict(fit, data.frame(aggregateInputRT=rt, fileName=fileName))

    return(rt.predicted)
}

.defaultSwathFragmentIon <- function (b, y) {
        Hydrogen <- 1.007825
        Oxygen <- 15.994915
        Nitrogen <- 14.003074
            
        #bn_ <- (b + (n-1) * Hydrogen) / n 

        b1_ <- (b )
        y1_ <- (y ) 

        b2_ <- (b + Hydrogen) / 2
        y2_ <- (y + Hydrogen) / 2 

        b3_ <- (b + 2 * Hydrogen) / 3
        y3_ <- (y + 2 * Hydrogen) / 3

        return( cbind(b1_, y1_, b2_, y2_, b3_, y3_) )
}

# TODO Check why is q3 and a1 not in-silico?
genSwathIonLib <- function(data, 
    mascotIonScoreCutOFF=20, 
    proteinIDPattern='', 
    file = "genSwathIonLib.txt", 
    max.mZ.Da.error = 0.1, 
    ignoreMascotIonScore = TRUE, 
    topN = 10,
    fragmentIonMzRange = c(200, 2000),
    fragmentIonRange = c(2,100), 
    fragmentIonFUN = .defaultSwathFragmentIon, 
    iRT = iRTpeptides,
    data.fit = data){

    if (fragmentIonRange[1] < 2){
        fragmentIonRange = c(2,100)
        warning("min fragmentIonRange should be at least set to 2. reset fragmentIonRange = c(2,100).")
    }

    x.peptideSeq<-unlist(lapply(data, function(x){x$peptideSequence}))

    x.varModMass<-(lapply(data, function(x){x$varModification}))

    x.AAmass<-.Call("aa2mass_main", x.peptideSeq , AA$Monoisotopic,  AA$letter1, PACKAGE="protViz")$output

    x.AAmodifiedMass <- mapply(function(x,y){x+y}, x.varModMass, x.AAmass, SIMPLIFY = FALSE)

    x.rt <- unlist(lapply(data, function(x){x$rt}))

    if (length(iRT) > 1 & length(data.fit) > 1){
            x.rt <- .normalize_rt(data, data.fit, iRT, plot=FALSE)
    }

    # determine b and y fragment ions while considering the var mods
    fi<-lapply(x.AAmodifiedMass, function(x){fragmentIon(x, fragmentIonFUN)[[1]]})
    fragmentIonTyp = names(fi[[1]]) 
    
    # find NN peak
    #findNN.idx<-mapply(function(x.fi,y.data){findNN_(c(x.fi$b,x.fi$y), y.data$mZ)}, fi, data, SIMPLIFY = FALSE)
    findNN.idx<-mapply(function(x.fi, y.data){findNN_(unlist(x.fi), y.data$mZ) }, fi, data, SIMPLIFY = FALSE)

    # determine mZ error
    mZ.error<-mapply(function(x, y.findNN.idx, z){
                abs(x$mZ[y.findNN.idx] - unlist(z))
            }, data, findNN.idx, fi, SIMPLIFY = FALSE)

    # prepare table for output
    output<-mapply (function(x, fi, findNN.idx, mZ.error_, rt){
        m <-length(2 * nchar(x$peptideSequence))
        q1 <-rep(x$pepmass, m)
        q3 <-x$mZ[findNN.idx]
        irt <- round(rep(rt, m), 2)
        decoy <- rep(0, m)
        prec_z <- rep(x$charge, m)

        frg_type <- gsub("(.*)([0-9]+)_([0-9]+)", "\\1",  names(unlist(fi)))
        frg_z <- gsub("(.*)([0-9]+)_([0-9]+)", "\\2",  names(unlist(fi)))
        frg_nr  <- rep(1:nrow(fi), ncol(fi))  

        group_id <- rep(paste(x$peptideSequence, ".", x$charge,";", x$pepmass, sep=''), m)

        intensity <- 100 * round(x$intensity[findNN.idx] / max(x$intensity[findNN.idx], na.rm=TRUE), 2)

        peptide_sequence <- rep(x$peptideSequence, m)

        peptideModSeq <- rep(x$peptideModSeq, m)

        res <- cbind(group_id, peptide_sequence, q1, q3, decoy, prec_z, frg_type, frg_nr, frg_z, relativeFragmentIntensity=intensity, irt, peptideModSeq, mZ.error_)
         
        massErrorFilter <- ( (mZ.error_ < max.mZ.Da.error) & (fragmentIonMzRange[1] < q3 & q3 < fragmentIonMzRange[2]) )

        res.massErrorFilter <- res[massErrorFilter, ]

        tryCatch(

        if (fragmentIonRange[1] <= sum(massErrorFilter) & sum(massErrorFilter) <= fragmentIonRange[2] & length(res.massErrorFilter[, 10]) > 0){

            intensity.idx <- rev( order(as.double(res.massErrorFilter[, 10] ) ) )

            if (length(intensity.idx) < topN){

                return( res.massErrorFilter[intensity.idx, ] )

            }else{

               return( res.massErrorFilter[intensity.idx[1:topN], ] )

            }

        }, error=function(e){
            # warning(e);
         })

    }, data, fi, findNN.idx, mZ.error, x.rt, SIMPLIFY = FALSE)


    xx<-cbind(group_id='', peptide_sequence='', q1='', q3='', decoy='', prec_z='', frg_type='', frg_nr='', frg_z='', relativeFragmentIntensity='', irt='', peptideModSeq='', mZ.error='')

    write.table(xx, file=file, col.names = TRUE, 
        row.names = FALSE, quote = FALSE, sep = "\t", append=FALSE)

    res<-lapply(output, function(x){
            write.table(x, file=file, row.names = FALSE, 
                col.names = FALSE, quote = FALSE, sep = "\t", append=TRUE)
            })

    return(list(output, rt=unlist(x.rt), rt.org=unlist(lapply(data, function(x){x$rt}))))
}
