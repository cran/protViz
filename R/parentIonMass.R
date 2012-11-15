#R
#fetuinPeptides<-c('LCPGR', 'GSVIQK', 'QYGFCK', 'AHYDLR', 'EVVDPTK', 'CNLLAEK', 'ALGGEDVR', 'HTLNQIDSVK', 'QDGQFSVLFTK', 'CDSSPDSAEDVR', 'TPIVGQPSIPGGPVR', 'LCPDCPLLAPLNDSR', 'QQTQHAVEGDCDIHVLK', 'HTFSGVASVESSSGEAFHVGK', 'EPACDDPDTEQAALAAVDYINK', 'AQFVPLPVSVSVEFAVAATDCIAK', 'SFVLLFCLAQLWGCHSIPLDPVAGYK', 'VVHAVEVALATFNAESNGSYLQLVEISR', 'RPTGEVYDIEIDTLETTCHVLDPTPLANCSVR', 'VTCTLFQTQPVIPQPQPDGAEAEAPSAVPDAAGPTPSAAGPPVASVVVGPSVVAVPLPLHR')

parentIonMass <- function(sequence) {
    if (!is.character(sequence)) {
        stop ("argument x must be a character")
    }else{
           out <- .C("computeParentIonMass",
                   n=as.integer(length(sequence)),
                   pepSeq=as.character(sequence),
                   pim=as.double(rep(0,length(sequence)))
                   )
    }
    return(out$pim)
}
