#R

test_fragmentIons <-
function(){

    m<-fragmentIon('HTLNQIDSVK')[[1]]

    b_assumed<-c(138.0662, 239.1139, 352.1979, 466.2409, 594.2994, 
        707.3835, 822.4104, 909.4425, 1008.5109, 1136.6058)

    lapply(1:length(m$b), function(i){
        checkEqualsNumeric(m$b[i], b_assumed[i], tolerance=1.0e-3)
    })


    checkEqualsNumeric(m$y[5], 561.3242, tolerance=1.0e-3)
}

test_fragmentIons()
