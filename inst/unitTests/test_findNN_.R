#R


test_findNN_ <-
function(){
    checkEqualsNumeric(findNN_(3.5, 1:5), findNN_(3.5, 1:6), tolerance=0.0)

    DB<-sort(rnorm(100, mean=100, sd=10))
    checkEqualsNumeric(unique(DB[findNN_(DB,DB)] - DB), 0, tolerance=0.0)
}

test_findNN_()
