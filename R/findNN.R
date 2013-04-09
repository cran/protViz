#R

# TODO 
# compute score by sum of error div. by number of hits

findNN<-function(q, vec) {

    if (is.unsorted(vec)){
        return (list(error="vec is not sorted"))
    }

    out <- .C("findNN",
        m=as.integer(length(q)),
        n=as.integer(length(vec)),
        q=as.double(q),
        vec=as.double(vec),
        NN=as.integer(rep(-1, length(q))),
        PACKAGE = "protViz")

    return (out$NN+1)
}
