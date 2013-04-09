# ProgenesisDataImporter from Hubi

#file = "FeatureDataMap_first100_TechReps.csv"
#file = "Top1000Lines.csv"
# result = readProgenesis(file)
#ProgenesisFeatureMap_136 <- readProgenesis("MasterMap_216min.csv")

# 20120430 -> adapted for version 4 of LCMS .. adding info on spectral counts in master map.. 

pgImporter <- function(file, sep=';'){
  
  header = as.matrix(read.table(file, nrows=3, sep=sep, quote='\"', stringsAsFactors=FALSE))
  
  ### extract the info about
  ### - the scores that are in the file
  ### - the grouping
  ##  - the number of replicates
  idx = which(header[1,] != "")
  nSamples = idx[2] - idx[1]
  scoreNames = header[1, idx]
  
  grouping = unlist(header[2, idx[1]:(idx[2]-1)])
  currentGroup = grouping[1]
  for (i in 2:length(grouping)){
    if (grouping[i] == ""){
      grouping[i] = currentGroup
    } else {
      currentGroup = grouping[i]
    }
  }
  dataNames = header[3, ]
  
  ## read the real data
  data = read.table(file, skip=3, sep=sep, quote="\"", stringsAsFactors=FALSE)
  
  
  ### build the result
  result = list()
  result$grouping = grouping
  result$scoreNames = scoreNames
  
  ## the ouput of the progenesis analysis
  ii = 1:(idx[1]-1)
  output = data[ , ii]
  colnames(output) = dataNames[ii]  
  result[["output"]] = output
  
  ## the peptide info
  ii = (idx[length(idx)-1]+nSamples):(ncol(data)-nSamples-1)
  peptideInfo = data[ , ii]
  colnames(peptideInfo) = dataNames[ii]
  result[["peptideInfo"]] = peptideInfo
  
  ## now all the data that there is for the individual samples  
  for (i in 1:length(idx)){
    ii = idx[i]:(idx[i]+nSamples-1)
    x = data[ , ii]
    colnames(x) = dataNames[ii]
    result[[scoreNames[i]]] = x
  }
  
  return(result)
}
