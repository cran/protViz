BUILD:
$R_HOME/bin/R CMD build --resave-data --compact-vignettes --md5 protViz

R CMD build --no-vignettes --resave-data --compact-vignettes --md5 protViz \
&& sudo R CMD INSTALL --build  protViz_0.1.27.tar.gz 

TODO:

#R
> combn(c('S1','T','Y','S2'),2)
     [,1] [,2] [,3] [,4] [,5] [,6]
[1,] "S1" "S1" "S1" "T"  "T"  "Y" 
[2,] "T"  "Y"  "S2" "Y"  "S2" "S2"
> combn(c('S1','T','Y','S2'),3)
     [,1] [,2] [,3] [,4]
[1,] "S1" "S1" "S1" "T" 
[2,] "T"  "T"  "Y"  "Y" 
[3,] "Y"  "S2" "S2" "S2"
> combn(c('S1','T','Y','S2'),3)
     [,1] [,2] [,3] [,4]
[1,] "S1" "S1" "S1" "T" 
[2,] "T"  "T"  "Y"  "Y" 
[3,] "Y"  "S2" "S2" "S2"

