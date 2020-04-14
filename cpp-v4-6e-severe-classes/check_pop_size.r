rm(list=ls())
setwd("/mnt/DDrive/Github/Bio/Projects/covid-19/covid19modeling/cpp-v4-6e-severe-classes")


dthrum = c(0,31,60,91,121,152,182,213,244,274,305,335)


raw_stdout = system("./odesim none -tf 300 -beta 1.0 1.0 1.0 0.8 2>&1", intern=TRUE)
B = read.table(text=raw_stdout, sep="\t", header=FALSE)

# number of rows
nrB=dim(B)[1]
ncB=dim(B)[2]


Pop_total = rowSums( B[,2:(ncB-9)] )


#par(mfrow=c(3,1))
#plot( M[,1], M[,2], ylim=c(0,60000), xlim=c(60,200),pch=16 )
#lines( B[2:nrB,1], DeltaJ, col='red' )
#plot( M[,1], M[,2], ylim=c(0,200), xlim=c(60,110),pch=16 )
#lines( B[2:nrB,1], DeltaJ, col='red' )
plot( B[,1], Pop_total, ylim=c(1058000,1060000), xlim=c(60,310),pch=16 )

