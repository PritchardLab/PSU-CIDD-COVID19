rm(list=ls())
setwd("/mnt/DDrive/Github/Bio/Projects/covid-19/covid19modeling/cpp-v4-6e-severe-classes")

library(curl)
curl_download("http://covidtracking.com/api/states/daily.csv","daily.csv")
daily=read.csv("daily.csv")
#str(daily)

A = daily[ which(daily$state=='RI'),]

dthrum = c(0,31,60,91,121,152,182,213,244,274,305,335)

nr = dim(A)[1]
M = matrix(, nrow = nr, ncol = 4) # columns are daynum, new cases, new negative tests, new deaths,

for (r in 1:nr)
{
  tmp = A$date[r] - 20200000
  thism = floor(tmp/100)
  thisd = tmp%%100
  daynum = dthrum[thism] + thisd
  
  M[nr-r+1,1] = daynum
  if(r < nr)
  {
    M[nr-r+1,2] = A$positive[r] - A$positive[r+1]
    M[nr-r+1,3] = A$negative[r] - A$negative[r+1]
    M[nr-r+1,4] = A$death[r] - A$death[r+1]
  }
  else
  {
    M[nr-r+1,2] = A$positive[r] - 0
    M[nr-r+1,3] = A$negative[r] - 0
    M[nr-r+1,4] = A$death[r] - 0
  }
  
}


#lines( M[,1], 100*M[,2]/(M[,2]+M[,3]))





raw_stdout = system("./odesim none -tf 300 -beta 1.0 2>&1", intern=TRUE)
B = read.table(text=raw_stdout, sep="\t", header=FALSE)

# number of rows
nrB=dim(B)[1]
ncB=dim(B)[2]

# according to the 0-based indices in C++, the infected individuals are in indices 99 to 134 inclusive
# this will be indices 101 to 136 inclusive here in R, since we start at 1, and since time is the first column

# the variable I_total below is now a 'vector' (but not an 'array', thanks R)
I_total = rowSums( B[,101:136] )
J_total = rowSums( B[,272:280] )

Pop_total = rowSums( B[,2:(ncB-9)] )

DeltaJ = J_total[2:nrB] - J_total[1:(nrB-1)]

par(mfrow=c(3,1))
plot( M[,1], M[,2], ylim=c(0,60000), xlim=c(60,200),pch=16 )
lines( B[2:nrB,1], DeltaJ, col='red' )
plot( M[,1], M[,2], ylim=c(-30,200), xlim=c(60,110),pch=16 )
lines( B[2:nrB,1], DeltaJ, col='red' )
plot( B[,1], Pop_total, ylim=c(1058000,1060000), xlim=c(60,310),pch=16 )

