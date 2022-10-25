## This code is used to show the result of running the MxT with parsimonious covariance structure $\Sigma_c$: CS and $\Sigma_r$: U
## In order to run, the packages R.matlab and MixMatrix need to be installed.

#install.packages('R.matlab')
library(R.matlab)
data = readMat('AUS_without.mat') 
X = data$X

#install.packages('MixMatrix')
library(MixMatrix)

result<-try(
  MLmatrixt(X, tol = 1e-8, df = 10, fixed = FALSE, max.iter = 1000, row.variance = "CS")
)
if("try-error" %in% class(result))
{
  print("error")
  writeMat('res_case4.mat',Sc = NaN,Sr = NaN)
}else{
  writeMat('res_case4.mat',Sc = result$U,Sr = result$V)
}