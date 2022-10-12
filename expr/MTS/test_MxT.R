## This code is used to show the cases that MxT fails to run on AUSLAN data set

#install.packages('R.matlab')
library(R.matlab)
data = readMat('AUS_without.mat') # AUS_caseIV.mat
X = data$X

#install.packages('MixMatrix')
library(MixMatrix)

result0 = MLmatrixt(X, tol = 1e-8, df = 10, fixed = FALSE, max.iter = 1000) 

result1 = MLmatrixt(X, tol = 1e-8, df = 10, fixed = FALSE, max.iter = 1000, row.variance = "AR(1)")

result2 = MLmatrixt(X, tol = 1e-8, df = 10, fixed = FALSE, max.iter = 1000, row.variance = "CS")
