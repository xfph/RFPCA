## This code is used to show whether MxT imposing parsimonious structure can run on the AUSLAN dataset
## In order to run, the packages R.matlab and MixMatrix need to be installed.

#install.packages('R.matlab')
library(R.matlab)
#install.packages('MixMatrix')
library(MixMatrix)

data = readMat('AUS_without.mat') # AUS_caseIV.mat
X = data$X

#result0 = MLmatrixt(X, tol = 1e-8, df = 10, fixed = FALSE, max.iter = 1000)

result1 = MLmatrixt(X, tol = 1e-8, df = 10, fixed = FALSE, max.iter = 1000, row.variance = "AR(1)")

result2 = MLmatrixt(X, tol = 1e-8, df = 10, fixed = FALSE, max.iter = 1000, col.variance = "AR(1)")

result3 = MLmatrixt(X, tol = 1e-8, df = 10, fixed = FALSE, max.iter = 1000, row.variance = "AR(1)", col.variance = "AR(1)")

result4 = MLmatrixt(X, tol = 1e-8, df = 10, fixed = FALSE, max.iter = 1000, row.variance = "CS")

result5 = MLmatrixt(X, tol = 1e-8, df = 10, fixed = FALSE, max.iter = 1000, col.variance = "CS")

result6 = MLmatrixt(X, tol = 1e-8, df = 10, fixed = FALSE, max.iter = 1000, row.variance = "CS", col.variance = "CS")