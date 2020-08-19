library(OpenMx)
require(OpenMx)

matrix1 <- matrix(1:9, nrow = 3, ncol = 3)
matrix2 <- matrix(1:9, nrow = 3, dimnames = list(c("X", "Y", "Z")))

A <- matrix(1:3)
B <- matrix(1:3)
C <- matrix(1:6, nrow = 3)
C[1,2] <- 3


#openMX commands 
natA <- mxMatrix(type = "Full", nrow = 3, ncol = 1, values = c(1,2,3), name = "A")
natB <- mxMatrix(type = "Full", nrow = 3, ncol = 1, values = c(1,2,3), name = "B")
natC <- mxMatrix(type = "Full", nrow = 3, ncol = 2, values = c(1,2,3, 4, 5, 6), name = "C")
natD <- mxMatrix(type = "Full", nrow = 3, ncol = 2, values = c(1,2,3, 3, 0.5, 6), name = "D")
natE <- mxMatrix(type = "Stand", nrow = 2, ncol = 2, values = c(0.2), name = "E") #standardized matrix
natF <- mxMatrix(type = "Iden", nrow = 2, ncol = 2, name = "F") #identity matrix
              
?mxMatrix()
?matrix()

#algebra1
algebra1 <- mxAlgebra(expression = A %*% t(A), name = "a1")
model1 <- mxModel("new", natA, natB, algebra1)
fit1 <- mxRun(model1)
fit1$a1

#algebra4
algebra4 <- mxAlgebra(expression = A + B, name = "a4")
model4 <- mxModel("new2", natA, natB, algebra4)
fit4 <- mxRun(model4)
fit4$a4

#algebra8
algebra8 <- mxAlgebra(expression = D %*% t(D), name = "a8")
model8 <- mxModel("new8", natD, algebra8)
fit8 <- mxRun(model8)
fit8$a8

#algebra12
algebra12 <- mxAlgebra(expression = E %*% F, name = "a12") #note inverse of identity matrix is identity matrix 
model12 <- mxModel("new12", natE, natF, algebra12)
fit12 <- mxRun(model12)
fit12$a12

result <- mxEval(list(algebra1, algebra2), fit1) ## try this out later, not working 
