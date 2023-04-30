##################################
###Computational Macroeconomics###
##################################

#Ananya Iyengar#

#Setting Working Directory

setwd("C:/Ananya Iyengar/Delhi School of Economics/202_Computational Macroeconomics/Comp Macro Assignment")

################################################################################

#Loading Packages

library(stats)

################################################################################

#Optimisation Exercise 1: Minimising Sum of Squares

#Globally concave function!

#Creating Data 

set.seed(134340)
x <- rnorm(500)
y <- rnorm(500) + 0.7*x
data1 <- data.frame(x,y)

#Creating function for residual sum of squares

sumsq <- function(data, par) {
  with(data, sum((par[1] + par[2]*x - y)^2))
}logl

#Using the optim function to minimise sum of squares

optim1 <- optim(par = c(0.5,0.5), fn = sumsq, data = data1)
optim1

#Running OLS with the same data

lm1 <- summary(lm(y~x, data = data1))
lm1


################################################################################

#Optimisation Exercise 2: Maximum Likelihood Estimation 

#Reparametrised to get a concave function!

#Let us use the same data set as above to assume a normal distribution and execute MLE

#Creating the log likelihood function:

loglik <- function(x, par) {
  m = par[1]
  s = par[2]
  n = length(x)
  ll = -(n/2)*(log(2*pi*s^2)) + (-1/(2*s^2)) * sum((x-m))
  return(-ll)
}

#Running the MLE function

mle <- optim(par = c(0.5,0.5), loglik, x = x)
mle

################################################################################

#Optimisation Exercise 3: Newton-Raphson for Logit Regression 

set.seed(1306)

#Simulating Data 

#Independent Variables
x1 = rnorm(50,3,2) + 0.3*c(1:50)
x2 = rbinom(50, 1, 0.6)
x3 = rpois(n = 50, lambda = 4)
x3[26:50] = x3[26:50] - rpois(n = 25, lambda = 3)

#Dependent Variable 
Y = c(rbinom(20, 1,0.1),rbinom(10, 1,0.25),rbinom(15, 1,0.75),rbinom(5, 1,0.9))


#Bias
x0 <- rep(1, 50)

X = cbind(x0, x1, x2, x3)

#Making the Newton Raphson Algorithm 

manual_logistic_regression = function(X,y,threshold = 1e-10, max_iter = 100) {
  calc_p = function(X,beta)
  {
    beta = as.vector(beta)
    return(exp(X%*%beta) / (1+ exp(X%*%beta)))
  }  
  #initial guess for beta
  beta = rep(0,ncol(X))
  
  #initial value bigger than threshold so that we can enter our while loop 
  diff = 10000 
  
  #counter to ensure we're not stuck in an infinite loop
  iter_count = 0
  
  while(diff > threshold ) #tests for convergence
  {
    #calculate probabilities using current estimate of beta
    p = as.vector(calc_p(X,beta))
    
    #calculate matrix of weights W
    W =  diag(p*(1-p)) 
    
    #calculate the change in beta
    beta_change = solve(t(X)%*%W%*%X) %*% t(X)%*%(y - p)
    
    #update beta
    beta = beta + beta_change
    
    #calculate how much we changed beta by in this iteration 
    #if this is less than threshold, we'll break the while loop 
    diff = sum(beta_change^2)
    
    #see if we've hit the maximum number of iterations
    iter_count = iter_count + 1
    if(iter_count > max_iter) {
      stop("No conv.")
    }
  }
  coef = c("(Intercept)" = beta[1], x1 = beta[2], x2 = beta[3], x3 = beta[4])
  return(coef)
  
}
  
#Running the algorithm

manual_logistic_regression(X,Y)


#Running the same model using the glm() command

output_logit <- glm(Y~x1+x2+x3, family = "binomial")
output_logit














