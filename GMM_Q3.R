################################################
###Generalised Method of Moments:Applications###
################################################

#Ananya Iyengar#

#Setting Working Directory

setwd("C:/Ananya Iyengar/Delhi School of Economics/202_Computational Macroeconomics/Comp Macro Assignment")

################################################################################

#Loading Packages

library(gmm)
library(sandwich)
library(AER)
library(lmtest)
library(dplyr)
library(ivreg)
library(momentfit)
library(readr)
library(tseries)
library(urca)
library(readxl)
library(stargazer)

################################################################################

#IV Estimates using GMM 

#We use the Schooling Returns data set as originally used by David Card.

data("SchoolingReturns", package = "ivreg")

attach(SchoolingReturns)


#Specification: log of wage as a function of education, experience, squared experience, ethnicity, residence in a city and region

#However, education and experience can be endogeneous i.e. correlated with the error term, as is selection into the labour market! 

#Therefore, we use an IV approach! 

#Geographical proximity to college is an instrument for education
#Age is an instrument for experience


#Creating formula 

model1 <- momentModel(log(wage) ~ education + poly(experience, 2) + ethnicity + smsa + south, ~ nearcollege + poly(age, 2) + ethnicity + smsa + south, data = SchoolingReturns, vcov = "MDS")

#Running two stage least squares estimation of the IV model using Base GMM with Heteroskedasticity correction

baseline_gmm_iv <- tsls(model1)
summary(baseline_gmm_iv)

#Running two stage least squares estimation of the IV model using Iterative GMM with Heteroskedasticity correction

it_gmm_iv <- gmm4(log(wage) ~ education + poly(experience, 2) + ethnicity + smsa + south,
                  ~ nearcollege + poly(age, 2) + ethnicity + smsa + south, type = "iter", vcov = "MDS", data = SchoolingReturns)
summary(it_gmm_iv)

#Running the inbuilt ivreg() command for tsls to compare results with that of GMM 

iv <- ivreg(log(wage) ~ education + poly(experience, 2) + ethnicity + smsa + south |
                nearcollege + poly(age, 2) + ethnicity + smsa + south,
              data = SchoolingReturns)
summary(iv)

detach(SchoolingReturns)

################################################################################
################################################################################
################################################################################


#Estimation of Hybrid New Keynesian Phillips Curve (Gali and Gertler) for India 


#Loading data sets

india <- read_excel("india.xlsx")
exch <- read_excel("exch.xls")
m3 <- read_excel("m3.xls")

################################################################################

#Converting the exch and m3 data sets into annual data 

exch_annual <- exch%>%group_by(year)%>%summarise_at(vars(real_ex_rate), list(name = mean))
colnames(exch_annual) <- c("year", "ex_rate")

m3_annual <- m3%>%group_by(year)%>%summarise_at(vars(m3_growth), list(name = mean))
colnames(m3_annual) <- c("year", "gm3")

################################################################################

#Merging Data Sets 

d1 <- inner_join(india, exch_annual, by = "year")

d2 <- inner_join(d1, m3_annual, by = "year")

data <- as.data.frame(d2)

#We have data from 1994 to 2018, corresponding to the start of full current account convertibility!

################################################################################

#Unit Root testing for Inflation


#We perform ADF, KPSS and ERS tests 

adf.test(data$inf) #does not reject null of unit root. since adf has low power in rejecting the null, move to kpss! 

kpss.test(data$inf, null = "Trend") #does not reject null of stationarity!

summary(ur.ers(data$inf, type = "DF-GLS", model = "trend")) #reject the null of random walk with drift! 


#Conclude that inflation is trend stationary! Coefficient of time trend is not significant. Matches results of Chowdhury and Sarkar (2014) 

#No structural break according to this paper! 


################################################################################

#Unit Root Testing for labour shares


#Performing ADF, KPSS and ERS tests

adf.test(data$lab_share) #does not reject the null of unit root

kpss.test(data$lab_share, null = "Trend") #rejects null of stationarity

summary(ur.ers(data$lab_share, type = "DF-GLS", model = "trend")) #does not reject null!

#thus we use the first lag of labour share as an instrument for labour share! This is one proxy for marginal cost! 

################################################################################

#creating the final data set

fd <- data %>% mutate_at(vars(lab_share), list(~ .x - lag(.x)))

fdd <- cbind(fd, data$year)

fdd <- fdd %>% select(year, lab_share)

colnames(fdd) <- c("year", "lab_diff")

data <- inner_join(fdd, data, by = "year")


#Another proxy for MC is detrended log GDP

data$lngdp <- log(data$real_gdp) 

data <- data %>% mutate_at(vars(lngdp), list(~ .x - lag(.x)))

#Creating variable for exchange rate growth

data <- data %>% mutate(exch_growth = (ex_rate - lag(ex_rate))/lag(ex_rate))

################################################################################

#Hybrid NKPC specification: pi_t = gammaf*pi_(t+1) + gammab*pi_(t-1) + lambda*xt + et

#Instruments: m3 growth, exchange rate growth, detrended log gdp

#Baseline GMM Models

nkpc1 <- gmm4(inf ~ lag(inf) + lead(inf) + lab_share, ~ lag(inf) + lead(inf) + lab_share + lab_diff + gm3 + exch_growth + lngdp, vcov = "MDS", data = data)
summary(nkpc1)

nkpc2 <- gmm4(inf ~ lag(inf) + lead(inf) + lngdp,  ~lag(inf) + lead(inf) + lngdp + gm3 + exch_growth, vcov = "MDS", data = data)
summary(nkpc2)


#Iterated GMM Model

nkpc3 <- gmm4(inf ~ lag(inf) + lead(inf) + lngdp,  ~lag(inf) + lead(inf) + lngdp + gm3 + exch_growth, type = "iter", vcov = "MDS", data = data)
summary(nkpc3)

#CU GMM Model

nkpc4 <- gmm4(inf ~ lag(inf) + lead(inf) + lngdp, ~lag(inf) + lead(inf) + lngdp + gm3 + exch_growth, type = "cu", vcov = "MDS", data = data)
summary(nkpc4)


################################################################################




