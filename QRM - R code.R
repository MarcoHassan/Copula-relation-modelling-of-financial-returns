setwd("C:/Users/Marco Hassan/Desktop/Quantitative Risk Management/Assignment")
 
#Import Data
   library(readxl)

Data1 <- read_excel("qrm17HSG_assignmentdata.xlsx")

#Select the last 500 observations
X <- tail(Data1$FIN_Return, 500)
Y <- tail(Data1$HLTH_Return, 500)
 
#Analyze the relationship between the two variables
library(psych)

 
#Create a Matrix for the last 500 observations
e.4 <- matrix(c(X,Y), ncol = 2)
 
colnames(e.4) <-c("Financial", "Health")
 
pairs.panels(e.4, main= "Distribution and Correlation Empirical Observations")
 
#MLE estimation bivariate normal dist.
   
#To check the correctness of the MLE estimation compute first the sample parametrs
mu1 <- mean(X)
mu2 <- mean(Y)
sig1 <- sd(X)
sig2 <- sd(Y)
rho <- cor(X, Y)
 
#From the analytical solution for the MLE of the parameters we know that the first 4 parammeters should be equal to the sample parameters. 
#The rho should also be equal to the sample rho, adjusted by the term n/(n-1)
   
#To estimate these parameters numerically via R (this useful for more complex functions where ana analytical solution for the MLE parameters may not be found)
   
#Use optim function of R
   
#Define a function for the log likelihood of the bivariatee normal dist.
   
#Step 1: specify the log-likelihood function.
wrap <- function(parms, dat){
   mymu1  = parms[1]
   mymu2  = parms[2]
   mysig1 = parms[3]
   mysig2 = parms[4]
   myrho  = parms[5]
   myx1 <- dat[,1]
   myx2 <- dat[,2]
   n = length(myx1)
 
   f <- function(x1=myx1, x2=myx2, mu1=mymu1, mu2=mymu2, sig1=mysig1, sig2=mysig2, rho=myrho){
     -n*(log(sig1) + log(sig2) + 0.5*log(1-rho^2)) - 0.5/(1-rho^2)*(
       sum((x1-mu1)^2)/sig1^2 + sum((x2-mu2)^2)/sig2^2 - 2*rho*sum((x1-mu1)*(x2-mu2))/(sig1*sig2)
     )
   }
   f(x1=myx1, x2=myx2, mu1=mymu1, mu2=mymu2, sig1=mysig1, sig2=mysig2, rho=myrho)
   
 }
 
#Step 2: get very small value for bounding the paramter space to avoid things such as log(0) and consequentely to get a solution.
 
eps <- .Machine$double.eps  

#Step 3: Run the optimation
fit.2 <- optim(rep(0.5,5), wrap, dat= e.4, 
                method="L-BFGS-B", 
                lower = c(-Inf, -Inf, eps, eps, -1+eps), 
                upper = c(Inf, Inf, 100, 100, 1-eps), 
                control = list(fnscale=-1))

#Step 4: Get the result
fit.2

 
#Notice these parameters are very close to the analytical solution. 
#Correct would now to use the analytical parameters. 
#Due to the very small difference I will use here the numerically computed parameters.
#This will have just a negligible effect on the results (statistically speaking these parameters are not even statistically different on very high confidence levels)
   
#Generate a Vector of Returns and a Correlation Matrix
estim.mean <- as.vector(fit.2$par[1:2])
 
#Generate estim. Var-Cov-Matrix
estim.var <- (cbind(as.vector(c(fit.2$par[3], fit.2$par[5])), as.vector(fit.2$par[4:5])))
 
   
#For Model 3
   
#To get the dependece structur in order to get the copula that best fit in the data generate two vectors of uniform distributed variables [0,1] apllying the CDF of the marginals, which are considered normal here
   
# Create column for cumulative dist of the standardized normal X
Data1$R1 <- NA

#Create a column for cumulative dist of the standardized normal y
Data1$R2 <- NA
 
#Generate the cumulative value for standard normal distributed X
for(i in 2698:nrow(Data1)) {
   Data1$R1[i] <- (pnorm((Data1$FIN_Return[i]-fit.2$par[1])/fit.2$par[3]))
 }
 
#Generate the cumulative value for standard normal dist. Y
for(i in 2698:nrow(Data1)) {
   Data1$R2[i] <- (pnorm((Data1$HLTH_Return[i]-fit.2$par[2])/fit.2$par[4]))
 }
 
#View data - notice atm first X-500 standardized with wrong parameters
View(Data1)

#Select the last 500 stand norm values for model 3
U1 <- tail(Data1$R1, 500)
U2 <- tail(Data1$R2, 500)
 
#Create a Matrix combining the two
e <- matrix(c(U1, U2), ncol = 2) 
 
   
#Libraries to work with copulas
library(copula)

#Test the structure of the two uniform dist. cumukative stand. norm.+
pairs.panels(e, main="Distribution and Correlation CDF on Z-Values")
 
#Model 2
# Other way to get the structure 
 
#create pseudo observation vectors trough the function pobs. This carry the info on the dependence structure and can hence also be used to search for the copula that best fits the empirical data.
# Pseudo observations are the observations in the [0,1] interval
u1 <- pobs(as.matrix(cbind(X,Y)))[,1]
u2 <- pobs(as.matrix(cbind(X,Y)))[,2]
 
#Create a new Matrix combining the two
e.2 <- matrix(c(u1, u2), ncol = 2)
 
#Test the structure of the two uniform dist. cumukative stand. norm.+
pairs.panels(e.2, main="Dist. and Correlation of generated pseudo-numbers")
 
#Create a gumbel copula of dimension 2
cop_model <- gumbelCopula(dim=2)
 
#Works for the pseudo - not for the CDF generated uniform dist [0,1] numbers
#Find the best fit for the gumbel copula
fit <- fitCopula(cop_model, e.2, method = 'ml')
fit.4 <- fitCopula(cop_model, e, method = 'ml')
 
#get the MLE alpha for the gumbel copula
fit
 
#To choose the model performing better copare the kendall's taus of the different estimated parameters to the empirical dist tau.
   
#Tau of the empirical dist.
cor(X, Y, method = "kendall")
 
#Tau normal dist model
tau(gumbelCopula(param=1.683))

#Tau for the pseudo observations model
tau(gumbelCopula(param=1.624))

 
#The pseudo observations model seems to work better
   
   
#To estimate best model fit for the M4 model
   
#Here the coupla is gaussian with correlation parameter rho.2
#Notice in M4 the epsilon i correspond to the standard normal distributed X and Y.
   
#generate them
U3 <- qnorm(U1)
U4 <- qnorm(U2)
 
#Method 1 - use pseudo observation matrix
#generate pseudo observation matrix for the two epsilon variables
u3 <- pobs(as.matrix(cbind(U3,U4)))[,1]
u4 <- pobs(as.matrix(cbind(U3,U4)))[,2]
 
#Create a new Matrix combining the two
e.3 <- matrix(c(u3, u4), ncol = 2)
 
#Analyze the distribution
pairs.panels(e.3)
 
#This is the same as before as the dependence and the copula do not change when apllying linear transformations.
   
#Define a gaussian copula of dimension 2
cop_model2 <- normalCopula(dim=2)
 
#find the correlation parameter that best fits the copula
fit.3 <- fitCopula(cop_model2, e.3, method='ml')
 
#Give the Copula back
fit.3

 
#Method 2
#Use cumulative t-distribution to generate uniformely [0,1] distributed numbers.

#First however required to estimate the degrees of freedom that best fits the marginals
#Originally we decided to use df=50 (here the normal dist is already well approximated), as this is a good tdist that fits the distribution of the epsilon as defined in the exercise.
   
#In a second moment we decided to use the df for the t-ditribution, which best fits the data
   
#For that: 
#Step 1: Derive the true number of degrees of freedom through MLE
   
#df for financial data
library(stats4)
 
#Define t-distribution function and fit it with the empirical epsilons
wrap.2 <- function(df){
   -log(prod(((gamma((df+1)/2))/((gamma(df/2))*sqrt(df*pi)))*(((1+(U3**2)/df))**(-(df+1)/2))))
 }
 
#Through ML find the df for which the assumed t-dist best fit the data
mle(wrap.2, start = list(df=30), method = "L-BFGS-B") 


#Same approach to find the df for the health care sector empirical data
   
#Step1: Define t-dist and fit the data
   wrap.3 <- function(df){
   -log(prod(((gamma((df+1)/2))/((gamma(df/2))*sqrt(df*pi)))*(((1+(U4**2)/df))**(-(df+1)/2))))
 }
 
#Through ML find the df for which the assumed t-dist best fit the data
mle(wrap.3, start = list(df=30), method = "L-BFGS-B") 


#Step 1: Generate CDF for the two vectors - margins assumed - t-dist with 4 degrees of freedom
U5 <- pt(U3, df=13.3082)
U6 <- pt(U4, df=16.93849)
 
   
#Step 2: Combine the uniformely distributed random numbers into a matrix
e.5 <- matrix(c(U5,U6), ncol=2)
 
#Analyze the comovement
pairs.panels(e.5)
 
#What happens is the degree of freedom is reduced to 1 - extreme fat tails
U7 <- pt(U3, df =1)
U8 <- pt(U4, df=1)
 
#Matrix creation
e.6 <- matrix(c(U7,U8), ncol=2)
 
#plot relation
pairs.panels(e.6)
 
#Notice how the higher the df the higher the value given to the observations in the tails
#This follows from the lower probability of such events in the case of the normal compared to the tdist with low degree of freedom
#What follows is a higher CDF for the z-value in the normal dist. case.
   
#Estimate best rho for normalCoupula in the df=1 case
fit.5 <- fitCopula(cop_model2, e.5, method="ml")
 
#Output
fit.5

#Estimate best rho for normalCoupula in teh df=30 case
fit.6 <- fitCopula(cop_model2, e.6 , method="ml")
 
#Output
fit.6

 
#Compare models
   
#Higher df - less estimated correlation for the normalCopula
#Hence possible to model higher correlation among the defaults by choosing a lower df for the marginal of teh epsilon
   
#How do thee estimated model fit the data
   
#For df=1 Kendall's tau
tau(normalCopula(0.8567))

 
#For empirical df
tau(normalCopula(0.6546))
 
#Both very far from the empirical dependence structure
   
#Best fit for dependence is again given by the model, which uses pseudo observations
tau(normalCopula(0.5849))

 
   
#For the simulation we decided therefore to use the estimated parameters of the copulas derived using pseudo observations
#In such a way we keep the dependence structure of the empirical data but applying different marginals we can model the full distribution of the returns in order to do stresstests and capture the VaR assuming different distributions of the returns.
   
#Simulation of Model 2
   
#Assumption - returns are bivaratiate distributed with MLE parameters estimated in fit.2
   
#Step 1: Get the Cholesky decomposition for the estimated Var-Cov Matrix
D <- chol(estim.var)
 
#Step 2: generate 10000 standard normal distributed paira
library(MASS)
 
Stand.Norm.Matrix <-matrix(c(rnorm(10000*10000,mean = 0,sd=1)),nrow= 10000, ncol=2)
 
colnames(Stand.Norm.Matrix)<-c("r1","r2")
 
#Compute the 10000 simulated returns
M2 <-t(estim.mean+D%*%t(Stand.Norm.Matrix))
 
#Simulation of Model 4
   
#Needed for the estimation: degrees of freedom of the epsilon that best fit the distribution of the historical data.
   
#We use for simualtion of Model 4 - epsilons are t-dist. with degree of freedom 50, and normal copula with rho=0.5849
   
#Step 1: Define a normalCopula with parameter=0.5849
nCop <- normalCopula(param=0.6546, dim=2)
 
   
#Step 2: generate the multivariate dist. for the epsilons, with coupla - normalCopula [rho=0.58489], marginals - 2 tdist with df=50 each
epsil.multi.dist. <- mvdc(copula = nCop,
                           margins = c("t", "t"),
                           paramMargins = list(list(df=13.3082),
                                                                              +                                               list(df =16.93849)))
 
#Step 3: run 10000 simulations for epsilons pairs
sim.epsilon <- rMvdc(10000, epsil.multi.dist.)
 
View(sim.epsilon)
 
#Check similarity in dependece structure empirical data vs simulated epsilons
   
cor(X, Y, method = "kendall") #empirical

cor(sim.epsilon[,1], sim.epsilon[,2], method ="kendall") #simulated

 
#Similar - simulation OK
   
# As a plus - what happens if we assume a palykurtic dist. of the marginals?
# - do if time left
   
#Export simulated epsilons to evaluate VaR and Expected shortfall in MatLab
write.csv(sim.epsilon, "new simulated epsilons for Matlab")
