# calculates B_{M,j}(t)
Btilde <- function(M, j, t) {
  #M - degree of bernstein polynomials, number of basis functions: M+1
  #j - indicates which of the M+1 bernstein polynomials are being evaluated 
  #t - the value the polynomial is being evaluated at, in [-1,1]
  fun <- 1/(2*(M+1))*dbeta((t+1)/2, j+1, M-j+1)
  return(fun)
}

#calculates D_alpha matrix
D_mat_fun <- function(M, t) {
  #t is a vector of length n, t in [-1,1]
  #M is the degree of the bernstein polynomials
  
  values <- sapply(M:0, Btilde, M = M, t = t) #rows are over t, columns are over M:0, values[i,j] is Btilde(M,(M-j),t[i])
  D_mat <- apply(values, 1, cumsum)
  
  return(t(D_mat[(M+1):1,]))
}

######################## Data Generation ########################
set.seed(8)
n <- 100  # number of total individuals
T_max <- 10  # maximum number of time points for an individual 
T_vec <- sample(1:T_max,n,replace = TRUE) # irregular no of time points for individuals
T <- sum(T_vec) # Total number of timepoints/observations
obs <- sum(T_vec) # total number of observations
p <- 4   # dimension of regression coefficients 
M <- 22  # Dimension of basis coefficients is M+1

A <- matrix(1,M+1,M+1)
A[upper.tri(A)] <- 0

#####################  regression  coefficients  ########################
beta <- c(2.75,0.85,2,1.25)   
beta_true <- beta/sqrt(sum(beta^2))

####################  Covariates  ######################
set.seed(6)
x <- matrix(0,obs,p)
x[,1] <- rep(rnorm(n),T_vec)  
x[,2] <- rnorm(obs)          
x[,3] <- rep(rnorm(n),T_vec) 
x[,4] <- rnorm(obs)
x <- x/max(sqrt(rowSums(x*x)))
t <- x%*%beta_true

################### Random Effects ####################
sigma_z <- 0.1
z <- rnorm(n,0,sigma_z)
z <- rep(z,T_vec)

################### Basis Coefficients ################
phi_true <- c(0.15,sample(c(0,0.05,0.1,0.2,0.3,0.4),M,replace = T))
phi <- phi_true
mu <- 1/( 1 + exp(- (D_mat_fun(M,t)%*%as.matrix(phi) + z)))

t_n <- sort(x%*%beta_true)
si_n <- D_mat_fun(M,t_n)%*%as.matrix(phi)

################## Precision Parameter of Beta Distribution ############
psi <- 3
psi <- rep(psi,T)

################# Missing at Random #############
miss_prob <- 0.2
a0 <- log(miss_prob/(1-miss_prob)) 
si_miss <- rep(a0,obs)
prob_miss <- 1/( 1 + exp(- si_miss))

################ Response Variable #################
Y <- rbeta(obs,mu*psi, (1-mu)*psi)

ind_miss <- sample(1:n,miss_prob*n,replace = FALSE)
for (i in ind_miss) {
  Y[(sum(T_vec[1:(i-1)])+1):(sum(T_vec[1:i]))] <- NA
}
pos_obs <- which(!is.na(Y))  # observed positions

fels_Data <- list(m=n,p=p,N=T,T_vec=T_vec,y=replace_na(Y,0.99),
                  y_obs=as.numeric(!is.na(Y)),x=x,M=M,A=A)

