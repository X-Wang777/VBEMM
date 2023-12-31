VBEMM
===============

This is a variational Bayesian Expectation Maximization-Maximization (VBEM-M) algorithm for log-linear cognitive diagnostic model (LCDM).

Install the package
---------------

    install.packages("devtools")
    library(devtools)
    install_github("X-Wang777/VBEMM")
    library(VBEMM)

Examples
---------------

#### 1. Example for DINA model

First, genetate a set of data,

    N <- 1000 ## number of students
    J <- 30   ## length of test
    K <- 5    ## number of attribute
    sigma=0.3 ## correlation between attribute
    Q0 <- Qgenerate(J,K,0.5) ## generate a Q-matrix
    s0 <- rep(0.1,J) ## true values of slipping parameter
    g0 <- rep(0.1,J) ## true values of guessing parameter
    generate_data <- DINAdatagenerate(N,J,K,s0,g0,Q0,sigma)
    y <- generate_data$y

Then, execute the VBEM-M algorithm,

    result <- VBEMM(data=y,K=K,Q=Q0,model='DINA')

Finally, get the estimation of slipping parameter $s_0$ and guessing parameter $g_0$,

    lambda_est <-r$para_est$lambda ## slope parameter in LCDM form
    eta_est <- r$para_est$eta      ## intercept parameter in LCDM form
    g_est <- 1/(1+exp(-(eta)))     ## convert to guessing parameter in DINA
    s_est <- 1-1/(1+exp(-(eta +lambda))) ## convert to slipping parameter in DINA
        
#### 2. Example for saturated LCDM 

When based on LCDM, the slope parameter **$\lambda$** is a $J \times 2^K-1$ dimention matrix. Additionally, we need a 0-1 matrix to indicate the active (non-zero) elements in 
 **$\lambda$**, which is denoted as 'lambdaindex' as follows,

    lambdaindex <- matrix(1,J,2^K-1) 
    lambdaindex <- Q_saturatedLCDM(Q0)*lambdaindex
    generate_data <- LCDMdatagenerate(N,J,K,lambda0,eta0,Q0,sigma)
    y <- generate_data$y
    r <- VBEMM(data=y,K=K,Q=Q0,lambdaindex =lambdaindex,model='LCDM',trunc_num=-Inf)













