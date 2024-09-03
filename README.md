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

#### 1. Example for DINA mode

    result <- VBEMM(data,K,Q,model='DINA')#data is a n_student by n_item binary response matrix, K is the attribute number, Q is a n_item by K Q-matrix.

Then, we get the estimation of slipping parameter $s_0$ and guessing parameter $g_0$,

    lambda_est <-result$para_est$lambda ## slope parameter in LCDM form
    eta_est <- result$para_est$eta      ## intercept parameter in LCDM form
    g_est <- 1/(1+exp(-(eta)))     ## convert to guessing parameter in DINA
    s_est <- 1-1/(1+exp(-(eta +lambda))) ## convert to slipping parameter in DINA
        
#### 2. Example for saturated LCDM 

When based on LCDM, the slope parameter **$\lambda$** is a $J \times 2^K-1$ dimention matrix. 

    r <- VBEMM(data,K,Q,model='LCDM')

#### 3. Example for constrained LCDM 

When estimate a constrained LCDM, such as LLM, we need a additional parameter 'lambdaindex' which is a 0-1 matrix to indicate the active (non-zero) elements in **$\lambda$**

    r <- VBEMM(data,K,Q,lambdaindex =lambdaindex,model='LCDM')










