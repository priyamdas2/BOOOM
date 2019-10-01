\name{BOOOMproc}
\alias{BOOOMproc}
\docType{package}
\title{
BOOOMproc
}
\description{
BOOOMproc
}

\usage{
BOOOMproc(X0, A, B, desired_min = -10^(60), s_init = 1,
          no_runs = 1000,max_iter = 10000, rho = 2, phi = 10^(-20),
          tol_fun = 10^(-6), tol_fun_2 = 10^(-15),
          desired_improv_rate = 10^(-20), improv_rate_period = 20,
          total_iter = 10000000, print_output = 0)
}

\examples{
library(pracma)
set.seed(1)
A_nrow <- 100
A_ncol <- 100
A <- rand(A_nrow,A_ncol)
X_true <- randortho(A_ncol)
B <- A\%*\%X_true
X_1 <- diag(A_ncol)

final_mat <- BOOOMproc(X_1, A, B, desired_min = 4986, print_output = 1)
# final_mat\%*\%t(final_mat)
(norm(A\%*\%final_mat - B,"F"))^2
# should be equal to final objective function value returned

}
