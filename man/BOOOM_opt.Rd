\name{BOOOM_opt}
\alias{BOOOM_opt}
\docType{package}
\title{
BOOOM_opt
}
\description{
BOOOM_opt
}

\usage{
BOOOM_opt(X0, Function func, desired_min = -10^(60), s_init = 1,
          no_runs = 1000,max_iter = 10000, rho = 2, phi = 10^(-20),
          tol_fun = 10^(-6), tol_fun_2 = 10^(-15),
          desired_improv_rate = 10^(-20), improv_rate_period = 20,
          total_iter = 10000000, print_output = 0)
}

\examples{
library(pracma)
set.seed(1)
A_nrow <- 4
A_ncol <- 3
A <- rand(A_nrow,A_ncol)
X_true <- randortho(A_ncol)
B <- A\%*\%X_true

func <- function(X){return(norm(A\%*\%X - B,"F"))}
func(X_true)
X_1 <- diag(A_ncol)
final_mat <- BOOOM_opt(X_1, func, print_output = 0)
func(final_mat)
final_mat <- BOOOM_opt(X_1, func, desired_min = 0.5, print_output = 1)
func(final_mat)
}
