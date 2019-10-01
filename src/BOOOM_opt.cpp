#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector BOOOM_opt(NumericMatrix X0, Function func, double desired_min = -10^(60),
                        double s_init = 1, double no_runs = 1000,
                        double max_iter = 10000, double rho = 2, double phi = 10^(-20),
                        double tol_fun = 10^(-6), double tol_fun_2 = 10^(-15),
                        double desired_improv_rate = 10^(-20), double improv_rate_period = 20,
                        double total_iter = 10000000, double print_output = 0){
  int n = X0(_,1).size();
  int p = X0(1,_).size();
  NumericMatrix prod_of_mats(n,n);
  double sum_here;
  double check_orthogonal = 1;
  double epsilon;
  double current_lh;
  double epsilon_temp = 0;
  NumericVector array_of_values(max_iter);
  double cos_value_now;
  double sin_value_now;

  double min_value;

  double best_row;
  double best_col;
  double current_new_lh;
  double temp_sum;




  NumericMatrix THETA(n,p);
  NumericMatrix THETA_possibility;
  NumericMatrix total_lh(p,p);
  NumericMatrix FINAL_MAT;

  double improv_rate = 0;
  double num_total_iter = 0;
  double start_value = 0;


  if(n < 2)
  {Rcout << "Minimum number of rows must be 2. " << std::endl;
  }

  if(n > p)
  {Rcout << "The number of rows must be smaller than or equal to number of columns! " << std::endl;
    return NULL;
  }
  else
  {
    // Checking orthonormality

    for(int ii = 0; ii < n; ++ii)
    {for(int jj = 0; jj < n; ++jj)
    {sum_here = 0;
      for(int kk = 0; kk < p; ++kk)
      {sum_here = sum_here + X0(ii,kk)*X0(kk,jj);
      }
      prod_of_mats(ii,jj) = sum_here;
      if(ii == jj && sum_here !=1)
      {check_orthogonal = 0;
      }
      if(ii != jj && sum_here !=0)
      {check_orthogonal = 0;
      }
    }
    }

    if(check_orthogonal == 0)
    {Rcout << "Rows are NOT orthonormal, taking DEFAULT starting point! " << std::endl;
      for(int ii = 0; ii < n; ++ii)
      {for(int jj = 0; jj < p; ++jj)
      {X0(ii,jj) = 0;
        if(ii == jj)
        {X0(ii,jj) = 1;
        }
      }
      }
    }
    else
    {Rcout << "Starting point is feasible, ROWS are orthonormal, now BOOOM-ing... " << std::endl;}



    // Starting runs

    start_value = Rcpp::as<double>(func(X0));
    NumericMatrix OLD_RUN_THETA = clone(X0);

    for(int iii = 0; iii < no_runs; ++iii)
    {epsilon = s_init*PI;
      if(iii == 0)
      {THETA = clone(X0);
      }
      else
      {THETA = clone(OLD_RUN_THETA);
      }


      // Starting iterations


      for(int i = 0; i < max_iter; ++i)
      {num_total_iter = num_total_iter + 1;
        current_lh = Rcpp::as<double>(func(THETA));


        if(print_output == 1)
        {Rcout << "Run number : " << iii+1 << std::endl;
          Rcout << "iteration number (within current run) : " << i+1 << std::endl;
          Rcout << "total number of iterations : " << num_total_iter << std::endl;
          Rcout << "Current Obj function value :" << current_lh << std::endl;
        }

        epsilon_temp = epsilon;
        cos_value_now = cos(epsilon_temp);
        sin_value_now = sin(epsilon_temp);

        THETA_possibility = clone(THETA);



        for(int rr = 0; rr < p; ++rr)
        {for(int cc = 0; cc < p; ++cc)
        {THETA_possibility = clone(THETA);
          if(rr == cc)
          {total_lh(rr,cc) = current_lh + 123456;
          }
          else
          {// Positive + Negative movement
            for(int vv = 0; vv < n; ++vv)
            {THETA_possibility(vv,rr) = cos_value_now*THETA(vv,rr) + sin_value_now*THETA(vv,cc);
              THETA_possibility(vv,cc) = -sin_value_now*THETA(vv,rr) + cos_value_now*THETA(vv,cc);
            }
            total_lh(rr,cc) = Rcpp::as<double>(func(THETA_possibility));

          }
        }
        }

        // Finding position of minimum values


        min_value = current_lh;
        for(int rr = 0; rr < p; ++rr)
        {for(int cc = 0; cc < p; ++cc)
        {if(rr != cc)
        {if(total_lh(rr,cc) < min_value)
        {min_value = total_lh(rr,cc);
          best_row = rr;
          best_col = cc;
        }
        }
        }
        }

        //Rcout << "min_value" << min_value << std::endl;
        //Rcout << "current_lh" << current_lh << std::endl;



        FINAL_MAT = clone(THETA);

        if(min_value < current_lh)
        {if(best_row != best_col)
        {for(int vv = 0; vv < n; ++vv)
        {FINAL_MAT(vv,best_row) = cos_value_now*THETA(vv,best_row) + sin_value_now*THETA(vv,best_col);
          FINAL_MAT(vv,best_col) = -sin_value_now*THETA(vv,best_row) + cos_value_now*THETA(vv,best_col);
        }
        }
        }




        //Rcout << "min_value" << min_value << std::endl;





        current_new_lh = Rcpp::as<double>(func(FINAL_MAT));


        if(current_new_lh != min_value)   // Just looking for error (if any)
        {Rcout << "Error in the movement part " << std::endl;
        }

        array_of_values[i] = current_new_lh;
        THETA = clone(FINAL_MAT);

        if(i > 0)
        {if(abs(current_new_lh - current_lh) < tol_fun)
        {if(epsilon > phi)
        {epsilon = epsilon/rho;
        }
        else
        {break;
        }
        }
        }
        if(current_new_lh < desired_min)
        {break;
        }


        if(i > improv_rate_period)
        {improv_rate = abs(array_of_values[i] - array_of_values[i-improv_rate_period])/improv_rate_period;
          if(improv_rate < desired_improv_rate)
          {break;
          }
        }

        if(num_total_iter > total_iter)
        {break;
        }

      }  // ....................... End of iteration loop


      temp_sum = 0;
      for(int rr = 0; rr < n; ++rr)
      {for(int cc = 0; cc < p; ++cc)
      {temp_sum = temp_sum + (THETA(rr,cc)-OLD_RUN_THETA(rr,cc))*(THETA(rr,cc)-OLD_RUN_THETA(rr,cc));
      }
      }
      if(temp_sum < tol_fun_2)
      {break;
      }
      OLD_RUN_THETA = clone(THETA);

      if(current_new_lh < desired_min)
      {break;
      }

      if(num_total_iter > total_iter)
      {break;
      }

    }  // ....................... End of run loop

    Rcout << "********** TASK COMPLETE ****************" << std::endl;
    Rcout << "Initial Obj function value :" << start_value << std::endl;
    Rcout << "Final Obj function value :" << current_new_lh << std::endl;


    NumericVector result_array(n*p);

    for(int cc = 0; cc < p; ++cc)    // By column because ".attr" transforms {1,2,3,4}-> [1,3; 2,4]
    {for(int rr = 0; rr < n; ++rr)
    {result_array(n*cc + rr) = THETA(rr,cc);
    }
    }

    result_array.attr("dim") = Dimension(n, p);
    return(result_array);



  }

}
