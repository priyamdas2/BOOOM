#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector BOOOMproc(NumericMatrix X0,NumericMatrix A,NumericMatrix B, double desired_min = -10^(60),
                        double s_init = 1,double no_runs = 1000, double max_iter = 10000, double rho = 2,
                        double phi = 10^(-20),double tol_fun = 10^(-6), double tol_fun_2 = 10^(-15),
                        double desired_improv_rate = 10^(-20), double improv_rate_period = 20,
                        double total_iter = 10000000, double print_output = 0){

  int nrowsA = A(_,1).size();
  int ncolsA = A(1,_).size();
  int nrowsB = B(_,1).size();
  int ncolsB = B(1,_).size();


  int n = X0(_,1).size();
  int p = X0(1,_).size();
  Rcout << "n = ncolsA =  " << n << std::endl;
  Rcout << "p = nrowsB =  " << p << std::endl;

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

  double overall_new_min;
  double best_row;
  double best_col;
  double current_new_lh;
  double temp_sum;
  double change_status = 0;




  NumericMatrix THETA(n,p);
  NumericMatrix THETA_possibility;
  NumericMatrix total_lh(n,p);

  NumericMatrix FINAL_MAT;



  NumericVector new_sum_rr(nrowsA);
  NumericVector old_sum_rr(nrowsA);
  double diffs_rr = 0;
  NumericVector new_sum_cc(nrowsA);
  NumericVector old_sum_cc(nrowsA);
  double diffs_cc = 0;

  double improv_rate = 0;
  double num_total_iter = 0;


  double diff;
  double diff_square;
  double SUM_AX_minus_B = 0;
  double start_value = 0;
  double arbit;

  if(n < 2)
  {Rcout << "Minimum number of rows must be 2. " << std::endl;
  }


  if(n > p)
  {Rcout << "The number of rows of X must be smaller than or equal to number of columns! " << std::endl;
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


    // *** Only calculating square of Frob norm from X0,A,B  ************

    SUM_AX_minus_B = 0;

    for(int ii = 0; ii < nrowsA; ++ii)
    {for(int jj = 0; jj < ncolsB; ++jj)
    {sum_here = 0;
      for(int kk = 0; kk < ncolsA; ++kk)
      {sum_here = sum_here + A(ii,kk)*X0(kk,jj) ;
      }
      diff = (sum_here - B(ii,jj));
      diff_square = diff*diff;
      SUM_AX_minus_B = SUM_AX_minus_B+diff_square;
    }
    }


    start_value = SUM_AX_minus_B;

    // ************* Frob norm square calculated **************************

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

        // Calculating objective function value ****************

        // *** Only calculating square of Frob norm from X0,A,B  ************

        SUM_AX_minus_B = 0;

        for(int ii = 0; ii < nrowsA; ++ii)
        {for(int jj = 0; jj < ncolsB; ++jj)
        {sum_here = 0;
          for(int kk = 0; kk < ncolsA; ++kk)
          {sum_here = sum_here + A(ii,kk)*THETA(kk,jj) ;
          }
          diff = (sum_here - B(ii,jj));
          diff_square = diff*diff;
          SUM_AX_minus_B = SUM_AX_minus_B+diff_square;
        }
        }
        current_lh = SUM_AX_minus_B;

        // *****************************************************************
        //current_lh = Rcpp::as<double>(func_proc(THETA,A,B));
        //Rcout << "SUM_AX_minus_B" << SUM_AX_minus_B << std::endl;
        //Rcout << "Rcpp::as<double> current_lh" << current_lh << std::endl;


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


            diffs_rr = 0;
            diffs_cc = 0;
            for(int hhh = 0; hhh <  nrowsA; ++hhh)
            {new_sum_rr[hhh] = 0;
              old_sum_rr[hhh] = 0;
              new_sum_cc[hhh] = 0;
              old_sum_cc[hhh] = 0;
              for(int kkk = 0; kkk < ncolsA; ++kkk)
              {new_sum_rr[hhh] = new_sum_rr[hhh] + A(hhh,kkk)*THETA_possibility(kkk,rr);
                old_sum_rr[hhh] = old_sum_rr[hhh] + A(hhh,kkk)*THETA(kkk,rr);

                new_sum_cc[hhh] = new_sum_cc[hhh] + A(hhh,kkk)*THETA_possibility(kkk,cc);
                old_sum_cc[hhh] = old_sum_cc[hhh] + A(hhh,kkk)*THETA(kkk,cc);
              }
              diffs_rr = diffs_rr + pow(new_sum_rr[hhh] - B(hhh,rr),2)  - pow(old_sum_rr[hhh] - B(hhh,rr),2);
              diffs_cc = diffs_cc + pow(new_sum_cc[hhh] - B(hhh,cc),2)  - pow(old_sum_cc[hhh] - B(hhh,cc),2);
            }
            arbit = current_lh + diffs_rr + diffs_cc;
            total_lh(rr,cc) = arbit;
            //total_lh(rr,cc) = Rcpp::as<double>(func_proc(THETA_possibility,A,B));


            //Rcout << "arbit positive" << arbit << std::endl;
            //Rcout << "total_lh(rr,cc)" << total_lh(rr,cc) << std::endl;

          }

        }
        }

        // Finding position of minimum values


        min_value = current_lh;

        for(int rr = 0; rr < p; ++rr)
        {for(int cc = 0; cc < p; ++cc)
        {if(rr != cc)
        {if(total_lh(rr,cc) <= min_value)
        {min_value = total_lh(rr,cc);
          best_row = rr;
          best_col = cc;
        }
        }
        }
        }


        FINAL_MAT = clone(THETA);

        change_status = 0;

        overall_new_min = current_lh;


        if(min_value < current_lh)
        {change_status = 1;
          overall_new_min = min_value;
          for(int vv = 0; vv < n; ++vv)
          {FINAL_MAT(vv,best_row) = cos_value_now*THETA(vv,best_row) + sin_value_now*THETA(vv,best_col);
            FINAL_MAT(vv,best_col) = -sin_value_now*THETA(vv,best_row) + cos_value_now*THETA(vv,best_col);
          }
        }


        if(min_value > current_lh)
        {overall_new_min = current_lh;
          FINAL_MAT = clone(THETA);
        }



        // *** Only calculating square of Frob norm from X0,A,B  ************

        SUM_AX_minus_B = 0;

        for(int ii = 0; ii < nrowsA; ++ii)
        {for(int jj = 0; jj < ncolsB; ++jj)
        {sum_here = 0;
          for(int kk = 0; kk < ncolsA; ++kk)
          {sum_here = sum_here + A(ii,kk)*FINAL_MAT(kk,jj) ;
          }
          diff = (sum_here - B(ii,jj));
          diff_square = diff*diff;
          SUM_AX_minus_B = SUM_AX_minus_B+diff_square;
        }
        }
        current_new_lh = SUM_AX_minus_B;

        // *****************************************************************







        //if(abs(current_new_lh - overall_new_min) > pow(10,-3))  // Just looking for error (if any)
        //{Rcout << "Error in the movement part " << std::endl;
        // Rcout << "current_new_lh " << current_new_lh <<  std::endl;
        // Rcout << "overall_new_min " << overall_new_min <<  std::endl;
        // Rcout << "Change status " << change_status <<  std::endl;
        // Rcout << "***********************************************" <<  std::endl;
        // }

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

      if(iii>1)
      {temp_sum = 0;
        for(int rr = 0; rr < n; ++rr)
        {for(int cc = 0; cc < p; ++cc)
        {temp_sum = temp_sum + (THETA(rr,cc)-OLD_RUN_THETA(rr,cc))*(THETA(rr,cc)-OLD_RUN_THETA(rr,cc));
        }
        }
        if(temp_sum < tol_fun_2)
        {break;
        }
      }
      if(current_new_lh < desired_min)
      {break;
      }

      OLD_RUN_THETA = clone(THETA);


      if(num_total_iter > total_iter)
      {break;
      }

    }  // ....................... End of run loop

    Rcout << "********** TASK COMPLETE ****************" << std::endl;
    //Rcout << "temp_sum :" << temp_sum << std::endl;
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
