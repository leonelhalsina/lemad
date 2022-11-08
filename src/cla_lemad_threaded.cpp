#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <tuple>

#include <thread>
#include <chrono>

#include "odeint.h"
#include "util.h"

#include <cmath>

#include <RcppParallel.h>


#include "threaded_ll.h"

template< typename OD_TYPE>
struct combine_states_cla {
  
  combine_states_cla(int d, const OD_TYPE& od) : d_(d), od_(od) {}
  
  state_vec operator()(const std::tuple< state_vec, state_vec >& input_states) {
    
    state_vec nodeN =  std::get<0>(input_states);
    state_vec nodeM =  std::get<1>(input_states);
    
    double ll1 = nodeN.back(); nodeN.pop_back();
    double ll2 = nodeM.back(); nodeM.pop_back();
    
    state_vec mergeBranch = std::vector<double>(d_, 0.0);
    
    for (size_t i = 0; i < d_; ++i) {
      for (size_t j = 0; j < d_; ++j) {
        for (size_t k = 0; k < d_; ++k) {
          
          double a = od_.get_l(i, j, k);
          
          if (a != 0.0) {
            double mult = (nodeN[j + d_] * nodeM[k + d_] +
                           nodeM[j + d_] * nodeN[k + d_]);
            
            mergeBranch[i] += a * mult;
          }
        }
      }
      mergeBranch[i] *= 0.5;
    }
    
    long double loglik = ll1 + ll2;
    
    normalize_loglik(mergeBranch, loglik);

    state_vec newstate(d_);
    for (int i = 0; i < d_; ++i) {
      newstate[i] = nodeM[i];
    }
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());
    newstate.push_back(loglik);
    
    return newstate;
  }
  
  size_t d_;
  OD_TYPE od_;
};

// [[Rcpp::export]]
Rcpp::List calc_cla_ll_threaded(const Rcpp::NumericVector& ances,
                                const Rcpp::NumericMatrix& states_R,
                                const Rcpp::NumericMatrix& forTime_R,
                                const Rcpp::List& lambdas_R,
                                const Rcpp::NumericVector& mus_R,
                                const Rcpp::NumericMatrix& Q,
                                int num_threads = 1,
                                std::string method = "odeint::bulirsch_stoer",
                                bool is_complete_tree = false) {
  try {
    
    std::vector< std::vector< double >> states_cpp, for_time_cpp, Q_cpp;
    numericmatrix_to_vector(states_R, states_cpp);
    numericmatrix_to_vector(forTime_R, for_time_cpp);
    numericmatrix_to_vector(Q, Q_cpp);
    
    std::vector< int > ances_cpp(ances.begin(), ances.end());
    
    std::vector<double> mus_cpp(mus_R.begin(), mus_R.end());
    
    std::vector< std::vector< std::vector< double > >> ll_cpp;
    for (size_t i = 0; i < lambdas_R.size(); ++i) {
      Rcpp::NumericMatrix temp = lambdas_R[i];
      std::vector< std::vector< double >> temp2;
      for (size_t j = 0; j < temp.nrow(); ++j) {
        std::vector<double> row;
        for (size_t k = 0; k < temp.ncol(); ++k) {
          row.push_back(temp(j, k));  
        }
        temp2.push_back(row);
      }
      ll_cpp.push_back(temp2);
    }
    
    if (is_complete_tree) {
      ode_cla_d od_(ll_cpp, mus_cpp, Q_cpp);
      
      threaded_ll<ode_cla_d, combine_states_cla<ode_cla_d>  > ll_calc(od_, ances_cpp, 
                                                       for_time_cpp, states_cpp, 
                                                       num_threads, method);
      return ll_calc.calc_ll();
    } else {
      
      ode_cla od_(ll_cpp, mus_cpp, Q_cpp);
      
      threaded_ll<ode_cla, combine_states_cla<ode_cla> > ll_calc(od_, ances_cpp, 
                                                                 for_time_cpp, states_cpp, 
                                                                 num_threads, method);
      return ll_calc.calc_ll();
    } 
    
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}