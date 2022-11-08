#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <tuple>

#include <thread>
#include <chrono>

#include "odeint.h"
#include "util.h"

#include <RcppParallel.h>

#include "threaded_ll.h"

template<typename OD_TYPE>
struct combine_states {
  
  combine_states(int d, const OD_TYPE& od) : d_(d), od_(od) {}
  
  state_vec operator()(const std::tuple< state_vec, state_vec >& input_states) {
    state_vec vec1 =  std::get<0>(input_states);
    state_vec vec2 =  std::get<1>(input_states);
    
    double ll1 = vec1.back(); vec1.pop_back();
    double ll2 = vec2.back(); vec2.pop_back();
    
    state_vec mergeBranch(d_);
    for (int i = 0; i < d_; ++i) {
      mergeBranch[i] = vec1[i + d_] * vec2[i + d_] * od_.get_l(i);
    }
    
    long double loglik = ll1 + ll2;
    normalize_loglik(mergeBranch, loglik);
    
    state_vec newstate(d_);
    for (int i = 0; i < d_; ++i) {
      newstate[i] = vec2[i];
    }
    newstate.insert(newstate.end(), mergeBranch.begin(), mergeBranch.end());
    newstate.push_back(loglik);

    return newstate;
  }
  
  size_t d_;
  OD_TYPE od_;
};



// [[Rcpp::export]]
Rcpp::List calc_ll_threaded(const Rcpp::NumericVector& ll,
                            const Rcpp::NumericVector& mm,
                            const Rcpp::NumericMatrix& Q,
                            const Rcpp::NumericVector& ances,
                            const Rcpp::NumericMatrix& for_time,
                            const Rcpp::NumericMatrix& states,
                            int num_threads,
                            std::string method = "odeint::bulirsch_stoer",
                            bool is_complete_tree = false) {
  try {
    std::vector< int > ances_cpp(ances.begin(), ances.end());
    
    std::vector< std::vector< double >> for_time_cpp;
    numericmatrix_to_vector(for_time, for_time_cpp);
    
    std::vector< std::vector< double >> states_cpp;
    numericmatrix_to_vector(states, states_cpp);
    
    if (is_complete_tree) {
      ode_standard_ct od_(ll, mm, Q);
      
      threaded_ll<ode_standard_ct, combine_states<ode_standard_ct>> 
        ll_calc(od_, ances_cpp, 
                for_time_cpp, states_cpp, 
                num_threads, method);
      return ll_calc.calc_ll();
    
    } else {
      ode_standard od_(ll, mm, Q);
      threaded_ll<ode_standard, combine_states< ode_standard>> 
        ll_calc(od_, ances_cpp, 
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