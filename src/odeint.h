//
//  odeint.h
//  calc_ll_secsse
//
//  Created by Thijs Janzen on 08/03/2021.
//  Copyright © 2021 Thijs Janzen. All rights reserved.
//

#ifndef odeint_h
#define odeint_h

#include <iostream>

// [[Rcpp::depends(BH)]]
#include "Rcpp.h"
#include "boost/numeric/odeint.hpp"

#include "util.h"

class ode_standard {
public:
  
  ode_standard(const std::vector<double>& l,
                      const std::vector<double>& m,
                      const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q) {
    d = l.size();
  }
  
  ode_standard(const Rcpp::NumericVector& l,
               const Rcpp::NumericVector& m,
               const Rcpp::NumericMatrix& q) {
    
    l_ = std::vector<double>(l.begin(), l.end());
    m_ = std::vector<double>(m.begin(), m.end());
    numericmatrix_to_vector(q, q_);
    
    d = l_.size();
  }
  
  void operator()( const std::vector< double > &x ,
                std::vector<  double > &dxdt,
                const double /* t */ ) {
    for (size_t i = 0; i < d; ++i) {
      
      if (l_[i] != 0.0) {
      
        dxdt[i] = m_[i] - (l_[i] + m_[i]) * x [i]  +
                  l_[i] * x[i] * x[i];
        
        long double FF3 = -1.0 * l_[i] - m_[i] + 2 * l_[i] * x[i];
        dxdt[i + d] = FF3 * x[ i + d];
      } else {
        dxdt[i] = - 1.0 * m_[i] * x [i] + m_[i];
        
        dxdt[i + d] = -1.0 * m_[i] * x[ i + d];
      }
      
      for (size_t j = 0; j < d; ++j) {
        long double diff_e = x[j] - x[i];
        dxdt[i] += diff_e * q_[i][j];
        
        long double diff_d = x[j + d] - x[i + d];
        dxdt[i + d] += diff_d * q_[i][j];
      }
    }
    return;
  }
  
  
  double get_l(int index) const {
    return l_[index];
  } 
  
  size_t get_d() const {
    return d;
  }
  
private:
  std::vector< double > l_;
  std::vector< double > m_;
  std::vector< std::vector< double >> q_;
  size_t d;
};

class ode_standard_ct {
public:
  
  ode_standard_ct(const std::vector<double>& l,
               const std::vector<double>& m,
               const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q) {
    d = l.size();
  }
  
  ode_standard_ct(const Rcpp::NumericVector& l,
                  const Rcpp::NumericVector& m,
                  const Rcpp::NumericMatrix& q) {
    
    l_ = std::vector<double>(l.begin(), l.end());
    m_ = std::vector<double>(m.begin(), m.end());
    numericmatrix_to_vector(q, q_);
    
    d = l_.size();
  }
  
  void operator()( const std::vector< double > &x ,
                std::vector<  double > &dxdt,
                const double /* t */ ) {
    
    for (int i = 0; i < d; ++i) {
      long double diff_1 = (m_[i] - (l_[i] * x[i]));
      dxdt[i] =  diff_1 * (1 - x[i]);
      
      dxdt[i + d] = -1.0 * (l_[i] + m_[i]) * x[i + d];
    }
      
    for (int j = 0; j < d; ++j) {
      for (int k = 0; k < d; ++k) {
        long double diff_e = x[k] - x[j];
        dxdt[j] +=  q_[j][k] * diff_e;
        
        long double diff_d = x[k + d] - x[j + d];
        dxdt[j + d] += q_[j][k] * diff_d;
      }
    }
    
    
    return;
  }
  
  
  double get_l(int index) const {
    return l_[index];
  } 
  
  size_t get_d() const {
    return d;
  }
  
private:
  std::vector< double > l_;
  std::vector< double > m_;
  std::vector< std::vector< double >> q_;
  size_t d;
};



class ode_cla {
public:
  
  ode_cla(const std::vector<std::vector<std::vector<double>>>& l,
          const std::vector<double>& m,
          const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q), d(m.size()) {
    
    lambda_sum = std::vector<long double>(d, 0.0);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          lambda_sum[i] += l_[i][j][k];
        }
      }
    }
  }
  
  // rename to operator() to use this version of the operator, withouth kahan sum.
  void operator()(const std::vector< double > &x ,
                std::vector< double > &dxdt,
                const double /* t */ ) const {
    
    for (int i = 0; i < d; ++i) {
      double Df = 0.0;
      double Ef = 0.0;
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          if (l_[i][j][k] != 0.0) { // slightly safer.
            
            long double ff1 = (x[j] * x[k + d] + x[j + d] * x[k]);
            long double ff2 = (x[j] * x[k]);
            
            Df += l_[i][j][k] * ff1;
            Ef += l_[i][j][k] * ff2;
          }
        }
      }

    //  dxdt[i]     = Ef + m_[i] + x[i] * (m_[i] - lambda_sum[i]);
    //  dxdt[i + d] = Df - x[i + d] * (lambda_sum[i] + m_[i]); //remaining: Q_ds - x[i+d] * Q_one
      
      dxdt[i]     = Ef + m_[i] - (lambda_sum[i] + m_[i]) * x[i];
      dxdt[i + d] = Df + (-lambda_sum[i] - m_[i]) * x[i + d];
      
      for (size_t j = 0; j < d; ++j) {
        // q_[i][j] is always non-zero.
        long double temp1 = (x[j]     - x[i]);
        dxdt[i]     += q_[i][j] * temp1;
        
        long double temp2 = (x[j + d] - x[i + d]);
        dxdt[i + d] += q_[i][j] * temp2;
      }
    }
    return;
  }
  
  double get_l(size_t i, size_t j, size_t k) const {
    return l_[i][j][k];
  } 

  size_t get_d() const {
    return d;
  }
  
private:
  const std::vector< std::vector< std::vector< double > > > l_;
  const std::vector< double > m_;
  const std::vector< std::vector< double >> q_;
  const size_t d;
  std::vector< long double > lambda_sum;
};

class ode_cla_backup {
public:
  
  ode_cla_backup(const std::vector<std::vector<std::vector<double>>>& l,
          const std::vector<double>& m,
          const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q), d(m.size()) {
    
    lambda_sum = std::vector<long double>(d, 0.0);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          lambda_sum[i] += l_[i][j][k];
        }
      }
    }
  }
  
  void operator()(const std::vector< double > &x ,
                std::vector< double > &dxdt,
                const double /* t */ ) const {
    
    for (int i = 0; i < d; ++i) {
      
      long double lamEE = 0.0;
      long double lamDE = 0.0;
      
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          
          if (l_[i][j][k] != 0) {
            long double FF1 = x[j] * x[k];
            lamEE += l_[i][j][k] * FF1;
            
            long double FF3 = x[ d + j ] * x[k];
            long double FF2 = x[d + k] * x[k];
            lamDE += l_[i][j][k] * (FF3 + FF2);  
          }
        }
      }
    
      long double FF1 = m_[i] - (lambda_sum[i] + m_[i]) * x[i];
      dxdt[i] = FF1 + lamEE;
      
      long double FF2 = (-lambda_sum[i] - m_[i]) * x[i + d];
      dxdt[i + d] = FF2 + lamDE;
      
      for (int j = 0; j < d; ++j) {
        long double diff = x[j] - x[i];
        dxdt[i] += diff * q_[i][j];
        
        long double diff2 = x[j + d] - x[i + d];
        dxdt[i + d] += diff2 * q_[i][j];
      }
    }

    return;
  }
  
  double get_l(size_t i, size_t j, size_t k) const {
    return l_[i][j][k];
  } 
  
  size_t get_d() const {
    return d;
  }
  
private:
  const std::vector< std::vector< std::vector< double > > > l_;
  const std::vector< double > m_;
  const std::vector< std::vector< double >> q_;
  const size_t d;
  std::vector< long double > lambda_sum;
};

class ode_cla_d {
public:
  
  ode_cla_d(const std::vector<std::vector<std::vector<double>>>& l,
          const std::vector<double>& m,
          const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q), d(m.size()) { 
   
    lambda_sum = std::vector<long double>(d, 0.0);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < d; ++k) {
          lambda_sum[i] += l_[i][j][k];
        }
      }
    }
  }
  

  void operator()(const std::vector< double > &x ,
                std::vector< double > &dxdt,
                const double /* t */ ) const {
    
    for (int i = 0; i < d; ++i) {
      dxdt[i + d] = -1.0 * (lambda_sum[i] + m_[i]) * x[i + d];
      
      for (int j = 0; j < d; ++j) {
        long double dx = x[j + d] - x[i + d];
        dxdt[i + d] +=  q_[i][j] * dx;
      }
    }
  }
  
  double get_l(size_t i, size_t j, size_t k) const {
    return l_[i][j][k];
  } 
  
  size_t get_d() const {
    return d;
  }
  
private:
  const std::vector< std::vector< std::vector< double > > > l_;
  const std::vector< double > m_;
  const std::vector< std::vector< double >> q_;
  const size_t d;
  std::vector<long double> lambda_sum;
};

class ode_cla_e {
public:
  
  ode_cla_e(const std::vector<std::vector<std::vector<double>>>& l,
            const std::vector<double>& m,
            const std::vector<std::vector<double>>& q) :
  l_(l), m_(m), q_(q), d(m.size()) { 
  }
  
  void operator()(const std::vector< double > &x ,
                std::vector< double > &dxdt,
                const double /* t */ ) const {
    
    
    for (int i = 0; i < d; ++i) {
      dxdt[i] = 0.0;
      if (m_[i] != 0.0) {
        dxdt[i] = m_[i] * (1.0 - x[i]);
      }
      
      for (int j = 0; j < d; ++j) {
        long double diff = (x[j] - x[i]);
        dxdt[i] += q_[i][j] * diff;
        
        for (int k = 0; k < d; ++k) {
          if (l_[i][j][k] != 0.0) {
            long double diff2 = (x[j] * x[k] - x[i]);
            dxdt[i] += l_[i][j][k] * diff2;
          }
        }
      }
    }
  }
  
  double get_l(size_t i, size_t j, size_t k) const {
    return l_[i][j][k];
  } 
  
  size_t get_d() const {
    return d;
  }
  
private:
  const std::vector< std::vector< std::vector< double > > > l_;
  const std::vector< double > m_;
  const std::vector< std::vector< double >> q_;
  const size_t d;
};


namespace odeintcpp {

namespace bno = boost::numeric::odeint;

constexpr char default_stepper_name[] = "odeint::runge_kutta_fehlberg78";
constexpr double default_atol = 1E-10;
constexpr double default_rtol = 1E-10;
constexpr double default_init_steps = 10;

template <
  typename STEPPER,
  typename STATE,
  typename ODE
>
void integrate(STEPPER stepper, ODE ode, STATE& y, double t0, double t1, double dt)
{
  bno::integrate_adaptive(stepper, ode, y, t0, t1, dt);
}


template <
  typename STATE,
  typename ODE
>
void integrate(const std::string& stepper_name, std::unique_ptr<ODE> ode, STATE& y, double t0, double t1, double dt, double atol, double rtol)
{
  if ("odeint::runge_kutta_cash_karp54" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_cash_karp54<STATE>>(atol, rtol),
              std::ref(*ode), y, t0, t1, dt);
  }
  else if ("odeint::runge_kutta_fehlberg78" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_fehlberg78<STATE>>(atol, rtol), std::ref(*ode), y, t0, t1, dt);
  }
  else if ("odeint::runge_kutta_dopri5" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_dopri5<STATE>>(atol, rtol), std::ref(*ode), y, t0, t1, dt);
  }
  else if ("odeint::bulirsch_stoer" == stepper_name) {
    integrate(bno::bulirsch_stoer<STATE>(atol, rtol), std::ref(*ode), y, t0, t1, dt);
  }
  else if ("odeint::runge_kutta4" == stepper_name) {
    integrate(bno::runge_kutta4<STATE>(), std::ref(*ode), y, t0, t1, dt);
  }
 
  else {
    throw std::runtime_error("odeintcpp::integrate: unknown stepper");
  }
}


template <
  typename STATE,
  typename ODE
>
void integrate(const std::string& stepper_name, std::unique_ptr<ODE> ode, STATE& y, double t0, double t1)
{
  integrate(stepper_name, std::move(ode), y, t0, t1, (t1 - t0) / default_init_steps, default_atol, default_rtol);
}


template <
  typename STATE,
  typename ODE
>
void integrate(std::unique_ptr<ODE> ode, STATE& y, double t0, double t1, double dt, double atol, double rtol)
{
  integrate(default_stepper_name, std::move(ode), y, t0, t1, dt, atol, rtol);
}


template <
  typename STATE,
  typename ODE
>
void integrate(std::unique_ptr<ODE> ode, STATE& y, double t0, double t1)
{
  integrate(default_stepper_name, std::move(ode), y, t0, t1);
}
}

#endif /* odeint_h */