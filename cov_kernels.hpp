#ifndef COVARIANCE_KERNELS_HPP
#define COVARIANCE_KERNELS_HPP

#include <string>
#include <map>

#include "nonLinearPars.hpp"

class BaseCovKernel {
public:
  std::string type;
  double rmax;
  
  BaseCovKernel(){};
  ~BaseCovKernel(){};
  
  virtual double getCovariance(double r) = 0;
  virtual void setParameters(std::map<std::string,BaseNlpar*> pars) = 0;
};


class GaussKernel: public BaseCovKernel {
public:
  double sdev;

  GaussKernel(double rmax,double sdev);
  double getCovariance(double r);
  void setParameters(std::map<std::string,BaseNlpar*> pars);
private:
  double fac = 0.39894228; // 1/sqrt(2*pi)
};

class PowerLawKernel: public BaseCovKernel {
public:
  double ampl;
  double slope;

  PowerLawKernel(double rmax,double ampl,double slope);
  double getCovariance(double r);
  void setParameters(std::map<std::string,BaseNlpar*> pars);
};

class FactoryCovKernel {//This is a singleton class.
public:
  FactoryCovKernel(FactoryCovKernel const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryCovKernel const&) = delete;

  static FactoryCovKernel* getInstance(){
    static FactoryCovKernel dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseCovKernel* createCovKernel(std::string kernel_type,std::map<std::string,BaseNlpar*> pars){
    if( kernel_type == "gauss" ){
      return new GaussKernel(pars["rmax"]->val,pars["sdev"]->val);
    } else if( kernel_type == "power_law" ){
      return new PowerLawKernel(pars["rmax"]->val,pars["ampl"]->val,pars["slope"]->val);
    } else {
      return NULL;
    }
  }

private:
  FactoryCovKernel(){};
};


#endif /* COVARIANCE_KERNELS_HPP */
