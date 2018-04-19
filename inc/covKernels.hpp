#ifndef COVARIANCE_KERNELS_HPP
#define COVARIANCE_KERNELS_HPP

#include <string>
#include <vector>
#include <iostream>

class Nlpar;

class BaseCovKernel {
public:
  std::string type;
  double rmax;
  
  BaseCovKernel(){};
  ~BaseCovKernel(){};
  
  virtual double getCovariance(double r) = 0;
  virtual void setParameters(std::vector<Nlpar*> pars) = 0;
};


class GaussKernel: public BaseCovKernel {
public:
  double sdev;

  GaussKernel(std::vector<Nlpar*> pars);
  double getCovariance(double r);
  void setParameters(std::vector<Nlpar*> pars);
private:
  double fac = 0.39894228; // 1/sqrt(2*pi)
};

class PowerLawKernel: public BaseCovKernel {
public:
  double ampl;
  double slope;

  PowerLawKernel(std::vector<Nlpar*> pars);
  double getCovariance(double r);
  void setParameters(std::vector<Nlpar*> pars);
};

class FactoryCovKernel {//This is a singleton class.
public:
  FactoryCovKernel(FactoryCovKernel const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryCovKernel const&) = delete;

  static FactoryCovKernel* getInstance(){
    static FactoryCovKernel dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseCovKernel* createCovKernel(const std::string kernel_type,std::vector<Nlpar*> pars){
    if( kernel_type == "gauss" ){
      return new GaussKernel(pars);
    } else if( kernel_type == "power_law" ){
      return new PowerLawKernel(pars);
    } else {
      return NULL;
    }
  }

private:
  FactoryCovKernel(){};
};


#endif /* COVARIANCE_KERNELS_HPP */
