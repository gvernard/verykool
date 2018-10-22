#ifndef COVARIANCE_KERNELS_HPP
#define COVARIANCE_KERNELS_HPP

#include <string>
#include <vector>
#include <iostream>

class Nlpar;

class BaseCovKernel {
public:
  std::string type;
  double cmax;
  
  BaseCovKernel(){};
  ~BaseCovKernel(){};
  
  virtual double getCovariance(double r) = 0;
  virtual double getCovarianceSelf() = 0;
  virtual void setParameters(std::vector<Nlpar*> pars) = 0;
};


class GaussKernel: public BaseCovKernel {
public:
  double sdev;

  GaussKernel(std::vector<Nlpar*> pars);
  double getCovariance(double r);
  double getCovarianceSelf();
  void setParameters(std::vector<Nlpar*> pars);
private:
 double fac;
};

class ModGaussKernel: public BaseCovKernel {
public:
  double sdev;

  ModGaussKernel(std::vector<Nlpar*> pars);
  double getCovariance(double r);
  double getCovarianceSelf();
  void setParameters(std::vector<Nlpar*> pars);
};

class ExpGaussKernel: public BaseCovKernel {
public:
  double expo;
  double sdev;

  ExpGaussKernel(std::vector<Nlpar*> pars);
  double getCovariance(double r);
  double getCovarianceSelf();
  void setParameters(std::vector<Nlpar*> pars);
private:
  double fac;
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
    if( kernel_type == "modgauss" ){
      return new ModGaussKernel(pars);
    } else if( kernel_type == "gauss" ){
      return new GaussKernel(pars);
    } else if( kernel_type == "expgauss" ){
      return new ExpGaussKernel(pars);
    } else {
      return NULL;
    }
  }

private:
  FactoryCovKernel(){};
};


#endif /* COVARIANCE_KERNELS_HPP */
