#include "covKernels.hpp"

#include <cmath>

#include "nonLinearPars.hpp"



ModGaussKernel::ModGaussKernel(std::vector<Nlpar*> pars){
  this->type = "modgauss";
  this->setParameters(pars);
}
double ModGaussKernel::getCovariance(double r){
  double cov = exp(-r/this->sdev);
  return cov;
}
double ModGaussKernel::getCovarianceSelf(){
  double cov = 1.0;
  return cov;
}
void ModGaussKernel::setParameters(std::vector<Nlpar*> pars){
  this->rmax = Nlpar::getValueByName("rmax",pars);
  this->sdev = Nlpar::getValueByName("sdev",pars);
}



GaussKernel::GaussKernel(std::vector<Nlpar*> pars){
  this->type = "gauss";
  this->setParameters(pars);
}
double GaussKernel::getCovariance(double r){
  double cov = this->fac*exp(-r*r/(2*this->sdev*this->sdev));
  return cov;
}
double GaussKernel::getCovarianceSelf(){
  double cov = this->fac;
  return cov;
}
void GaussKernel::setParameters(std::vector<Nlpar*> pars){
  this->rmax = Nlpar::getValueByName("rmax",pars);
  this->sdev = Nlpar::getValueByName("sdev",pars);
  this->fac  = 1.0/(this->sdev*sqrt(2*M_PI));
}



ExpGaussKernel::ExpGaussKernel(std::vector<Nlpar*> pars){
  this->type  = "expgauss";
  this->setParameters(pars);
}
double ExpGaussKernel::getCovariance(double r){
  double cov = this->fac*exp(-pow(r,this->expo)/(2*this->sdev*this->sdev));
  return cov;
}
double ExpGaussKernel::getCovarianceSelf(){
  double cov = this->fac;
  return cov;
}
void ExpGaussKernel::setParameters(std::vector<Nlpar*> pars){
  this->rmax = Nlpar::getValueByName("rmax",pars);
  this->expo = Nlpar::getValueByName("expo",pars);
  this->sdev = Nlpar::getValueByName("sdev",pars);
  this->fac  = 1.0/(this->sdev*sqrt(2*M_PI));
}
