#include "covKernels.hpp"

#include <cmath>

#include "nonLinearPars.hpp"



ModGaussKernel::ModGaussKernel(std::vector<Nlpar*> pars){
  this->type = "modgauss";
  this->setParameters(pars);
}
ModGaussKernel::ModGaussKernel(const ModGaussKernel& other){
  this->sdev = other.sdev;
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
  this->cmax = Nlpar::getValueByName("cmax",pars);
  this->sdev = Nlpar::getValueByName("sdev",pars);
}



GaussKernel::GaussKernel(std::vector<Nlpar*> pars){
  this->type = "gauss";
  this->setParameters(pars);
}
GaussKernel::GaussKernel(const GaussKernel& other){
  this->sdev = other.sdev;
}
double GaussKernel::getCovariance(double r){
  //  double cov = this->fac*exp(-r*r/(2*this->sdev*this->sdev));
  double cov = exp(-r*r/(2*this->sdev*this->sdev));
  return cov;
}
double GaussKernel::getCovarianceSelf(){
  //  double cov = this->fac;
  double cov = 1.1;
  return cov;
}
void GaussKernel::setParameters(std::vector<Nlpar*> pars){
  this->cmax = Nlpar::getValueByName("cmax",pars);
  this->sdev = Nlpar::getValueByName("sdev",pars);
  //  this->fac  = 1.0/(this->sdev*sqrt(2*M_PI));
}



ExpGaussKernel::ExpGaussKernel(std::vector<Nlpar*> pars){
  this->type  = "expgauss";
  this->setParameters(pars);
}
ExpGaussKernel::ExpGaussKernel(const ExpGaussKernel& other){
  this->sdev = other.sdev;
  this->expo = other.expo;
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
  this->cmax = Nlpar::getValueByName("cmax",pars);
  this->expo = Nlpar::getValueByName("expo",pars);
  this->sdev = Nlpar::getValueByName("sdev",pars);
  this->fac  = 1.0/(this->sdev*sqrt(2*M_PI));
}
