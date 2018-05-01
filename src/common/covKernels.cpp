#include "covKernels.hpp"

#include <cmath>

#include "nonLinearPars.hpp"


ModGaussKernel::ModGaussKernel(std::vector<Nlpar*> pars){
  this->type = "modgauss";
  this->setParameters(pars);
}
double ModGaussKernel::getCovariance(double r){
  //  double cov = this->ampl*exp(-r*r/(2*this->sdev*this->sdev));
  double cov = this->ampl*exp(-r/this->sdev);
  return cov;
}
double ModGaussKernel::getCovarianceSelf(){
  double cov = this->ampl;
  return cov;
}
void ModGaussKernel::setParameters(std::vector<Nlpar*> pars){
  this->rmax = Nlpar::getValueByName("rmax",pars);
  this->ampl = Nlpar::getValueByName("ampl",pars);
  this->sdev = Nlpar::getValueByName("sdev",pars);
}


GaussKernel::GaussKernel(std::vector<Nlpar*> pars){
  this->type = "gauss";
  this->setParameters(pars);
}
double GaussKernel::getCovariance(double r){
  double cov = this->fac*exp(-r*r/(2*this->sdev*this->sdev))/this->sdev;
  return cov;
}
double GaussKernel::getCovarianceSelf(){
  double cov = this->fac;
  return cov;
}
void GaussKernel::setParameters(std::vector<Nlpar*> pars){
  this->rmax = Nlpar::getValueByName("rmax",pars);
  this->sdev = Nlpar::getValueByName("sdev",pars);
}


PowerLawKernel::PowerLawKernel(std::vector<Nlpar*> pars){
  this->type  = "power_law";
  this->setParameters(pars);
}
double PowerLawKernel::getCovariance(double r){
  double cov = this->ampl*pow(r,this->slope);
  return cov;
}
double PowerLawKernel::getCovarianceSelf(){
  double cov = this->ampl;
  return cov;
}
void PowerLawKernel::setParameters(std::vector<Nlpar*> pars){
  this->rmax  = Nlpar::getValueByName("rmax",pars);
  this->ampl  = Nlpar::getValueByName("ampl",pars);
  this->slope = Nlpar::getValueByName("slope",pars);
}
