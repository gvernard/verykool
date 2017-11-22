#include "cov_kernels.hpp"

#include <cmath>

GaussKernel::GaussKernel(double rmax,double sdev){
  this->type = "gauss";
  this->rmax = rmax;
  this->sdev = sdev;
}
double GaussKernel::getCovariance(double r){
  double cov = this->fac*exp(-r*r/(2*this->sdev*this->sdev))/this->sdev;
  return cov;
}
void GaussKernel::setParameters(std::map<std::string,BaseNlpar*> pars){
  this->rmax = pars["rmax"]->val;
  this->sdev = pars["sdev"]->val;
}

PowerLawKernel::PowerLawKernel(double rmax,double ampl,double slope){
  this->type  = "power_law";
  this->rmax  = rmax;
  this->ampl  = ampl;
  this->slope = slope;
}
double PowerLawKernel::getCovariance(double r){
  double cov = this->ampl*pow(r,this->slope);
  return cov;
}
void PowerLawKernel::setParameters(std::map<std::string,BaseNlpar*> pars){
  this->rmax  = pars["rmax"]->val;
  this->ampl  = pars["ampl"]->val;
  this->slope = pars["slope"]->val;
}

