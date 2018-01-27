#include "nonLinearPars.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

//Abstract class: BaseNlpar
//===============================================================================================================
Nlpar::Nlpar(std::string a,int b,int c,double d,double e,double f,double g){
  nam = a;
  fix = b;
  per = c;
  val = d;
  err = e;
  min = f;
  max = g;
  ran = max - min;
}

std::string Nlpar::getName(){
  return this->nam;
}
double Nlpar::getValue(){
  return this->val;
}
int Nlpar::getActive(){
  if( this->fix ){
    return 0;
  } else {
    return 1;
  }
}
void Nlpar::setNewPrior(BasePrior* prior){
  this->pri = prior;
}

std::vector<std::string> Nlpar::getVectorNames(std::vector<Nlpar*> pars){
  std::vector<std::string> names;
  for(int i=0;i<pars.size();i++){
    names.push_back( pars[i]->getName() );
  }
  return names;
}

std::vector<double> Nlpar::getVectorValues(std::vector<Nlpar*> pars){
  std::vector<double> values;
  for(int i=0;i<pars.size();i++){
    values.push_back( pars[i]->getValue() );
  }
  return values;
}

std::vector<int> Nlpar::getVectorActive(std::vector<Nlpar*> pars){
  std::vector<int> active;
  for(int i=0;i<pars.size();i++){
    if( pars[i]->getActive() ){
      active.push_back(i);
    }
  }
  return active;
}

double Nlpar::getValueByName(std::string nam,std::vector<Nlpar*> pars){
  for(int i=0;i<pars.size();i++){
    if( nam == pars[i]->getName() ){
      return pars[i]->getValue();
    }
  }
}
//===============================================================================================================


//Derived class from BasePrior: Uni
//===============================================================================================================
//virtual
Uni::Uni(Nlpar* p){
  type = "uni";
  mother = p;
  p0  = 1.0/mother->ran;
}

//virtual
double Uni::prior(double x){
  return this->p0;
}

//virtual
double Uni::fromUnitCube(double u){
  return this->mother->ran*u + this->mother->min;
}

//virtual
void Uni::printPars(){
  std::cout << "Uniform prior in the range: " << this->mother->min << "  -  " << this->mother->max <<std::endl;
}

//virtual
std::map<std::string,double> Uni::getPars(){
  std::map<std::string,double> mymap;
  return mymap;
}


//Derived class from BasePrior: Gauss
//===============================================================================================================
//virtual
Gauss::Gauss(Nlpar* p,double a,double b){
  type = "gauss";
  mother = p;
  mean = a;
  sdev = b;

  const double pi = 3.14159265358979323846;
  fac       = 1.0/(sdev*sqrt(2*pi));
  two_sdev2 = 2*pow(sdev,2);
  den       = this->F( (this->mother->max-mean)/sdev ) - this->F( (this->mother->min-mean)/sdev );
}

//virtual
double Gauss::prior(double x){
  double num = this->fac*exp( -pow(x-this->mean,2)/this->two_sdev2 );
  double p = num/this->den;
  return p;
}

//virtual
double Gauss::fromUnitCube(double u){
  return this->mother->ran*u + this->mother->min;
}

//virtual
void Gauss::printPars(){
  std::cout << "Gaussian prior in the range: " << this->mother->min << "  -  " << this->mother->max <<std::endl;
  std::cout << " with mean: " << this->mean << std::endl;
  std::cout << "  and sdev: " << this->sdev << std::endl;
}

//virtual
std::map<std::string,double> Gauss::getPars(){
  std::map<std::string,double> mymap;
  mymap["mean"] = this->mean;
  mymap["sdev"] = this->sdev;
  return mymap;
}

//private
double Gauss::F(double x){
  const double pi = 3.14159265358979323846;
  double signx = 0;
  if( x>=0 ){ signx =  1; }
  if( x<0  ){ signx = -1; }
  
  double f = 1.0 + signx*sqrt( 1.0-exp(-2*x*x/pi) );
  f /= 2.0;
  return f;
}

//Derived class from BasePrior: Exp
//===============================================================================================================
//virtual
Exp::Exp(Nlpar* p,double a){
  type = "exp";
  mother = p;
  beta = a;
  fac = beta * (exp(-this->mother->min/beta) - exp(-this->mother->max/beta));
}

//virtual
double Exp::prior(double x){
  double p = this->fac*exp(-x/this->beta);
  return p;
}

//virtual
double Exp::fromUnitCube(double u){
  return this->mother->ran*u + this->mother->min;
}

//virtual
void Exp::printPars(){
  std::cout << "Exponential prior exp(-x/b) in the range: " << this->mother->min << "  -  " << this->mother->max <<std::endl;
  std::cout << " with b: " << this->beta << std::endl;
}

//virtual
std::map<std::string,double> Exp::getPars(){
  std::map<std::string,double> mymap;
  mymap["beta"] = this->beta;
  return mymap;
}


//Derived class from BasePrior: Log
//===============================================================================================================
//virtual
Log::Log(Nlpar* p){
  type = "log";
  mother = p;
  fac = 1.0/log(this->mother->max/this->mother->min);
}

//virtual
double Log::prior(double x){
  double p = this->fac/x;
  return p;
}

//virtual
double Log::fromUnitCube(double u){
  return this->mother->ran*u + this->mother->min;
}

//virtual
void Log::printPars(){
  std::cout << "Logarithmic prior 1/x in the range: " << this->mother->min << "  -  " << this->mother->max <<std::endl;
}

//virtual
std::map<std::string,double> Log::getPars(){
  std::map<std::string,double> mymap;
  return mymap;
}
