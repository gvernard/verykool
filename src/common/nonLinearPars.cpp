#include "nonLinearPars.hpp"

#include <fstream>
#include <cmath>

//Abstract class: BaseNlpar
//===============================================================================================================
std::vector<std::string> BaseNlpar::getActive(std::map<std::string,BaseNlpar*> nlpars){
  std::vector<std::string> active;
  typedef std::map<std::string,BaseNlpar*>::iterator it_type;
  for(it_type iterator=nlpars.begin();iterator!=nlpars.end();iterator++){
    if( iterator->second->fix == 0 ){
      active.push_back(iterator->first);
    }
  }
  return active;
}

void getPriorPars(){
}

void setPriorToUni(){
}
void setPriorToGauss(){
}
void setPriorToLog(){
}
//===============================================================================================================


//Derived class from BaseNlpar: Uni
//===============================================================================================================

//virtual
Uni::Uni(std::string a,int b,int c,double d,double e,double f,double g,std::map<std::string,double> h){
  nam = a;
  fix = b;
  per = c;
  val = d;
  err = e;
  min = f;
  max = g;
  ran = max-min;
  pri = h;
  
  p0  = 1./ran;
}

//virtual
Uni::~Uni(){}

//virtual
double Uni::prior(double x){
  return this->p0;
}

//virtual
void Uni::fromUnitCube(double u){
  this->val = this->ran*u + this->min;
}



//Derived class from BaseNlpar: Gauss
//===============================================================================================================

//virtual
Gauss::Gauss(std::string a,int b,int c,double d,double e,double f,double g,std::map<std::string,double> h){
  nam = a;
  fix = b;
  per = c;
  val = d;
  err = e;
  min = f;
  max = g;
  ran = max-min;
  pri = h;
    
  const double pi = 3.14159265358979323846;
  fac       = 1./(pri["sdev"]*sqrt(2*pi));
  two_sdev2 = 2*pow(pri["sdev"],2);
  den       = this->F( (max-pri["mean"])/pri["sdev"] ) - this->F( (min-pri["mean"])/pri["sdev"] );
}

//virtual
Gauss::~Gauss(){}

//virtual
double Gauss::prior(double x){
  double num = this->fac*exp( -pow(x-this->pri["mean"],2)/this->two_sdev2 );
  double p = num/this->den;
  return p;
}

//virtual
void Gauss::fromUnitCube(double u){
  this->val = this->ran*u + this->min;
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

//Derived class from BaseNlpar: Log
//===============================================================================================================

//virtual
Log::Log(std::string a,int b,int c,double d,double e,double f,double g,std::map<std::string,double> h){
  nam = a;
  fix = b;
  per = c;
  val = d;
  err = e;
  min = f; // logarithm
  max = g; // logarithm
  ran = max-min; // in logarithmic space
  pri = h;

  p0  = 1./ran;
  fac = 1.0/log(10);
}

//virtual
Log::~Log(){}

//virtual
double Log::prior(double x){
  double p = this->fac*this->p0/x;
  return p;
}

//virtual
void Log::fromUnitCube(double u){
  double exp = this->ran*u + this->min;
  this->val = pow(10,exp);
}


