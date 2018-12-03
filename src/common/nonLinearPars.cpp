#include "nonLinearPars.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

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

std::map<std::string,double> Nlpar::getSigmaIntervals(const std::vector<double>& x,const std::vector<double>& p,int sigma_interval){
  double Xsigma;
  if( sigma_interval == 3 ){
    Xsigma = 0.9973;
  } else if( sigma_interval == 2 ){
    Xsigma = 0.9545;
  } else {
    Xsigma = 0.6827;
  }
  double low_limit  = 0.5 - Xsigma/2.0;
  double high_limit = 0.5 + Xsigma/2.0;

  int Nsamples = x.size();

  // sort the probability and one parameter simultaneously
  std::vector< std::pair<double,double> > combi;
  for(int i=0;i<Nsamples;i++){
    combi.push_back(std::make_pair(p[i],x[i]));
  }
  std::sort(combi.begin(),combi.end(),compareFunc);
  
  // reduce the sorted arrays by the parameter value
  std::vector<double> redu_prob;
  std::vector<double> redu_par;
  std::vector<double> redu_count;
  
  double v0 = combi[0].second;
  int i=1;
  while( i < Nsamples ){
    double v1 = combi[i].second;
    int count = 1;
    double tmp_prob = combi[i-1].first;
    while( v0 == v1 && i < Nsamples ){
      count++;
      tmp_prob += combi[i].first;
      i++;
      v1 = combi[i].second;
    }
    redu_par.push_back(v0);
    redu_count.push_back(count);
    redu_prob.push_back(tmp_prob);
    v0 = v1;
    i++;
  }
  std::vector< std::pair<double,double> >().swap(combi); // clear combi vector
  
  // convert from probability density to probability
  std::vector<double> pcum(redu_prob.size()-1);
  std::vector<double> xbin(redu_prob.size()-1);
  for(int i=1;i<redu_prob.size();i++){
    xbin[i-1] = (redu_par[i] + redu_par[i-1])/2.0;
    pcum[i-1] = (redu_prob[i] + redu_prob[i-1])*(redu_par[i] - redu_par[i-1])/2.0;
  }
  std::vector< double>().swap(redu_par);   // clear vector
  std::vector< double>().swap(redu_count); // clear vector
  std::vector< double>().swap(redu_prob);  // clear vector

  // convert to cumulative probability
  for(int i=1;i<pcum.size();i++){
    pcum[i] = pcum[i-1] + pcum[i];
  }

  // normalize
  for(int i=0;i<pcum.size();i++){
    pcum[i] /= pcum[pcum.size()-1];
  }
  
  // calculate mean value and limits
  std::map<std::string,double> output;
  output["mean"] = Nlpar::calculateInterval(xbin,pcum,0.5);
  output["low"]  = 0.0;
  output["high"] = 0.0;
  if( pcum[0] > low_limit ){ // only an upper limit can be computed
    output["high"] = Nlpar::calculateInterval(xbin,pcum,Xsigma);
  } else if ( pcum[pcum.size()-1] < high_limit ){ // only a lower limit can be computed
    output["low"]  = Nlpar::calculateInterval(xbin,pcum,1.0-Xsigma);
  } else { // both limits can be computed
    output["low"]  = Nlpar::calculateInterval(xbin,pcum,low_limit);
    output["high"] = Nlpar::calculateInterval(xbin,pcum,high_limit);    
  }

  return output;
}

bool Nlpar::compareFunc(const std::pair<double,double> &a,const std::pair<double,double> &b){
  return a.second < b.second;
}

double Nlpar::calculateInterval(const std::vector<double>& x,const std::vector<double>& p,double limit){
  double x0 = 0.0;
  for(int i=0;i<x.size();i++){
    if( p[i] >= limit ){
      x0 = x[i-1] + (x[i] - x[i-1])*(limit - p[i-1])/(p[i] - p[i-1]); // linear interpolation
      break;
    }
  }
  return x0;
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

Nlpar* Nlpar::getParByName(std::string name,std::vector<Nlpar*> pars){
  for(int i=0;i<pars.size();i++){
    if( pars[i]->nam == name ){
      return pars[i];
    }
  }
  return NULL;
}

bool Nlpar::getSampleReg(std::vector<Nlpar*> pars){
  for(int i=0;i<pars.size();i++){
    if( pars[i]->fix == 0 && pars[i]->nam.substr(0,6) != "lambda" ){
      return true;
    }
  }
  return false;
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
