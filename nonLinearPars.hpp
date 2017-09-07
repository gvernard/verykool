#ifndef NON_LINEAR_PARS_HPP
#define NON_LINEAR_PARS_HPP

#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include <vector>

//used in minimizers
struct myactive {
  int index;
  std::string nam;
  int per;
};

class BaseNlpar{
public:
  std::string nam;
  int fix;
  int per;
  double val;
  double err;
  double min;
  double max;
  double ran;
  std::map<std::string,double> pri;

  BaseNlpar(){};
  ~BaseNlpar(){
    this->pri.clear();
  };

  virtual double prior(double x) = 0;
  virtual void fromUnitCube(double u) = 0;
  //  virtual void getPriorPars() = 0; maybe no need to be virtual (loop over the 'pri' map)
  //  virtual void setPriorPars() = 0;

  static std::vector<std::string> getActive(std::map<std::string,BaseNlpar*> nlpars);

  /*
  static void printNlpars(std::map<std::string,BaseNlpar*> nlpars){
    typedef std::map<std::string,BaseNlpar*>::iterator it_type;
    for(it_type iterator=nlpars.begin();iterator!=nlpars.end();iterator++){
      printf("nam: %s, val: %f\n",iterator->second->nam,iterator->second->val);
    }
  };
  */
};


class Uni: public BaseNlpar{
public:
  Uni(std::string a,int b,int c,double d,double e,double f,double g,std::map<std::string,double> h);
  ~Uni();
  
  double prior(double x);
  void fromUnitCube(double u);

private:
  double p0; // the constant probability
};

class Gauss: public BaseNlpar{
public:
  Gauss(std::string a,int b,int c,double d,double e,double f,double g,std::map<std::string,double> h);
  ~Gauss();
  
  double prior(double x);
  void fromUnitCube(double u);

private:
  double two_sdev2; // 2 times the sdev squared
  double fac;       // 1/(s*sqrt(2*pi)) the normalization factor of a normal distribution
  double den;       // the denominator of a truncated normal distribution (difference between two cumulative distributions)

  double F(double x);
};

class Log: public BaseNlpar{
public:
  Log(std::string a,int b,int c,double d,double e,double f,double g,std::map<std::string,double> h);
  ~Log();

  double prior(double x);
  void fromUnitCube(double u);

private:
  double fac; // 1/ln(10)
  double p0;  // the constant probability in log space  
};


//This is a singleton class.
class FactoryNlpar{
public:
  FactoryNlpar(FactoryNlpar const&) = delete;//Stop the compiler generating methods of copying the object.
  void operator=(FactoryNlpar const&) = delete;

  static FactoryNlpar* getInstance(){
    static FactoryNlpar dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseNlpar* createNlpar(std::string nam,int fix,int per,double val,double err,double min,double max,std::string pri_nam,std::map<std::string,double> pri_par){
    if( pri_nam == "uni" ){ 
      return new Uni(nam,fix,per,val,err,min,max,pri_par);
    } else if ( pri_nam == "gauss" ){
      return new Gauss(nam,fix,per,val,err,min,max,pri_par);
    } else if ( pri_nam == "log" ){
      return new Log(nam,fix,per,val,err,min,max,pri_par);
    } else {
      return NULL;
    }
  }

  
private:
  FactoryNlpar(){};
};




#endif /* NON_LINEAR_PARS_HPP */
