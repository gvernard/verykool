#ifndef MINIMIZERS_HPP
#define MINIMIZERS_HPP

#include <string>
#include <vector>
#include <map>

//#include "simplex.hpp"


class ImagePlane;
class BaseSourcePlane;
class CollectionMassModels;
class BaseLikelihoodModel;


// BaseMinimizer parent/base class
//================================================================================================================================================
class BaseMinimizer {
public:
  std::string name;
  std::string type;
  BaseMinimizer(std::string name){
    this->name = name;
  };
  ~BaseMinimizer(){};

  virtual void minimize(std::map<std::string,std::string> minimizer,BaseLikelihoodModel* mypars,const std::string output) = 0;
  virtual void output(std::string output) = 0;
  virtual void finalizeMinimizer(std::string output) = 0;
};



// Nothing class (does nothing)
//================================================================================================================================================
class Nothing: public BaseMinimizer {
public:
  Nothing(std::string name) : BaseMinimizer(name) {
    this->type = "nothing";
  };
  ~Nothing(){};
  
  void minimize(std::map<std::string,std::string> minimizer,BaseLikelihoodModel* mypars,const std::string output);
  void output(std::string output);
  void finalizeMinimizer(std::string output){};
};


// MultiNest class
//================================================================================================================================================
class MultiNest: public BaseMinimizer {
public:
  MultiNest(std::string name) : BaseMinimizer(name) {
    this->type = "multinest";
  };
  ~MultiNest(){};
  int counter;
  int total_samples;
  int replacements;

  void minimize(std::map<std::string,std::string> opt,BaseLikelihoodModel* pars,const std::string output);
  void output(std::string output);
  void finalizeMinimizer(std::string output);
};

struct extras{
  BaseLikelihoodModel* pars;
  std::string output;
  MultiNest* minimizer;
};

// These two functions have to be declared outside the class (declaring them static does not work) because of the definition of the nested::run function (requires function pointer, not method pointer).
void MultiNestLogLike(double* Cube,int& ndim,int& npars,double& lnew,void* e);
void MultiNestDumper(int& nSamples,int& nlive,int& nPar,double** physLive,double** posterior,double** paramConstr,double& maxLogLike,double& logZ,double& INSlogZ,double& logZerr,void* e);


// Iterator class
//================================================================================================================================================
class Iterator: public BaseMinimizer {
public:
  Iterator(std::string name) : BaseMinimizer(name) {
    this->type = "iterator";
  };
  ~Iterator(){};
  int iterations;
  int output_counter;
  
  void minimize(std::map<std::string,std::string> opt,BaseLikelihoodModel* mypars,const std::string output);
  void output(std::string output);
  void finalizeMinimizer(std::string output);
};





/*
// The implementation of mySimplex derived class
//================================================================================================================================================
class mySimplex: public BaseMinimizer,public Simplex {
public:
  mySimplex(ImagePlane* a,BaseSourcePlane* b,CollectionMassModels* c,std::vector<std::map<std::string,BaseNlpar*> > d,mymatrices* e,precomp* f);
  ~mySimplex();
  
  //  double operator()(const std::vector<double> pars){};

  void minimize(std::map<std::string,std::string> minimizer,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,std::vector<std::map<std::string,BaseNlpar*> > nlpars,mymatrices* mat,precomp* pcomp,const std::string output);
  void output();


private:
  ImagePlane* image;
  BaseSourcePlane* source;
  CollectionMassModels* mycollection;
  std::vector<std::map<std::string,BaseNlpar*> > nlpars;
  mymatrices* mat;
  precomp* pcomp;
  std::vector<myactive> active;

  double LogLike(double* pars);
  std::vector<double> initPars(std::vector<std::map<std::string,BaseNlpar*> > nlpars);
  void printPars(const std::vector<double>& c);
};
*/








//class Minuit
//================================================================================================================================================
//class Emcee
//================================================================================================================================================
//class Dnest4
//================================================================================================================================================





// BaseMinimizer class factory
//================================================================================================================================================
class FactoryMinimizer {//This is a singleton class.
public:
  FactoryMinimizer(FactoryMinimizer const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryMinimizer const&) = delete;

  static FactoryMinimizer* getInstance(){
    static FactoryMinimizer dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseMinimizer* createMinimizer(std::string name,std::map<std::string,std::string> minimizer,const std::string output){

    if( minimizer["type"] == "test" ){
      fflush(stdout);
      return new Nothing(name);
    } else if( minimizer["type"] == "multinest" ){
      fflush(stdout);
      return new MultiNest(name);
      //    } else if( minimizer["type"] == "simplex" ){
      //      printf("%-25s","using Simplex");
      //      fflush(stdout);
      //      return new mySimplex(image,source,collection,nlpars,matrices,pcomp);
    } else if( minimizer["type"] == "iterator" ){
      fflush(stdout);
      return new Iterator(name);
    } else {
      return NULL;
    }
  }

private:
  FactoryMinimizer(){};
};



#endif /* MINIMIZERS_HPP */
