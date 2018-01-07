#ifndef MINIMIZERS_HPP
#define MINIMIZERS_HPP

#include <string>
#include <vector>
#include <map>

#include "multinest.h"
//#include "simplex.hpp"

class BaseParameterModel;
class ImagePlane;
class BaseSourcePlane;
class CollectionMassModels;
class mymatrices;
class precomp;



// BaseMinimizer parent/base class
//================================================================================================================================================
class BaseMinimizer {
public:
  BaseMinimizer(){};
  ~BaseMinimizer(){};

  virtual void minimize(std::map<std::string,std::string> minimizer,BaseParameterModel* mypars,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,mymatrices* mat,precomp* pcomp,const std::string output) = 0;
  virtual void output() = 0;
};



// Nothing class (does nothing)
//================================================================================================================================================
class Nothing: public BaseMinimizer {
public:
  Nothing(){};
  ~Nothing(){};
  
  void minimize(std::map<std::string,std::string> minimizer,BaseParameterModel* mypars,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,mymatrices* mat,precomp* pcomp,const std::string output){};
  void output(){};
};


// MultiNest class and related stuff
//================================================================================================================================================
struct extras{
  BaseParameterModel* pars;
  ImagePlane* image;
  BaseSourcePlane* source;
  CollectionMassModels* mycollection;
  mymatrices* mat;
  std::string output;
  precomp* pcomp;
  int counter;
};

// These two functions have to be declared outside the class (declaring them static does not work) because of the definition of the nested::run function (requires function pointer, not method pointer).
void MultiNestLogLike(double* Cube,int& ndim,int& npars,double& lnew,void* e);
void MultiNestDumper(int& nSamples,int& nlive,int& nPar,double** physLive,double** posterior,double** paramConstr,double& maxLogLike,double& logZ,double& INSlogZ,double& logZerr,void* e);

class MultiNest: public BaseMinimizer {
public:
  MultiNest(){};
  ~MultiNest(){};
  
  void minimize(std::map<std::string,std::string> minimizer,BaseParameterModel* pars,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,mymatrices* mat,precomp* pcomp,const std::string output);
  void output();
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

  BaseMinimizer* createMinimizer(std::map<std::string,std::string> minimizer,BaseParameterModel* pars,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* collection,mymatrices* matrices,precomp* pcomp,const std::string output){

    if( minimizer["type"] == "test" ){
      printf("%-25s","using given parameters");
      fflush(stdout);
      return new Nothing();
    } else if( minimizer["type"] == "multinest" ){
      printf("%-25s","using MultiNest");
      fflush(stdout);
      return new MultiNest();
      //    } else if( minimizer["type"] == "simplex" ){
      //      printf("%-25s","using Simplex");
      //      fflush(stdout);
      //      return new mySimplex(image,source,collection,nlpars,matrices,pcomp);
    } else {
      return NULL;
    }
  }

private:
  FactoryMinimizer(){};
};



#endif /* MINIMIZERS_HPP */
