#ifndef MINIMIZERS_HPP
#define MINIMIZERS_HPP

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>
#include <functional>

#include "multinest.h"

#include "inputOutput.hpp"
#include "massModels.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "tableAlgebra.hpp"
#include "nonLinearPars.hpp"




static double mainLogLike(ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,std::vector<std::map<std::string,BaseNlpar*> > nlpars,mymatrices* mat,precomp* pcomp){
  for(int i=0;i<mycollection->models.size();i++){
    mycollection->models[i]->setMassPars(nlpars[2+i]);
  }
  mycollection->setPhysicalPars(nlpars[0]);
  
  if( source->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(source);
    ada->createAdaGrid(image,mycollection);
    ada->createDelaunay();
    source->constructH(&mat->H);
  }


  
  if( source->reg == "covariance_kernel" ){
    getInverseCovarianceMatrix(source,nlpars[1],pcomp);
  }


  source->constructL(image,mycollection,&mat->L);
  setAlgebraRuntime(image,source,nlpars[1]["lambda"]->val,mat,pcomp);
  solveLinearSparseS(image,source,pcomp);
  double loglike = getLogLike(image,source,nlpars[1]["lambda"]->val,pcomp,nlpars);
  return loglike;



  //  return pow(pars[0]-6,2)+pow(pars[1]-6,2)+pow(pars[2]-6,2)+pow(pars[3]-6,2)+pow(pars[4]-6,2)+pow(pars[5]-6,2);
}







// BaseMinimizer parent/base class
//================================================================================================================================================
class BaseMinimizer {
public:
  BaseMinimizer(){};
  ~BaseMinimizer(){};

  virtual void minimize(std::map<std::string,std::string> minimizer,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,std::vector<std::map<std::string,BaseNlpar*> > nlpars,mymatrices* mat,precomp* pcomp,const std::string output) = 0;
  virtual void output() = 0;
};



// Nothing class (does nothing)
//================================================================================================================================================
class Nothing: public BaseMinimizer {
public:
  Nothing(){};
  ~Nothing(){};
  
  void minimize(std::map<std::string,std::string> minimizer,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,std::vector<std::map<std::string,BaseNlpar*> > nlpars,mymatrices* mat,precomp* pcomp,const std::string output){};
  void output(){};
};


// MultiNest class and related stuff
//================================================================================================================================================
struct extras{
  ImagePlane* image;
  BaseSourcePlane* source;
  CollectionMassModels* mycollection;
  std::vector<std::map<std::string,BaseNlpar*> > nlpars;
  std::vector<myactive> active;
  mymatrices* mat;
  std::string output;
  std::vector<std::map<std::string,double> >* map_pars;
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
  
  void minimize(std::map<std::string,std::string> minimizer,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,std::vector<std::map<std::string,BaseNlpar*> > nlpars,mymatrices* mat,precomp* pcomp,const std::string output);
  void output();
};




// Simplex class and related stuff
//================================================================================================================================================
class Simplex: public BaseMinimizer {
public:
  Simplex(ImagePlane* a,BaseSourcePlane* b,CollectionMassModels* c,std::vector<std::map<std::string,BaseNlpar*> > d,mymatrices* e,precomp* f);
  ~Simplex();
  
  double operator()(const std::vector<double> pars){};

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

  double LogLike(const std::vector<double> pars);
  std::vector<double> initPars(std::vector<std::map<std::string,BaseNlpar*> > nlpars);
  void printPars(const std::vector<double>& c);
  double myrandom();
  std::vector<double> SimplexMinimizer(std::vector<double> init,double tol,int iterations);
  //  template<class UC> std::vector<double> SimplexAnnealing(UC f,std::vector<double> init,double tol,int iterations,std::vector<std::vector<double> > x);
};



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

  BaseMinimizer* createMinimizer(std::map<std::string,std::string> minimizer,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* collection,std::vector<std::map<std::string,BaseNlpar*> > nlpars,mymatrices* matrices,precomp* pcomp,const std::string output){
    if( minimizer["type"] == "0" ){
      printf("%-25s","using given parameters");
      fflush(stdout);
      return new Nothing();
    } else if( minimizer["type"] == "multinest" ){
      printf("%-25s","using MultiNest");
      fflush(stdout);
      return new MultiNest();
    } else if( minimizer["type"] == "simplex" ){
      printf("%-25s","using Simplex");
      fflush(stdout);
      return new Simplex(image,source,collection,nlpars,matrices,pcomp);
    } else {
      return NULL;
    }
  }

private:
  FactoryMinimizer(){};
};



#endif /* MINIMIZERS_HPP */
