#ifndef MASS_MODELS_HPP
#define MASS_MODELS_HPP

#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <iostream>

#include "nonLinearPars.hpp"


extern "C"{
  void fastelldefl_(double* x1,double* x2,double* b,double* g,double* q,double* s2,double* defl);
}


class BaseMassModel{
public:
  int n;
  std::map<std::string,double> mpars;

  BaseMassModel(){};
  ~BaseMassModel(){
    mpars.clear();
  };
  
  void setMassPars(std::vector<Nlpar*> nlpars);
  void printMassPars();
  virtual void defl(double xin,double yin,double& xout,double& yout) = 0;
};


class Sie: public BaseMassModel{
public:
  Sie(std::vector<Nlpar*> nlpars);
  void defl(double xin,double yin,double& xout,double& yout);
};

class Spemd: public BaseMassModel{
public:
  Spemd(std::vector<Nlpar*> nlpars);
  void defl(double xin,double yin,double& xout,double& yout);
};



class FactoryMassModel{//This is a singleton class.
public:
  FactoryMassModel(FactoryMassModel const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryMassModel const&) = delete;

  static FactoryMassModel* getInstance(){
    static FactoryMassModel dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseMassModel* createMassModel(const std::string &modelname,std::vector<Nlpar*> nlpars){
    if( modelname == "sie" ){
      return new Sie(nlpars);
    } else if ( modelname == "spemd" ){
      return new Spemd(nlpars);
    } else {
      return NULL;
    }
  }

private:
  FactoryMassModel(){};
};



class CollectionMassModels{
public:
  int n;
  std::vector<BaseMassModel*> models;
  std::map<std::string,double> mpars;
  
  CollectionMassModels(){
    this->mpars["g1"] = 0.0;
    this->mpars["g2"] = 0.0;
  }
  CollectionMassModels(std::vector<Nlpar*> nlpars){
    this->setPhysicalPars(nlpars);
  };
  ~CollectionMassModels(){
    for(int i=0;i<this->models.size();i++){
      delete(this->models[i]);
    }
    mpars.clear();
  };

  void setPhysicalPars(std::vector<Nlpar*> nlpars){
    for(int i=0;i<nlpars.size();i++){
      this->mpars[nlpars[i]->nam] = nlpars[i]->val;
    }
    this->mpars["phi"] *= 0.01745329251;
      this->mpars["g1"]  = this->mpars["g"]*cos(2*this->mpars["phi"]);
    this->mpars["g2"]  = this->mpars["g"]*sin(2*this->mpars["phi"]);
  }

  void all_defl(double xin,double yin,double& xout,double& yout){
    double ax   = 0.;
    double ay   = 0.;
    double dumx = 0.;
    double dumy = 0.;
    for(int i=0;i<this->models.size();i++){
      this->models[i]->defl(xin,yin,dumx,dumy);
      ax += dumx;
      ay += dumy;
    }
    xout = (1.0-this->mpars["g1"])*xin - this->mpars["g2"]*yin - ax;
    yout = (1.0+this->mpars["g1"])*yin - this->mpars["g2"]*xin - ay;
  }
};

#endif /* MASS_MODELS_HPP */
