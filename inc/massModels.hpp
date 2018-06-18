#ifndef MASS_MODELS_HPP
#define MASS_MODELS_HPP

#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <iostream>

#include "nonLinearPars.hpp"
#include "imagePlane.hpp"
#include "tableDefinition.hpp"


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

class Pert: public BaseMassModel{
public:
  ImagePlane* dpsi;
  double* dpsi_dx;
  double* dpsi_dy;
  mytable Aint;
  mytable Bdev;

  Pert(int Ni,int Nj,double width,double height);
  Pert(std::string filename,int Ni,int Nj,double width,double height);
  Pert(ImagePlane* new_dpsi);
  ~Pert(){
    delete(dpsi);
    free(dpsi_dx);
    free(dpsi_dy);
  }
  void defl(double xin,double yin,double& xout,double& yout);
  void updateDpsi(double* new_dpsi);
  void createAint(ImagePlane* data);
  void createBdev();

private:
  double di;
  double dj;
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

  BaseMassModel* createMassModel(const std::string &modelname,std::map<std::string,std::string> pars){
    if( modelname == "pert" ){
      return new Pert(pars["filename"],std::stoi(pars["Ni"]),std::stoi(pars["Nj"]),std::stof(pars["width"]),std::stof(pars["height"]));
    } else {
      return NULL;
    }
  }

private:
  FactoryMassModel(){};
};



class CollectionMassModels {
public:
  int n;
  std::vector<BaseMassModel*> models;
  std::map<std::string,double> mpars;
  
  CollectionMassModels();
  CollectionMassModels(std::vector<Nlpar*> nlpars);
  ~CollectionMassModels();
  void setPhysicalPars(std::vector<Nlpar*> nlpars);
  void all_defl(double xin,double yin,double& xout,double& yout);
  void all_defl(ImagePlane* image);
};

#endif /* MASS_MODELS_HPP */
