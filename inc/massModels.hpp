#ifndef MASS_MODELS_HPP
#define MASS_MODELS_HPP

#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "nonLinearPars.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "tableDefinition.hpp"

extern "C"{
  void fastelldefl_(double* x1,double* x2,double* b,double* g,double* q,double* s2,double* defl);
}


class BaseMassModel{
public:
  int n;
  std::string type;
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
  FixedSource* dpsi;
  double* dpsi_dx;
  double* dpsi_dy;
  mytable Aint;
  mytable Bdev;
  std::string interp = "bicubic";
  gsl_spline2d* x_spline = NULL;
  gsl_spline2d* y_spline = NULL;
  gsl_interp_accel* xacc = gsl_interp_accel_alloc();
  gsl_interp_accel* yacc = gsl_interp_accel_alloc();
  double* x_tmp;
  double* y_tmp;
  double* za_x;
  double* za_y;
  
  Pert(int Ni,int Nj,double width,double height,std::string reg);
  Pert(int Ni,int Nj,ImagePlane* image,std::string reg);
  Pert(std::string filename,int Ni,int Nj,double width,double height,std::string reg);
  Pert(FixedSource* new_dpsi);
  ~Pert(){
    delete(dpsi);
    free(dpsi_dx);
    free(dpsi_dy);
    gsl_spline2d_free(x_spline);
    gsl_spline2d_free(y_spline);
    gsl_interp_accel_free(xacc);
    gsl_interp_accel_free(yacc);
    free(x_tmp);
    free(y_tmp);
    free(za_x);
    free(za_y);
  }
  void defl(double xin,double yin,double& xout,double& yout);
  void replaceDpsi(double* new_dpsi);
  void addDpsi(double* corrections);
  void updatePert();
  void createAint(ImagePlane* data);
  void createBdev();
  //  void tableDefl(int Nm,double* xdefl,double* ydefl);
  void createCrosses(ImagePlane* image);
  void getConvergence(ImagePlane* kappa);
  
private:
  double di;
  double dj;
  void derivativeDirection(int q,int qmax,double den,int* rel_ind,double* coeff);
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
      return new Pert(pars["filename"],std::stoi(pars["Ni"]),std::stoi(pars["Nj"]),std::stof(pars["width"]),std::stof(pars["height"]),pars["reg"]);
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
