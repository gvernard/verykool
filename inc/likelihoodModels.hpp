#ifndef LIKELIHOOD_MODELS_HPP
#define LIKELIHOOD_MODELS_HPP

#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <unordered_map>
#include <iostream>

#include "json/json.h"

class Nlpar;
class ImagePlane;
class BaseSourcePlane;
class CollectionMassModels;
class SmoothAlgebra;
class PertAlgebra;
class Pert;
class BaseMinimizer;


class BaseLikelihoodModel {
public:
  std::string name;
  std::vector<Nlpar*> active;
  std::vector<std::string> active_names;
  std::unordered_map<std::string,double> terms;
  double* maps    = NULL;
  double* means   = NULL;
  double* s1_low  = NULL;
  double* s1_high = NULL;
  ImagePlane* image;
  ImagePlane* model;
  ImagePlane* res;
  CollectionMassModels* collection;
  BaseSourcePlane* source;

  BaseLikelihoodModel(){};
  ~BaseLikelihoodModel();

  std::vector<double> getActiveValues();
  std::vector<std::string> getActiveNlparNames();
  void updateActive(std::vector<double> values);
  void printActive();
  //  void printTerms();
  double getLogLike();

  virtual void printTerms() = 0;
  virtual void getAllNamesValues(std::vector<std::string>& names,std::vector<double>& values) = 0;
  virtual std::vector<Nlpar*> getMassModelPars(int i) = 0;
  virtual std::vector<Nlpar*> getPhysicalPars() = 0;
  virtual Json::Value getActiveJson() = 0;

  virtual void initializeAlgebra() = 0;
  virtual void updateLikelihoodModel() = 0;
  virtual void initialOutputLikelihoodModel(std::string output) = 0;
  virtual void outputLikelihoodModel(std::string output) = 0;
  void finalizeLikelihoodModel(std::string output);
};


class SmoothLikelihood : public BaseLikelihoodModel {
public:
  std::vector<Nlpar*> physical;
  std::vector<std::string> lens_names;
  std::vector< std::vector<Nlpar*> > lenses;
  std::vector<Nlpar*> reg_s;
  SmoothAlgebra* algebra;

  SmoothLikelihood(std::vector<Nlpar*> a,std::vector<Nlpar*> b,std::vector< std::vector<Nlpar*> > d,std::vector<std::string> e,ImagePlane* f,BaseSourcePlane* g,CollectionMassModels* h);
  ~SmoothLikelihood();

  //non-virtual
  void updateCovKernelLimits();
  std::vector<Nlpar*> getRegPars();
  void deriveLinearDpsi(Pert* pert_mass_model,ImagePlane* img_grid);

  //virtual
  void printTerms();
  void getAllNamesValues(std::vector<std::string>& names,std::vector<double>& values);
  std::vector<Nlpar*> getPhysicalPars();
  std::vector<Nlpar*> getMassModelPars(int i);
  Json::Value getActiveJson();
  void getModel();
  void getResidual();

  void initializeAlgebra();
  void updateLikelihoodModel();
  void initialOutputLikelihoodModel(std::string output);
  void outputLikelihoodModel(std::string output);
};


class PertLikelihood : public BaseLikelihoodModel {
public:
  std::vector<Nlpar*> reg_s;
  std::vector<Nlpar*> reg_dpsi;
  BaseSourcePlane* source0;
  Pert* pert_mass_model;
  PertAlgebra* algebra;

  PertLikelihood(std::vector<Nlpar*> reg_s,std::vector<Nlpar*> reg_dpsi,ImagePlane* c,BaseSourcePlane* d,BaseSourcePlane* e,CollectionMassModels* f,Pert* g);
  ~PertLikelihood();

  //non-virtual
  void pertResiduals(std::string output,ImagePlane* image,BaseSourcePlane* source0,Pert* pert_mass_model);

  //virtual
  void printTerms();
  void getAllNamesValues(std::vector<std::string>& names,std::vector<double>& values);
  std::vector<Nlpar*> getPhysicalPars(){};
  std::vector<Nlpar*> getMassModelPars(int i){};
  Json::Value getActiveJson();
  void getModel();
  void getResidual();

  void initializeAlgebra();
  void updateLikelihoodModel();
  void initialOutputLikelihoodModel(std::string output);
  void outputLikelihoodModel(std::string output);
};


/*
class PertIterationLikelihood : public PertLikelihood {
public:
  CollectionMassModels* collection;

  PertIterationLikelihood(std::vector<Nlpar*> reg_s,std::vector<Nlpar*> reg_dpsi,std::string a,std::string b,BaseSourcePlane* c,Pert* d,CollectionMassModels* e);

  //non-virtual
  void initializePert(SmoothLikelihood* smooth_like);

  //virtual
  void initializeAlgebra();
  void updateLikelihoodModel();
  void outputLikelihoodModel(std::string output);
};
*/

#endif /* LIKELIHOOD_MODELS_HPP */
