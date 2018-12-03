#ifndef LIKELIHOOD_MODELS_HPP
#define LIKELIHOOD_MODELS_HPP

#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <unordered_map>

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

  BaseLikelihoodModel(){};
  ~BaseLikelihoodModel(){
    free(maps);
    free(means);
    free(s1_low);
    free(s1_high);
  };

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
};



class SmoothLikelihood : public BaseLikelihoodModel {
public:
  std::vector<Nlpar*> physical;
  std::vector<Nlpar*> reg;
  std::vector< std::vector<Nlpar*> > lenses;
  std::vector<std::string> lens_names;
  ImagePlane* image;
  BaseSourcePlane* source;
  CollectionMassModels* collection;
  SmoothAlgebra* algebra;
  ImagePlane* model;
  ImagePlane* res;

  SmoothLikelihood(std::vector<Nlpar*> a,std::vector<Nlpar*> b,std::vector< std::vector<Nlpar*> > d,std::vector<std::string> e,ImagePlane* f,BaseSourcePlane* g,CollectionMassModels* h);
  ~SmoothLikelihood();

  //non-virtual
  std::vector<Nlpar*> getRegPars();
  void getModel();
  void getResidual();
  void deriveLinearDpsi(Pert* pert_mass_model,ImagePlane* img_grid);
  
  //virtual
  void printTerms();
  void getAllNamesValues(std::vector<std::string>& names,std::vector<double>& values);
  std::vector<Nlpar*> getPhysicalPars();
  std::vector<Nlpar*> getMassModelPars(int i);
  Json::Value getActiveJson();

  void initializeAlgebra();
  void updateLikelihoodModel();
  void initialOutputLikelihoodModel(std::string output);
  void outputLikelihoodModel(std::string output);
};



class PertLikelihood : public BaseLikelihoodModel {
public:
  std::vector<Nlpar*> reg_s;
  std::vector<Nlpar*> reg_dpsi;
  std::string reg_s_type;
  std::string reg_dpsi_type;
  SmoothLikelihood* smooth_like;
  BaseSourcePlane* source;
  Pert* pert_mass_model;
  PertAlgebra* algebra;

  PertLikelihood(std::vector<Nlpar*> reg_s,std::vector<Nlpar*> reg_dpsi,std::string a,std::string b,BaseSourcePlane* c,Pert* d);
  ~PertLikelihood();

  //non-virtual
  void initializePert(SmoothLikelihood* smooth_like);
  void pertResiduals(std::string output,ImagePlane* image,BaseSourcePlane* source,Pert* pert_mass_model);

  //virtual
  void printTerms();
  void getAllNamesValues(std::vector<std::string>& names,std::vector<double>& values);
  std::vector<Nlpar*> getPhysicalPars(){};
  std::vector<Nlpar*> getMassModelPars(int i){};
  Json::Value getActiveJson();

  void initializeAlgebra();
  void updateLikelihoodModel();
  void initialOutputLikelihoodModel(std::string output);
  void outputLikelihoodModel(std::string output);
};



class PertIterationLikelihood : public PertLikelihood {
public:
  CollectionMassModels* collection;

  PertIterationLikelihood(std::vector<Nlpar*> reg_s,std::vector<Nlpar*> reg_dpsi,std::string a,std::string b,BaseSourcePlane* c,Pert* d,CollectionMassModels* e);

  //non-virtual
  void initializePert(SmoothLikelihood* smooth_like);

  //virtual
  void initializeAlgebra();
  void updateLikelihoodModel();
  void initialOutputLikelihoodModel(std::string output);
  void outputLikelihoodModel(std::string output);
};


/*
class PertIterationLikelihood : public BaseLikelihoodModel {
public:
  std::vector<Nlpar*> reg_s;
  std::vector<Nlpar*> reg_dpsi;
  std::string reg_s_type;
  std::string reg_dpsi_type;
  SmoothLikelihood* smooth_like;
  BaseSourcePlane* source;
  Pert* pert_mass_model;
  PertAlgebra* algebra;
  CollectionMassModels* collection;

  PertIterationLikelihood(std::vector<Nlpar*> reg_s,std::vector<Nlpar*> reg_dpsi,std::string a,std::string b,BaseSourcePlane* c,Pert* d,CollectionMassModels* e);
  ~PertIterationLikelihood();

  //non-virtual
  void initializePert(SmoothLikelihood* smooth_like);

  //virtual
  std::vector<std::string> getFullNames(){};
  std::vector<std::string> getActiveFullNames(){};
  std::vector<double> getValues(){};
  std::vector<Nlpar*> getPhysicalPars(){};
  std::vector<Nlpar*> getMassModelPars(int i){};
  Json::Value getActiveNamesValues(){};
  void initializeAlgebra();
  void updateLikelihoodModel();
  double getLogLike();
  void initialOutputLikelihoodModel(std::string output);
  void outputLikelihoodModel(std::string output);
};
*/


class FactoryLikelihoodModel {//This is a singleton class.
public:
  FactoryLikelihoodModel(FactoryLikelihoodModel const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryLikelihoodModel const&) = delete;

  static FactoryLikelihoodModel* getInstance(){
    static FactoryLikelihoodModel dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseLikelihoodModel* createLikelihoodModel(std::string path,std::string run,std::string like_model,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* collection,Pert* pert_mass_model){
    Json::Value root;
    Json::Value::Members jmembers;
    std::string runname = path + run;
    std::ifstream fin(runname+"vkl_input.json");
    fin >> root;
    
    if( like_model == "standard" ){
      std::vector<Nlpar*> physical = nlparsFromJsonVector(root["physical"]["nlpars"]);
      std::vector<Nlpar*> reg = nlparsFromJsonVector(root["reg"]["nlpars"]);
      std::vector< std::vector<Nlpar*> > lenses;
      std::vector<std::string> lens_names;
      jmembers = root["lenses"].getMemberNames();
      for(int i=0;i<jmembers.size();i++){
	lenses.push_back( nlparsFromJsonVector(root["lenses"][jmembers[i]]["nlpars"]) );
	lens_names.push_back( jmembers[i] );
      }
      return new SmoothLikelihood(physical,reg,lenses,lens_names,image,source,collection);
    } else if( like_model == "perturbations_standard" ){
      std::vector<Nlpar*> reg_s    = nlparsFromJsonVector(root["perturbations"]["reg_s"]["nlpars"]);
      std::vector<Nlpar*> reg_dpsi = nlparsFromJsonVector(root["perturbations"]["reg_dpsi"]["nlpars"]);
      std::string reg_s_type       = root["perturbations"]["reg_s"]["type"].asString();
      std::string reg_dpsi_type    = root["perturbations"]["reg_dpsi"]["type"].asString();
      return new PertLikelihood(reg_s,reg_dpsi,reg_s_type,reg_dpsi_type,source,pert_mass_model);
    } else if( like_model == "perturbations_iter" ){
      std::vector<Nlpar*> reg_s    = nlparsFromJsonVector(root["perturbations"]["reg_s"]["nlpars"]);
      std::vector<Nlpar*> reg_dpsi = nlparsFromJsonVector(root["perturbations"]["reg_dpsi"]["nlpars"]);
      std::string reg_s_type       = root["perturbations"]["reg_s"]["type"].asString();
      std::string reg_dpsi_type    = root["perturbations"]["reg_dpsi"]["type"].asString();
      return new PertIterationLikelihood(reg_s,reg_dpsi,reg_s_type,reg_dpsi_type,source,pert_mass_model,collection);
    } else {
      return NULL;
    }
  }

private:
  FactoryLikelihoodModel(){};
  std::vector<Nlpar*> nlparsFromJsonVector(const Json::Value myjson);
};



#endif /* LIKELIHOOD_MODELS_HPP */
