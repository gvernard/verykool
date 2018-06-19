#ifndef LIKELIHOOD_MODELS_HPP
#define LIKELIHOOD_MODELS_HPP

#include <string>
#include <vector>
#include <fstream>
#include <map>

#include "json/json.h"



class Nlpar;
class ImagePlane;
class BaseSourcePlane;
class CollectionMassModels;
class StandardAlgebra;
class PerturbationsAlgebra;
class Pert;

class BaseLikelihoodModel {
public:
  std::vector<Nlpar*> active;
  std::map<std::string,double> terms;

  BaseLikelihoodModel(){};
  ~BaseLikelihoodModel(){};

  std::vector<double> getActiveValues();
  std::vector<std::string> getActiveNames();
  void updateActive(std::vector<double> values);
  void printActive();
  void printTerms();

  virtual std::vector<std::string> getFullNames() = 0;
  virtual std::vector<std::string> getActiveFullNames() = 0;
  virtual std::vector<double> getValues() = 0;
  virtual std::vector<Nlpar*> getPhysicalPars() = 0;
  virtual std::vector<Nlpar*> getMassModelPars(int i) = 0;
  virtual Json::Value getActiveNamesValues() = 0;
  virtual void initializeAlgebra() = 0;
  virtual void updateLikelihoodModel() = 0;
  virtual double getLogLike() = 0;
  virtual void outputLikelihoodModel(std::string output) = 0;
};


class StandardLikelihood : public BaseLikelihoodModel {
public:
  std::vector<Nlpar*> physical;
  std::vector<Nlpar*> reg;
  std::vector< std::vector<Nlpar*> > lenses;
  std::vector<std::string> lens_names;
  ImagePlane* image;
  BaseSourcePlane* source;
  CollectionMassModels* collection;
  StandardAlgebra* algebra;
  ImagePlane* model;
  ImagePlane* res;

  StandardLikelihood(std::vector<Nlpar*> a,std::vector<Nlpar*> b,std::vector< std::vector<Nlpar*> > c,std::vector<std::string> d,ImagePlane* e,BaseSourcePlane* f,CollectionMassModels* g);
  ~StandardLikelihood();

  //non-virtual
  std::vector<Nlpar*> getRegPars();
  void getModel();
  void getResidual();
  
  //virtual
  std::vector<std::string> getFullNames();
  std::vector<std::string> getActiveFullNames();
  std::vector<Nlpar*> getPhysicalPars();
  std::vector<Nlpar*> getMassModelPars(int i);
  std::vector<double> getValues();
  Json::Value getActiveNamesValues();
  void initializeAlgebra();
  void updateLikelihoodModel();
  double getLogLike();
  void outputLikelihoodModel(std::string output);
};


class SourceCovarianceKernel : public StandardLikelihood {
public:
  using StandardLikelihood::StandardLikelihood;

  double getLogLike();
};


class PerturbationsLikelihood : public BaseLikelihoodModel {
public:
  std::vector<Nlpar*> reg;
  ImagePlane* image;
  BaseSourcePlane* source;
  CollectionMassModels* collection;
  BaseLikelihoodModel* smooth_like;
  Pert* pert_mass_model;
  PerturbationsAlgebra* algebra;

  PerturbationsLikelihood(std::vector<Nlpar*> reg,ImagePlane* a,BaseSourcePlane* b,CollectionMassModels* c,Pert* d);
  ~PerturbationsLikelihood();

  //non-virtual
  void initializePert(BaseLikelihoodModel* smooth_like);

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
  void outputLikelihoodModel(std::string output);
};






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
      return new StandardLikelihood(physical,reg,lenses,lens_names,image,source,collection);

    } else if( like_model == "source_covariance_kernel" ){

      std::vector<Nlpar*> physical = nlparsFromJsonVector(root["physical"]["nlpars"]);
      std::vector<Nlpar*> reg = nlparsFromJsonVector(root["reg"]["nlpars"]);
      std::vector< std::vector<Nlpar*> > lenses;
      std::vector<std::string> lens_names;
      jmembers = root["lenses"].getMemberNames();
      for(int i=0;i<jmembers.size();i++){
	lenses.push_back( nlparsFromJsonVector(root["lenses"][jmembers[i]]["nlpars"]) );
	lens_names.push_back( jmembers[i] );
      }
      return new SourceCovarianceKernel(physical,reg,lenses,lens_names,image,source,collection);

    } else if( like_model == "perturbations_standard" ){

      std::vector<Nlpar*> reg = nlparsFromJsonVector(root["perturbations"]["reg"]["nlpars"]);
      return new PerturbationsLikelihood(reg,image,source,collection,pert_mass_model);

    } else {
      return NULL;
    }
  }

private:
  FactoryLikelihoodModel(){};
  std::vector<Nlpar*> nlparsFromJsonVector(const Json::Value myjson);
};






#endif /* LIKELIHOOD_MODELS_HPP */
