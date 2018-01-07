#ifndef PARAMETER_MODELS_HPP
#define PARAMETER_MODELS_HPP

#include <string>
#include <vector>
#include <fstream>

#include "json/json.h"



class Nlpar;
class ImagePlane;
class BaseSourcePlane;
class CollectionMassModels;
class mymatrices;
class precomp;


class BaseParameterModel {
public:
  std::vector<Nlpar*> active;
  
  BaseParameterModel(){};
  ~BaseParameterModel(){};

  std::vector<double> getActiveValues();
  void updateActive(std::vector<double> values);

  virtual std::vector<Nlpar*> getActiveNlpars() = 0;
  virtual std::vector<std::string> getFullNames() = 0;
  virtual std::vector<std::string> getActiveFullNames() = 0;
  virtual std::vector<double> getValues() = 0;
  virtual std::vector<Nlpar*> getPhysicalPars() = 0;
  virtual std::vector<Nlpar*> getMassModelPars(int i) = 0;
  virtual Json::Value getActiveNamesValues() = 0;

  virtual void updateParameterModel(ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,mymatrices* mat,precomp* pcomp) = 0;
  virtual double getLogLike(ImagePlane* image,BaseSourcePlane* source,precomp* pcomp) = 0;
};


class Standard : public BaseParameterModel {
public:
  std::vector<Nlpar*> physical;
  std::vector<Nlpar*> reg;
  std::vector< std::vector<Nlpar*> > lenses;

  Standard(std::vector<Nlpar*> a,std::vector<Nlpar*> b,std::vector< std::vector<Nlpar*> > c);
  ~Standard();
  
  std::vector<Nlpar*> getRegPars();

  //virtual
  std::vector<Nlpar*> getActiveNlpars();
  std::vector<std::string> getFullNames();
  std::vector<std::string> getActiveFullNames();
  std::vector<Nlpar*> getPhysicalPars();
  std::vector<Nlpar*> getMassModelPars(int i);
  std::vector<double> getValues();
  Json::Value getActiveNamesValues();

  void updateParameterModel(ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,mymatrices* mat,precomp* pcomp);
  double getLogLike(ImagePlane* image,BaseSourcePlane* source,precomp* pcomp);
};


class FactoryParameterModel {//This is a singleton class.
public:
  FactoryParameterModel(FactoryParameterModel const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryParameterModel const&) = delete;

  static FactoryParameterModel* getInstance(){
    static FactoryParameterModel dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseParameterModel* createParameterModel(std::string path,std::string run){
    Json::Value root;
    Json::Value::Members jmembers;
    
    std::string runname = path + run;
    std::ifstream fin(runname+"vkl_input.json");
    fin >> root;
    
    if( root["parameter_model"] == "standard" ){
      std::vector<Nlpar*> physical = nlparsFromJsonVector(root["physical"]);
      std::vector<Nlpar*> reg = nlparsFromJsonVector(root["reg"]["nlpars"]);
      std::vector< std::vector<Nlpar*> > lenses;
      for(int i=0;i<root["lenses"].size();i++){
	lenses.push_back( nlparsFromJsonVector(root["lenses"][i]["nlpars"]) );
      }

      return new Standard(physical,reg,lenses);
    } else {
      return NULL;
    }
  }

private:
  FactoryParameterModel(){};
  std::vector<Nlpar*> nlparsFromJsonVector(const Json::Value myjson);
};






#endif /* PARAMETER_MODELS_HPP */
