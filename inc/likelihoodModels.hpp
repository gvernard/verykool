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

  virtual void initializeAlgebra(ImagePlane* image,BaseSourcePlane* source) = 0;
  virtual void updateLikelihoodModel(ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection) = 0;
  virtual double getLogLike(ImagePlane* image,BaseSourcePlane* source) = 0;
  virtual void outputLikelihoodModel(ImagePlane* image,BaseSourcePlane* source,std::string output) = 0;
};


class StandardLikelihood : public BaseLikelihoodModel {
public:
  std::vector<Nlpar*> physical;
  std::vector<Nlpar*> reg;
  std::vector< std::vector<Nlpar*> > lenses;

  StandardLikelihood(std::vector<Nlpar*> a,std::vector<Nlpar*> b,std::vector< std::vector<Nlpar*> > c,std::vector<std::string> d);
  ~StandardLikelihood();
  
  std::vector<Nlpar*> getRegPars();
  void getSourceErrors(int i,double* errors);
  void getMockData(ImagePlane* mockdata,BaseSourcePlane* source);

  //virtual
  std::vector<std::string> getFullNames();
  std::vector<std::string> getActiveFullNames();
  std::vector<Nlpar*> getPhysicalPars();
  std::vector<Nlpar*> getMassModelPars(int i);
  std::vector<double> getValues();
  Json::Value getActiveNamesValues();

  void initializeAlgebra(ImagePlane* image,BaseSourcePlane* source);
  void updateLikelihoodModel(ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection);
  double getLogLike(ImagePlane* image,BaseSourcePlane* source);
  void outputLikelihoodModel(ImagePlane* image,BaseSourcePlane* source,std::string output);

private:
  StandardAlgebra* algebra;
  std::vector<std::string> lens_names;
};



class SourceCovarianceKernel : public StandardLikelihood {
public:
  using StandardLikelihood::StandardLikelihood;

  double getLogLike(ImagePlane* image,BaseSourcePlane* source);

};




class FactoryLikelihoodModel {//This is a singleton class.
public:
  FactoryLikelihoodModel(FactoryLikelihoodModel const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryLikelihoodModel const&) = delete;

  static FactoryLikelihoodModel* getInstance(){
    static FactoryLikelihoodModel dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseLikelihoodModel* createLikelihoodModel(std::string path,std::string run){
    Json::Value root;
    Json::Value::Members jmembers;
    
    std::string runname = path + run;
    std::ifstream fin(runname+"vkl_input.json");
    fin >> root;
    
    if( root["parameter_model"] == "standard" ){

      std::vector<Nlpar*> physical = nlparsFromJsonVector(root["physical"]["nlpars"]);
      std::vector<Nlpar*> reg = nlparsFromJsonVector(root["reg"]["nlpars"]);
      std::vector< std::vector<Nlpar*> > lenses;
      std::vector<std::string> lens_names;
      jmembers = root["lenses"].getMemberNames();
      for(int i=0;i<jmembers.size();i++){
	lenses.push_back( nlparsFromJsonVector(root["lenses"][jmembers[i]]["nlpars"]) );
	lens_names.push_back( jmembers[i] );
      }
      return new StandardLikelihood(physical,reg,lenses,lens_names);

    } else if( root["parameter_model"] == "source_covariance_kernel" ){

      std::vector<Nlpar*> physical = nlparsFromJsonVector(root["physical"]["nlpars"]);
      std::vector<Nlpar*> reg = nlparsFromJsonVector(root["reg"]["nlpars"]);
      std::vector< std::vector<Nlpar*> > lenses;
      std::vector<std::string> lens_names;
      jmembers = root["lenses"].getMemberNames();
      for(int i=0;i<jmembers.size();i++){
	lenses.push_back( nlparsFromJsonVector(root["lenses"][jmembers[i]]["nlpars"]) );
	lens_names.push_back( jmembers[i] );
      }
      return new SourceCovarianceKernel(physical,reg,lenses,lens_names);

    } else {
      return NULL;
    }
  }

private:
  FactoryLikelihoodModel(){};
  std::vector<Nlpar*> nlparsFromJsonVector(const Json::Value myjson);
};






#endif /* LIKELIHOOD_MODELS_HPP */
