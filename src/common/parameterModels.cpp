#include "parameterModels.hpp"

#include <string>
#include <vector>
#include <map>

#include "json/json.h"

#include "nonLinearPars.hpp"
#include "massModels.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "tableAlgebra.hpp"


//Abstract class: BaseParameterModel
//===============================================================================================================
std::vector<double> BaseParameterModel::getActiveValues(){
  std::vector<double> values;
  for(int i=0;i<this->active.size();i++){
    values.push_back( this->active[i]->val );
  }
  return values;
}

void BaseParameterModel::updateActive(std::vector<double> values){
  for(int i=0;i<this->active.size();i++){
    this->active[i]->val = values[i];
  }
}





//Derived class from BaseParameterModel: Standard
//===============================================================================================================
Standard::Standard(std::vector<Nlpar*> a,std::vector<Nlpar*> b,std::vector< std::vector<Nlpar*> > c){
  physical = a;
  reg = b;
  lenses = c;
  active = getActiveNlpars();
}

std::vector<Nlpar*> Standard::getActiveNlpars(){
  std::vector<Nlpar*> active;

  for(int i=0;i<this->physical.size();i++){
    if( this->physical[i]->getActive() ){
      active.push_back( this->physical[i] );
    }
  }

  for(int i=0;i<this->reg.size();i++){
    if( this->reg[i]->getActive() ){
      active.push_back( this->reg[i] );
    }
  }

  for(int j=0;j<this->lenses.size();j++){
    for(int i=0;i<this->lenses[j].size();i++){
      if( this->lenses[j][i]->getActive() ){
	active.push_back( this->lenses[j][i] );
      }
    }
  }

  return active;
}


Standard::~Standard(){
  for(int i=0;i<physical.size();i++){
    delete(physical[i]);
  }
  for(int i=0;i<reg.size();i++){
    delete(reg[i]);
  }
  for(int i=0;i<lenses.size();i++){
    for(int j=0;j<lenses[i].size();j++){
      delete(lenses[i][j]);
    }
  }
}



std::vector<Nlpar*> Standard::getRegPars(){
  return this->reg;
}


//    virtual members
std::vector<Nlpar*> Standard::getPhysicalPars(){
  return this->physical;
}

std::vector<Nlpar*> Standard::getMassModelPars(int i){
  return this->lenses[i];
}

std::vector<std::string> Standard::getFullNames(){
  std::vector<std::string> full_names;

  for(int i=0;i<this->physical.size();i++){
    full_names.push_back( this->physical[i]->getName() );
  }

  for(int i=0;i<this->reg.size();i++){
    full_names.push_back( this->reg[i]->getName() );
  }
  
  for(int j=0;j<this->lenses.size();j++){
    for(int i=0;i<this->lenses[j].size();i++){
      full_names.push_back( std::to_string(j) + "_" + this->lenses[j][i]->getName() );
    }
  }
  
  return full_names; 
}

std::vector<std::string> Standard::getActiveFullNames(){
  std::vector<std::string> full_names;

  for(int i=0;i<this->physical.size();i++){
    if( this->physical[i]->getActive() ){
      full_names.push_back( this->physical[i]->getName() );
    }
  }

  for(int i=0;i<this->reg.size();i++){
    if( this->reg[i]->getActive() ){
      full_names.push_back( this->reg[i]->getName() );
    }
  }

  for(int j=0;j<this->lenses.size();j++){
    for(int i=0;i<this->lenses[j].size();i++){
      if( this->lenses[j][i]->getActive() ){
	full_names.push_back( std::to_string(j) + "_" + this->lenses[j][i]->getName() );
      }
    }
  }

  return full_names; 
}

std::vector<double> Standard::getValues(){
  std::vector<double> values;

  for(int i=0;i<this->physical.size();i++){
    values.push_back( this->physical[i]->getValue() );
  }

  for(int i=0;i<this->reg.size();i++){
    values.push_back( this->reg[i]->getValue() );
  }
  
  for(int j=0;j<this->lenses.size();j++){
    for(int i=0;i<this->lenses[j].size();i++){
      values.push_back( this->lenses[j][i]->getValue() );
    }
  }
  
  return values;
}

Json::Value Standard::getActiveNamesValues(){
  Json::Value all;

  Json::Value physical = Json::Value(Json::arrayValue);
  for(int i=0;i<this->physical.size();i++){
    if( this->physical[i]->getActive() ){
      Json::Value dum;
      dum["nam"] = this->physical[i]->getName();
      dum["val"] = this->physical[i]->getValue();
      physical.append(dum);
    }
  }
  all["physical"] = physical;

  Json::Value reg = Json::Value(Json::arrayValue);
  for(int i=0;i<this->reg.size();i++){
    if( this->reg[i]->getActive() ){
      Json::Value dum;
      dum["nam"] = this->reg[i]->getName();
      dum["val"] = this->reg[i]->getValue();
      reg.append(dum);
    }
  }
  all["reg"] = reg;

  Json::Value lenses = Json::Value(Json::arrayValue);
  for(int j=0;j<this->lenses.size();j++){
    Json::Value lens = Json::Value(Json::arrayValue);
    for(int i=0;i<this->lenses[j].size();i++){
      if( this->lenses[j][i]->getActive() ){
	Json::Value dum;
	dum["nam"] = this->lenses[j][i]->getName();
	dum["val"] = this->lenses[j][i]->getValue();
	lens.append(dum);
      }
    }
    lenses.append(lens);
  }
  all["lenses"] = lenses;

  return all;
}

void Standard::updateParameterModel(ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,mymatrices* mat,precomp* pcomp){
  for(int i=0;i<mycollection->models.size();i++){
    mycollection->models[i]->setMassPars(this->lenses[i]);
  }
  mycollection->setPhysicalPars(this->physical);

  if( source->sample_reg ){
    source->kernel->setParameters(this->reg);
  }

  if( source->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(source);
    ada->createAdaGrid(image,mycollection);
    ada->createDelaunay();
  }

  if( source->sample_reg || source->type == "adaptive" ){
    source->constructH(&mat->H);
  }

  //  std::cout << "1" << std::endl;
  source->constructL(image,mycollection,&mat->L);
  //  std::cout << "2" << std::endl;
  setAlgebraRuntime(image,source,this->reg,mat,pcomp);
  //  std::cout << "3" << std::endl;
  solveLinearSparseS(image,source,pcomp);
  //  std::cout << "4" << std::endl;
}

double Standard::getLogLike(ImagePlane* image,BaseSourcePlane* source,precomp* pcomp){
  double pi  = 3.14159265358979323846;

  if( source->reg == "covariance_kernel" ){
    getMin(image,source,pcomp);
    double g   = pcomp->chi2/2. + pcomp->reg/2.;
    double f1  = image->lookup.size()*log10(2*pi)/2.;
    //double f2  = source->Sm*log10(nlpars[1]["lambda"]->val)/2.;
    double val = -g -f1 +(pcomp->detC)/2. +(pcomp->detHtH)/2. -(pcomp->detA)/2.;
    
    std::vector<double> par_vals = this->getValues();
    for(int i=0;i<par_vals.size();i++){
      printf(" %12.8f",par_vals[i]);
    }
    printf("\n");
    printf("%16.7f  %16.7f  %16.7f  %16.7f  %16.7f  %16.7f  %16.7f\n",pcomp->chi2,pcomp->reg,-g,-f1,pcomp->detC/2.,pcomp->detHtH/2.,-pcomp->detA/2.);
    printf("%16.7f\n\n",val);

    return val;
  } else {
    getMin(image,source,pcomp);
    double g   = pcomp->chi2/2.0 + Nlpar::getValueByName("lambda",this->reg)*pcomp->reg/2.0;
    double f1  = image->lookup.size()*log10(2*pi)/2.0;
    double f2  = source->Sm*log10(Nlpar::getValueByName("lambda",this->reg))/2.0;
    double val = -g -f1 +f2 +(pcomp->detC)/2.0 +(pcomp->detHtH)/2.0 -(pcomp->detA)/2.0;
    
    std::vector<double> par_vals = this->getValues();
    for(int i=0;i<par_vals.size();i++){
      printf(" %12.8f",par_vals[i]);
    }
    printf("\n");
    printf("%16.7f  %16.7f  %16.7f  %16.7f  %16.7f  %16.7f  %16.7f  %16.7f  %16.7f\n",pcomp->chi2,pcomp->reg,Nlpar::getValueByName("lambda",this->reg)*pcomp->reg,-g,-f1,f2,pcomp->detC/2.0,pcomp->detHtH/2.0,-pcomp->detA/2.0);
    printf("%16.7f\n\n",val);

    return val;
  }

  //  return pow(nlpars[2]["b"]->val-6,2) + pow(nlpars[2]["q"]->val-6,2) + pow(nlpars[2]["pa"]->val-6,2) + pow(nlpars[2]["x0"]->val-6,2) + pow(nlpars[2]["y0"]->val-6,2);
}

















//Factory of BaseParameterModel 
//===============================================================================================================
std::vector<Nlpar*> FactoryParameterModel::nlparsFromJsonVector(const Json::Value myjson){
  std::vector<Nlpar*> pars;

  for(int i=0;i<myjson.size();i++){
    const Json::Value nlpar = myjson[i];

    std::string nam = nlpar["nam"].asString();
    int fix         = nlpar["fix"].asInt();
    int per         = nlpar["per"].asInt();
    double val      = nlpar["val"].asDouble();
    double err      = nlpar["err"].asDouble();
    double min      = nlpar["min"].asDouble();
    double max      = nlpar["max"].asDouble();
    Nlpar* par = new Nlpar(nam,fix,per,val,err,min,max);
    
    std::map<std::string,std::string> prior_pars;
    Json::Value::Members jmembers = nlpar["pri"].getMemberNames();
    for(int j=0;j<jmembers.size();j++){
      prior_pars[jmembers[j]] = nlpar["pri"][jmembers[j]].asString();
    }
    par->pri = FactoryPrior::getInstance()->createPrior(par,prior_pars);
    //par->setNewPrior(prior);
    
    pars.push_back(par);
  }
  
  return pars;
}


