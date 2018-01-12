#include "likelihoodModels.hpp"

#include <string>
#include <vector>
#include <map>

#include "json/json.h"

#include "nonLinearPars.hpp"
#include "massModels.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "eigenAlgebra.hpp"


//Abstract class: BaseLikelihoodModel
//===============================================================================================================
std::vector<double> BaseLikelihoodModel::getActiveValues(){
  std::vector<double> values;
  for(int i=0;i<this->active.size();i++){
    values.push_back( this->active[i]->val );
  }
  return values;
}

std::vector<std::string> BaseLikelihoodModel::getActiveNames(){
  std::vector<std::string> names;
  for(int i=0;i<this->active.size();i++){
    names.push_back( this->active[i]->nam );
  }
  return names;
}

void BaseLikelihoodModel::updateActive(std::vector<double> values){
  for(int i=0;i<this->active.size();i++){
    this->active[i]->val = values[i];
  }
}

void BaseLikelihoodModel::printActive(){
  for(int i=0;i<this->active.size();i++){
    printf(" %12s",(this->active[i]->nam).c_str());
  }
  printf("\n");
  for(int i=0;i<this->active.size();i++){
    printf(" %12.5f",this->active[i]->val);
  }
  printf("\n");
}

void BaseLikelihoodModel::printTerms(){
  for(std::map<std::string,double>::iterator it=this->terms.begin();it!=this->terms.end();it++){
    printf(" %16s",(it->first).c_str());
  }
  printf("\n");
  for(std::map<std::string,double>::iterator it=this->terms.begin();it!=this->terms.end();it++){
    printf(" %16.5f",it->second);
  }
  printf("\n");
}



//Derived class from BaseLikelihoodModel: StandardLikelihood
//===============================================================================================================
StandardLikelihood::StandardLikelihood(std::vector<Nlpar*> a,std::vector<Nlpar*> b,std::vector< std::vector<Nlpar*> > c,std::vector<std::string> d){
  physical = a;
  reg = b;
  lenses = c;
  lens_names = d;

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

  this->algebra = new StandardAlgebra(this);

  terms["chi2"]    = 0.0;
  terms["reg"]     = 0.0;
  terms["Nilog2p"] = 0.0;
  terms["Nslogl"]  = 0.0;
  terms["detC"]    = 0.0;
  terms["detHtH"]  = 0.0;
  terms["detA"]    = 0.0;
  terms["like"]    = 0.0;
}


StandardLikelihood::~StandardLikelihood(){
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

  delete(algebra);
}



std::vector<Nlpar*> StandardLikelihood::getRegPars(){
  return this->reg;
}

//    virtual members
std::vector<Nlpar*> StandardLikelihood::getPhysicalPars(){
  return this->physical;
}

std::vector<Nlpar*> StandardLikelihood::getMassModelPars(int i){
  return this->lenses[i];
}

std::vector<std::string> StandardLikelihood::getFullNames(){
  std::vector<std::string> full_names;

  for(int i=0;i<this->physical.size();i++){
    full_names.push_back( this->physical[i]->getName() );
  }

  for(int i=0;i<this->reg.size();i++){
    full_names.push_back( this->reg[i]->getName() );
  }
  
  for(int j=0;j<this->lenses.size();j++){
    for(int i=0;i<this->lenses[j].size();i++){
      full_names.push_back( lens_names[j] + "_" + this->lenses[j][i]->getName() );
    }
  }
  
  return full_names; 
}

std::vector<std::string> StandardLikelihood::getActiveFullNames(){
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
	full_names.push_back( lens_names[j] + "_" + this->lenses[j][i]->getName() );
      }
    }
  }

  return full_names; 
}

std::vector<double> StandardLikelihood::getValues(){
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

Json::Value StandardLikelihood::getActiveNamesValues(){
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
  if( physical.size() != 0 ){
    all["physical"] = physical;
  }

  Json::Value reg = Json::Value(Json::arrayValue);
  for(int i=0;i<this->reg.size();i++){
    if( this->reg[i]->getActive() ){
      Json::Value dum;
      dum["nam"] = this->reg[i]->getName();
      dum["val"] = this->reg[i]->getValue();
      reg.append(dum);
    }
  }
  if( reg.size() != 0 ){
    all["reg"] = reg;
  }

  Json::Value lenses;
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
    lenses[lens_names[j]] = lens;
  }
  if( lenses.size() != 0 ){
    all["lenses"] = lenses;
  }

  return all;
}


void StandardLikelihood::initializeAlgebra(ImagePlane* image,BaseSourcePlane* source){
  this->algebra->setAlgebraInit(image,source);
}



void StandardLikelihood::updateLikelihoodModel(ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection){
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
    source->constructH();
  }

  source->constructL(image,mycollection);
  this->algebra->setAlgebraRuntime(image,source,Nlpar::getValueByName("lambda",this->reg));
  this->algebra->solveLinearSparseS(image,source);
}

double StandardLikelihood::getLogLike(ImagePlane* image,BaseSourcePlane* source){
  double pi  = 3.14159265358979323846;

  this->terms["Nilog2p"] = -(image->lookup.size()*log10(2*pi)/2.0); // for some reason i need the outer parenthsis here, otherwise there is amemory problem
  this->terms["Nslogl"]  = source->Sm*log10(Nlpar::getValueByName("lambda",this->reg))/2.0;
  double val = terms["chi2"] + Nlpar::getValueByName("lambda",this->reg)*terms["reg"] + terms["Nilog2p"] + terms["Nslogl"] + terms["detC"] + terms["detHtH"] + terms["detA"];
  this->terms["like"] = val;

  return val;
  //  return pow(nlpars[2]["b"]->val-6,2) + pow(nlpars[2]["q"]->val-6,2) + pow(nlpars[2]["pa"]->val-6,2) + pow(nlpars[2]["x0"]->val-6,2) + pow(nlpars[2]["y0"]->val-6,2);
}


void StandardLikelihood::outputLikelihoodModel(ImagePlane* image,BaseSourcePlane* source,std::string output){
  // Output reconstructed source
  source->outputSource(output + "vkl_voronoi.dat");
  
  // Output errors of reconstructed source
  double* errors = (double*) calloc(source->Sm,sizeof(double));
  this->algebra->getSourceErrors(source->Sm,errors);
  source->outputSourceErrors(errors,output + "vkl_voronoi_errors.dat");
  free(errors);

  // Create mock data (lensed MAP source)
  ImagePlane mock_data(image->Nj,image->Ni,image->width,image->height);
  this->algebra->getMockData(&mock_data,source);
  
  // Output image of the model (model)
  mock_data.writeImage(output + "vkl_image.fits");
  
  // Output residual image (diference between mydata and mock_data)
  ImagePlane residual(image->Nj,image->Ni,image->width,image->height);
  for(int i=0;i<image->Ni;i++){
    for(int j=0;j<image->Nj;j++){
      residual.img[i*image->Nj+j] = image->img[i*image->Nj+j] - mock_data.img[i*image->Nj+j];
    }
  }
  residual.writeImage(output + "vkl_residual.fits");
  



  // Output various quantities in json format
  Json::Value json_output;
  
  // Print and write to json the MAP parameters from the parameter model
  Json::Value pars;
  std::vector<std::string> all_names = this->getFullNames();
  std::vector<double> all_values = this->getValues();
  for(int i=0;i<all_names.size();i++){
    std::cout << all_names[i] << " " << all_values[i] << std::endl;    
    pars[all_names[i]] = all_values[i];
  }
  json_output["full_pars"] = pars;
  
  Json::Value full_active;
  std::vector<std::string> names = this->getActiveFullNames();
  std::vector<double> values = this->getActiveValues();
  for(int i=0;i<names.size();i++){
    full_active[names[i]] = values[i];
  }
  json_output["full_active"] = full_active;
  
  // Write structured active parameter names and values
  Json::Value json_active = this->getActiveNamesValues();
  json_output["json_active"] = json_active;




  // Write generic parameters
  Json::Value other;
  other["Nsource"] = source->Sm;
  other["Ndata"]   = image->Nm;
  other["Nmask"]   = static_cast<int>( image->lookup.size() );
  other["Psize"]   = image->width/image->Ni;
  json_output["generic"] = other;

  /*
  // Write likelihood terms
  Json::Value terms;
  terms["ED"]        =  pcomp->chi2/2.0;
  terms["ES"]        =  pcomp->reg/2.0;
  terms["logdetA"]   = -pcomp->detA/2.0;
  terms["logdetHtH"] =  pcomp->detHtH/2.0;
  terms["logdetC"]   =  pcomp->detC/2.0;
  if( source->reg == "covariance_kernel" ){
    // do nothing
  } else {
    terms["logl"]      =  source->Sm*log10(nlpars[1]["lambda"]->val)/2.0;
  }
  double pi          =  3.14159265358979323846;
  terms["log2pi"]    = -image->Nm*log10(2*pi)/2.0;
  json_output["terms"] = terms;
  */


  std::ofstream jsonfile(output+"vkl_output.json");
  jsonfile << json_output;
  jsonfile.close();
}



//Derived class from BaseLikelihoodModel: SourceCovarianceKernel
//===============================================================================================================
double SourceCovarianceKernel::getLogLike(ImagePlane* image,BaseSourcePlane* source){
  /*
  double pi  = 3.14159265358979323846;

  getMin(image,source,pcomp);
  double g   = pcomp->chi2/2. + pcomp->reg/2.;
  double f1  = image->lookup.size()*log10(2*pi)/2.;
  //double f2  = source->Sm*log10(nlpars[1]["lambda"]->val)/2.;
  double val = -g -f1 +(pcomp->detC)/2. +(pcomp->detHtH)/2. -(pcomp->detA)/2.;
  
  return val;
  */
  //  return pow(nlpars[2]["b"]->val-6,2) + pow(nlpars[2]["q"]->val-6,2) + pow(nlpars[2]["pa"]->val-6,2) + pow(nlpars[2]["x0"]->val-6,2) + pow(nlpars[2]["y0"]->val-6,2);
  return 0;
}




//Factory of BaseLikelihoodModel 
//===============================================================================================================
std::vector<Nlpar*> FactoryLikelihoodModel::nlparsFromJsonVector(const Json::Value myjson){
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


