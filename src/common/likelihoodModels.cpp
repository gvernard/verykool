#include "likelihoodModels.hpp"

#include <string>
#include <vector>
#include <map>

#include "json/json.h"

#include "nonLinearPars.hpp"
#include "imagePlane.hpp"
#include "massModels.hpp"
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
StandardLikelihood::StandardLikelihood(std::vector<Nlpar*> a,std::vector<Nlpar*> b,std::vector< std::vector<Nlpar*> > c,std::vector<std::string> d,ImagePlane* e,BaseSourcePlane* f,CollectionMassModels* g){
  physical = a;
  reg = b;
  lenses = c;
  lens_names = d;
  image = e;
  source = f;
  collection = g;
  this->algebra = new StandardAlgebra(this);

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


void StandardLikelihood::initializeAlgebra(){
  this->source->constructH();
  this->algebra->setAlgebraInit(this->image,this->source);
}



void StandardLikelihood::updateLikelihoodModel(){
  for(int i=0;i<this->collection->models.size();i++){
    if( this->collection->models[i]->type != "pert" ){
      this->collection->models[i]->setMassPars(this->lenses[i]);
    }
  }
  this->collection->setPhysicalPars(this->physical);

  if( this->source->sample_reg ){
    this->source->kernel->setParameters(this->reg);
  }

  this->collection->all_defl(this->image);

  if( this->source->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(this->source);
    ada->createAdaGrid(this->image,this->collection);
    ada->createDelaunay();
  }

  if( this->source->sample_reg || this->source->type == "adaptive" ){
    this->source->constructH();
  }

  this->source->createInterpolationWeights(this->image);
  this->source->constructL(this->image);
  this->algebra->setAlgebraRuntime(this->image,this->source,Nlpar::getValueByName("lambda",this->reg));
  this->algebra->solveLinearSparseS(this->image,this->source);
}

double StandardLikelihood::getLogLike(){
  double pi  = 3.14159265358979323846;

  this->terms["Nilog2p"] = -(this->image->lookup.size()*log10(2*pi)/2.0); // for some reason i need the outer parenthsis here, otherwise there is a memory problem
  this->terms["Nslogl"]  = this->source->Sm*log10(Nlpar::getValueByName("lambda",this->reg))/2.0;
  double val = terms["chi2"] + Nlpar::getValueByName("lambda",this->reg)*terms["reg"] + terms["Nilog2p"] + terms["Nslogl"] + terms["detC"] + terms["detHtH"] + terms["detA"];
  this->terms["like"] = val;

  return val;
}


void StandardLikelihood::outputLikelihoodModel(std::string output){
  // Output reconstructed source
  this->source->outputSource(output);

  // Output errors of reconstructed source
  double* errors = (double*) calloc(this->source->Sm,sizeof(double));
  this->algebra->getSourceErrors(this->source->Sm,errors);
  this->source->outputSourceErrors(errors,output);
  free(errors);

  // Create mock data (lensed MAP source)
  ImagePlane mock_data(this->image->Nj,this->image->Ni,this->image->width,this->image->height);
  this->algebra->getMockData(&mock_data,this->source);
  
  // Output image of the model (model)
  mock_data.writeImage(output + "vkl_image.fits");
  
  // Output residual image (diference between mydata and mock_data)
  ImagePlane residual(this->image->Nj,this->image->Ni,this->image->width,this->image->height);
  for(int i=0;i<this->image->Ni;i++){
    for(int j=0;j<this->image->Nj;j++){
      residual.img[i*this->image->Nj+j] = this->image->img[i*this->image->Nj+j] - mock_data.img[i*this->image->Nj+j];
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
  other["Nsource"] = this->source->Sm;
  other["Ndata"]   = this->image->Nm;
  other["Nmask"]   = static_cast<int>( this->image->lookup.size() );
  other["Psize"]   = this->image->width/this->image->Ni;
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

// void  SourceCovarianceKernel::updateLikelihoodModel(ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection){
//
//   SAME AS StandardLikelihood
//
//}

double SourceCovarianceKernel::getLogLike(){
  double pi  = 3.14159265358979323846;

  this->terms["Nilog2p"] = -(this->image->lookup.size()*log10(2*pi)/2.0); // for some reason i need the outer parenthsis here, otherwise there is a memory problem
  double val = terms["chi2"] + terms["reg"] + terms["Nilog2p"] + terms["detC"] + terms["detHtH"] + terms["detA"];
  this->terms["like"] = val;

  return val;
}



//Derived class from BaseLikelihoodModel: PerturbationsLikelihood
//===============================================================================================================
PerturbationsLikelihood::PerturbationsLikelihood(std::vector<Nlpar*> reg,ImagePlane* a,BaseSourcePlane* b,CollectionMassModels* c,Pert* d){
  this->reg = reg;
  this->image = a;
  this->source = b;
  this->collection = c;
  this->pert_mass_model = d;
  this->algebra = new PerturbationsAlgebra(this);

  for(int i=0;i<this->reg.size();i++){
    if( this->reg[i]->getActive() ){
      active.push_back( this->reg[i] );
    }
  }

  terms["A"] = 0.0;
  terms["B"] = 0.0;
  terms["C"] = 0.0;
}

PerturbationsLikelihood::~PerturbationsLikelihood(){
  delete(algebra);
}

//non-virtual
void PerturbationsLikelihood::initializePert(BaseLikelihoodModel* smooth_like){
  // Add pointer to smooth likelihood model
  this->smooth_like = smooth_like;

  // Add additive perturbations to the mass collection
  Pert* additive_pert = new Pert(*this->pert_mass_model);
  for(int i=0;i<this->pert_mass_model->dpsi->Sm;i++){
    additive_pert->dpsi->src[i] = 0.0;
  }
  StandardLikelihood* specific_pointer = dynamic_cast<StandardLikelihood*>(this->smooth_like);
  specific_pointer->collection->models.push_back(additive_pert);
  
  this->initializeAlgebra();
}

//virtual
void PerturbationsLikelihood::initializeAlgebra(){
  this->algebra->setAlgebraInit(this->image,this->pert_mass_model);
}

//virtual
void PerturbationsLikelihood::updateLikelihoodModel(){
  this->smooth_like->updateLikelihoodModel();

  this->source->constructDs(this->image);
  if( this->pert_mass_model->dpsi->sample_reg ){
    // update regularization parameters
    this->pert_mass_model->dpsi->constructH();
  }


  this->algebra->setAlgebraRuntime(this->image,this->source,this->pert_mass_model,this->smooth_like,Nlpar::getValueByName("lambda",this->reg));
  //  this->algebra->setAlgebraRuntime(,Nlpar::getValueByName("lambda",this->reg));
  //  this->algebra->solvePert();

  // Update perturbations in collection
}

//virtual
double PerturbationsLikelihood::getLogLike(){
  double like = 0.0;
  return like;
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


