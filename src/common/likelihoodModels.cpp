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
//non-virtual
std::vector<double> BaseLikelihoodModel::getActiveValues(){
  std::vector<double> values;
  for(int i=0;i<this->active.size();i++){
    values.push_back( this->active[i]->val );
  }
  return values;
}

//non-virtual
std::vector<std::string> BaseLikelihoodModel::getActiveNames(){
  std::vector<std::string> names;
  for(int i=0;i<this->active.size();i++){
    names.push_back( this->active[i]->nam );
  }
  return names;
}

//non-virtual
void BaseLikelihoodModel::updateActive(std::vector<double> values){
  for(int i=0;i<this->active.size();i++){
    this->active[i]->val = values[i];
  }
}

//non-virtual
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

//non-virtual
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

  model = new ImagePlane(*this->image);
  res = new ImagePlane(*this->image);

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
  delete(model);
  delete(res);
}


//non-virtual
std::vector<Nlpar*> StandardLikelihood::getRegPars(){
  return this->reg;
}

//non-virtual
void StandardLikelihood::getModel(){
  this->algebra->getMockData(this->model,this->source);
}

//non-virtual
void StandardLikelihood::getResidual(){
  for(int i=0;i<this->image->Nm;i++){
    this->res->img[i] = this->image->img[i] - this->model->img[i];
  }
}


//virtual
std::vector<Nlpar*> StandardLikelihood::getPhysicalPars(){
  return this->physical;
}

//virtual
std::vector<Nlpar*> StandardLikelihood::getMassModelPars(int i){
  return this->lenses[i];
}

//virtual
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

//virtual
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

//virtual
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

//virtual
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



//virtual
void StandardLikelihood::initializeAlgebra(){
  this->algebra->setAlgebraInit(this->image,this->source);
}

//virtual
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
  this->source->createInterpolationWeights(this->image);
  this->source->constructL(this->image);
  if( this->source->sample_reg || this->source->type == "adaptive" ){
    this->source->constructH();
  }

  this->algebra->setAlgebraRuntime(this->image,this->source,Nlpar::getValueByName("lambda",this->reg));
  this->algebra->solveLinearSparseS(this->image,this->source);
}

//virtual
double StandardLikelihood::getLogLike(){
  double pi  = 3.14159265358979323846;

  this->terms["Nilog2p"] = -(this->image->lookup.size()*log10(2*pi)/2.0); // for some reason i need the outer parenthsis here, otherwise there is a memory problem
  this->terms["Nslogl"]  = this->source->Sm*log10(Nlpar::getValueByName("lambda",this->reg))/2.0;
  double val = terms["chi2"] + Nlpar::getValueByName("lambda",this->reg)*terms["reg"] + terms["Nilog2p"] + terms["Nslogl"] + terms["detC"] + terms["detHtH"] + terms["detA"];
  this->terms["like"] = val;

  return val;
}

//virtual
void StandardLikelihood::initialOutputLikelihoodModel(std::string output){
  // File with the parameter names
  std::ofstream f_names(output+"plt_corner.paramnames",std::ofstream::out);
  std::vector<std::string> active_full_names = this->getActiveFullNames();
  for(int i=0;i<active_full_names.size();i++){
    f_names << active_full_names[i];
    f_names << " ";
    f_names << active_full_names[i];
    f_names << std::endl;
  }
  f_names.close();

  // File with the parameter ranges
  std::ofstream f_ranges(output+"plt_corner.ranges",std::ofstream::out);
  for(int i=0;i<this->active.size();i++){
    f_ranges << active_full_names[i];
    f_ranges << " ";
    f_ranges << this->active[i]->min;
    f_ranges << " ";
    f_ranges << this->active[i]->max;
    f_ranges << std::endl;
  }
  f_ranges.close();
}



//virtual
void StandardLikelihood::outputLikelihoodModel(std::string output){
  // Output reconstructed source
  this->source->outputSource(output);

  // Output errors of reconstructed source
  double* errors = (double*) calloc(this->source->Sm,sizeof(double));
  this->algebra->getSourceErrors(this->source->Sm,errors);
  this->source->outputSourceErrors(errors,output);
  free(errors);

  // Create mock data (lensed MAP source)
  this->getModel();
  this->model->writeImage(output + "vkl_image.fits");

  // Output residual image (diference between mydata and mock_data)
  this->getResidual();
  this->res->writeImage(output + "vkl_residual.fits");
  



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

//virtual
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

  // Initialize Pert
  this->pert_mass_model->createAint(this->image);
  this->pert_mass_model->dpsi->constructH();

  // Add additive perturbations to the mass collection
  Pert* additive_pert = new Pert(*this->pert_mass_model);
  double* zeros = (double*) calloc(this->pert_mass_model->dpsi->Sm,sizeof(double));
  additive_pert->replaceDpsi(zeros);
  free(zeros);
  additive_pert->updatePert();
  this->collection->models.push_back(additive_pert);
  
  // Get residuals of the smooth model
  StandardLikelihood* specific_pointer = dynamic_cast<StandardLikelihood*>(this->smooth_like);
  specific_pointer->getModel();
  specific_pointer->getResidual();  

  this->initializeAlgebra();
}

//virtual
void PerturbationsLikelihood::initializeAlgebra(){
  StandardLikelihood* specific_pointer = dynamic_cast<StandardLikelihood*>(this->smooth_like);
  this->algebra->setAlgebraInit(this->image,this->pert_mass_model,specific_pointer->res->img);
}

//virtual
void PerturbationsLikelihood::updateLikelihoodModel(){
  this->smooth_like->updateLikelihoodModel();

  this->source->constructDs(this->image);
  if( this->pert_mass_model->dpsi->sample_reg ){
    this->pert_mass_model->dpsi->kernel->setParameters(this->reg);
    this->pert_mass_model->dpsi->constructH();
  }

  this->algebra->setAlgebraRuntime(this->image,this->source,this->pert_mass_model,this->smooth_like,Nlpar::getValueByName("lambda",this->reg));
  this->algebra->solvePert(this->image,this->pert_mass_model,this->smooth_like);

  // Update perturbations in mass model collection
  Pert* pert_pointer = dynamic_cast<Pert*>(this->collection->models.back());
  pert_pointer->addDpsi(this->pert_mass_model->dpsi->src);
  pert_pointer->updatePert();
}

//virtual
double PerturbationsLikelihood::getLogLike(){
  double like = 0.0;
  return like;
}

//virtual
void PerturbationsLikelihood::initialOutputLikelihoodModel(std::string output){
}

//virtual
void PerturbationsLikelihood::outputLikelihoodModel(std::string output){
  // Output reconstructed source
  this->source->outputSource(output + "pert_");

  // Update smooth model to include the latest potential corrections
  StandardLikelihood* like_pointer = dynamic_cast<StandardLikelihood*>(this->smooth_like);
  like_pointer->updateLikelihoodModel();

  // Output errors of reconstructed source
  double* errors = (double*) calloc(this->source->Sm,sizeof(double));
  like_pointer->algebra->getSourceErrors(this->source->Sm,errors);
  this->source->outputSourceErrors(errors,output + "pert_");
  free(errors);
  
  // Create mock data (lensed MAP source)
  like_pointer->getModel();
  like_pointer->model->writeImage(output + "pert_vkl_image.fits");
  
  // Output residual image (diference between mydata and mock_data)
  like_pointer->getResidual();
  like_pointer->res->writeImage(output + "pert_vkl_residual.fits");

  // Output perturbations
  Pert* pert_pointer = dynamic_cast<Pert*>(this->collection->models.back());
  pert_pointer->dpsi->outputSource(output + "perturbations_");

  // Ouput latest corrections
  //  pert_mass_model->dpsi->outputSource(output + "perturbations0_");
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


