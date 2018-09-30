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


//Derived class from BaseLikelihoodModel: SmoothLikelihood
//===============================================================================================================
SmoothLikelihood::SmoothLikelihood(std::vector<Nlpar*> a,std::vector<Nlpar*> b,std::vector< std::vector<Nlpar*> > c,std::vector<std::string> d,ImagePlane* e,BaseSourcePlane* f,CollectionMassModels* g){
  physical = a;
  reg = b;
  lenses = c;
  lens_names = d;
  image = e;
  source = f;
  collection = g;
  this->algebra = new SmoothAlgebra(this);

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


SmoothLikelihood::~SmoothLikelihood(){
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
std::vector<Nlpar*> SmoothLikelihood::getRegPars(){
  return this->reg;
}

//non-virtual
void SmoothLikelihood::getModel(){
  this->algebra->getMockData(this->model,this->source);
}

//non-virtual
void SmoothLikelihood::getResidual(){
  for(int i=0;i<this->image->Nm;i++){
    this->res->img[i] = this->image->img[i] - this->model->img[i];
  }
}


//virtual
std::vector<Nlpar*> SmoothLikelihood::getPhysicalPars(){
  return this->physical;
}

//virtual
std::vector<Nlpar*> SmoothLikelihood::getMassModelPars(int i){
  return this->lenses[i];
}

//virtual
std::vector<std::string> SmoothLikelihood::getFullNames(){
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
std::vector<std::string> SmoothLikelihood::getActiveFullNames(){
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
std::vector<double> SmoothLikelihood::getValues(){
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
Json::Value SmoothLikelihood::getActiveNamesValues(){
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
void SmoothLikelihood::initializeAlgebra(){
  this->algebra->setAlgebraInit(this->image,this->source);
}

//virtual
void SmoothLikelihood::updateLikelihoodModel(){
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
double SmoothLikelihood::getLogLike(){
  double pi  = 3.14159265358979323846;

  this->terms["Nilog2p"] = -(this->image->lookup.size()*log10(2*pi)/2.0); // for some reason i need the outer parenthsis here, otherwise there is a memory problem
  this->terms["Nslogl"]  = this->source->Sm*log10(Nlpar::getValueByName("lambda",this->reg))/2.0;
  double val = terms["chi2"] + Nlpar::getValueByName("lambda",this->reg)*terms["reg"] + terms["Nilog2p"] + terms["Nslogl"] + terms["detC"] + terms["detHtH"] + terms["detA"];
  this->terms["like"] = val;

  return val;
}

//virtual
void SmoothLikelihood::initialOutputLikelihoodModel(std::string output){
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
void SmoothLikelihood::outputLikelihoodModel(std::string output){
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


  // Write likelihood terms
  Json::Value terms;
  std::map<std::string,double>::iterator it;
  for(it=this->terms.begin();it!=this->terms.end();it++){
    terms[it->first] = it->second;
  }
  json_output["terms"] = terms;


  std::ofstream jsonfile(output+"vkl_output.json");
  jsonfile << json_output;
  jsonfile.close();
}



//Derived class from BaseLikelihoodModel: SourceCovarianceKernel
//===============================================================================================================

// void  SourceCovarianceKernel::updateLikelihoodModel(ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection){
//
//   SAME AS SmoothLikelihood
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



//Derived class from BaseLikelihoodModel: PertLikelihood
//===============================================================================================================
PertLikelihood::PertLikelihood(std::vector<Nlpar*> reg_s,std::vector<Nlpar*> reg_dpsi,std::string a,std::string b,BaseSourcePlane* c,Pert* d){
  this->reg_s = reg_s;
  this->reg_dpsi = reg_dpsi;
  this->reg_s_type = a;
  this->reg_dpsi_type = b;
  this->source = c; // copy the source here
  this->source->reg = reg_s_type;
  this->pert_mass_model = d;
  this->algebra = new PertAlgebra(this);

  for(int i=0;i<this->reg_s.size();i++){
    if( this->reg_s[i]->getActive() ){
      active.push_back( this->reg_s[i] );
    }
  }
  for(int i=0;i<this->reg_dpsi.size();i++){
    if( this->reg_dpsi[i]->getActive() ){
      active.push_back( this->reg_dpsi[i] );
    }
  }

  // dum likelihood terms for the moment
  terms["A"] = 0.0;
  terms["B"] = 0.0;
  terms["C"] = 0.0;
}

PertLikelihood::~PertLikelihood(){
  delete(algebra);
}

//non-virtual
void PertLikelihood::initializePert(SmoothLikelihood* smooth_like){
  // Add pointer to smooth likelihood model
  this->smooth_like = smooth_like;

  // Use source from smooth likelihood model
  this->smooth_like->source->constructDs(this->smooth_like->image);
  this->smooth_like->source->constructH();

  // Initialize Pert
  this->pert_mass_model->createCrosses(this->smooth_like->image);
  this->pert_mass_model->dpsi->constructH();

  this->initializeAlgebra();  
}

//virtual
void PertLikelihood::initializeAlgebra(){
  this->algebra->setAlgebraInit(this->source,this->pert_mass_model);
}

//virtual
void PertLikelihood::updateLikelihoodModel(){
  
  // Update covariance kernel for the source.
  // For identity,gradient, and curvature, H stays the same at every call and lambda values are accessed in setAlgebraRuntime
  if( this->source->sample_reg && this->source->reg == "covariance" ){
    this->source->kernel->setParameters(this->reg_s);
    this->source->constructH();
  }

  // Update covariance kernel for the perturbations
  // For identity,gradient, and curvature, H stays the same at every call and lambda values are accessed in setAlgebraRuntime
  if( this->pert_mass_model->dpsi->sample_reg && this->pert_mass_model->dpsi->reg == "covariance" ){
    this->pert_mass_model->dpsi->kernel->setParameters(this->reg_dpsi);
    this->pert_mass_model->dpsi->constructH();
  }
  
  this->algebra->setAlgebraRuntime(this->source,this->pert_mass_model,this->smooth_like);
  this->algebra->solveSourcePert(this->source,this->pert_mass_model,this->smooth_like);
}

//virtual
double PertLikelihood::getLogLike(){
  double like = 0.0;
  return like;
}

//virtual
void PertLikelihood::initialOutputLikelihoodModel(std::string output){
}

//virtual
void PertLikelihood::outputLikelihoodModel(std::string output){
  // Output reconstructed source
  this->source->outputSource(output + "pert_");


  // Replace source pointer in the smooth model with the new perturbed source
  BaseSourcePlane* dum_source = this->smooth_like->source;
  this->smooth_like->source = this->source;

  // Create mock data (lensed MAP source)
  this->smooth_like->getModel();
  this->smooth_like->model->writeImage(output + "pert_vkl_image.fits");
  
  // Output residual image (diference between mydata and mock_data)
  this->smooth_like->getResidual();
  this->smooth_like->res->writeImage(output + "pert_vkl_residual.fits");
  
  // Replace source pointer in the smooth model with the previous source
  this->smooth_like->source = dum_source;

  
  // Output perturbations
  this->pert_mass_model->dpsi->outputSource(output + "perturbations_");
  ImagePlane* kappa = new ImagePlane(this->pert_mass_model->dpsi->Si,this->pert_mass_model->dpsi->Sj,this->pert_mass_model->dpsi->width,this->pert_mass_model->dpsi->height);
  for(int i=0;i<kappa->Nm;i++){
    kappa->cells[i] = NULL;
    kappa->crosses[i] = NULL;
  }
  pert_mass_model->getConvergence(kappa);
  kappa->writeImage(output + "convergence.fits");
  delete(kappa);
  
  // Ouput latest corrections
  pert_mass_model->dpsi->outputSource(output + "perturbations0_");
}











//Derived class from BaseLikelihoodModel: PertIterationLikelihood
//===============================================================================================================
PertIterationLikelihood::PertIterationLikelihood(std::vector<Nlpar*> reg_s,std::vector<Nlpar*> reg_dpsi,std::string a,std::string b,BaseSourcePlane* c,Pert* d,CollectionMassModels* e){
  this->reg_s = reg_s;
  this->reg_dpsi = reg_dpsi;
  this->reg_s_type = a;
  this->reg_dpsi_type = b;
  this->source = c; // copy the source here
  this->source->reg = reg_s_type;
  this->pert_mass_model = d;
  this->collection = e;
  this->algebra = new PertIterationAlgebra(this);

  for(int i=0;i<this->reg_s.size();i++){
    if( this->reg_s[i]->getActive() ){
      active.push_back( this->reg_s[i] );
    }
  }
  for(int i=0;i<this->reg_dpsi.size();i++){
    if( this->reg_dpsi[i]->getActive() ){
      active.push_back( this->reg_dpsi[i] );
    }
  }

  terms["A"] = 0.0;
  terms["B"] = 0.0;
  terms["C"] = 0.0;
}

PertIterationLikelihood::~PertIterationLikelihood(){
  delete(algebra);
}

//non-virtual
void PertIterationLikelihood::initializePert(SmoothLikelihood* smooth_like){
  // Add pointer to smooth likelihood model
  this->smooth_like = smooth_like;

  // Use source from smooth likelihood model
  this->smooth_like->source->constructDs(this->smooth_like->image);
  this->smooth_like->source->constructH();

  // Initialize Pert
  //  this->pert_mass_model->createAint(this->image);
  this->pert_mass_model->createCrosses(this->smooth_like->image);
  this->pert_mass_model->dpsi->constructH();

  
  // Add additive perturbations to the mass collection
  Pert* additive_pert = new Pert(this->pert_mass_model->dpsi->Si,this->pert_mass_model->dpsi->Sj,this->pert_mass_model->dpsi->width,this->pert_mass_model->dpsi->height,"identity");
  // 'additive_pert' is just a container, regularization is not used
  double* zeros = (double*) calloc(this->pert_mass_model->dpsi->Sm,sizeof(double));
  additive_pert->replaceDpsi(zeros);
  free(zeros);
  additive_pert->updatePert();
  this->collection->models.push_back(additive_pert);

  this->initializeAlgebra();
}

//virtual
void PertIterationLikelihood::initializeAlgebra(){
  this->algebra->setAlgebraInit(this->source,this->pert_mass_model);
}

//virtual
void PertIterationLikelihood::updateLikelihoodModel(){

  // Mass model has been updated in the previous step, so I need to deflect the rays and calculcate a new grid for an adaptive source
  this->collection->all_defl(this->smooth_like->image);
  if( this->source->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(this->source);
    ada->createAdaGrid(this->smooth_like->image,this->collection);
    ada->createDelaunay();
  }
  this->source->createInterpolationWeights(this->smooth_like->image);
  this->source->constructL(this->smooth_like->image);
  
  // Update covariance matrix for the source if needed
  if( this->source->sample_reg && this->source->reg == "covariance" ){
    this->source->kernel->setParameters(this->reg_s);
  }
  // Construct a new H for the source if needed
  if( this->source->type == "adaptive" || (this->source->sample_reg && this->source->reg == "covariance") ){
    this->source->constructH();
  }
  
  // Update covariance kernel for the perturbations
  // For identity,gradient, and curvature, H stays the same at every call and lambda values are accessed in setAlgebraRuntime
  if( this->pert_mass_model->dpsi->sample_reg && this->pert_mass_model->dpsi->reg == "covariance" ){
    this->pert_mass_model->dpsi->kernel->setParameters(this->reg_dpsi);
    this->pert_mass_model->dpsi->constructH();
  }
  
  this->algebra->setAlgebraRuntime(this->source,this->pert_mass_model,this->smooth_like);
  this->algebra->solveSourcePert(this->source,this->pert_mass_model,this->smooth_like);
  
  // Update perturbations in mass model collection
  Pert* pert_pointer = dynamic_cast<Pert*>(this->collection->models.back());
  pert_pointer->addDpsi(this->pert_mass_model->dpsi->src);
  pert_pointer->updatePert();
}

//virtual
double PertIterationLikelihood::getLogLike(){
  double like = 0.0;
  return like;
}

//virtual
void PertIterationLikelihood::initialOutputLikelihoodModel(std::string output){
}

//virtual
void PertIterationLikelihood::outputLikelihoodModel(std::string output){
  // Output reconstructed source
  this->source->outputSource(output + "pert_");

  
  // Replace source pointer in the smooth model with the new perturbed source
  BaseSourcePlane* dum_source = this->smooth_like->source;
  this->smooth_like->source = this->source;

  // Create mock data (lensed MAP source)
  this->smooth_like->getModel();
  this->smooth_like->model->writeImage(output + "pert_vkl_image.fits");
  
  // Output residual image (diference between mydata and mock_data)
  this->smooth_like->getResidual();
  this->smooth_like->res->writeImage(output + "pert_vkl_residual.fits");
  
  // Replace source pointer in the smooth model with the previous source
  this->smooth_like->source = dum_source;


  // Output perturbations
  Pert* pert_pointer = dynamic_cast<Pert*>(this->collection->models.back());
  //  Pert* pert_pointer = this->pert_mass_pointer;
  pert_pointer->dpsi->outputSource(output + "perturbations_");
  ImagePlane* kappa = new ImagePlane(pert_pointer->dpsi->Si,pert_pointer->dpsi->Sj,pert_pointer->dpsi->width,pert_pointer->dpsi->height);
  for(int i=0;i<kappa->Nm;i++){
    kappa->cells[i] = NULL;
    kappa->crosses[i] = NULL;
  }
  pert_pointer->getConvergence(kappa);
  kappa->writeImage(output + "convergence.fits");
  delete(kappa);

  // Output additive corrections
  pert_pointer->dpsi->outputSource(output + "perturbations0_");
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


