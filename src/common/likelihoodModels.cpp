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
SmoothLikelihood::SmoothLikelihood(std::vector<Nlpar*> a,std::vector<Nlpar*> b,std::vector< std::vector<Nlpar*> > d,std::vector<std::string> e,ImagePlane* f,BaseSourcePlane* g,CollectionMassModels* h){
  name = "smooth";
  physical = a;
  reg = b;
  lenses = d;
  lens_names = e;
  image = f;
  source = g;
  collection = h;
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

  this->means = (double*) malloc(this->active.size()*sizeof(double));
  this->sdevs = (double*) malloc(this->active.size()*sizeof(double));
  this->bests = (double*) malloc(this->active.size()*sizeof(double));
  this->maps  = (double*) malloc(this->active.size()*sizeof(double));

  terms["chi2"]    = 0.0;
  terms["reg"]     = 0.0;
  terms["Nilog2p"] = 0.0;
  terms["Nslogl"]  = 0.0;
  terms["detC"]    = 0.0;
  terms["detCs"]   = 0.0;
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


//virtual,
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
  // Set the new smooth lens, physical, and source covariance parameters (lambda is set in SmoothLikelihood->reg)
  for(int i=0;i<this->collection->models.size();i++){
    if( this->collection->models[i]->type != "pert" ){
      this->collection->models[i]->setMassPars(this->lenses[i]);
    }
  }
  this->collection->setPhysicalPars(this->physical);
  if( this->source->sample_reg ){
    this->source->kernel->setParameters(this->reg);
  }

  // Deflect the rays in the updated mass model
  this->collection->all_defl(this->image);

  // Set the adaptive source in the updated mass model
  if( this->source->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(this->source);
    ada->createAdaGrid(this->image,this->collection);
    ada->createDelaunay();
  }

  // Get the new interpolation weights of the rays in the source grid and construct the lensing matrix
  this->source->createInterpolationWeights(this->image);
  this->source->constructL(this->image);

  // Re-construct the source regularization matrix
  if( this->source->sample_reg || this->source->type == "adaptive" ){
    this->source->constructH();
  }

  // If lambda is allowed to vary then update the lambda term
  Nlpar* lambda = Nlpar::getParByName("lambda",this->reg);
  if( lambda->fix == 0 ){
    this->terms["Nslogl"] = source->Sm*log10(lambda->val)/2.0;
  }

  // Update all the needed algebraic tables, e.g. M, Mt, Mt*Cd*M + l*Cs, Cs and detCs (if needed)
  this->algebra->setAlgebraRuntime(this->image,this->source);

  // Solve for the source, calculate the chi2, reg, and detA terms
  this->algebra->solveSource(this->source);
}

//virtual
double SmoothLikelihood::getLogLike(){
  double val = terms["chi2"] + terms["reg"] + terms["Nilog2p"] + terms["Nslogl"] + terms["detC"] + terms["detCs"] + terms["detA"];
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



  /*  
  // Write linear Dpsi
  Pert* pert_mass_model = new Pert(20,20,image->width,image->height,"gradient");
  //  Pert* pert_mass_model = new Pert(image->Ni,image->Nj,image,"gradient");
  ImagePlane* img_grid = new ImagePlane(pert_mass_model->dpsi->Si,pert_mass_model->dpsi->Sj,pert_mass_model->dpsi->width,pert_mass_model->dpsi->height);
  for(int i=0;i<img_grid->Nm;i++){
    img_grid->active[i] = -1;
  }
  this->collection->all_defl(img_grid);
  this->source->createInterpolationWeights(img_grid);
  this->source->constructL(img_grid);

  this->deriveLinearDpsi(pert_mass_model,img_grid);
  //pert_mass_model->dpsi->outputSource(output + "linear_perturbations_");


  FILE* fh = fopen((output + "linear_perturbations_derivatives.dat").c_str(),"w");
  for(int i=0;i<pert_mass_model->dpsi->Sm;i++){
    fprintf(fh,"%10.5f\n",pert_mass_model->dpsi_dx[i]);
  }
  for(int i=0;i<pert_mass_model->dpsi->Sm;i++){
    fprintf(fh,"%10.5f\n",pert_mass_model->dpsi_dy[i]);
  }
  fclose(fh);

  delete(pert_mass_model);
  delete(img_grid);
  */
}


void SmoothLikelihood::deriveLinearDpsi(Pert* pert_mass_model,ImagePlane* img_grid){
  // calculate the residuals (smooth modelling residuals)
  this->getModel();
  this->getResidual();

  // calculate the derivative of the source
  this->source->constructDs(img_grid);

  for(int i=0;i<pert_mass_model->dpsi->Sm;i++){
    if( source->Ds.tri[2*i].v == 0.0 ){
      pert_mass_model->dpsi_dx[i] = 0.0;
    } else {
      pert_mass_model->dpsi_dx[i] = -0.5*this->res->img[i]/source->Ds.tri[2*i].v;
    }

    if( source->Ds.tri[2*i+1].v == 0.0 ){
      pert_mass_model->dpsi_dy[i] = 0.0;
    } else {
      pert_mass_model->dpsi_dy[i] = -0.5*this->res->img[i]/source->Ds.tri[2*i+1].v;
    }
  }


  //  this->algebra->solveDpsi(pert_mass_model,this->image,this->source);
}




//Derived class from BaseLikelihoodModel: PertLikelihood
//===============================================================================================================
PertLikelihood::PertLikelihood(std::vector<Nlpar*> reg_s,std::vector<Nlpar*> reg_dpsi,std::string a,std::string b,BaseSourcePlane* c,Pert* d){
  this->name = "pert";
  this->reg_s = reg_s;
  this->reg_dpsi = reg_dpsi;
  this->reg_s_type = a;
  this->reg_dpsi_type = b;
  this->source = c->clone(); // This is the right way to copy a derived class, e.g. "this->source = new BaseSourcePlane(*c)" does not work, and "this->source = c" only copies the pointer.
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

  this->means = (double*) malloc(this->active.size()*sizeof(double));
  this->sdevs = (double*) malloc(this->active.size()*sizeof(double));
  this->bests = (double*) malloc(this->active.size()*sizeof(double));
  this->maps  = (double*) malloc(this->active.size()*sizeof(double));

  terms["chi2"]    = 0.0; // (Mr*r-d)^T C_d^-1 (Mr*r-d)/2
  terms["reg"]     = 0.0; // r^T (R^T R) r /2
  terms["Nilog2p"] = 0.0; // -Nd*log(2pi)/2
  terms["Nslogls"] = 0.0; // Ns*log(l_s)/2
  terms["Nploglp"] = 0.0; // Npsi*log(l_psi)/2
  terms["detCd"]   = 0.0; // log(det C_d^-1)/2
  terms["detCs"]   = 0.0; // log(det C_s^-1)/2
  terms["detCp"]   = 0.0; // log(det C_dpsi^-1)/2
  terms["detA"]    = 0.0; // -log(detA)/2
  terms["like"]    = 0.0; // evidence
}

PertLikelihood::~PertLikelihood(){
  delete(algebra);
}

//non-virtual
void PertLikelihood::initializePert(SmoothLikelihood* smooth_like){
  // This function needs to be called after a SmoothLikelihood model has been run,
  // otherwise the "smooth_like->source" will be empty, and this class needs a s_0 to initialize properly.

  // Add pointer to smooth likelihood model
  this->smooth_like = smooth_like;


  // Use source from smooth likelihood model
  if( this->source->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(this->source);
    ada->createAdaGrid(this->smooth_like->image,this->smooth_like->collection);
    ada->createDelaunay();
  }
  // Create lensing matrix for the new source used here
  this->smooth_like->collection->all_defl(this->smooth_like->image);
  this->source->createInterpolationWeights(this->smooth_like->image);
  this->source->constructL(this->smooth_like->image);
  // Copy the actual smooth source values, otherwise "source->src" is empty
  for(int i=0;i<this->source->Sm;i++){
    this->source->src[i] = this->smooth_like->source->src[i];
  }
  // Get derivative and regularization matrix
  this->source->constructDs(this->smooth_like->image);
  this->source->constructH();


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
  // Update covariance kernel for the source (f needed)
  // For identity,gradient, and curvature, (no sampling of the regularization by definition) H stays the same at every call and lambda values are accessed in setAlgebraRuntime
  if( this->source->sample_reg ){
    this->source->kernel->setParameters(this->reg_s);
    this->source->constructH();
  }

  // Update covariance kernel for the perturbations (if needed)
  // For identity,gradient, and curvature, (no sampling of the regularization by definition) H stays the same at every call and lambda values are accessed in setAlgebraRuntime
  if( this->pert_mass_model->dpsi->sample_reg ){
    this->pert_mass_model->dpsi->kernel->setParameters(this->reg_dpsi);
    this->pert_mass_model->dpsi->constructH();
  }

  // If lambda is allowed to vary then update the lambda term for the source
  Nlpar* lambda_s = Nlpar::getParByName("lambda_s",this->reg_s);
  if( lambda_s->fix == 0 ){
    this->terms["Nslogls"] = this->source->Sm*log10(lambda_s->val)/2.0;
  }  
  
  // If lambda is allowed to vary then update the lambda term for the potential corrections
  Nlpar* lambda_dpsi = Nlpar::getParByName("lambda_dpsi",this->reg_dpsi);
  if( lambda_dpsi->fix == 0 ){
    this->terms["Nploglp"] = this->pert_mass_model->dpsi->Sm*log10(lambda_dpsi->val)/2.0;
  }

  // Update all the needed algebraic tables, e.g. M_r, Mt_r, Mt_r*Cd*M_r + RtR, Cs and detCs (if needed), Cp and detCp (if needed)
  this->algebra->setAlgebraRuntime(this->source,this->pert_mass_model);

  // Solve for the source and perturbations (r), calculate the chi2, reg, and detA terms
  this->algebra->solveSourcePert(this->source,this->pert_mass_model);
}

//virtual
double PertLikelihood::getLogLike(){
  double val = terms["chi2"] + terms["reg"] + terms["Nilog2p"] + terms["Nslogls"] + terms["Nploglp"] + terms["detCd"] + terms["detCs"] + terms["detCp"] + terms["detA"];
  this->terms["like"] = val;
  return val;
}

//virtual
void PertLikelihood::initialOutputLikelihoodModel(std::string output){
  std::cout << "Sample_reg for source: " << this->source->sample_reg << " " << this->source->reg <<  std::endl;
  std::cout << "Sample_reg for perturbations: " << this->pert_mass_model->dpsi->sample_reg << " " << this->pert_mass_model->dpsi->reg << std::endl;
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
  
  // Output errors of reconstructed source
  double* errors = (double*) calloc(this->source->Sm,sizeof(double));
  this->smooth_like->algebra->getSourceErrors(this->source->Sm,errors);
  this->smooth_like->source->outputSourceErrors(errors,output);
  free(errors);

  // Replace source pointer in the smooth model with the previous source
  this->smooth_like->source = dum_source;
  
  // Output current perturbations correction
  this->pert_mass_model->dpsi->outputSource(output + "perturbations_");

  // Output convergence
  ImagePlane* kappa = new ImagePlane(this->pert_mass_model->dpsi->Si,this->pert_mass_model->dpsi->Sj,this->pert_mass_model->dpsi->width,this->pert_mass_model->dpsi->height);
  pert_mass_model->getConvergence(kappa);
  kappa->writeImage(output + "convergence.fits");
  delete(kappa);
}


void PertLikelihood::pertResiduals(std::string output,ImagePlane* image,BaseSourcePlane* source,Pert* pert_mass_model){
  this->algebra->constructDsDpsi(image,source,pert_mass_model);
  this->algebra->solvePerturbationsResiduals(output,pert_mass_model);
}








//Derived class from BaseLikelihoodModel: PertIterationLikelihood
//===============================================================================================================
PertIterationLikelihood::PertIterationLikelihood(std::vector<Nlpar*> reg_s,std::vector<Nlpar*> reg_dpsi,std::string a,std::string b,BaseSourcePlane* c,Pert* d,CollectionMassModels* e) : PertLikelihood(reg_s,reg_dpsi,a,b,c,d) {
  this->name = "pert_iter";
  this->collection = e;
}

//non-virtual
void PertIterationLikelihood::initializePert(SmoothLikelihood* smooth_like){
  // Add pointer to smooth likelihood model
  this->smooth_like = smooth_like;


  // Use source from smooth likelihood model
  if( this->source->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(this->source);
    ada->createAdaGrid(this->smooth_like->image,this->smooth_like->collection);
    ada->createDelaunay();
  }
  // Create lensing matrix for the new source used here
  this->smooth_like->collection->all_defl(this->smooth_like->image);
  this->source->createInterpolationWeights(this->smooth_like->image);
  this->source->constructL(this->smooth_like->image);
  // Copy the actual smooth source values, otherwise "source->src" is empty
  for(int i=0;i<this->source->Sm;i++){
    this->source->src[i] = this->smooth_like->source->src[i];
  }
  // Get derivative and regularization matrix
  this->source->constructDs(this->smooth_like->image);
  this->source->constructH();


  // Initialize Pert
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
  // Update covariance kernel for the source (f needed)
  if( this->source->sample_reg ){
    this->source->kernel->setParameters(this->reg_s);
  }
  // Construct a new H for the source if needed
  if( this->source->type == "adaptive" || this->source->sample_reg ){
    this->source->constructH();
  }
  this->source->constructDs(this->smooth_like->image);


  // Update covariance kernel for the perturbations (if needed)
  // For identity,gradient, and curvature, (no sampling of the regularization by definition) H stays the same at every call and lambda values are accessed in setAlgebraRuntime
  if( this->pert_mass_model->dpsi->sample_reg ){
    this->pert_mass_model->dpsi->kernel->setParameters(this->reg_dpsi);
    this->pert_mass_model->dpsi->constructH();
  }


  // Update all the needed algebraic tables, e.g. M_r, Mt_r, Mt_r*Cd*M_r + RtR, Cs and detCs (if needed), Cp and detCp (if needed)
  this->algebra->setAlgebraRuntime(this->source,this->pert_mass_model);

  // Solve for the source and perturbations (r), calculate the chi2, reg, and detA terms
  this->algebra->solveSourcePert(this->source,this->pert_mass_model);

  
  // Update perturbations in mass model collection
  Pert* pert_pointer = dynamic_cast<Pert*>(this->collection->models.back());
  pert_pointer->addDpsi(this->pert_mass_model->dpsi->src);
  pert_pointer->updatePert();
}

//virtual
double PertIterationLikelihood::getLogLike(){
  double val = terms["chi2"] + terms["reg"] + terms["Nilog2p"] + terms["Nslogls"] + terms["Nploglp"] + terms["detCd"] + terms["detCs"] + terms["detCp"] + terms["detA"];
  this->terms["like"] = val;
  return val;
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
  
  // Output errors of reconstructed source
  double* errors = (double*) calloc(this->source->Sm,sizeof(double));
  this->smooth_like->algebra->getSourceErrors(this->source->Sm,errors);
  this->smooth_like->source->outputSourceErrors(errors,output);
  free(errors);

  // Replace source pointer in the smooth model with the previous source
  this->smooth_like->source = dum_source;

  // Output last perturbations (corrections)
  this->pert_mass_model->dpsi->outputSource(output + "perturbations_");

  // Output added perturbations
  Pert* pert_pointer = dynamic_cast<Pert*>(this->collection->models.back());
  pert_pointer->dpsi->outputSource(output + "added_perturbations_"); // output additive perturbations

  // Output convergence
  ImagePlane* kappa = new ImagePlane(pert_pointer->dpsi->Si,pert_pointer->dpsi->Sj,pert_pointer->dpsi->width,pert_pointer->dpsi->height);
  pert_pointer->getConvergence(kappa);
  kappa->writeImage(output + "added_convergence.fits");
  delete(kappa);
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


