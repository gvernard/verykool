#include "likelihoodModels.hpp"

#include <string>
#include <vector>
#include <map>
#include <unordered_map>

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
std::vector<std::string> BaseLikelihoodModel::getActiveNlparNames(){
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
    printf(" %12s",(this->active_names[i]).c_str());
    //printf(" %12s",(this->active[i]->nam).c_str());
  }
  printf("\n");
  for(int i=0;i<this->active.size();i++){
    printf(" %12.5f",this->active[i]->val);
  }
  printf("\n");
}

/*
//non-virtual
void BaseLikelihoodModel::printTerms(){
  for(std::unordered_map<std::string,double>::iterator it=this->terms.begin();it!=this->terms.end();it++){
    printf(" %16s",(it->first).c_str());
  }
  printf("\n");
  for(std::unordered_map<std::string,double>::iterator it=this->terms.begin();it!=this->terms.end();it++){
    printf(" %16.5f",it->second);
  }
  printf("\n");
}
*/

//non-virtual
double BaseLikelihoodModel::getLogLike(){
  this->terms["evidence"] = 0.0;
  double sum = 0.0;
  for(std::unordered_map<std::string,double>::iterator it=this->terms.begin();it!=this->terms.end();it++){
    sum += it->second;
  }
  this->terms["evidence"] = sum;
  return sum;
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
      this->active.push_back( this->physical[i] );
      this->active_names.push_back( this->physical[i]->nam );
    }
  }

  for(int i=0;i<this->reg.size();i++){
    if( this->reg[i]->getActive() ){
      this->active.push_back( this->reg[i] );
      this->active_names.push_back( this->reg[i]->nam );
    }
  }

  for(int j=0;j<this->lenses.size();j++){
    for(int i=0;i<this->lenses[j].size();i++){
      if( this->lenses[j][i]->getActive() ){
	this->active.push_back( this->lenses[j][i] );
	this->active_names.push_back( this->lens_names[j] + "_" + this->lenses[j][i]->nam );
      }
    }
  }

  this->maps    = (double*) malloc(this->active.size()*sizeof(double));
  this->means   = (double*) malloc(this->active.size()*sizeof(double));
  this->s1_low  = (double*) malloc(this->active.size()*sizeof(double));
  this->s1_high = (double*) malloc(this->active.size()*sizeof(double));

  terms["Nilog2pi"] = 0.0;
  terms["Nslogl"]   = 0.0;
  if( this->source->reg == "curvature_in_identity_out" || this->source->reg == "covariance_kernel_in_identity_out" ){
    terms["Noutlog2pi"] = 0.0;
  }
  terms["detCd"]    = 0.0;
  terms["detA"]     = 0.0;
  terms["detCs"]    = 0.0; 
  if( this->source->reg == "curvature_in_identity_out" || this->source->reg == "covariance_kernel_in_identity_out" ){
    terms["detLout"] = 0.0;
  }
  terms["chi2"]     = 0.0;
  terms["reg"]      = 0.0;
  if( this->source->reg == "curvature_in_identity_out" || this->source->reg == "covariance_kernel_in_identity_out" ){
    terms["reg_out"] = 0.0;
  }
  terms["evidence"] = 0.0;
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
void SmoothLikelihood::updateCovKernelLimits(){
  /*
  // Update limits of source covariance kernel, if applicable
  if( this->source->type == "fixed" ){
    if( this->source->kernel->type == "modgauss" || this->source->kernel->type == "gauss" ){
      Nlpar* nlpar_sdev = Nlpar::getParByName("sdev",this->active);
      FixedSource* fixed = dynamic_cast<FixedSource*>(this->source); // needed to access width
      nlpar_sdev->min = fixed->width/(this->source->Sj*log(this->source->kernel->cmax));
      nlpar_sdev->max = 50.0*this->image->width;
    }
  } else if( this->source->type == "adaptive" ){
    if( this->source->kernel->type == "modgauss" || this->source->kernel->type == "gauss" ){
      Nlpar* nlpar_sdev = Nlpar::getParByName("sdev",this->active);
      nlpar_sdev->max = 50.0*this->image->width;
    }
  }
  */
}

//non-virtual
std::vector<Nlpar*> SmoothLikelihood::getRegPars(){
  return this->reg;
}

//virtual
void SmoothLikelihood::getModel(){
  this->algebra->getMockData(this->model,this->source);
}

//virtual
void SmoothLikelihood::getResidual(){
  for(int i=0;i<this->image->Nm;i++){
    this->res->img[i] = this->image->img[i] - this->model->img[i];
  }
}

//virtual
void SmoothLikelihood::printTerms(){
  if( this->source->reg == "curvature_in_identity_out" || this->source->reg == "covariance_kernel_in_identity_out" ){
    printf(" %16s","Nilog2pi");
    printf(" %16s","Nslogl");
    printf(" %16s","Noutlog2pi");
    printf(" %16s","detCd");
    printf(" %16s","detA");
    printf(" %16s","detCs");
    printf(" %16s","detLout");
    printf(" %16s","chi2");
    printf(" %16s","reg");
    printf(" %16s","reg_out");
    printf(" %16s","evidence");
    printf("\n");
    printf(" %16f",terms["Nilog2pi"]);
    printf(" %16f",terms["Nslogl"]);
    printf(" %16f",terms["Noutlog2pi"]);
    printf(" %16f",terms["detCd"]);
    printf(" %16f",terms["detA"]);
    printf(" %16f",terms["detCs"]);
    printf(" %16f",terms["detLout"]);
    printf(" %16f",terms["chi2"]);
    printf(" %16f",terms["reg"]);
    printf(" %16f",terms["reg_out"]);
    printf(" %16f",terms["evidence"]);
    printf("\n");
  } else {
    printf(" %16s","Nilog2pi");
    printf(" %16s","Nslogl");
    printf(" %16s","detCd");
    printf(" %16s","detA");
    printf(" %16s","detCs");
    printf(" %16s","chi2");
    printf(" %16s","reg");
    printf(" %16s","evidence");
    printf("\n");
    printf(" %16f",terms["Nilog2pi"]);
    printf(" %16f",terms["Nslogl"]);
    printf(" %16f",terms["detCd"]);
    printf(" %16f",terms["detA"]);
    printf(" %16f",terms["detCs"]);
    printf(" %16f",terms["chi2"]);
    printf(" %16f",terms["reg"]);
    printf(" %16f",terms["evidence"]);
    printf("\n");
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
void SmoothLikelihood::getAllNamesValues(std::vector<std::string>& names,std::vector<double>& values){
  std::map<std::string,int> active_ind;
  for(int i=0;i<this->active.size();i++){
    active_ind[this->active[i]->nam] = i;
  }
  double value;
  std::string name;

  for(int i=0;i<this->physical.size();i++){
    name = this->physical[i]->getName();
    if( this->physical[i]->getActive() ){
      value = this->maps[active_ind[name]];
    } else {
      value = this->physical[i]->getValue();
    }
    names.push_back(name);
    values.push_back(value);
  }
  for(int i=0;i<this->reg.size();i++){
    name = this->reg[i]->getName();
    if( this->reg[i]->getActive() ){
      value = this->maps[active_ind[name]];
    } else {
      value = this->reg[i]->getValue();
    }
    names.push_back(name);
    values.push_back(value);
  }
  for(int j=0;j<this->lenses.size();j++){
    for(int i=0;i<this->lenses[j].size();i++){
      name = this->lenses[j][i]->getName();
      if( this->lenses[j][i]->getActive() ){
	value = this->maps[active_ind[name]];
      } else {
	value = this->lenses[j][i]->getValue();
      }
      names.push_back(this->lens_names[j] + "_" + name);
      values.push_back(value);
    }
  }
}


//virtual
Json::Value SmoothLikelihood::getActiveJson(){
  Json::Value all;

  std::map<std::string,int> active_ind;
  for(int i=0;i<this->active.size();i++){
    active_ind[this->active[i]->nam] = i;
  }

  Json::Value physical = Json::Value(Json::arrayValue);
  for(int i=0;i<this->physical.size();i++){
    if( this->physical[i]->getActive() ){
      std::string name = this->physical[i]->getName();
      Json::Value dum;
      dum["nam"]     = name;
      dum["map"]     = this->maps[active_ind[name]];
      dum["mean"]    = this->means[active_ind[name]];
      dum["s1_low"]  = this->s1_low[active_ind[name]];
      dum["s1_high"] = this->s1_high[active_ind[name]];
      physical.append(dum);
    }
  }
  if( physical.size() != 0 ){
    all["physical"] = physical;
  }

  Json::Value reg = Json::Value(Json::arrayValue);
  for(int i=0;i<this->reg.size();i++){
    if( this->reg[i]->getActive() ){
      std::string name = this->reg[i]->getName();
      Json::Value dum;
      dum["nam"]     = name;
      dum["map"]     = this->maps[active_ind[name]];
      dum["mean"]    = this->means[active_ind[name]];
      dum["s1_low"]  = this->s1_low[active_ind[name]];
      dum["s1_high"] = this->s1_high[active_ind[name]];
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
	std::string name = this->lenses[j][i]->getName();
	Json::Value dum;
	dum["nam"]     = name;
	dum["map"]     = this->maps[active_ind[name]];
	dum["mean"]    = this->means[active_ind[name]];
	dum["s1_low"]  = this->s1_low[active_ind[name]];
	dum["s1_high"] = this->s1_high[active_ind[name]];
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
    std::string suffix = Nlpar::removeSuffix(this->reg);
    this->source->kernel->setParameters(this->reg);
    Nlpar::addSuffix(this->reg,suffix);
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
void SmoothLikelihood::initialOutputLikelihoodModel(std::string output){
  // generic quantities computed in the code (smooth part)
  Json::Value json_initial;
  json_initial["Ndata"]   = this->image->Nm;
  json_initial["Nmask"]   = this->image->Nmask;
  json_initial["Psize"]   = this->image->width/this->image->Ni;
  json_initial["Nsource"] = this->source->Sm;
  json_initial["Nsource_mask"] = this->source->Smask;
  json_initial["reg_source"]        = this->source->reg;
  json_initial["sample_reg_source"] = this->source->sample_reg;

  std::ofstream jsonfile(output+"smooth_initial_output.json");
  jsonfile << json_initial;
  jsonfile.close();
}

//virtual
void SmoothLikelihood::outputLikelihoodModel(std::string output){
  // Output reconstructed source
  this->source->outputSource(output + "smooth_");

  // Output errors of reconstructed source
  double* errors = (double*) calloc(this->source->Sm,sizeof(double));
  this->algebra->getSourceErrors(this->source->Sm,errors);
  this->source->outputSourceErrors(errors,output + "smooth_");
  free(errors);

  /*  
  // Write the regularization matrix of the source as an image
  ImagePlane* matrix_s = new ImagePlane(this->source->Sm,this->source->Sm,1.0,1.0);
  for(int i=0;i<this->source->H.tri.size();i++){
    int nx = this->source->H.tri[i].i;
    int ny = this->source->H.tri[i].j;
    matrix_s->img[nx*this->source->Sm + ny] = this->source->H.tri[i].v;
  }
  matrix_s->writeImage(output+"Hmatrix_source.fits");
  delete(matrix_s);
  */

  // Create mock data (lensed MAP source)
  this->getModel();
  this->model->writeImage(output + "smooth_model.fits");

  // Output residual image (diference between mydata and mock_data)
  this->getResidual();
  this->res->writeImage(output + "smooth_residual.fits");
  
  // Output various quantities in json format
  Json::Value json_output;
  
  // Print and write to json the parameters and their associated uncertainties from the parameter model
  // 1: collapsed list of ALL the parameter names and MAP values
  Json::Value pars;
  std::vector<double> values;
  std::vector<std::string> names;
  this->getAllNamesValues(names,values);
  for(int i=0;i<names.size();i++){
    //    std::cout << it->first << " " << it->second << std::endl;
    printf("%10s %10.5f\n",names[i].c_str(),values[i]);
    pars[names[i]] = values[i];
  }
  json_output["collapsed_all"] = pars;
  
  // 2: collapsed list of the ACTIVE parameter names and MAP, mean, 1-sigma lower and 1-sigma upper bounds
  Json::Value collapsed_active;
  for(int i=0;i<this->active_names.size();i++){
    Json::Value active_par;
    active_par["map"]     = this->maps[i];
    active_par["mean"]    = this->means[i];
    active_par["s1_low"]  = this->s1_low[i];
    active_par["s1_high"] = this->s1_high[i];
    collapsed_active[this->active_names[i]] = active_par;
  }
  json_output["collapsed_active"] = collapsed_active;
  
  // 3: same as above, but the non-linear parameter json structure is preserved
  Json::Value json_active = this->getActiveJson();
  json_output["json_active"] = json_active;

  // 4: likelihood terms
  Json::Value terms;
  for(std::unordered_map<std::string,double>::iterator it=this->terms.begin();it!=this->terms.end();it++){
    terms[it->first] = it->second;
  }
  json_output["terms"] = terms;

  std::ofstream jsonfile(output+"smooth_output.json");
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
      this->active.push_back( this->reg_s[i] );
      this->active_names.push_back( this->reg_s[i]->nam );
    }
  }
  for(int i=0;i<this->reg_dpsi.size();i++){
    if( this->reg_dpsi[i]->getActive() ){
      this->active.push_back( this->reg_dpsi[i] );
      this->active_names.push_back( this->reg_dpsi[i]->nam );
    }
  }

  this->maps    = (double*) malloc(this->active.size()*sizeof(double));
  this->means   = (double*) malloc(this->active.size()*sizeof(double));
  this->s1_low  = (double*) malloc(this->active.size()*sizeof(double));
  this->s1_high = (double*) malloc(this->active.size()*sizeof(double));

  terms["Nilog2pi"] = 0.0;
  terms["Nslogl_s"] = 0.0;
  if( this->source->reg == "curvature_in_identity_out" || this->source->reg == "covariance_kernel_in_identity_out" ){
    terms["Noutlog2pi_s"] = 0.0;
  }
  terms["Nplogl_p"] = 0.0;
  if( this->pert_mass_model->dpsi->reg == "curvature_in_identity_out" || this->pert_mass_model->dpsi->reg == "covariance_kernel_in_identity_out" ){
    terms["Noutlog2pi_p"] = 0.0;
  }
  terms["detCd"]    = 0.0;
  terms["detA"]     = 0.0;
  terms["detCs"]    = 0.0;
  if( this->source->reg == "curvature_in_identity_out" || this->source->reg == "covariance_kernel_in_identity_out" ){
    terms["detLout_s"] = 0.0;
  }
  terms["detCp"]    = 0.0;
  if( this->pert_mass_model->dpsi->reg == "curvature_in_identity_out" || this->pert_mass_model->dpsi->reg == "covariance_kernel_in_identity_out" ){
    terms["detLout_p"] = 0.0;
  }
  terms["chi2"]     = 0.0;
  terms["reg"]      = 0.0; // the regularization term already containts an Lout contribution if required
  terms["evidence"] = 0.0;
}

PertLikelihood::~PertLikelihood(){
  delete(algebra);
}

//virtual
void PertLikelihood::printTerms(){
  printf(" %16s","Nilog2pi");
  printf(" %16s","detCd");
  printf(" %16s","detA");
  printf(" %16s","Nslogl_s");
  if( this->source->reg == "curvature_in_identity_out" || this->source->reg == "covariance_kernel_in_identity_out" ){
    printf(" %16s","Noutlog2pi_s");
    printf(" %16s","detLout_s");
  }
  printf(" %16s","detCs");
  printf(" %16s","Nplogl_p");
  if( this->pert_mass_model->dpsi->reg == "curvature_in_identity_out" || this->pert_mass_model->dpsi->reg == "covariance_kernel_in_identity_out" ){
    printf(" %16s","Noutlog2pi_p");
    printf(" %16s","detLout_p");
  }
  printf(" %16s","detCp");
  printf(" %16s","chi2");
  printf(" %16s","reg");
  printf(" %16s","evidence");
  printf("\n");
  printf(" %16f",terms["Nilog2pi"]);
  printf(" %16f",terms["detCd"]);
  printf(" %16f",terms["detA"]);
  printf(" %16f",terms["Nslogl_s"]);
  if( this->source->reg == "curvature_in_identity_out" || this->source->reg == "covariance_kernel_in_identity_out" ){
    printf(" %16f",terms["Noutlog2pi_s"]);
    printf(" %16f",terms["detLout_s"]);
  }
  printf(" %16f",terms["detCs"]);
  printf(" %16f",terms["Nplogl_p"]);
  if( this->pert_mass_model->dpsi->reg == "curvature_in_identity_out" || this->pert_mass_model->dpsi->reg == "covariance_kernel_in_identity_out" ){
    printf(" %16f",terms["Noutlog2pi_p"]);
    printf(" %16f",terms["detLout_p"]);
  }
  printf(" %16f",terms["detCp"]);
  printf(" %16f",terms["chi2"]);
  printf(" %16f",terms["reg"]);
  printf(" %16f",terms["evidence"]);
  printf("\n");
}

//virtual
void PertLikelihood::getAllNamesValues(std::vector<std::string>& names,std::vector<double>& values){
  std::map<std::string,int> active_ind;
  for(int i=0;i<this->active.size();i++){
    active_ind[this->active[i]->nam] = i;
  }
  double value;
  std::string name;

  for(int i=0;i<this->reg_s.size();i++){
    name = this->reg_s[i]->getName();
    if( this->reg_s[i]->getActive() ){
      value = this->maps[active_ind[name]];
    } else {
      value = this->reg_s[i]->getValue();
    }
    names.push_back(name);
    values.push_back(value);
  }
  for(int i=0;i<this->reg_dpsi.size();i++){
    name = this->reg_dpsi[i]->getName();
    if( this->reg_dpsi[i]->getActive() ){
      value = this->maps[active_ind[name]];
    } else {
      value = this->reg_dpsi[i]->getValue();
    }
    names.push_back(name);
    values.push_back(value);
  }
}

//virtual
Json::Value PertLikelihood::getActiveJson(){
  Json::Value all;

  std::map<std::string,int> active_ind;
  for(int i=0;i<this->active.size();i++){
    active_ind[this->active[i]->nam] = i;
  }

  Json::Value reg_s = Json::Value(Json::arrayValue);
  for(int i=0;i<this->reg_s.size();i++){
    if( this->reg_s[i]->getActive() ){
      std::string name = this->reg_s[i]->getName();
      Json::Value dum;
      dum["nam"]     = name;
      dum["map"]     = this->maps[active_ind[name]];
      dum["mean"]    = this->means[active_ind[name]];
      dum["s1_low"]  = this->s1_low[active_ind[name]];
      dum["s1_high"] = this->s1_high[active_ind[name]];
      reg_s.append(dum);
    }
  }
  if( reg_s.size() != 0 ){
    all["reg_s"] = reg_s;
  }

  Json::Value reg_dpsi = Json::Value(Json::arrayValue);
  for(int i=0;i<this->reg_dpsi.size();i++){
    if( this->reg_dpsi[i]->getActive() ){
      std::string name = this->reg_dpsi[i]->getName();
      Json::Value dum;
      dum["nam"]     = name;
      dum["map"]     = this->maps[active_ind[name]];
      dum["mean"]    = this->means[active_ind[name]];
      dum["s1_low"]  = this->s1_low[active_ind[name]];
      dum["s1_high"] = this->s1_high[active_ind[name]];
      reg_dpsi.append(dum);
    }
  }
  if( reg_dpsi.size() != 0 ){
    all["reg_dpsi"] = reg_dpsi;
  }

  return all;
}

//virtual
void PertLikelihood::getModel(){
  BaseSourcePlane* tmp_source = this->source->clone(); // This is the right way to copy a derived class, e.g. "this->source = new BaseSourcePlane(*c)" does not work, and "this->source = c" only copies the pointer.
  for(int i=0;i<this->source->Sm;i++){
    tmp_source->src[i] = this->source->src[i];
  }

  // Add perturbations to the mass collection
  this->pert_mass_model->updatePert();
  this->smooth_like->collection->models.push_back(this->pert_mass_model);
  this->smooth_like->collection->all_defl(this->smooth_like->image);

  tmp_source->createInterpolationWeights(this->smooth_like->image);
  tmp_source->constructL(this->smooth_like->image);

  this->algebra->getMockData(this->model,tmp_source,this->smooth_like->algebra->B);

  // Remove perturbations from the mass collection
  this->smooth_like->collection->models.pop_back();
  this->smooth_like->collection->all_defl(this->smooth_like->image);
}

//virtual
void PertLikelihood::getResidual(){
  for(int i=0;i<this->smooth_like->image->Nm;i++){
    this->res->img[i] = this->smooth_like->image->img[i] - this->model->img[i];
  }
}

//non-virtual
void PertLikelihood::initializePert(SmoothLikelihood* smooth_like){
  // This function needs to be called after a SmoothLikelihood model has been run,
  // otherwise the "smooth_like->source" will be empty, and this class needs a s_0 to initialize properly.

  // Add pointer to smooth likelihood model
  this->smooth_like = smooth_like;
  this->model = new ImagePlane(*this->smooth_like->image);
  this->res   = new ImagePlane(*this->smooth_like->image);

  // Update limits of source covariance kernel, if applicable
  /*
  if( this->source->reg == "covariance_kernel" || this->source->reg == "covariance_kernel_in_identity_out" ){
    if( this->source->type == "fixed" ){
      if( this->source->kernel->type == "modgauss" ){
	Nlpar* nlpar_sdev = Nlpar::getParByName("sdev_s",this->active);
	FixedSource* fixed = dynamic_cast<FixedSource*>(this->source); // needed to access width
	nlpar_sdev->min = fixed->width/(this->source->Sj*log(this->source->kernel->cmax));
	nlpar_sdev->max = 50.0*smooth_like->image->width;
      }
    } else if( this->source->type == "adaptive" ){
      if( this->source->kernel->type == "modgauss" ){
	Nlpar* nlpar_sdev = Nlpar::getParByName("sdev_s",this->active);
	nlpar_sdev->max = 50.0*smooth_like->image->width;
      }
    }
  }

  // Update limits of perturbations covariance kernel, if applicable
  if( this->pert_mass_model->dpsi->reg == "covariance_kernel" || this->pert_mass_model->dpsi->reg == "covariance_kernel_in_identity_out" ){
    if( this->pert_mass_model->dpsi->kernel->type == "modgauss" ){
      Nlpar* nlpar_sdev = Nlpar::getParByName("sdev_dpsi",this->active);
      nlpar_sdev->min = this->pert_mass_model->dpsi->width/(this->pert_mass_model->dpsi->Sj*log(this->pert_mass_model->dpsi->kernel->cmax));
      nlpar_sdev->max = 50.0*smooth_like->image->width;
    }
  }
  */

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
    std::string suffix = Nlpar::removeSuffix(this->reg_s);
    this->source->kernel->setParameters(this->reg_s);
    Nlpar::addSuffix(this->reg_s,suffix);
    this->source->constructH();
  }

  // Update covariance kernel for the perturbations (if needed)
  // For identity,gradient, and curvature, (no sampling of the regularization by definition) H stays the same at every call and lambda values are accessed in setAlgebraRuntime
  if( this->pert_mass_model->dpsi->sample_reg ){
    std::string suffix = Nlpar::removeSuffix(this->reg_dpsi);
    this->pert_mass_model->dpsi->kernel->setParameters(this->reg_dpsi);
    Nlpar::addSuffix(this->reg_dpsi,suffix);
    this->pert_mass_model->dpsi->constructH();
  }

  // If lambda is allowed to vary then update the lambda term for the source
  Nlpar* lambda_s = Nlpar::getParByName("lambda_s",this->reg_s);
  if( lambda_s->fix == 0 ){
    this->terms["Nslogl_s"] = this->source->Sm*log10(lambda_s->val)/2.0;
  }  
  
  // If lambda is allowed to vary then update the lambda term for the potential corrections
  Nlpar* lambda_dpsi = Nlpar::getParByName("lambda_dpsi",this->reg_dpsi);
  if( lambda_dpsi->fix == 0 ){
    this->terms["Nplogl_p"] = this->pert_mass_model->dpsi->Sm*log10(lambda_dpsi->val)/2.0;
  }

  // Update all the needed algebraic tables, e.g. M_r, Mt_r, Mt_r*Cd*M_r + RtR, Cs and detCs (if needed), Cp and detCp (if needed)
  this->algebra->setAlgebraRuntime(this->source,this->pert_mass_model);

  // Solve for the source and perturbations (r), calculate the chi2, reg, and detA terms
  this->algebra->solveSourcePert(this->source,this->pert_mass_model);
}

//virtual
void PertLikelihood::initialOutputLikelihoodModel(std::string output){
  this->pert_mass_model->dpsi->outputMask(output + "pert_dpsi_mask.fits");

  // generic quantities computed in the code (perturbations part)
  Json::Value json_initial;
  json_initial["Ndata"]        = this->smooth_like->image->Nm;
  json_initial["Nmask"]        = this->smooth_like->image->Nmask;
  json_initial["Psize"]        = this->smooth_like->image->width/this->smooth_like->image->Ni;
  json_initial["Npert"]        = this->pert_mass_model->dpsi->Sm;
  json_initial["Npert_mask"]   = this->pert_mass_model->dpsi->Smask;
  json_initial["Ppert_size"]   = this->pert_mass_model->dpsi->width/this->pert_mass_model->dpsi->Sj;
  json_initial["Nsource"]      = this->source->Sm;
  json_initial["Nsource_mask"] = this->source->Smask;
  json_initial["reg_source"]        = this->source->reg;
  json_initial["sample_reg_source"] = this->source->sample_reg;
  json_initial["reg_dpsi"]          = this->pert_mass_model->dpsi->reg;
  json_initial["sample_reg_dpsi"]   = this->pert_mass_model->dpsi->sample_reg;

  std::ofstream jsonfile(output+"pert_initial_output.json");
  jsonfile << json_initial;
  jsonfile.close();
}

//virtual
void PertLikelihood::outputLikelihoodModel(std::string output){
  // Output reconstructed source
  this->source->outputSource(output + "pert_");

  // *** For the following I need to replace the source from the smooth model with the one from this likelihood model
  // Replace source pointer in the smooth model with the new perturbed source
  BaseSourcePlane* dum_source = this->smooth_like->source;
  this->smooth_like->source = this->source;

  // Output errors of reconstructed source
  double* errors = (double*) calloc(this->source->Sm,sizeof(double));
  this->smooth_like->algebra->getSourceErrors(this->source->Sm,errors);
  this->smooth_like->source->outputSourceErrors(errors,output + "pert_");
  free(errors);

  // Create mock data (lensed MAP source)
  this->getModel();
  this->model->writeImage(output + "pert_model.fits");
  
  // Output residual image (diference between mydata and mock_data)
  this->getResidual();
  this->res->writeImage(output + "pert_residual.fits");
  
  // Replace source pointer in the smooth model with the previous source
  this->smooth_like->source = dum_source;
  // *** I don't need to use the functions in the smooth model anymore.

  /*
  // Write the regularization matrix of the source as an image
  ImagePlane* matrix_s = new ImagePlane(this->source->Sm,this->source->Sm,1.0,1.0);
  for(int i=0;i<this->source->H.tri.size();i++){
    int nx = this->source->H.tri[i].i;
    int ny = this->source->H.tri[i].j;
    matrix_s->img[nx*this->source->Sm + ny] = this->source->H.tri[i].v;
  }
  matrix_s->writeImage(output+"Hmatrix_source.fits");
  delete(matrix_s);
  */

  /*
  // Write the regularization matrix of the perturbations as an image
  ImagePlane* matrix_p = new ImagePlane(this->pert_mass_model->dpsi->Sm,this->pert_mass_model->dpsi->Sm,1.0,1.0);
  for(int i=0;i<this->pert_mass_model->dpsi->H.tri.size();i++){
    int nx = this->pert_mass_model->dpsi->H.tri[i].i;
    int ny = this->pert_mass_model->dpsi->H.tri[i].j;
    matrix_p->img[nx*this->pert_mass_model->dpsi->Sm + ny] = this->pert_mass_model->dpsi->H.tri[i].v;
  }
  matrix_p->writeImage(output+"Hmatrix_pert.fits");
  delete(matrix_p);
  */
  
  // Output current perturbations correction
  this->pert_mass_model->dpsi->outputSource(output + "pert_dpsi.fits");

  // Output errors on perturbations
  //  double* errors = (double*) calloc(this->pert_mass_model->dpsi->Sm,sizeof(double));
  //  this->pert_mass_model->dpsi->outputSourceErrors(errors,output + "pert_dpsi_errors.fits");
  //  free(errors);

  // Output convergence
  ImagePlane* kappa = new ImagePlane(this->pert_mass_model->dpsi->Si,this->pert_mass_model->dpsi->Sj,this->pert_mass_model->dpsi->width,this->pert_mass_model->dpsi->height);
  pert_mass_model->getConvergence(kappa);
  kappa->writeImage(output + "pert_convergence.fits");
  delete(kappa);


  // Output various quantities in json format
  Json::Value json_output;
  
  // Print and write to json the parameters and their associated uncertainties from the parameter model
  // 1: collapsed list of ALL the parameter names and MAP values
  Json::Value pars;
  std::vector<double> values;
  std::vector<std::string> names;
  this->getAllNamesValues(names,values);
  for(int i=0;i<names.size();i++){
    std::cout << names[i] << " " << values[i] << std::endl;
    pars[names[i]] = values[i];
  }
  json_output["collapsed_all"] = pars;
  
  // 2: collapsed list of the ACTIVE parameter names and MAP, mean, 1-sigma lower and 1-sigma upper bounds
  Json::Value collapsed_active;
  for(int i=0;i<this->active_names.size();i++){
    Json::Value active_par;
    active_par["map"]     = this->maps[i];
    active_par["mean"]    = this->means[i];
    active_par["s1_low"]  = this->s1_low[i];
    active_par["s1_high"] = this->s1_high[i];
    collapsed_active[this->active_names[i]] = active_par;
  }
  json_output["collapsed_active"] = collapsed_active;
  
  // 3: same as above, but the non-linear parameter json structure is preserved
  Json::Value json_active = this->getActiveJson();
  json_output["json_active"] = json_active;

  // 4: likelihood terms
  Json::Value terms;
  for(std::unordered_map<std::string,double>::iterator it=this->terms.begin();it!=this->terms.end();it++){
    terms[it->first] = it->second;
  }
  json_output["terms"] = terms;

  std::ofstream jsonfile(output+"pert_output.json");
  jsonfile << json_output;
  jsonfile.close();
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
    std::string suffix = Nlpar::removeSuffix(this->reg_s);
    this->source->kernel->setParameters(this->reg_s);
    Nlpar::addSuffix(this->reg_s,suffix);
  }
  // Construct a new H for the source if needed
  if( this->source->type == "adaptive" || this->source->sample_reg ){
    this->source->constructH();
  }
  this->source->constructDs(this->smooth_like->image);


  // Update covariance kernel for the perturbations (if needed)
  // For identity,gradient, and curvature, (no sampling of the regularization by definition) H stays the same at every call and lambda values are accessed in setAlgebraRuntime
  if( this->pert_mass_model->dpsi->sample_reg ){
    std::string suffix = Nlpar::removeSuffix(this->reg_dpsi);
    this->pert_mass_model->dpsi->kernel->setParameters(this->reg_dpsi);
    Nlpar::addSuffix(this->reg_dpsi,suffix);
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
void PertIterationLikelihood::outputLikelihoodModel(std::string output){
  // Output reconstructed source
  this->source->outputSource(output + "pert_");
  
  // Replace source pointer in the smooth model with the new perturbed source
  BaseSourcePlane* dum_source = this->smooth_like->source;
  this->smooth_like->source = this->source;

  // Output errors of reconstructed source
  double* errors = (double*) calloc(this->source->Sm,sizeof(double));
  this->smooth_like->algebra->getSourceErrors(this->source->Sm,errors);
  this->smooth_like->source->outputSourceErrors(errors,output + "pert_");
  free(errors);

  // Create mock data (lensed MAP source)
  this->smooth_like->getModel();
  this->smooth_like->model->writeImage(output + "pert_model.fits");
  
  // Output residual image (diference between mydata and mock_data)
  this->smooth_like->getResidual();
  this->smooth_like->res->writeImage(output + "pert_residual.fits");
  
  // Replace source pointer in the smooth model with the previous source
  this->smooth_like->source = dum_source;

  // Output last perturbations (corrections)
  this->pert_mass_model->dpsi->outputSource(output + "pert_dpsi.fits");

  // Output added perturbations
  Pert* pert_pointer = dynamic_cast<Pert*>(this->collection->models.back());
  pert_pointer->dpsi->outputSource(output + "pert_added_dpsi.fits"); // output additive perturbations

  // Output convergence
  ImagePlane* kappa = new ImagePlane(pert_pointer->dpsi->Si,pert_pointer->dpsi->Sj,pert_pointer->dpsi->width,pert_pointer->dpsi->height);
  pert_pointer->getConvergence(kappa);
  kappa->writeImage(output + "pert_added_convergence.fits");
  delete(kappa);

  // Output various quantities in json format
  Json::Value json_output;
  
  // Print and write to json the parameters and their associated uncertainties from the parameter model
  // 1: collapsed list of ALL the parameter names and MAP values
  Json::Value pars;
  std::vector<double> values;
  std::vector<std::string> names;
  this->getAllNamesValues(names,values);
  for(int i=0;i<names.size();i++){
    pars[names[i]] = values[i];
  }
  json_output["collapsed_all"] = pars;
  
  // 2: collapsed list of the ACTIVE parameter names and MAP, mean, 1-sigma lower and 1-sigma upper bounds
  Json::Value collapsed_active;
  for(int i=0;i<this->active_names.size();i++){
    Json::Value active_par;
    active_par["map"]     = this->maps[i];
    active_par["mean"]    = this->means[i];
    active_par["s1_low"]  = this->s1_low[i];
    active_par["s1_high"] = this->s1_high[i];
    collapsed_active[this->active_names[i]] = active_par;
  }
  json_output["collapsed_active"] = collapsed_active;
  
  // 3: same as above, but the non-linear parameter json structure is preserved
  Json::Value json_active = this->getActiveJson();
  json_output["json_active"] = json_active;

  // 4: likelihood terms
  Json::Value terms;
  for(std::unordered_map<std::string,double>::iterator it=this->terms.begin();it!=this->terms.end();it++){
    terms[it->first] = it->second;
  }
  json_output["terms"] = terms;

  std::ofstream jsonfile(output+"pert_output.json");
  jsonfile << json_output;
  jsonfile.close();
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
    //std::cout << par->nam << " " << par->pri->type << std::endl;
    //std::cout << par->nam << " " << par->pri->prior(par->min) << " " << par->pri->prior(par->max) << std::endl;
    //par->setNewPrior(prior);
    
    pars.push_back(par);
  }
  
  return pars;
}


