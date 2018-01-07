#include "inputOutput.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "json/json.h"

#include "nonLinearPars.hpp"
#include "parameterModels.hpp"
#include "massModels.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "tableAlgebra.hpp"




void Initialization::initialize_program(std::string path,std::string run,Initialization*& init,BaseParameterModel*& mypars,ImagePlane*& mydata,CollectionMassModels*& mycollection,BaseSourcePlane*& mysource,mymatrices*& matrices,precomp*& pcomp){
  printf("%-50s","Starting initialization");
  fflush(stdout);
  


  // Instantiate init object ---------------------------------------------------------------------------------------------------------------------------------------
  init = new Initialization();
  init->parseInputJSON(path,run);



  // Initialize parameterModel -------------------------------------------------------------------------------------------------------------------------------------
  mypars = FactoryParameterModel::getInstance()->createParameterModel(path,run);



  // Read image data -----------------------------------------------------------------------------------------------------------------------------------------------
  mydata = new ImagePlane(init->imgpath,stoi(init->image["pix_x"]),stoi(init->image["pix_y"]),stof(init->image["siz_x"]),stof(init->image["siz_y"]));



  // Initialize mass model parameters ------------------------------------------------------------------------------------------------------------------------------
  mycollection = new CollectionMassModels(mypars->getPhysicalPars());
  mycollection->models.resize(init->mmodel.size());
  for(int k=0;k<mycollection->models.size();k++){
    std::vector<Nlpar*> lens = mypars->getMassModelPars(k);
    mycollection->models[k] = FactoryMassModel::getInstance()->createMassModel(init->mmodel[k],lens);
  }



  // Initialize source ---------------------------------------------------------------------------------------------------------------------------------------------
  //initialize here, but also in each iteration for an adaptive grid                                                  [<---iteration dependent for adaptive source]
  mysource = FactorySourcePlane::getInstance()->createSourcePlane(init->source);
  if( mysource->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(mysource);
    ada->createAdaGrid(mydata,mycollection);
    ada->createDelaunay();
  }
  if( init->source["sample_reg"] == "true" ){
    mysource->sample_reg = true;
  }
  if( mysource->reg == "covariance_kernel" ){
    Standard* my_standard_pars = dynamic_cast<Standard*>(mypars);
    mysource->kernel = FactoryCovKernel::getInstance()->createCovKernel(init->source["kernel"],my_standard_pars->getRegPars());
  }



  // Initialize matrices: S,C,B,H,L --------------------------------------------------------------------------------------------------------------------------------
  //This is independent of the specific linear algebra package used.
  //S: mask,
  //C: covariance matrix,
  //B: blurring matrix,
  //L: lensing matrix,
  //H: source regularization matrix
  matrices = new mymatrices();
  matrices->initMatrices(mydata,mysource,init->maskpath,init->noise_flag,init->covpath,init->psfpath,init->psf);



  // Initialize precomputed algebraic quantities -------------------------------------------------------------------------------------------------------------------
  //Precompute a number of algebraic quantities (table products etc)                                                                [<---algebra package dependent]
  pcomp = new precomp(mydata,mysource);
  setAlgebraInit(mydata,mysource,matrices,pcomp);




  // Initial output ------------------------------------------------------------------------------------------------------------------------------------------------
  init->outputInitial(mypars);


  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);
}


void Initialization::parseInputJSON(std::string path,std::string run){
  Json::Value root;
  Json::Value::Members jmembers;

  std::string runname = path + run;
  std::ifstream fin(runname+"vkl_input.json");
  fin >> root;

  //General parameters
  this->output = runname+root["output"].asString();
  if( root["maskpath"].asString() == "0" ){
    this->maskpath = "0";
  } else {
    this->maskpath = path+root["maskpath"].asString();
  }
  

  //Parameters for the image
  this->imgpath = path+root["imgpath"].asString();
  jmembers = root["iplane"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->image[jmembers[i]] = root["iplane"][jmembers[i]].asString();
  }

  //Parameters of the source
  this->interp = root["interp"].asString();
  jmembers = root["splane"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->source[jmembers[i]] = root["splane"][jmembers[i]].asString();
  }
  this->source["reg"] = root["reg"]["type"].asString();
  if( this->source["reg"] == "covariance_kernel" ){
    this->source["kernel"] = root["reg"]["subtype"].asString();
  }

  //Calculate the number of adaptive source pixels based on the image data
  if( this->source["type"] == "adaptive" && ( this->source["mode"] == "image" || this->source["mode"] == "grid" ) ){
    int spacing = std::stoi(this->source["spacing"]);
    int Ni      = std::stoi(this->image["pix_y"]);
    int Nj      = std::stoi(this->image["pix_x"]);

    int count = 0;
    int i0    = (int) floor( (spacing-1)/2. );
    int j0    = (int) floor( (spacing-1)/2. );
    for(int i=i0;i<Ni;i=i+spacing){
      for(int j=j0;j<Nj;j=j+spacing){
	count++;
      }
    }

    this->source["sm"] = std::to_string( count );
    //    std::cout << count << std::endl;
    //    std::cout << (int) floor(Ni/spacing)*floor(Nj/spacing) << std::endl;
  }

  //Parameters for the psf
  this->psfpath = path+root["psfpath"].asString();
  jmembers = root["psf"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->psf[jmembers[i]] = root["psf"][jmembers[i]].asString();
  }

  //Parameters for the noise
  this->noise_flag = root["noise_flag"].asString();
  this->covpath    = path+root["covpath"].asString();

  //Parameters for the minimizer/sampler
  this->minimizer["type"] = root["minimizer"]["type"].asString();
  jmembers = root[ this->minimizer["type"] ].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->minimizer[jmembers[i]] = root["minimizer"][jmembers[i]].asString();
    //    std::cout << root[ this->minimizer["type"] ][jmembers[i]].asString() << std::endl;
  }

  const Json::Value lenses = root["lenses"];
  for(int i=0;i<lenses.size();i++){
    this->mmodel.push_back( lenses[i]["subtype"].asString() );
  }

  /*  
  // Sampling parameters of the regularization
  int c = 0;
  typedef std::map<std::string,BaseNlpar*>::iterator it_nlpar;
  for(it_nlpar myiterator=nlpars[1].begin();myiterator!=nlpars[1].end();myiterator++){
    if( myiterator->second->fix == 0 ){
      c++;
    }
  }
  if( c > 0 && this->source["reg"] == "covariance_kernel" ){
    this->source["sample_reg"] = "true";
  }
  */

}

void Initialization::outputInitial(BaseParameterModel* mypars){

  // File with the parameter names
  std::ofstream f_names(this->output+"corner.paramnames",std::ofstream::out);
  std::vector<std::string> active_full_names = mypars->getActiveFullNames();
  for(int i=0;i<active_full_names.size();i++){
    f_names << active_full_names[i];
    f_names << " ";
    f_names << active_full_names[i];
    f_names << std::endl;
  }
  f_names.close();

  // File with the parameter ranges
  std::ofstream f_ranges(this->output+"corner.ranges",std::ofstream::out);
  for(int i=0;i<mypars->active.size();i++){
    f_ranges << active_full_names[i];
    f_ranges << " ";
    f_ranges << mypars->active[i]->min;
    f_ranges << " ";
    f_ranges << mypars->active[i]->max;
    f_ranges << std::endl;
  }
  f_ranges.close();

}



void Initialization::finalize_program(Initialization* init,BaseParameterModel* mypars,ImagePlane* mydata,CollectionMassModels* mycollection,BaseSourcePlane* mysource,mymatrices* matrices,precomp* pcomp){
  printf("%-25s\n","Starting output");
  fflush(stdout);
  
  mypars->updateParameterModel(mydata,mysource,mycollection,matrices,pcomp);
  mypars->getLogLike(mydata,mysource,pcomp);
  Initialization::outputGeneric(mypars,mydata,mysource,pcomp,init->output);
  
  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);  
}

void Initialization::outputGeneric(BaseParameterModel* mypars,ImagePlane* image,BaseSourcePlane* source,precomp* pcomp,std::string output){
  // Output reconstructed source
  source->outputSource(output + "vkl_voronoi.dat");
  
  // Output errors of reconstructed source
  double* errors = (double*) calloc(source->Sm,sizeof(double));
  getSourceErrors(source->Sm,errors,pcomp);
  source->outputSourceErrors(errors,output + "vkl_voronoi_errors.dat");
  
  // Create mock data (lensed MAP source)
  ImagePlane mock_data(image->Nj,image->Ni,image->width,image->height);
  getMockData(&mock_data,source,pcomp);
  
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
  std::vector<std::string> all_names = mypars->getFullNames();
  std::vector<double> all_values = mypars->getValues();
  for(int i=0;i<all_names.size();i++){
    std::cout << all_names[i] << " " << all_values[i] << std::endl;    
    pars[all_names[i]] = all_values[i];
  }
  json_output["full_pars"] = pars;
  
  Json::Value full_active;
  std::vector<std::string> names = mypars->getActiveFullNames();
  std::vector<double> values = mypars->getActiveValues();
  for(int i=0;i<names.size();i++){
    full_active[names[i]] = values[i];
  }
  json_output["full_active"] = full_active;
  
  // Write structured active parameter names and values
  Json::Value json_active = mypars->getActiveNamesValues();
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




