#include "inputOutput.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "json/json.h"

#include "nonLinearPars.hpp"
#include "likelihoodModels.hpp"
#include "imagePlane.hpp"
#include "massModels.hpp"
#include "sourcePlane.hpp"




void Initialization::initialize_program(std::string path,std::string run,Initialization*& init,BaseLikelihoodModel*& smooth_like,ImagePlane*& mydata,CollectionMassModels*& mycollection,BaseSourcePlane*& mysource,BaseLikelihoodModel*& pert_like,Pert*& pert_mass_model){
  printf("%-50s","Starting initialization");
  fflush(stdout);



  // Instantiate init object ---------------------------------------------------------------------------------------------------------------------------------------
  init = new Initialization();
  init->parseInputJSON(path,run);



  // Initializing components: image plane --------------------------------------------------------------------------------------------------------------------------
  mydata = new ImagePlane(init->imgpath,stoi(init->image["pix_x"]),stoi(init->image["pix_y"]),stof(init->image["siz_x"]),stof(init->image["siz_y"]));
  mydata->readB(init->psfpath,stoi(init->psf["pix_x"]),stoi(init->psf["pix_y"]),stoi(init->psf["crop_x"]),stoi(init->psf["crop_y"]));
  mydata->readC(init->noise_flag,init->covpath);
  mydata->readS(init->maskpath);



  // Initializing components: mass model collection ----------------------------------------------------------------------------------------------------------------
  mycollection = new CollectionMassModels();



  // Initializing components: source plane -------------------------------------------------------------------------------------------------------------------------
  mysource = FactorySourcePlane::getInstance()->createSourcePlane(init->source);



  // Initialize smooth likelihood model (requires the initialized components from above) ---------------------------------------------------------------------------
  smooth_like = FactoryLikelihoodModel::getInstance()->createLikelihoodModel(path,run,init->smooth_like,mydata,mysource,mycollection,pert_mass_model);



  // Update/initialize mass model parameters -----------------------------------------------------------------------------------------------------------------------
  mycollection->setPhysicalPars(smooth_like->getPhysicalPars());
  mycollection->models.resize(init->mmodel.size());
  for(int k=0;k<mycollection->models.size();k++){
    std::vector<Nlpar*> lens = smooth_like->getMassModelPars(k);
    mycollection->models[k] = FactoryMassModel::getInstance()->createMassModel(init->mmodel[k],lens);
  }
  mycollection->all_defl(mydata); // I need to deflect the image plane once I know the mass model parameters. Needed to create the adaptive grid.



  // Update/initialize source --------------------------------------------------------------------------------------------------------------------------------------
  //initialize here, but also in each iteration for an adaptive grid                                                  [<---iteration dependent for adaptive source]
  if( mysource->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(mysource);
    ada->createAdaGrid(mydata,mycollection);
    ada->createDelaunay();
  }
  if( init->source["reg"] == "covariance_kernel" ){
    SourceCovarianceKernel* mycovpars = dynamic_cast<SourceCovarianceKernel*>(smooth_like);
    std::vector<Nlpar*> covreg = mycovpars->getRegPars();
    mysource->kernel = FactoryCovKernel::getInstance()->createCovKernel(init->source["kernel"],covreg);
    for(int i=0;i<covreg.size();i++){
      if( covreg[i]->fix == 0 ){
	mysource->sample_reg = true;
      }
    }
  }
  mysource->constructH();



  // Initialize precomputed algebraic quantities for smooth model likelihood ---------------------------------------------------------------------------------------
  //Precompute a number of algebraic quantities (table products etc)                                                                [<---algebra package dependent]
  smooth_like->initializeAlgebra();



  // Initialize perturbations --------------------------------------------------------------------------------------------------------------------------------------
  if( init->perturbations.size() > 0 ){
    pert_mass_model = new Pert(std::stoi(init->perturbations["pix_x"]),std::stoi(init->perturbations["pix_y"]),mydata,init->perturbations["reg_dpsi"]);
    pert_like = FactoryLikelihoodModel::getInstance()->createLikelihoodModel(path,run,init->pert_like,mydata,mysource,mycollection,pert_mass_model);
  }



  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);
}




void Initialization::parseInputJSON(std::string path,std::string run){
  Json::Value root;
  Json::Value::Members jmembers;

  std::string runname = path + run;
  std::ifstream fin((runname+"vkl_input.json").c_str());
  fin >> root;


  //General parameters
  this->output = runname+root["output"].asString();
  if( root["maskpath"].asString() == "0" ){
    this->maskpath = "0";
  } else {
    this->maskpath = path+root["maskpath"].asString();
  }
  this->smooth_like = root["parameter_model"].asString();



  //Parameters for the perturbations
  if( root.isMember("perturbations") ){
    this->pert_like = root["perturbations"]["parameter_model"].asString();
    this->perturbations["pix_x"] = root["perturbations"]["pix_x"].asString();
    this->perturbations["pix_y"] = root["perturbations"]["pix_y"].asString();

    jmembers = root["perturbations"]["minimizer"].getMemberNames();
    for(int i=0;i<jmembers.size();i++){
      this->pert_minimizer[jmembers[i]] = root["perturbations"]["minimizer"][jmembers[i]].asString();
    }

    this->perturbations["reg_s"] = root["perturbations"]["reg_s"]["type"].asString();
    if( this->perturbations["reg_s"] == "covariance_kernel" ){
      this->perturbations["kernel"] = root["perturbations"]["reg_s"]["subtype"].asString();
    }

    this->perturbations["reg_dpsi"] = root["perturbations"]["reg_dpsi"]["type"].asString();
    if( this->perturbations["reg_dpsi"] == "covariance_kernel" ){
      this->perturbations["kernel"] = root["perturbations"]["reg_dpsi"]["subtype"].asString();
    }
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
  jmembers = root["minimizer"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->smooth_minimizer[jmembers[i]] = root["minimizer"][jmembers[i]].asString();
  }


  //Mass model name for the lenses
  jmembers = root["lenses"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->mmodel.push_back( root["lenses"][jmembers[i]]["subtype"].asString() );
  }
}



void Initialization::finalize_smooth(std::string output,BaseLikelihoodModel* smooth_like){
  printf("%-25s\n","Starting output smooth");
  fflush(stdout);
  
  smooth_like->updateLikelihoodModel();
  smooth_like->getLogLike();
  smooth_like->printActive();
  smooth_like->printTerms();
  printf("\n");
  smooth_like->outputLikelihoodModel(output);
  
  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);  
}


void Initialization::finalize_pert(std::string output,BaseLikelihoodModel* pert_like){
  printf("%-25s\n","Starting output perturbations");
  fflush(stdout);
  
  pert_like->updateLikelihoodModel();
  pert_like->getLogLike();
  pert_like->printActive();
  pert_like->printTerms();
  printf("\n");
  pert_like->outputLikelihoodModel(output);
  
  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);
}

