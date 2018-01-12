#include "inputOutput.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "json/json.h"

#include "nonLinearPars.hpp"
#include "likelihoodModels.hpp"
#include "massModels.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"




void Initialization::initialize_program(std::string path,std::string run,Initialization*& init,BaseLikelihoodModel*& mypars,ImagePlane*& mydata,CollectionMassModels*& mycollection,BaseSourcePlane*& mysource){
  printf("%-50s","Starting initialization");
  fflush(stdout);
  


  // Instantiate init object ---------------------------------------------------------------------------------------------------------------------------------------
  init = new Initialization();
  init->parseInputJSON(path,run);



  // Initialize likelihood model -----------------------------------------------------------------------------------------------------------------------------------
  mypars = FactoryLikelihoodModel::getInstance()->createLikelihoodModel(path,run);



  // Read image data -----------------------------------------------------------------------------------------------------------------------------------------------
  mydata = new ImagePlane(init->imgpath,stoi(init->image["pix_x"]),stoi(init->image["pix_y"]),stof(init->image["siz_x"]),stof(init->image["siz_y"]));
  mydata->readB(init->psfpath,stoi(init->psf["pix_x"]),stoi(init->psf["pix_y"]),stoi(init->psf["crop_x"]),stoi(init->psf["crop_y"]));
  mydata->readC(init->noise_flag,init->covpath);
  mydata->readS(init->maskpath);




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
  //  if( mysource->reg == "covariance_kernel" ){
  //    Standard* my_standard_pars = dynamic_cast<Standard*>(mypars);
  //    mysource->kernel = FactoryCovKernel::getInstance()->createCovKernel(init->source["kernel"],my_standard_pars->getRegPars());
  //  }



  // Initialize precomputed algebraic quantities -------------------------------------------------------------------------------------------------------------------
  //Precompute a number of algebraic quantities (table products etc)                                                                [<---algebra package dependent]
  mypars->initializeAlgebra(mydata,mysource);




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
  std::ifstream fin((runname+"vkl_input.json").c_str());
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
  jmembers = root["minimizer"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->minimizer[jmembers[i]] = root["minimizer"][jmembers[i]].asString();
  }


  //Mass model name for the lenses
  jmembers = root["lenses"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->mmodel.push_back( root["lenses"][jmembers[i]]["subtype"].asString() );
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

void Initialization::outputInitial(BaseLikelihoodModel* mypars){

  // File with the parameter names
  std::ofstream f_names(this->output+"plt_corner.paramnames",std::ofstream::out);
  std::vector<std::string> active_full_names = mypars->getActiveFullNames();
  for(int i=0;i<active_full_names.size();i++){
    f_names << active_full_names[i];
    f_names << " ";
    f_names << active_full_names[i];
    f_names << std::endl;
  }
  f_names.close();

  // File with the parameter ranges
  std::ofstream f_ranges(this->output+"plt_corner.ranges",std::ofstream::out);
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



void Initialization::finalize_program(Initialization* init,BaseLikelihoodModel* mypars,ImagePlane* mydata,CollectionMassModels* mycollection,BaseSourcePlane* mysource){
  printf("%-25s\n","Starting output");
  fflush(stdout);
  
  mypars->updateLikelihoodModel(mydata,mysource,mycollection);
  mypars->getLogLike(mydata,mysource);
  mypars->printActive();
  mypars->printTerms();
  printf("\n");
  mypars->outputLikelihoodModel(mydata,mysource,init->output);
  
  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);  
}




