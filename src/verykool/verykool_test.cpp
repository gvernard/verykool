#include <iostream>
#include <string>
#include <vector>

#include "inputOutput.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "massModels.hpp"
#include "minimizers.hpp"
#include "likelihoodModels.hpp"

// This version runs only on one CPU, withouth being capable to use MPI.
// This is to make it easier to use debugging tools like GDM and Valgrind.
// Only the 'Nothing' minimizer is supposed to work in this code.

int main(int argc,char* argv[]){
  // Check command line arguments
  if( argc != 3 ){
    std::cout << "2 command line arguments are required!!!" << std::endl;
    return 0;
  }


  //=============== BEGIN:INITIALIZATION =======================
  // Initialize variables
  Initialization* init = 0;
  ImagePlane*           vkl_data = 0;
  BaseSourcePlane*      vkl_source = 0;
  BaseSourcePlane*      vkl_source0 = 0;
  CollectionMassModels* vkl_collection = 0;
  Pert*                 vkl_pert_mass_model = 0;
  BaseMinimizer*        vkl_minimizer = 0;
  BaseLikelihoodModel*  vkl_likeModel = 0;

  printf("%-50s","Starting initialization");
  fflush(stdout);

  Initialization::initialize_program(argv[1],argv[2],init,vkl_data,vkl_collection,vkl_source,vkl_source0,vkl_pert_mass_model,vkl_minimizer);

  if( init->likeModel == "standard" ){
    vkl_minimizer->name = "smooth";
    vkl_likeModel = new SmoothLikelihood(init->nlpars_physical,init->nlpars_reg_s,init->nlpars_lenses,init->lens_names,vkl_data,vkl_source,vkl_collection);
  } else if( init->likeModel == "perturbations_standard" ){
    vkl_minimizer->name = "pert";
    vkl_likeModel = new PertLikelihood(init->nlpars_reg_s,init->nlpars_reg_dpsi,vkl_data,vkl_source,vkl_source0,vkl_collection,vkl_pert_mass_model);
  } else if( init->likeModel == "both" ){
    vkl_minimizer->name = "both";
    vkl_likeModel = new BothLikelihood(init->nlpars_physical,init->nlpars_lenses,init->lens_names,init->nlpars_reg_s,init->nlpars_reg_dpsi,vkl_data,vkl_source,vkl_collection,vkl_pert_mass_model);
  } else if( init->likeModel == "perturbations_iter" ){
    // Something
  }

  vkl_likeModel->initializeAlgebra();

  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);
  //================= END:INITIALIZATION =======================


  //=============== BEGIN:MINIMIZATION =========================
  // Initial output
  vkl_likeModel->initialOutputLikelihoodModel(init->output);
  printf("Starting %s minimization using %s",vkl_minimizer->name.c_str(),vkl_minimizer->type.c_str());
  fflush(stdout);
  
  vkl_minimizer->minimize(init->minimizer,vkl_likeModel,init->output);
  
  // Finalize output etc
  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);
  vkl_likeModel->finalizeLikelihoodModel(init->output);
  vkl_minimizer->finalizeMinimizer(init->output);
  //================= END:MINIMIZATION =========================

  
  // Cleanup pointers
  delete(init);
  delete(vkl_data);
  delete(vkl_source);
  delete(vkl_source0);
  delete(vkl_collection);
  delete(vkl_pert_mass_model);
  delete(vkl_minimizer);
  delete(vkl_likeModel);

  return 0;
}
