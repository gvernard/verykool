#include <iostream>
#include <string>
#include <vector>


#include "inputOutput.hpp"
#include "imagePlane.hpp"
#include "massModels.hpp"
#include "sourcePlane.hpp"
#include "likelihoodModels.hpp"
#include "minimizers.hpp"


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
  ImagePlane* mydata = 0;
  BaseSourcePlane* mysource = 0;
  CollectionMassModels* mycollection = 0;
  BaseLikelihoodModel* smooth_like = 0;
  BaseLikelihoodModel* pert_like = 0;
  Pert* pert_mass_model = 0;
  BaseMinimizer* smooth_minimizer = 0;

  Initialization::initialize_program(argv[1],argv[2],init,smooth_like,mydata,mycollection,mysource,pert_like,pert_mass_model);
  //================= END:INITIALIZATION =======================





  //=============== BEGIN:SMOOTH MODEL =========================
  // Initial output
  smooth_like->initialOutputLikelihoodModel(init->output);

  printf("%-25s","Starting smooth minimization ");
  fflush(stdout);
 
  smooth_minimizer = FactoryMinimizer::getInstance()->createMinimizer("smooth",init->smooth_minimizer,smooth_like,init->output);
  smooth_minimizer->minimize(init->smooth_minimizer,smooth_like,init->output);

  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);

  // Finalize output etc
  Initialization::finalizeLikelihoodModel(init->output,smooth_like);
  smooth_minimizer->finalizeMinimizer(init->output);

  delete(smooth_minimizer);
  //================= END:SMOOTH MODEL =========================





  /*    
  Pert* tmp_pert_mass_model = new Pert("/net/argo/data/users/gvernard/RESULTS/VKL_MODELS/test_modelling_perturbations/test_4/data/perturbations.fits",mydata->Ni,mydata->Nj,mydata->width,mydata->height,"identity");
  tmp_pert_mass_model->createCrosses(mydata);
  mysource->constructDs(mydata);

  BaseLikelihoodModel* tmp_pert_like = FactoryLikelihoodModel::getInstance()->createLikelihoodModel(argv[1],argv[2],"perturbations_standard",mydata,mysource,mycollection,tmp_pert_mass_model);
  PertLikelihood* pert_pointer = dynamic_cast<PertLikelihood*>(tmp_pert_like);
  pert_pointer->initializePert(dynamic_cast<SmoothLikelihood*>(smooth_like));
  pert_pointer->pertResiduals(init->output,mydata,mysource,tmp_pert_mass_model);
  delete(tmp_pert_like);
  delete(tmp_pert_mass_model);
  */  
  
  


  //=============== BEGIN:PERTURBATIONS =========================
  if( init->perturbations.size() > 0 ){
    pert_like->initialOutputLikelihoodModel(init->output);

    //Initialize perturbations
    SmoothLikelihood* smooth_pointer = dynamic_cast<SmoothLikelihood*>(smooth_like);

    // I need to cast either to 'PetLikelihood' or to 'PertIterationLikelihood'
    if( init->pert_like == "perturbations_standard" ){
      PertLikelihood* pert_pointer   = dynamic_cast<PertLikelihood*>(pert_like);
      pert_pointer->initializePert(smooth_pointer);
    } else if( init->pert_like == "perturbations_iter" ){
      PertIterationLikelihood* pert_pointer   = dynamic_cast<PertIterationLikelihood*>(pert_like);
      pert_pointer->initializePert(smooth_pointer);      
    }


    printf("%-25s","Starting perturbation minimization ");
    fflush(stdout);
    
    BaseMinimizer* pert_minimizer = FactoryMinimizer::getInstance()->createMinimizer("pert",init->pert_minimizer,pert_like,init->output);
    pert_minimizer->minimize(init->pert_minimizer,pert_like,init->output);

    printf("%+7s\n","...done");
    std::cout << std::string(200,'=') << std::endl;
    fflush(stdout);

    Initialization::finalizeLikelihoodModel(init->output,pert_like);
    pert_minimizer->finalizeMinimizer(init->output);
    
    delete(pert_minimizer);
  }
  //================= END:PERTURBATIONS =========================




  // Cleanup pointers
  delete(init);
  delete(mydata);
  delete(mysource);
  delete(mycollection);
  delete(smooth_like);
  delete(pert_like);
  delete(pert_mass_model);
  

  return 0;
}
