#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "inputOutput.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "massModels.hpp"
#include "minimizers.hpp"
#include "likelihoodModels.hpp"
#include "eigenAlgebra.hpp"
#include "nonLinearPars.hpp"

int main(int argc,char* argv[]){

  Initialization* init = 0;
  ImagePlane*           vkl_data = 0;
  BaseSourcePlane*      vkl_source = 0;
  BaseSourcePlane*      vkl_source0 = 0;
  CollectionMassModels* vkl_collection = 0;
  Pert*                 vkl_pert_mass_model = 0;
  BaseMinimizer*        vkl_minimizer = 0;
  BaseLikelihoodModel*  vkl_likeModel = 0;

  Initialization::initialize_program(argv[1],argv[2],init,vkl_data,vkl_collection,vkl_source,vkl_source0,vkl_pert_mass_model,vkl_minimizer);

  //vkl_likeModel = new BothLikelihood(init->nlpars_physical,init->nlpars_lenses,init->lens_names,init->nlpars_reg_s,init->nlpars_reg_dpsi,vkl_data,vkl_source,vkl_collection,vkl_pert_mass_model);
  vkl_likeModel = new PertLikelihood(init->nlpars_reg_s,init->nlpars_reg_dpsi,vkl_data,vkl_source,vkl_source0,vkl_collection,vkl_pert_mass_model);

  vkl_likeModel->initializeAlgebra();
  //vkl_likeModel->algebra->constructDsDpsi(vkl_data,vkl_source0,vkl_pert_mass_model);

  // And here I just need to output the operator into some file

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
