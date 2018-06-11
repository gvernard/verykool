#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <stdlib.h> 

#include "inputOutput.hpp"
#include "nonLinearPars.hpp"
#include "likelihoodModels.hpp"



int main(int argc,char* argv[]){
  std::string path(argv[1]);
  std::string run(argv[2]);
  ImagePlane* dum_image = 0;
  BaseSourcePlane* dum_source = 0;
  CollectionMassModels* dum_collection = 0;
  Pert* dum_pert = 0;
  BaseLikelihoodModel* dum_smooth_like;
  BaseLikelihoodModel* dum_pert_like;
  std::string like_model_name = "standard";

  Initialization* init = new Initialization();
  Initialization::initialize_program(path,run,init,dum_smooth_like,dum_image,dum_collection,dum_source,dum_pert_like,dum_pert);
  BaseLikelihoodModel* mypars = FactoryLikelihoodModel::getInstance()->createLikelihoodModel(path,run,like_model_name,dum_image,dum_source,dum_collection,dum_pert);

  if( init->smooth_minimizer["type"] != "cosmosis_emcee" ){
    std::cout << "wrong minimizer, it has to be cosmosis_emcee" << std::endl;
    return 0;
  } 


  int ncols = mypars->active.size();
  int nrows = std::stoi(init->smooth_minimizer["walkers"]);


  FILE* fh = fopen((path+run+"emcee_start.txt").c_str(),"w");
  for(int i=0;i<nrows;i++){
    for(int j=0;j<ncols;j++){
      double dum = drand48();
      fprintf(fh," %12.5f",mypars->active[j]->pri->fromUnitCube(dum));
    }
    fprintf(fh,"\n");
  }
  fclose(fh);


  delete(init);
  delete(mypars);
  return 0;
}
