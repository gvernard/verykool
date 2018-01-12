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


  Initialization* init = new Initialization();
  init->parseInputJSON(path,run);
  BaseLikelihoodModel* mypars = FactoryLikelihoodModel::getInstance()->createLikelihoodModel(path,run);

  if( init->minimizer["type"] != "cosmosis_emcee" ){
    std::cout << "wrong minimizer, it has to be cosmosis_emcee" << std::endl;
    return 0;
  } 


  int ncols = mypars->active.size();
  int nrows = stoi(init->minimizer["walkers"]);


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
