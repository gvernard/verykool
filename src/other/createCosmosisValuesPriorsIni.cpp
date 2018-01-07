#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>

#include "nonLinearPars.hpp"
#include "parameterModels.hpp"

int main(int argc,char* argv[]){
  std::string path(argv[1]);
  std::string run(argv[2]);


  BaseParameterModel* mypars = FactoryParameterModel::getInstance()->createParameterModel(path,run);

  std::vector<std::string> names = mypars->getActiveFullNames();
  
  FILE* fh = fopen((path+run+"cosmosis_values.ini").c_str(),"w");
  fprintf(fh,"[nlpars]\n");
  for(int i=0;i<names.size();i++){
    fprintf(fh,"%s = %f %f %f\n",names[i].c_str(),mypars->active[i]->min,mypars->active[i]->val,mypars->active[i]->max);
  }
  fclose(fh);
  

  FILE* fh2 = fopen((path+run+"cosmosis_priors.ini").c_str(),"w");
  fprintf(fh2,"[nlpars]\n");
  for(int i=0;i<names.size();i++){
    std::map<std::string,double> pri_pars = mypars->active[i]->pri->getPars();
    if( mypars->active[i]->pri->type == "uni" ){
      fprintf(fh2,"%s = %s %f %f\n",names[i].c_str(),"uniform",mypars->active[i]->min,mypars->active[i]->max);
    } else if( mypars->active[i]->pri->type == "exp" ){
      fprintf(fh2,"%s = %s %f\n",names[i].c_str(),"exponential",pri_pars["beta"]);
    } else if( mypars->active[i]->pri->type == "gauss" ){
      fprintf(fh2,"%s = %s %f %f\n",names[i].c_str(),"gaussian",pri_pars["mean"],pri_pars["sdev"]);
    }
  }
  fclose(fh2);


  delete(mypars);
  return 0;
}
