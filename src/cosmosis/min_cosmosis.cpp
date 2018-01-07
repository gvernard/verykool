#include <vector>
#include <cstdlib>

#include "cosmosis/datablock/datablock.hh"

#include "parameterModels.hpp"
#include "inputOutput.hpp"
#include "imagePlane.hpp"
#include "massModels.hpp"
#include "sourcePlane.hpp"
#include "tableAlgebra.hpp"




class vkl_vars {
public:
  Initialization* init;
  BaseParameterModel* mypars;
  ImagePlane* mydata;
  CollectionMassModels* mycollection;
  BaseSourcePlane* mysource;
  mymatrices* matrices;
  precomp* pcomp;

  vkl_vars(){};
  ~vkl_vars(){};
};




extern "C" {


  
  void* setup(cosmosis::DataBlock* options){
    //DATABLOCK_STATUS status = DBS_SUCCESS;
    //DATABLOCK_STATUS status = 0;
    int status = 0;

    std::string path;
    std::string run;
    status |= options->get_val("verykool","path",path);
    status |= options->get_val("verykool","run",run);
    
    //    vkl_vars* dum = (vkl_vars*) malloc(sizeof(vkl_vars));
    vkl_vars* dum = new vkl_vars();

    Initialization::initialize_program(path,run,dum->init,dum->mypars,dum->mydata,dum->mycollection,dum->mysource,dum->matrices,dum->pcomp);


    printf("%-25s","Starting minimization");
    printf("%6s%9s\n","using ",dum->init->minimizer["type"].c_str());
    fflush(stdout);

    return (void*) dum;
  }


  
  DATABLOCK_STATUS execute(cosmosis::DataBlock* block,void* config){
    
    //    DATABLOCK_STATUS status = DBS_SUCCESS;
    int status = 0;

    vkl_vars* dum = (vkl_vars*) config;


    std::vector<std::string> full_names = dum->mypars->getActiveFullNames();
    std::vector<double> values;
    for(int i=0;i<full_names.size();i++){
      double val;
      status |= block->get_val("nlpars",full_names[i].c_str(),val);
      values.push_back(val);
    }

    dum->mypars->updateActive(values);
    dum->mypars->updateParameterModel(dum->mydata,dum->mysource,dum->mycollection,dum->matrices,dum->pcomp);
    double like = dum->mypars->getLogLike(dum->mydata,dum->mysource,dum->pcomp);
    
    
    //double like = 99.99;
    
    status |= block->put_val("LIKELIHOODS","verykool_like",like);
    
    return (DATABLOCK_STATUS) status;
  }

  
  
  int cleanup(void* config){
    printf("%+7s\n","...done");
    std::cout << std::string(200,'=') << std::endl;
    fflush(stdout);
    
    vkl_vars* dum = (vkl_vars*) config;
    Initialization::finalize_program(dum->init,dum->mypars,dum->mydata,dum->mycollection,dum->mysource,dum->matrices,dum->pcomp);
    delete(dum);
  }


  
}
