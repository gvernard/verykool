#include <vector>
#include <cstdlib>

#include "cosmosis/datablock/datablock.hh"

#include "likelihoodModels.hpp"
#include "inputOutput.hpp"
#include "imagePlane.hpp"
#include "massModels.hpp"
#include "sourcePlane.hpp"




class vkl_vars {
public:
  Initialization* init;
  ImagePlane* image;
  BaseSourcePlane* source;
  CollectionMassModels* collection;
  BaseLikelihoodModel* smooth_like;
  BaseLikelihoodModel* pert_like;
  Pert* pert_mass_model;

  std::string cosmosis_output;
  std::string path;
  std::string run;

  vkl_vars(){};
  ~vkl_vars(){};
};




extern "C" {


  
  void* setup(cosmosis::DataBlock* options){
    //DATABLOCK_STATUS status = DBS_SUCCESS;
    //DATABLOCK_STATUS status = 0;
    int status = 0;

    //    vkl_vars* dum = (vkl_vars*) malloc(sizeof(vkl_vars));
    vkl_vars* dum = new vkl_vars();

    std::string path;
    std::string run;
    status |= options->get_val("verykool","path",dum->path);
    status |= options->get_val("verykool","run",dum->run);
    status |= options->get_val("output","filename",dum->cosmosis_output);
    

    Initialization::initialize_program(dum->path,dum->run,dum->init,dum->smooth_like,dum->image,dum->collection,dum->source,dum->pert_like,dum->pert_mass_model);
    

    printf("%-25s","Starting minimization");
    printf("%6s%9s\n","using ",dum->init->smooth_minimizer["type"].c_str());
    fflush(stdout);

    return (void*) dum;
  }


  
  DATABLOCK_STATUS execute(cosmosis::DataBlock* block,void* config){
    
    //    DATABLOCK_STATUS status = DBS_SUCCESS;
    int status = 0;

    vkl_vars* dum = (vkl_vars*) config;


    std::vector<std::string> full_names = dum->smooth_like->getActiveFullNames();
    std::vector<double> values;
    for(int i=0;i<full_names.size();i++){
      double val;
      status |= block->get_val("nlpars",full_names[i].c_str(),val);
      values.push_back(val);
    }

    dum->smooth_like->updateActive(values);
    dum->smooth_like->updateLikelihoodModel();
    double like = dum->smooth_like->getLogLike();
    
    
    //double like = 99.99;
    
    status |= block->put_val("LIKELIHOODS","verykool_like",like);
    
    return (DATABLOCK_STATUS) status;
  }

  
  
  int cleanup(void* config){
    printf("%+7s\n","...done");
    std::cout << std::string(200,'=') << std::endl;
    fflush(stdout);
    
    vkl_vars* dum = (vkl_vars*) config;
    Initialization::finalize_smooth(dum->init,dum->smooth_like,dum->image,dum->collection,dum->source);


    // Rewrite the cosmosis EMCEE output to the correct format for a corner plot
    std::cout << dum->init->smooth_minimizer["type"] << std::endl;
    if( dum->init->smooth_minimizer["type"] == "cosmosis_emcee" ){
      std::ifstream mystream(dum->cosmosis_output+".txt");
      std::string line;
      FILE* fh     = fopen((dum->path+dum->run+"output/corner.txt").c_str(),"w");
      int ncols    = dum->smooth_like->active.size()+1;
      double* vals = (double*) malloc(ncols*sizeof(double));
      while( std::getline(mystream,line) ){
	//	std::cout <<  line.substr(0,20) << std::endl;
	if( line.substr(0,1).compare("#") != 0 ){
	  std::istringstream ss(line);
	  for(int i=0;i<ncols;i++){
	    ss >> vals[i];
	  }
	  fprintf(fh," 1 %12.4f",vals[ncols-1]);
	  for(int i=0;i<ncols-1;i++){
	    fprintf(fh," %12.4f",vals[i]);
	  }
	  fprintf(fh,"\n");
	}
      }
      fclose(fh);
      free(vals);
    }

    delete(dum);
  }


  
}
