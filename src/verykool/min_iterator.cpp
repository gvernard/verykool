#include "minimizers.hpp"
#include "massModels.hpp"
#include "likelihoodModels.hpp"

#include <iostream>
#include <string>

void Iterator::minimize(std::map<std::string,std::string> opt,BaseLikelihoodModel* like_model,const std::string output){
  int maxiter = std::stoi(opt["maxiter"]);
  std::cout << "Doing " << maxiter << " iterations" << std::endl;

  // This is just to fill the containers of the 'best' parameters
  for(int i=0;i<like_model->active.size();i++){
    like_model->means[i] = like_model->active[i]->val;
    like_model->sdevs[i] = 0;
    like_model->bests[i] = like_model->active[i]->val;
    like_model->maps[i]  = like_model->active[i]->val;
  }

  for(int i=0;i<maxiter;i++){
    /*
    if( i == (int) floor(maxiter/2.0) ){
      double l_old,l_new;
      for(int k=0;k<like_model->active.size();k++){
	if( like_model->active[k]->nam == "lambda_dpsi" ){
	  l_old = like_model->active[k]->val;
	  l_new = l_old/10.0;
	  like_model->active[k]->val = l_new;
	  break;
	}
      }
      std::cout << "Changing the value of lambda_dpsi from " << l_old << " to " << l_new << std::endl;
    }
    */    

    std::cout << i << std::endl;
    like_model->updateLikelihoodModel();
    //    like_model->printActive();
    like_model->getLogLike();
    like_model->printTerms();
    like_model->outputLikelihoodModel(output + std::to_string(i) + "_");    
  }

  /*
  std::cout << "final update" << std::endl;
  like_model->updateLikelihoodModel();
  //    like_model->printActive();
  like_model->getLogLike();
  like_model->printTerms();
  like_model->outputLikelihoodModel(output);
  */
}

