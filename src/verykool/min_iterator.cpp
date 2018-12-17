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
    like_model->maps[i]    = like_model->active[i]->val;
    like_model->means[i]   = like_model->active[i]->val;
    like_model->s1_low[i]  = 0;
    like_model->s1_high[i] = 0;
  }

  this->output_counter = 0;
  for(int i=0;i<maxiter;i++){
    this->iterations = i;
    std::cout << i << std::endl;

    /*    
    if( ((i+1) % 60) == 0 ){
      Nlpar* lambda = Nlpar::getParByName("lambda_dpsi",like_model->active);
      double l_old = lambda->val;
      double l_new = l_old*0.1;
      lambda->val = l_new;
      std::cout << "Changing the value of lambda_dpsi from " << l_old << " to " << l_new << std::endl;
    }
    */

    like_model->updateLikelihoodModel();
    like_model->printActive();
    like_model->getLogLike();
    like_model->printTerms();

    if( ((i+1) % 10) == 0 ){
      std::cout << "Writing output" << std::endl;
      like_model->outputLikelihoodModel(output + std::to_string(this->output_counter) + "_");
      this->output_counter++;
    }
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

void Iterator::output(std::string output){
  Json::Value min_out;
  min_out["iterations"] = this->iterations;

  std::ofstream jsonfile(output + std::to_string(this->output_counter) + "_" + this->name + "_minimizer_output.json");
  jsonfile << min_out;
  jsonfile.close();
}

void Iterator::finalizeMinimizer(std::string output){
  std::ifstream src2((output + std::to_string(this->output_counter) + this->name + "_minimizer_output.json").c_str(),std::ios::binary);
  std::ofstream dst2((output + this->name + "_minimizer_output.json").c_str(), std::ios::binary);
  dst2 << src2.rdbuf();
}
