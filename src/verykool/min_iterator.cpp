#include "minimizers.hpp"
#include "massModels.hpp"
#include "likelihoodModels.hpp"

#include <iostream>
#include <string>

void Iterator::minimize(std::map<std::string,std::string> opt,BaseLikelihoodModel* like_model,const std::string output){
  int maxiter = std::stoi(opt["maxiter"]);
  std::cout << "Doing " << maxiter << " iterations" << std::endl;

  for(int i=0;i<maxiter;i++){
    std::cout << i << std::endl;
    like_model->updateLikelihoodModel();

    // Output the potential corrections at each step
    //    PerturbationsLikelihood* specific_pointer = dynamic_cast<PerturbationsLikelihood*>(like_model);
    //    specific_pointer->pert_mass_model->dpsi->outputSource(output + std::to_string(i) + "_perturbations_");

    // Output the cumulative total perturbations up to each step
    PerturbationsLikelihood* like_pointer = dynamic_cast<PerturbationsLikelihood*>(like_model);
    Pert* pert_pointer = dynamic_cast<Pert*>(like_pointer->collection->models.back());
    pert_pointer->dpsi->outputSource(output + std::to_string(i) + "_perturbations_");
    
  }

}

