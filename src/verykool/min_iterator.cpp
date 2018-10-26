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
    std::cout << i << std::endl;
    like_model->updateLikelihoodModel();
    //    like_model->printActive();
    like_model->getLogLike();
    like_model->printTerms();

    // Output the potential corrections at each step
    PertIterationLikelihood* specific_pointer = dynamic_cast<PertIterationLikelihood*>(like_model);
    specific_pointer->pert_mass_model->dpsi->outputSource(output + std::to_string(i) + "_perturbations_");
    

    // Output the cumulative total perturbations up to each step
    /*
    PertLikelihood* like_pointer = dynamic_cast<PertLikelihood*>(like_model);
    Pert* pert_pointer = dynamic_cast<Pert*>(like_pointer->collection->models.back());
    pert_pointer->dpsi->outputSource(output + std::to_string(i) + "_perturbations_");
    
    ImagePlane* kappa = new ImagePlane(pert_pointer->dpsi->Si,pert_pointer->dpsi->Sj,pert_pointer->dpsi->width,pert_pointer->dpsi->height);
    for(int j=0;j<kappa->Nm;j++){
      kappa->cells[j] = NULL;
      kappa->crosses[j] = NULL;
    }
    pert_pointer->getConvergence(kappa);
    kappa->writeImage(output + std::to_string(i) + "_convergence.fits");
    delete(kappa);
    */
  }

  /*
  std::cout << "final update" << std::endl;
  like_model->updateLikelihoodModel();
  //    like_model->printActive();
  like_model->getLogLike();
  like_model->printTerms();
  */
}

