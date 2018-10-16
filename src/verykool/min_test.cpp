#include "minimizers.hpp"
#include "likelihoodModels.hpp"
#include "nonLinearPars.hpp"

void Nothing::minimize(std::map<std::string,std::string> minimizer,BaseLikelihoodModel* mypars,const std::string output){
  mypars->updateLikelihoodModel();
  double dum = mypars->getLogLike();
  
  for(int i=0;i<mypars->active.size();i++){
    mypars->means[i] = mypars->active[i]->val;
    mypars->sdevs[i] = 0;
    mypars->bests[i] = mypars->active[i]->val;
    mypars->maps[i]  = mypars->active[i]->val;
  }
}
