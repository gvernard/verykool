#include "minimizers.hpp"
#include "likelihoodModels.hpp"
#include "nonLinearPars.hpp"

void Nothing::minimize(std::map<std::string,std::string> minimizer,BaseLikelihoodModel* mypars,const std::string output){
  mypars->updateLikelihoodModel();
  double dum = mypars->getLogLike();
  
  for(int i=0;i<mypars->active.size();i++){
    mypars->maps[i]    = mypars->active[i]->val;
    mypars->means[i]   = mypars->active[i]->val;
    mypars->s1_low[i]  = 0;
    mypars->s1_high[i] = 0;
  }
}

void Nothing::output(std::string output){
}
