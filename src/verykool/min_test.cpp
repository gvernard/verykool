#include "minimizers.hpp"

#include "likelihoodModels.hpp"

void Nothing::minimize(std::map<std::string,std::string> minimizer,BaseLikelihoodModel* mypars,const std::string output){
  mypars->updateLikelihoodModel();
  double dum = mypars->getLogLike();
}
