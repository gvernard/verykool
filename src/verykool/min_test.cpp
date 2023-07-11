#include "minimizers.hpp"
#include "likelihoodModels.hpp"
#include "nonLinearPars.hpp"

void Nothing::minimize(std::map<std::string,std::string> minimizer,BaseLikelihoodModel* lmodel,const std::string output){
}

void Nothing::outputMinimizer(std::string output){
}

void Nothing::finalizeMinimizer(std::string output,BaseLikelihoodModel* lmodel){
  printf("%-16s%20s\n","Starting output ",this->name.c_str());
  fflush(stdout);


  Json::Value min_out;
  min_out["minimizer"] = this->name;
  
  // Write the mean, sdev, maxlike, and MAP values for each parameter
  // Calculating the mean and the lower and upper 1-sigma bounds
  Json::Value parameters,par_json;
  double mean,sdev,maxlike,map;
  std::string name;
  std::vector<double> map_pars;
  for(int i=0;i<lmodel->active.size();i++){
    name    = lmodel->active_names[i];
    mean    = lmodel->active[i]->val;
    sdev    = 0.0;
    maxlike = lmodel->active[i]->val;
    map     = lmodel->active[i]->val;
    par_json["mean"] = mean;
    par_json["sdev"] = sdev;
    par_json["maxlike"] = maxlike;
    par_json["map"] = map;
    parameters[name] = par_json;
    map_pars.push_back( map );
  }
  min_out["parameters"] = parameters;
  
  std::ofstream jsonfile(output + lmodel->name + "_minimizer_output.json");
  jsonfile << min_out;
  jsonfile.close();

  
  lmodel->updateLikelihoodModel();
  lmodel->getLogLike();
  lmodel->printActive();
  lmodel->printTerms();
  printf("\n");
  lmodel->outputLikelihoodModel(output);

  printf("%7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);
}
