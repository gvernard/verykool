#ifndef MINUIT_FUNCS_HPP
#define MINUIT_FUNCS_HPP

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "Minuit2/FCNBase.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnPrint.h"

#include "massModels.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "tableAlgebra.hpp"
#include "nonLinearPars.hpp"

class myMinuitFCN : public ROOT::Minuit2::FCNBase {
public:
  myMinuitFCN(const std::vector<double> pars) {};
  ~myMinuitFCN() {}
  virtual double Up() const {return 0.1;}
  virtual double operator()(const std::vector<double>&) const;
  
  ImagePlane* image;
  BaseSourcePlane* source;
  CollectionMassModels* collection;
  std::vector<std::map<std::string,BaseNlpar*> > nlpars;
  mymatrices* matrices;
  precomp* pcomp;
};
void myMinuit(ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,std::vector<std::map<std::string,BaseNlpar*> > nlpars,mymatrices* mat,precomp* pcomp,const std::string output,std::map<std::string,std::string> opt);


#endif /* MINUIT_FUNCS_HPP */
