#ifndef MAIN_LOG_LIKE_HPP
#define MAIN_LOG_LIKE_HPP

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>


#include "inputOutput.hpp"
#include "massModels.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "tableAlgebra.hpp"
#include "nonLinearPars.hpp"




static double mainLogLike(ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,std::vector<std::map<std::string,BaseNlpar*> > nlpars,mymatrices* mat,precomp* pcomp){
  for(int i=0;i<mycollection->models.size();i++){
    mycollection->models[i]->setMassPars(nlpars[2+i]);
  }
  mycollection->setPhysicalPars(nlpars[0]);

  

  if( source->sample_reg ){
    source->kernel->setParameters(nlpars[1]);
  }

  if( source->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(source);
    ada->createAdaGrid(image,mycollection);
    ada->createDelaunay();
  }

  if( source->sample_reg || source->type == "adaptive" ){
    source->constructH(&mat->H);
  }



  //  std::cout << "1" << std::endl;
  source->constructL(image,mycollection,&mat->L);
  //  std::cout << "2" << std::endl;
  setAlgebraRuntime(image,source,nlpars[1],mat,pcomp);
  //  std::cout << "3" << std::endl;
  solveLinearSparseS(image,source,pcomp);
  //  std::cout << "4" << std::endl;
  double loglike = getLogLike(image,source,pcomp,nlpars);
  //  std::cout << "5" << std::endl;
  return loglike;


  //  return pow(nlpars[2]["b"]->val-6,2) + pow(nlpars[2]["q"]->val-6,2) + pow(nlpars[2]["pa"]->val-6,2) + pow(nlpars[2]["x0"]->val-6,2) + pow(nlpars[2]["y0"]->val-6,2);
}

#endif /* MAIN_LOG_LIKE_HPP */
