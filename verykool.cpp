#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
#include <mpi.h>

#include "inputOutput.hpp"
#include "nonLinearPars.hpp"
#include "imagePlane.hpp"
#include "massModels.hpp"
#include "sourcePlane.hpp"
#include "tableAlgebra.hpp"
#include "minimizers.hpp"




int main(int argc,char* argv[]){
  if( argc != 3 ){
    std::cout << "2 command line arguments are required!!!" << std::endl;
    return 0;
  }
  int nprocs,myrank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);




  
  //=============== BEGIN:INITIALIZATION =======================
  printf("%-50s","Starting initialization");
  fflush(stdout);
  


  Initialization init;
  init.parseInputJSON(argv[1],argv[2]);



  //Read image data -----------------------------------------------------------------------------------------------------------------------------------------------
  ImagePlane my_data(init.imgpath,stoi(init.image["pix_x"]),stoi(init.image["pix_y"]),stof(init.image["siz_x"]),stof(init.image["siz_y"]));



  //Initialize mass model parameters ------------------------------------------------------------------------------------------------------------------------------
  CollectionMassModels* mycollection = new CollectionMassModels(init.nlpars[0]);
  mycollection->models.resize(init.mmodel.size());
  for(int k=0;k<mycollection->models.size();k++){
    mycollection->models[k] = FactoryMassModel::getInstance()->createMassModel(init.mmodel[k],init.nlpars[2+k]);
  }



  //Initialize source ---------------------------------------------------------------------------------------------------------------------------------------------
  //initialize here, but also in each iteration for an adaptive grid                                                  [<---iteration dependent for adaptive source]
  BaseSourcePlane* mysource = FactorySourcePlane::getInstance()->createSourcePlane(init.source);

  if( mysource->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(mysource);
    ada->createAdaGrid(&my_data,mycollection);
    ada->createDelaunay();
  }



  //Initialize matrices: S,C,B,H,L --------------------------------------------------------------------------------------------------------------------------------
  //This is independent of the specific linear algebra package used.
  //S: mask,
  //C: covariance matrix,
  //B: blurring matrix,
  //L: lensing matrix,
  //H: source regularization matrix
  mymatrices matrices;
  matrices.initMatrices(&my_data,mysource,init.maskpath,init.noise_flag,init.covpath,init.psfpath,init.psf);



  //Initialize precomputed algebraic quantities -------------------------------------------------------------------------------------------------------------------
  //Precompute a number of algebraic quantities (table products etc)                                                                [<---algebra package dependent]
  precomp pcomp(&my_data,mysource);
  setAlgebraInit(&my_data,mysource,&matrices,&pcomp);



  outputInitial(init.nlpars,init.output);



  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);
  //================= END:INITIALIZATION =======================





  //=============== BEGIN:MINIMIZATION =========================
  printf("%-25s","Starting minimization");
  fflush(stdout);
 
  BaseMinimizer* myminimizer = FactoryMinimizer::getInstance()->createMinimizer(init.minimizer,&my_data,mysource,mycollection,init.nlpars,&matrices,&pcomp,init.output);
  myminimizer->minimize(init.minimizer,&my_data,mysource,mycollection,init.nlpars,&matrices,&pcomp,init.output);

  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);
  //================= END:MINIMIZATION =========================
  




  //=================== BEGIN:OUTPUT ===========================
  //==> beginning of code executed only by MPI process with id 0
  if( myrank == 0 ){
    printf("%-25s\n","Starting output");
    fflush(stdout);
    
    // call mainLogLike to set the values of all the variables,matrices, etc, before the output
    double logL= mainLogLike(&my_data,mysource,mycollection,init.nlpars,&matrices,&pcomp);
    outputGeneric(&my_data,mysource,init.nlpars,&pcomp,init.output);
    myminimizer->output();
        
    printf("%+7s\n","...done");
    std::cout << std::string(200,'=') << std::endl;
    fflush(stdout);  
  }
  //==> end of code executed only by process with id 0
  //===================== END:OUTPUT ===========================
  

  


  delete(mycollection);
  delete(mysource);
  delete(myminimizer);
  MPI_Finalize();
  return 0;
}
