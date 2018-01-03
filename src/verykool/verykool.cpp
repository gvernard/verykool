#include <iostream>
#include <string>
#include <vector>
#include <mpi.h>

#include "inputOutput.hpp"
#include "imagePlane.hpp"
#include "massModels.hpp"
#include "sourcePlane.hpp"
#include "tableAlgebra.hpp"
#include "minimizers.hpp"




int main(int argc,char* argv[]){
  // Check command line arguments
  if( argc != 3 ){
    std::cout << "2 command line arguments are required!!!" << std::endl;
    return 0;
  }


  // Initialize MPI
  int nprocs,myrank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);



  // Initialize variables
  Initialization* init;
  ImagePlane* mydata;
  CollectionMassModels* mycollection;
  BaseSourcePlane* mysource;
  std::vector<myactive> active;
  mymatrices* matrices;
  precomp* pcomp;
  initialize_program(argv[1],argv[2],init,mydata,mycollection,mysource,active,matrices,pcomp);


  //=============== BEGIN:MINIMIZATION =========================
  printf("%-25s","Starting minimization");
  fflush(stdout);
 
  BaseMinimizer* myminimizer = FactoryMinimizer::getInstance()->createMinimizer(init->minimizer,mydata,mysource,mycollection,init->nlpars,matrices,pcomp,init->output);
  myminimizer->minimize(init->minimizer,mydata,mysource,mycollection,init->nlpars,matrices,pcomp,init->output);

  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);
  //================= END:MINIMIZATION =========================


  // Finalize output etc
  if( myrank == 0 ){
    finalize_program(init,mydata,mycollection,mysource,matrices,pcomp);
    myminimizer->output();
  }
  


  // Cleanup pointers
  delete(init);
  delete(mydata);
  delete(mycollection);
  delete(mysource);
  delete(matrices);
  delete(pcomp);
  delete(myminimizer);


  // Finalize MPI
  MPI_Finalize();


  return 0;
}
