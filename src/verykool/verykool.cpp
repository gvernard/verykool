#include <iostream>
#include <string>
#include <vector>
#include <mpi.h>

#include "inputOutput.hpp"
#include "imagePlane.hpp"
#include "massModels.hpp"
#include "sourcePlane.hpp"
#include "likelihoodModels.hpp"
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



  //=============== BEGIN:INITIALIZATION =======================
  // Initialize variables
  Initialization* init = 0;
  ImagePlane* mydata = 0;
  BaseSourcePlane* mysource = 0;
  CollectionMassModels* mycollection = 0;
  BaseLikelihoodModel* smooth_like = 0;
  BaseLikelihoodModel* pert_like = 0;
  Pert* pert_mass_model = 0;
  BaseMinimizer* smooth_minimizer = 0;

  Initialization::initialize_program(argv[1],argv[2],init,smooth_like,mydata,mycollection,mysource,pert_like,pert_mass_model);
  //================= END:INITIALIZATION =======================






  /*
  Pert* test = new Pert(40,40,4.05,4.05,"curvature");
  double* dpsi = (double*) malloc(test->dpsi->Sm*sizeof(double));
  for(int i=0;i<test->dpsi->Sm;i++){
    dpsi[i] = (static_cast <double> (rand()) / static_cast <double> (RAND_MAX))*0.2 - 0.1;
  }
  test->replaceDpsi(dpsi);
  test->updatePert();
  test->createAint(mydata);
  free(dpsi);

  double xdefl,ydefl;
  for(int i=0;i<10;i++){
    test->defl(mydata->x[i],mydata->y[i],xdefl,ydefl);
    printf("%8.3f%8.3f%8.3f%8.3f\n",mydata->x[i],mydata->y[i],xdefl,ydefl);
  }
  

  double* xdeflarr = (double*) malloc(2*mydata->Nm*sizeof(double));
  double* ydeflarr = (double*) malloc(2*mydata->Nm*sizeof(double));
  test->tableDefl(mydata->Nm,xdeflarr,ydeflarr);
  for(int i=0;i<10;i++){
    printf("%8.3f%8.3f%8.3f%8.3f\n",mydata->x[i],mydata->y[i],xdeflarr[i],ydeflarr[i]);
  }
  free(xdeflarr);
  free(ydeflarr);


  MPI_Finalize();
  return 0;
  */



  //=============== BEGIN:SMOOTH MODEL =========================
  // Initial output
  if( myrank == 0 ){
    smooth_like->initialOutputLikelihoodModel(init->output);
    //myminimizer->output();
  }

  printf("%-25s","Starting smooth minimization ");
  fflush(stdout);
 
  smooth_minimizer = FactoryMinimizer::getInstance()->createMinimizer(init->smooth_minimizer,smooth_like,init->output);
  smooth_minimizer->minimize(init->smooth_minimizer,smooth_like,init->output);

  printf("%+7s\n","...done");
  std::cout << std::string(200,'=') << std::endl;
  fflush(stdout);

  // Finalize output etc
  if( myrank == 0 ){
    Initialization::finalize_smooth(init->output,smooth_like);
    //myminimizer->output();
  }

  delete(smooth_minimizer);
  //================= END:SMOOTH MODEL =========================





  /*
  //  mycollection->all_defl(mydata);
  //  mysource->createInterpolationWeights(mydata);
  mysource->constructDs(mydata);

  for(int i=0;i<mysource->Sm;i++){
    //  std::cout << mysource->s_dx[i] << std::endl;
    mysource->src[i] = hypot(mysource->s_dx[i],mysource->s_dy[i]);
  }
  mysource->outputSource(init->output);

  MPI_Finalize();
  return 0;
  */





  //=============== BEGIN:PERTURBATIONS =========================
  if( init->perturbations.size() > 0 ){
    if( myrank == 0 ){
      pert_like->initialOutputLikelihoodModel(init->output);
    }

    //Initialize perturbations
    PerturbationsLikelihood* specific_pointer = dynamic_cast<PerturbationsLikelihood*>(pert_like);
    specific_pointer->initializePert(smooth_like);

    printf("%-25s","Starting perturbation minimization ");
    fflush(stdout);
    
    BaseMinimizer* pert_minimizer = FactoryMinimizer::getInstance()->createMinimizer(init->pert_minimizer,pert_like,init->output);
    pert_minimizer->minimize(init->pert_minimizer,pert_like,init->output);

    printf("%+7s\n","...done");
    std::cout << std::string(200,'=') << std::endl;
    fflush(stdout);

    
    if( myrank == 0 ){
      Initialization::finalize_pert(init->output,pert_like);
      //myminimizer->output();
    }
    
    delete(pert_minimizer);
  }
  //================= END:PERTURBATIONS =========================




  // Cleanup pointers
  delete(init);
  delete(mydata);
  delete(mysource);
  delete(mycollection);
  delete(smooth_like);
  delete(pert_like);
  delete(pert_mass_model);
  

  // Finalize MPI
  MPI_Finalize();


  return 0;
}
