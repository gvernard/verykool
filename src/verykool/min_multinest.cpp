#include "minimizers.hpp"

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>

#include "nonLinearPars.hpp"
#include "likelihoodModels.hpp"




//Functions related to the calculation of the Bayesian evidence and parameter estimation using MultiNest
//==============================================================================================================================
void MultiNest::minimize(std::map<std::string,std::string> opt,BaseLikelihoodModel* pars,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,const std::string output){

  int myndims = pars->active.size();
  
  int IS         = 0;				// do Nested Importance Sampling?
  int mmodal     = 0;				// do mode separation?
  int ceff       = 0;				// run in constant efficiency mode?
  int nlive      = stoi(opt["nlive"]);		// number of live points
  double efr     = stof(opt["efr"]);		// set the required efficiency
  double tol     = stof(opt["tol"]);		// tol, defines the stopping criteria
  int ndims      = myndims;			// dimensionality (no. of free parameters)
  int nPar       = myndims;			// total no. of parameters including free & derived parameters
  int nClsPar    = myndims;			// no. of parameters to do mode separation on
  int updInt     = 10;  			// after how many iterations feedback is required & the output files should be updated
	                                        // note: posterior files are updated & dumper routine is called after every updInt*10 iterations
  double Ztol    = -1e90;			// all the modes with logZ < Ztol are ignored
  int maxModes   = 100;				// expected max no. of modes (used only for memory allocation)

  int pWrap[ndims];				// which parameters to have periodic boundary conditions?
  for(int i=0;i<ndims;i++){
    pWrap[i] = pars->active[i]->per;
  }

  std::string strroot = output + "mn-";         // root for output files
  char* root = new char[strroot.size() + 1];
  std::copy(strroot.begin(),strroot.end(),root);
  root[strroot.size()] = '\0';


  int seed       = stoi(opt["seed"]);		// random no. generator seed, if < 0 then take the seed from system clock
  int fb         = 1;				// need feedback on standard output?
  int resume     = 0;				// resume from a previous job?
  int outfile    = 1;				// write output files?
  int initMPI    = 0;				// initialize MPI routines? Set it to F if you want your main program to handle MPI initialization
  double logZero = -1e90;			// points with loglike < logZero will be ignored by MultiNest
  int maxiter    = stoi(opt["maxiter"]);	// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
                                                // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied

  //  void* context  = 0;			// not required by MultiNest, any additional information user wants to pass

  int counter = 0;
  extras myextras = {pars,image,source,mycollection,output,counter};

  nested::run(IS,mmodal,ceff,nlive,tol,efr,ndims,nPar,nClsPar,maxModes,updInt,Ztol,root,seed,pWrap,fb,resume,outfile,initMPI,logZero,maxiter,MultiNestLogLike,MultiNestDumper,&myextras);
}


/*
void LogLike(double* Cube,int& ndim,int& npars,double& lnew,void* myextras){
  extras* e = (extras*) myextras;
  double pi = 3.14159265359;

  std::map<std::string,double> means;
  means["b"]  = 4.87;
  means["q"]  = 0.6667;
  means["pa"] = 164;
  means["x0"] = 0;
  means["y0"] = 0;

  std::map<std::string,double> sdevs;
  sdevs["b"]  = 0.5;
  sdevs["q"]  = 0.05;
  sdevs["pa"] = 5;
  sdevs["x0"] = 0.1;
  sdevs["y0"] = 0.1;

  double val = 0;

  lnew = 0;
  for(int i=0;i<e->active.size();i++){
    e->nlpars[e->active[i].index][e->active[i].nam]->fromUnitCube(Cube[i]);
    val = e->nlpars[e->active[i].index][e->active[i].nam]->val;
    //    std::cout << means[e->active[i].nam] << " " ;
    std::cout << val << " " ;

    lnew += -pow(val-means[e->active[i].nam],2)/(2*pow(sdevs[e->active[i].nam],2)) - 0.5*log(2*pi*sdevs[e->active[i].nam]);
    //    lnew += -pow(val-means[e->active[i].nam],2);
  }
  std::cout << std::endl;
}
*/


// Input arguments:
//    ndim 						= dimensionality (total number of free parameters) of the problem
//    npars 						= total number of free plus derived parameters
//    context						= void pointer, any additional information
//
// Input/Output arguments:
//    Cube[npars] 					= on entry has the ndim parameters in unit-hypercube. On exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//
// Output arguments:
//    lnew 						= log(likelihood)
void MultiNestLogLike(double* Cube,int& ndim,int& npars,double& lnew,void* myextras){
  extras* e = (extras*) myextras; // need to cast back void pointer

  //Update values for active parameters
  std::vector<double> new_pars;
  for(int i=0;i<e->pars->active.size();i++){
    new_pars.push_back( e->pars->active[i]->pri->fromUnitCube(Cube[i]) );
  }
  e->pars->updateActive(new_pars);
  e->pars->updateLikelihoodModel(e->image,e->source,e->mycollection);
  lnew = e->pars->getLogLike(e->image,e->source);
  e->pars->printActive();
  e->pars->printTerms();
}





// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
// Arguments:
//
//    nSamples 						   = total number of samples in posterior distribution
//    nlive 						   = total number of live points
//    nPar 						   = total number of parameters (free + derived)
//    physLive[1][nlive * (nPar + 1)] 			   = 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
//    posterior[1][nSamples * (nPar + 2)] 		   = posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
//
//    paramConstr[1][4*nPar]:
//       paramConstr[0][0] to paramConstr[0][nPar - 1] 	      = mean values of the parameters
//       paramConstr[0][nPar] to paramConstr[0][2*nPar - 1]   = standard deviation of the parameters
//       paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
//       paramConstr[0][nPar*3] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
//
//    maxLogLike					   = maximum loglikelihood value
//    logZ						   = log evidence value from the default (non-INS) mode
//    INSlogZ						   = log evidence value from the INS mode
//    logZerr						   = error on log evidence value
//    context						   = void pointer, any additional information
void MultiNestDumper(int& nSamples,int& nlive,int& nPar,double** physLive,double** posterior,double** paramConstr,double& maxLogLike,double& logZ,double& INSlogZ,double& logZerr,void* myextras){
  // the posterior distribution postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
  extras* e = (extras*) myextras; // need to cast back void pointer
  FILE* fh;

  e->counter++;
  e->pars->outputLikelihoodModel(e->image,e->source,e->output + std::to_string(e->counter) + "_");


  // lastlive holds the parameter values for the last set of live points, and has the same structure as nlpars, with an array of nlive points for each parameter
  std::vector<double*> lastlive(nPar);
  for(int i=0;i<nPar;i++){
    lastlive[i] = (double*) calloc(nlive,sizeof(double));
    for(int j=0;j<nlive;j++){
      lastlive[i][j] = e->pars->active[i]->pri->fromUnitCube( physLive[0][i*nlive+j] );
    }
  }

  // File with the last live points
  fh = fopen( (e->output+"lastlive.txt").c_str() ,"w");
  for(int j=0;j<nlive;j++){
    for(int i=0;i<nPar;i++){
      fprintf(fh," %25.18e",lastlive[i][j]);
    }
    fprintf(fh,"\n");
  }
  fclose(fh);



  // postdist holds the posterior probability, loglike, and the parameters
  double postdist[nSamples][nPar+2];

  for(int j=0;j<nSamples;j++){
    if( posterior[0][(nPar+1)*nSamples+j] > 1.e-99 ){
      postdist[j][0] = posterior[0][(nPar+1)*nSamples+j];//post goes first
      //    postdist[j][1] = -2*log10(posterior[0][0*nSamples+j]);//and then -2*loglikelihood
      postdist[j][1] = -posterior[0][(nPar)*nSamples+j];//and then -2*loglikelihood
      
      for(int i=0;i<nPar;i++){
	postdist[j][2+i] = e->pars->active[i]->pri->fromUnitCube( posterior[0][i*nSamples+j] );
      }
    }
  }

  // File for the corner plot
  fh = fopen( (e->output + std::to_string(e->counter) + "_corner.txt").c_str() ,"w");
  for(int j=0;j<nSamples;j++){
    if( postdist[j][0] > 1.e-99 ){
      double dum = 1.;
      fprintf(fh," %25.18e",dum);
      fprintf(fh," %25.18e",postdist[j][0]);
      for(int i=2;i<nPar+2;i++){
	fprintf(fh," %25.18e",postdist[j][i]);
      }
      fprintf(fh,"\n");
    }
  }
  fclose(fh);

  std::ifstream  src((e->output + std::to_string(e->counter) + "_corner.txt").c_str(),std::ios::binary);
  std::ofstream  dst((e->output + "corner.txt").c_str(), std::ios::binary);
  dst << src.rdbuf();



  // map_pars holds the maximum a-posteriori parameters and has the same structure as nlpars
  double* means = (double*) calloc(nPar,sizeof(double));
  double* sdevs = (double*) calloc(nPar,sizeof(double));
  double* bests = (double*) calloc(nPar,sizeof(double));
  double* maps  = (double*) calloc(nPar,sizeof(double));
  for(int i=0;i<nPar;i++){
    means[i] = e->pars->active[i]->pri->fromUnitCube( paramConstr[0][i] );
    sdevs[i] = paramConstr[0][nPar+i] * e->pars->active[i]->ran;
    bests[i] = e->pars->active[i]->pri->fromUnitCube( paramConstr[0][i] );
    maps[i]  = e->pars->active[i]->pri->fromUnitCube( paramConstr[0][i] );
    //    (*(e->map_pars))[e->active[i].index][e->active[i].nam] = e->nlpars[e->active[i].index][e->active[i].nam]->val;  
  }

  /*
  for(int i=2;i<nPar+2;i++){
    double mean = 0.0;
    double sdev = 0.0;
    for(int j=0;j<nSamples;j++){
      mean += postdist[j][0]*postdist[j][i];
    }
    for(int j=0;j<nSamples;j++){
      sdev += postdist[j][0]*pow(postdist[j][i]-mean,2);
    }
    sdev = sqrt(sdev);
    std::cout << mean << " " << sdev << std::endl;
  }
  */

  // File with the mean, sdev, best-fit, and MAP parameters
  fh = fopen( (e->output+"vkl_stats.txt").c_str() ,"w");
  for(int i=0;i<nPar;i++){
    fprintf(fh,"%25s %25.18e %25.18e %25.18e %25.18e\n",e->pars->active[i]->nam.c_str(),means[i],sdevs[i],bests[i],maps[i]);
  }
  fclose(fh);






  free(means);
  free(sdevs);
  free(bests);
  free(maps);
}


void MultiNest::output(){
  // MultiNest completed
}
