#include "minimizers.hpp"

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>

#include "multinest.h"
#include "json/json.h"

#include "nonLinearPars.hpp"
#include "likelihoodModels.hpp"




//Functions related to the calculation of the Bayesian evidence and parameter estimation using MultiNest
//==============================================================================================================================
void MultiNest::minimize(std::map<std::string,std::string> opt,BaseLikelihoodModel* pars,const std::string output){
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

  //  std::cout << strroot << std::endl;
  //  std::cout << root << std::endl;

  int seed       = stoi(opt["seed"]);		// random no. generator seed, if < 0 then take the seed from system clock
  int fb         = 1;				// need feedback on standard output?
  int resume     = 0;				// resume from a previous job?
  int outfile    = 1;				// write output files?
  int initMPI    = 0;				// initialize MPI routines? Set it to F if you want your main program to handle MPI initialization
  double logZero = -1e90;			// points with loglike < logZero will be ignored by MultiNest
  int maxiter    = stoi(opt["maxiter"]);	// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
                                                // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied

  //  void* context  = 0;			// not required by MultiNest, any additional information user wants to pass

  this->counter = 0;
  extras myextras = {pars,output,this};


  // Initial output from the minimizer:
  // File with the parameter names
  std::ofstream f_names(output + this->name + "_postdist.paramnames",std::ofstream::out);
  for(int i=0;i<pars->active_names.size();i++){
    f_names << pars->active_names[i];
    f_names << " ";
    f_names << pars->active_names[i];
    f_names << std::endl;
  }
  f_names.close();

  // File with the parameter ranges
  std::ofstream f_ranges(output + this->name + "_postdist.ranges",std::ofstream::out);
  for(int i=0;i<pars->active.size();i++){
    f_ranges << pars->active_names[i];
    f_ranges << " ";
    f_ranges << pars->active[i]->min;
    f_ranges << " ";
    f_ranges << pars->active[i]->max;
    f_ranges << std::endl;
  }
  f_ranges.close();

  nested::run(IS,mmodal,ceff,nlive,tol,efr,ndims,nPar,nClsPar,maxModes,updInt,Ztol,root,seed,pWrap,fb,resume,outfile,initMPI,logZero,maxiter,MultiNestLogLike,MultiNestDumper,&myextras);
}



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
    //    printf(" %25.18e",e->pars->active[i]->pri->fromUnitCube(Cube[i]));
  }
  //  std::cout << std::endl;
      
  e->pars->updateActive(new_pars);
  e->pars->updateLikelihoodModel();
  lnew = e->pars->getLogLike();
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

  e->minimizer->counter++;
  e->pars->outputLikelihoodModel(e->output + std::to_string(e->minimizer->counter) + "_");

  /*
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
  */

  // postdist holds the posterior probability, loglike, and the parameters, i.e. it is the restructured internal 'posterior' array
  int Ncols = nPar + 2;
  double postdist[nSamples][Ncols];
  for(int j=0;j<nSamples;j++){
    postdist[j][0] = posterior[0][nSamples*(Ncols-2)+j];
    postdist[j][1] = posterior[0][nSamples*(Ncols-1)+j];
    for(int i=0;i<Ncols-2;i++){
      postdist[j][i+2] = e->pars->active[i]->pri->fromUnitCube( posterior[0][nSamples*i+j] );
    }
  }

  // File for the corner plot
  fh = fopen( (e->output + std::to_string(e->minimizer->counter) + "_" + e->minimizer->name + "_postdist.txt").c_str() ,"w");
  for(int j=0;j<nSamples;j++){
    fprintf(fh," %25.18e",1.0);
    fprintf(fh," %25.18e",postdist[j][1]);
    for(int i=2;i<Ncols;i++){
      fprintf(fh," %25.18e",postdist[j][i]);
    }
    fprintf(fh,"\n");
  }
  fclose(fh);

  // Calculating the mean and the lower and upper 1-sigma bounds
  std::vector<double> prob(nSamples);
  for(int j=0;j<nSamples;j++){
    prob[j] = postdist[j][1];
  }
  for(int i=0;i<Ncols-2;i++){
    std::vector<double> par(nSamples);
    for(int j=0;j<nSamples;j++){
      par[j] = postdist[j][i+2];
    }
    std::map<std::string,double> stats = Nlpar::getSigmaIntervals(par,prob,1);
    e->pars->maps[i]    = e->pars->active[i]->pri->fromUnitCube( paramConstr[0][i] );
    e->pars->means[i]   = stats["mean"];
    e->pars->s1_low[i]  = stats["low"];
    e->pars->s1_high[i] = stats["high"];
  }

  e->minimizer->output(e->output);
}

//virtual
void MultiNest::output(std::string output){
  // Read the mn-resume file and write the number of total samples and replacements
  std::string dum;
  std::ifstream infile(output + "mn-resume.dat");
  infile >> dum >> this->replacements >> this->total_samples;
  infile.close();

  Json::Value min_out;
  min_out["total_samples"] = this->total_samples;
  min_out["replacements"]  = this->replacements;

  std::ofstream jsonfile(output + std::to_string(this->counter) + "_" + this->name + "_minimizer_output.json");
  jsonfile << min_out;
  jsonfile.close();
}

//non-virtual
void MultiNest::finalizeMinimizer(std::string output){
  std::ifstream src((output + std::to_string(this->counter) + "_" + this->name + "_postdist.txt").c_str(),std::ios::binary);
  std::ofstream dst((output + this->name + "_postdist.txt").c_str(), std::ios::binary);
  dst << src.rdbuf();
  std::ifstream src2((output + std::to_string(this->counter) + "_" + this->name + "_minimizer_output.json").c_str(),std::ios::binary);
  std::ofstream dst2((output + this->name + "_minimizer_output.json").c_str(), std::ios::binary);
  dst2 << src2.rdbuf();
}
