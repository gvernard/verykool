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
void MultiNest::minimize(std::map<std::string,std::string> opt,BaseLikelihoodModel* lmodel,const std::string output){
  int myndims = lmodel->active.size();
  
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
    pWrap[i] = lmodel->active[i]->per;
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
  extras myextras = {lmodel,output,this};


  // Initial output from the minimizer:
  // File with the parameter names
  std::ofstream f_names(output + lmodel->name + "_postdist.paramnames",std::ofstream::out);
  for(int i=0;i<lmodel->active_names.size();i++){
    f_names << lmodel->active_names[i];
    f_names << " ";
    f_names << lmodel->active_names[i];
    f_names << std::endl;
  }
  f_names.close();

  // File with the parameter ranges
  std::ofstream f_ranges(output + lmodel->name + "_postdist.ranges",std::ofstream::out);
  for(int i=0;i<lmodel->active.size();i++){
    f_ranges << lmodel->active_names[i];
    f_ranges << " ";
    f_ranges << lmodel->active[i]->min;
    f_ranges << " ";
    f_ranges << lmodel->active[i]->max;
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
  for(int i=0;i<e->lmodel->active.size();i++){
    new_pars.push_back( e->lmodel->active[i]->pri->fromUnitCube(Cube[i]) );
    //    printf(" %25.18e",e->pars->active[i]->pri->fromUnitCube(Cube[i]));
  }
  //  std::cout << std::endl;
      
  e->lmodel->updateActive(new_pars);
  e->lmodel->updateLikelihoodModel();
  lnew = e->lmodel->getLogLike();
  e->lmodel->printActive();
  e->lmodel->printTerms();
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

  
  // Minimizer specific output
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
    postdist[j][0] = posterior[0][nSamples*(Ncols-2)+j]; // loglike
    postdist[j][1] = posterior[0][nSamples*(Ncols-1)+j]; // posterior probability
    for(int i=0;i<Ncols-2;i++){
      postdist[j][i+2] = e->lmodel->active[i]->pri->fromUnitCube( posterior[0][nSamples*i+j] );
    }
  }

  // File for the corner plot
  fh = fopen( (e->output + std::to_string(e->minimizer->counter) + "_" + e->lmodel->name + "_postdist.txt").c_str() ,"w");
  for(int j=0;j<nSamples;j++){
    fprintf(fh," %25.18e",postdist[j][0]); // loglike
    fprintf(fh," %25.18e",postdist[j][1]); // posterior
    for(int i=2;i<Ncols;i++){
      fprintf(fh," %25.18e",postdist[j][i]);
    }
    fprintf(fh,"\n");
  }
  fclose(fh);


  
  Json::Value min_out;
  min_out["minimizer"] = e->minimizer->type;
  
  // Read the mn-resume file and write the number of total samples and replacements
  std::string dum;
  int total_samples;
  int replacements;
  std::ifstream infile(e->output + "mn-resume.dat");
  infile >> dum >> replacements >> total_samples;
  infile.close();
  min_out["total_samples"] = total_samples;
  min_out["replacements"]  = replacements;


  // Write the mean, sdev, maxlike, and MAP values for each parameter
  // Calculating the mean and the lower and upper 1-sigma bounds
  Json::Value parameters,par_json;
  double mean,sdev,maxlike,map;
  std::string name;
  std::vector<double> map_pars;
  printf("%10s %15s %15s %15s %15s %10s\n","","Mean (dubious)","Sdev (dubious)","Maxlike","MAP","");
  for(int i=0;i<nPar;i++){
    name    = e->lmodel->active_names[i]; // Need to have the right prefix, e.g. in the case of multiple lenses
    mean    = e->lmodel->active[i]->pri->fromUnitCube( paramConstr[0][i] );
    sdev    = e->lmodel->active[i]->pri->fromUnitCube( paramConstr[0][nPar+i] );
    maxlike = e->lmodel->active[i]->pri->fromUnitCube( paramConstr[0][2*nPar+i] );
    map     = e->lmodel->active[i]->pri->fromUnitCube( paramConstr[0][3*nPar+i] );
    par_json["mean"] = mean;
    par_json["sdev"] = sdev;
    par_json["maxlike"] = maxlike;
    par_json["map"] = map;
    parameters[name] = par_json;
    map_pars.push_back( map );
    printf("%10s %15.10f %15.10f %15.10f %15.10f %-10s\n",name.c_str(),mean,sdev,maxlike,map,name.c_str());
  }
  min_out["parameters"] = parameters;
  
  std::ofstream jsonfile(e->output + std::to_string(e->minimizer->counter) + "_" + e->lmodel->name + "_minimizer_output.json");
  jsonfile << min_out;
  jsonfile.close();



  //Update values for active parameters
  e->lmodel->updateActive(map_pars);
  e->lmodel->updateLikelihoodModel();
  e->lmodel->outputLikelihoodModel(e->output + std::to_string(e->minimizer->counter) + "_");  
}

//non-virtual
void MultiNest::finalizeMinimizer(std::string output,BaseLikelihoodModel* lmodel){
  printf("%-16s%20s\n","Starting output ",this->name.c_str());
  fflush(stdout);

  
  // Copy files to the final ones (without the step prefix)
  std::ifstream src((output + std::to_string(this->counter) +  "_" + lmodel->name + "_postdist.txt").c_str(),std::ios::binary);
  std::ofstream dst((output + lmodel->name + "_postdist.txt").c_str(), std::ios::binary);
  dst << src.rdbuf();
  std::ifstream src2((output + std::to_string(this->counter) + "_" + lmodel->name + "_minimizer_output.json").c_str(),std::ios::binary);
  std::ofstream dst2((output + lmodel->name + "_minimizer_output.json").c_str(), std::ios::binary);
  dst2 << src2.rdbuf();

  
  // Update likelihood model with MAP parameters
  Json::Value out_json;
  std::ifstream fin((output + lmodel->name + "_minimizer_output.json").c_str());
  fin >> out_json;

  printf("%s\n","   Printing results for the Max. Likelihood solution:");
  std::unordered_map<std::string,double> maps;
  Json::Value::Members jmembers;
  std::string key,key1,key2;
  jmembers = out_json["parameters"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    Json::Value tmp = jmembers[i]; // This is needed to cast the key into a json value before converting it to a string.
    std::string str = tmp.asString();
    double value = out_json["parameters"][jmembers[i]]["maxlike"].asDouble();

    std::stringstream ss(str);
    if( str.length() > 2 ){
      std::string last_two = str.substr( str.length()-2,str.length() );

      if( last_two.compare("_s") == 0 ){
	key = str;
      } else {
	std::getline(ss,key,'_');
	std::getline(ss,key,'_');
      }
    } else {
      key = str;
    }

    std::pair<std::string,double> mypair (key,value);
    maps.insert( mypair );
  }

  
  lmodel->updateActive(maps);
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
