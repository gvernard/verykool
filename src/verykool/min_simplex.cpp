#include "minimizers.hpp"

#include <iostream>
#include <iterator>

#include "nonLinearPars.hpp"
#include "mainLogLike.hpp"

// mySimplex derived class
//================================================================================================================================================
mySimplex::mySimplex(ImagePlane* a,BaseSourcePlane* b,CollectionMassModels* c,std::vector<std::map<std::string,BaseNlpar*> > d,mymatrices* e,precomp* f){
  image        = a;
  source       = b;
  mycollection = c;
  nlpars       = d;
  mat          = e;
  pcomp        = f;

  //set active parameters (in public variable 'active') in order: physical, other, lens parameters for each lens (from nlpars; order matters!!!)
  myactive tmp;
  std::vector<std::string> names;

  for(int j=0;j<nlpars.size();j++){
    names = BaseNlpar::getActive(nlpars[j]);
    for(int i=0;i<names.size();i++){
      tmp = {j,names[i],nlpars[j][names[i]]->per};
      active.push_back(tmp);
    }
    names.clear();
  }

}

mySimplex::~mySimplex(){
  /*
  typedef std::map<std::string,BaseNlpar*>::iterator it_type;
  for(int i=0;i<this->nlpars.size();i++){
    for(it_type iterator=this->nlpars[i].begin();iterator!=this->nlpars[i].end();iterator++){
      delete this->nlpars[i][iterator->first];
    }
    this->nlpars[i].clear();
  }
  this->nlpars.clear(); 
  */

  this->active.clear();
}

void mySimplex::printPars(const std::vector<double>& c){
  std::cout.precision(12);
  copy(c.begin(),c.end(),std::ostream_iterator<double>(std::cout," "));
  std::cout << std::endl;
}

std::vector<double> mySimplex::initPars(std::vector<std::map<std::string,BaseNlpar*> > nlpars){
  std::vector<double> init_pars;
  std::vector<std::string> names;
 
  //set active parameters (and periodicity) in order: physical, other, lens parameters for each lens
  for(int j=0;j<nlpars.size();j++){
    names = BaseNlpar::getActive(nlpars[j]);
    for(int i=0;i<names.size();i++){
      //      std::cout << nlpars[j][names[i]]->val << std::endl;
      init_pars.push_back(nlpars[j][names[i]]->val*0.95+0.001);
    }
    names.clear();
  }

  return init_pars;
}

double mySimplex::LogLike(double* pars){
  //Update values for nlpars
  for(int i=0;i<this->active.size();i++){
    this->nlpars[this->active[i].index][this->active[i].nam]->val = pars[i];
  }

  return mainLogLike(this->image,this->source,this->mycollection,this->nlpars,this->mat,this->pcomp);
}




void mySimplex::minimize(std::map<std::string,std::string> opt,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,std::vector<std::map<std::string,BaseNlpar*> > nlpars,mymatrices* mat,precomp* pcomp,const std::string output){
  // output and opt and used only in this function, the rest are passed on

  // Initial parameter values
  std::vector<double> init_pars = initPars(nlpars);
  int ndim = init_pars.size();

  double* pars0     = (double*) malloc(ndim*sizeof(double)); // Initial point for minimization without simulated annealing
  double* stepsizes = (double*) malloc(ndim*sizeof(double)); // Initial displacements to define starting vertices of simplex
  for(int i=0;i<ndim;i++){
    pars0[i] = init_pars[i];
    stepsizes[i] = 0.1 * pars0[i];
  }
  double initial_temp = stof(opt["initial_temp"]); // initial temperature for simulated annealing. Set to zero for a regular simplex.minimizing without simulated annealing.
  double precision    = stof(opt["precision"]);  // Required precision for convergence to local minimum
  int random_seed     = stoi(opt["seed"]); // used for random number generator
  double f;

  this->initialize_simplex(pars0,ndim,stepsizes,precision,random_seed);
  this->funcptr = static_cast<double (Simplex::*)(double*)> (&mySimplex::LogLike);
  this->simplex_set_function(funcptr);

  if( initial_temp == 0 ){
    printf("%-30s\n","Downhill simplex");
    printf("\n");
    printf("%-30s%d\n","Number of active parameters: ",init_pars.size());
    printf("%-30s\n","Initial parameters: ");
    this->printPars(init_pars);

    int iterations; // this will keep track of how many iterations were required to converge
    this->downhill_simplex(iterations,stoi(opt["maxiter"]),0); // runs downhill simplex with a maximum allowed # of iterations = 1000, and temperature of zero (i.e. no annealing)
  } else {
    printf("%-30s\n","Downhill simplex with simulated annealing");
    printf("\n");
    printf("%-30s%d\n","Number of active parameters: ",init_pars.size());
    printf("%-30s\n","Initial parameters: ");
    this->printPars(init_pars);

    this->simplex_set_fmin(stof(opt["min_function_val"])); // optional feature
    this->set_annealing_schedule_parameters(initial_temp,stof(opt["final_temp"]),stof(opt["cooling_factor"]),stoi(opt["iter_final_temp"]),stoi(opt["iter_per_temp"]));
    this->downhill_simplex_anneal(true);
  }
  this->simplex_minval(pars0,f); // store best-fit point and function value


  std::vector<double> out_pars(ndim);
  for(int i=0;i<ndim;i++){
    out_pars[i] = pars0[i];
  }

  std::cout << "Final parameters: " << std::endl;
  printPars(out_pars);
}



void mySimplex::output(){
  // Simplex finished
}

