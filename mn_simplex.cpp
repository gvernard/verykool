#include "minimizers.hpp"

Simplex::Simplex(ImagePlane* a,BaseSourcePlane* b,CollectionMassModels* c,std::vector<std::map<std::string,BaseNlpar*> > d,mymatrices* e,precomp* f){
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

Simplex::~Simplex(){
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

//double Simplex::operator()(const std::vector<double> pars){
double Simplex::LogLike(const std::vector<double> pars){
  //Update values for nlpars
  for(int i=0;i<this->active.size();i++){
    this->nlpars[this->active[i].index][this->active[i].nam]->val = pars[i];
  }

  return mainLogLike(this->image,this->source,this->mycollection,this->nlpars,this->mat,this->pcomp);
}











void Simplex::printPars(const std::vector<double>& c){
  std::cout.precision(12);
  copy(c.begin(),c.end(),std::ostream_iterator<double>(std::cout," "));
  std::cout << std::endl;
}

std::vector<double> Simplex::initPars(std::vector<std::map<std::string,BaseNlpar*> > nlpars){
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


void Simplex::minimize(std::map<std::string,std::string> opt,ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* mycollection,std::vector<std::map<std::string,BaseNlpar*> > nlpars,mymatrices* mat,precomp* pcomp,const std::string output){
  // output and opt and used only in this function, the rest are passed on

  // Initial parameter values
  std::vector<double> init_pars = initPars(nlpars);

  printf("\n");
  printf("%-30s%.1e\n","Tolerance: ",stof(opt["tol"]));
  printf("%-30s%d\n","Max. iterations: ",stoi(opt["maxiter"]));
  printf("%-30s%d\n","Number of active parameters: ",init_pars.size());
  printf("%-30s\n","Initial parameters: ");
  printPars(init_pars);

  //  mySimplexClass papari(image,source,mycollection,nlpars,mat,pcomp);
  std::vector<double> out_pars = SimplexMinimizer(init_pars,stof(opt["tol"]),stoi(opt["maxiter"]));
  //  std::vector<double> out_pars = SimplexAnnealing(papari,init_pars,stof(opt["tol"]),stoi(opt["maxiter"]));

  std::cout << "Final parameters: " << std::endl;
  printPars(out_pars);
}



void Simplex::output(){
  // Simplex finished
}











/*
 From: https://github.com/VinGarcia/downhill-simplex/blob/master/demo.cpp

 Copyright (C) 2010 Botao Jia
 This file is an implementation of the downhill simplex optimization algorithm using C++.
 To use BT::Simplex correctly, the followings are needed, inclusively.
 1. f: a function object or a function which takes a vector<class Type> and returns a Type, inclusively.
       Signature, e.g. for double Type, 
            double f(vector<double> x);
 2. init: an inital guess of the fitted parameter values which minmizes the value of f. 
          init must be a vector<double>
          init must have the exactly same dimension as the vector taken by f.
          init must order the parameters, such that the order follows the vector taken by f.
          e.g. f takes vector<double> x, where x[0] represents parameter1; x[1] represents parameter2, etc.
          And init must follow this order exactly, init[0] is the initial guess for parameter1,
          init[1] is the initial guess for parameter2, etc.
3 to 5 are all optional:
3. tol: the termination criteria. 
        It measures the difference of the simplex centroid*(N+1) of consecutive iterations.
4. x:  an initial simplex, which is calculated according to the initial trial parameter values.
5. iterations: maximum iterations.
The return value of BT::Simplex is a vector<Type>, 
which represents the optimized parameter values that minimize f.
The order of the parameter is the same as in the vector<Type> init.
 You can redistribute it and/or modify it at will.
 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  
*/

std::vector<double> Simplex::SimplexMinimizer(std::vector<double> init,double tol,int iterations){

  int N = init.size();                       // space dimension
  const double a=1.0, b=1.0, g=0.5, h=0.5;   // coefficients: a->reflection->xr, b->expansion->xe, g->contraction->xc, h->full contraction to x1
  std::vector<double> xcentroid_old(N,0);    // simplex center * (N+1)
  std::vector<double> xcentroid_new(N,0);    // simplex center * (N+1)
  std::vector<double> vf(N+1,0);             // f evaluated at simplex vertexes       
  int x1=0, xn=0, xnp1=0;                    // x1: f(x1) = min { f(x1), f(x2)...f(x_{n+1} },  xnp1: f(xnp1) = max { f(x1), f(x2)...f(x_{n+1} },  xn: f(xn)<f(xnp1) && f(xn)> all other f(x_i)


  // If no initial simplex is specified construct the trial simplex based upon the initial guess parameters
  std::vector<std::vector<double> > x;
  if( x.size() == 0 ){
    std::vector<double> del( init );
    std::transform(del.begin(),del.end(),del.begin(),std::bind2nd(std::divides<double>(),20));// '20' is picked assuming initial trail close to true
    
    for(int i=0;i<N;i++){
      std::vector<double> tmp(init);
      tmp[i] += del[i];
      x.push_back(tmp);
    }
    x.push_back(init);// x.size()=N+1, x[i].size()=N
    
    // xcentroid
    std::transform(init.begin(),init.end(),xcentroid_old.begin(),std::bind2nd(std::multiplies<double>(),N+1));
  }
  // constructing the simplex finished
  



  // Begin optimization
  int cnt = 0; // counter of iterations
  for(cnt=0;cnt<iterations;cnt++){
    
    for(int i=0;i<N+1;i++){
      vf[i] = this->LogLike(x[i]);
    }
    
    // find index of max, second max, min of vf.
    x1   = 0;
    xn   = 0;
    xnp1 = 0;  
    for(unsigned i=0;i<vf.size();i++){
      if( vf[i] < vf[x1] ){
	x1 = i;
      }
      if( vf[i] > vf[xnp1] ){
	xnp1 = i;
      }
    }

    xn = x1;
    for(unsigned i=0;i<vf.size();i++){ 
      if( vf[i] < vf[xnp1] && vf[i] > vf[xn] ){
	xn = i;
      }
    }
    // x1, xn, xnp1 are found
    

    // xg: centroid of the N best vertexes
    std::vector<double> xg(N,0);
    for(unsigned i=0;i<x.size();i++){
      if( (int)i != xnp1 ){
	std::transform(xg.begin(),xg.end(),x[i].begin(),xg.begin(),std::plus<double>());
      }
    }
    std::transform(xg.begin(),xg.end(),x[xnp1].begin(),xcentroid_new.begin(),std::plus<double>());
    std::transform(xg.begin(),xg.end(),xg.begin(),std::bind2nd(std::divides<double>(),N));
    // xg found, xcentroid_new updated
    


    //termination condition
    double diff=0;  // calculate the difference of the simplex centers and see if the difference is less than the termination criteria
    for(int i=0;i<N;i++){
      diff += fabs(xcentroid_old[i]-xcentroid_new[i]);
    }
    
    if( diff/N < tol ){
      // terminate the optimizer
      std::cout << std::endl << "Number of iterations: " << cnt << std::endl;
      break;
    } else {
      // update simplex center
      xcentroid_old.swap(xcentroid_new);
    }

    
    // reflection:
    std::vector<double> xr(N,0); 
    for(int i=0;i<N;i++){
      xr[i] = xg[i]+a*(xg[i]-x[xnp1][i]);
    }
    // reflection, xr found

    
    double fxr = this->LogLike(xr); // record function at xr
    
    if( vf[x1] <= fxr && fxr <= vf[xn] ){
      
      std::copy(xr.begin(),xr.end(),x[xnp1].begin());
      
    } else if( fxr < vf[x1] ){
      
      // expansion:
      std::vector<double> xe(N,0);
      for(int i=0;i<N;i++){
	xe[i] = xr[i]+b*(xr[i]-xg[i]);
      }
      if( this->LogLike(xe) < fxr ){
	std::copy(xe.begin(),xe.end(),x[xnp1].begin());
      } else {
	std::copy(xr.begin(),xr.end(),x[xnp1].begin());
      }
      // expansion finished,  xe is not used outside the scope
      
    } else if( fxr > vf[xn] ){
      
      // contraction:
      std::vector<double> xc(N,0);
      for(int i=0;i<N;i++){
	xc[i] = xg[i] + g*(x[xnp1][i]-xg[i]);
      }
      if( this->LogLike(xc) < vf[xnp1] ){
	std::copy(xc.begin(),xc.end(),x[xnp1].begin());
      } else {
	for(unsigned i=0;i<x.size();i++){
	  if( (int)i != x1 ){
	    for(int j=0;j<N;j++){
	      x[i][j] = x[x1][j] + h * ( x[i][j]-x[x1][j] );
	    }
	  }
	}
      }
      // contraction finished, xc is not used outside the scope
      
    }
    
    //    printPars(x[x1]);
  }
  // optimization is finished
  

  if( cnt == iterations ){
    // max number of iteration is reached before tol is satisfied
    std::cout << "Iteration limit reached, result may not be optimal" << std::endl;
  }
  return x[x1];
  
}



double Simplex::myrandom(){
  double r = ((double) rand() / (RAND_MAX));
  return r;
}


/*
template<class OP> std::vector<double> Simplex::SimplexAnnealing(OP f,std::vector<double> init,double tol,int iterations,std::vector<std::vector<double> > x = std::vector<std::vector<double> >()){

  int T_iter  = 100;
  int iter    = 20;
  int tot_iter = T_iter*iter;

  double T0  = 10;
  int alpha  = 1.0;



  int N = init.size();                       // space dimension
  const double a=1.0, b=1.0, g=0.5, h=0.5;   // coefficients: a->reflection->xr, b->expansion->xe, g->contraction->xc, h->full contraction to x1
  std::vector<double> xcentroid_old(N,0);    // simplex center * (N+1)
  std::vector<double> xcentroid_new(N,0);    // simplex center * (N+1)
  std::vector<double> vf(N+1,0);             // f evaluated at simplex vertexes       
  int x1=0, xn=0, xnp1=0;                    // x1: f(x1) = min { f(x1), f(x2)...f(x_{n+1} },  xnp1: f(xnp1) = max { f(x1), f(x2)...f(x_{n+1} },  xn: f(xn)<f(xnp1) && f(xn)> all other f(x_i)


  // If no initial simplex is specified construct the trial simplex based upon the initial guess parameters
  if( x.size() == 0 ){
    std::vector<double> del( init );
    std::transform(del.begin(),del.end(),del.begin(),std::bind2nd(std::divides<double>(),20));// '20' is picked assuming initial trail close to true
    
    for(int i=0;i<N;i++){
      std::vector<double> tmp(init);
      tmp[i] += del[i];
      x.push_back(tmp);
    }
    x.push_back(init);// x.size()=N+1, x[i].size()=N
    
    // xcentroid
    std::transform(init.begin(),init.end(),xcentroid_old.begin(),std::bind2nd(std::multiplies<double>(),N+1));
  }
  // constructing the simplex finished
  



  // Begin optimization
  int cnt = 0; // counter of iterations
  for(cnt=0;cnt<T_iter;cnt++){


    double T = T0*pow(1.0 - (double)(cnt*iter)/double(tot_iter),alpha);
    std::cout << std::endl << "Temperature: " << T << std::endl << std::endl;

    for(int k=0;k<iter;k++){
      for(int i=0;i<N+1;i++){
	vf[i] = f(x[i]) + T*log(myrandom());
      }
      
      // find index of max, second max, min of vf.
      x1   = 0;
      xn   = 0;
      xnp1 = 0;  
      for(unsigned i=0;i<vf.size();i++){
	if( vf[i] < vf[x1] ){
	  x1 = i;
	}
	if( vf[i] > vf[xnp1] ){
	  xnp1 = i;
	}
      }
      
      xn = x1;
      for(unsigned i=0;i<vf.size();i++){ 
	if( vf[i] < vf[xnp1] && vf[i] > vf[xn] ){
	  xn = i;
	}
      }
      // x1, xn, xnp1 are found
      
      
      // xg: centroid of the N best vertexes
      std::vector<double> xg(N,0);
      for(unsigned i=0;i<x.size();i++){
	if( (int)i != xnp1 ){
	  std::transform(xg.begin(),xg.end(),x[i].begin(),xg.begin(),std::plus<double>());
	}
      }
      std::transform(xg.begin(),xg.end(),x[xnp1].begin(),xcentroid_new.begin(),std::plus<double>());
      std::transform(xg.begin(),xg.end(),xg.begin(),std::bind2nd(std::divides<double>(),N));
      // xg found, xcentroid_new updated
      
      
      
      //termination condition
      double diff=0;  // calculate the difference of the simplex centers and see if the difference is less than the termination criteria
      for(int i=0;i<N;i++){
	diff += fabs(xcentroid_old[i]-xcentroid_new[i]);
      }
      
      if( diff/N < tol ){
	// terminate the optimizer
	std::cout << std::endl << "Number of iterations: " << cnt << std::endl;
	break;
      } else {
	// update simplex center
	xcentroid_old.swap(xcentroid_new);
      }
      
      
      // reflection:
      std::vector<double> xr(N,0); 
      for(int i=0;i<N;i++){
	xr[i] = xg[i]+a*(xg[i]-x[xnp1][i]);
      }
      // reflection, xr found
      
      
      double fxr = f(xr) - T*log(myrandom()); // record function at xr
      
      if( vf[x1] <= fxr && fxr <= vf[xn] ){
	
	std::copy(xr.begin(),xr.end(),x[xnp1].begin());
	
      } else if( fxr < vf[x1] ){
	
	// expansion:
	std::vector<double> xe(N,0);
	for(int i=0;i<N;i++){
	  xe[i] = xr[i]+b*(xr[i]-xg[i]);
	}
	if( (f(xe) - T*log(myrandom())) < fxr ){
	  std::copy(xe.begin(),xe.end(),x[xnp1].begin());
	} else {
	  std::copy(xr.begin(),xr.end(),x[xnp1].begin());
	}
	// expansion finished,  xe is not used outside the scope
	
      } else if( fxr > vf[xn] ){
	
	// contraction:
	std::vector<double> xc(N,0);
	for(int i=0;i<N;i++){
	  xc[i] = xg[i] + g*(x[xnp1][i]-xg[i]);
	}
	if( (f(xc) - T*log(myrandom())) < vf[xnp1] ){
	  std::copy(xc.begin(),xc.end(),x[xnp1].begin());
	} else {
	  for(unsigned i=0;i<x.size();i++){
	    if( (int)i != x1 ){
	      for(int j=0;j<N;j++){
		x[i][j] = x[x1][j] + h * ( x[i][j]-x[x1][j] );
	      }
	    }
	  }
	}
	// contraction finished, xc is not used outside the scope
	
      }
      
      //    printPars(x[x1]);
    }
    



  }
  // optimization is finished
  

  if( cnt == iterations ){
    // max number of iteration is reached before tol is satisfied
    std::cout << "Iteration limit reached, result may not be optimal" << std::endl;
  }
  return x[x1];
  
}
*/
