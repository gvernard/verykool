#include "myminuit.hpp"

//Functions related to likelihood minimization using Minuit2
//==============================================================================================================================
double myMinuitFCN::operator()(const std::vector<double>& pars) const {

  // set struct of nlpars according to pars vector which is updated by Minuit (some pars are already fixed)
  typedef std::map<std::string,BaseNlpar*>::iterator it_type;
  // need to use a dum map because this is not working for some reason: this->nlpars[0]["g"]->val, but this is: dum["g"]->val
  std::map<std::string,BaseNlpar*> phys,other,lens;

  /*
  printf("\n");
  for(int i=0;i<pars.size();i++){
    printf(" %16.7f",pars[i]);
  }
  */


  int pcount = 0;
  // set physical parameters
  phys = this->nlpars[0];
  for(it_type it=phys.begin();it!=phys.end();it++){
    phys[it->second->nam]->val = pars[pcount];
    pcount++;
  }
  this->collection->setPhysicalPars(this->nlpars[0]);

  // set other parameters
  other = this->nlpars[1];
  for(it_type it=other.begin();it!=other.end();it++){
    other[it->second->nam]->val = pars[pcount];
    pcount++;
  }

  // set parameters of all the lenses
  for(int k=0;k<this->collection->models.size();k++){
    lens = this->nlpars[2+k];
    for(it_type it=lens.begin();it!=lens.end();it++){
      lens[it->second->nam]->val = pars[pcount];
      pcount++;
    }
    this->collection->models[k]->setMassPars(this->nlpars[2+k]);
  }



  std::cout << " 1 " << std::endl;
  if( this->source->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(this->source);
    ada->createAdaGrid(this->image,this->collection);
    ada->createDelaunay();
    this->source->constructH(&this->matrices->H);
  }
  
  std::cout << " 2 " << std::endl;
  this->source->constructL(this->image,this->collection,&this->matrices->L);
  std::cout << " 3 " << std::endl;
  setAlgebraRuntime(this->image,this->source,other["lambda"]->val,this->matrices,this->pcomp);
  std::cout << " 4 " << std::endl;
  solveSource(this->image,this->source,this->pcomp);
  std::cout << " 5 " << std::endl;
  double logL = myLogLike(this->image,this->source,other["lambda"]->val,this->pcomp,this->nlpars);
  


  printf("%16.7f\n\n",-logL);



  // we maximize the evidence, so need a minus sign for the minuit minimization routines
  return -logL;
}

void myMinuit(ImagePlane* image,BaseSourcePlane* source,CollectionMassModels* collection,std::vector<std::map<std::string,BaseNlpar*> > nlpars,mymatrices* mat,precomp* pcomp,const std::string output,std::map<std::string,std::string> opt){

  ROOT::Minuit2::MnUserParameters upar;

  typedef std::map<std::string,BaseNlpar*>::iterator it_type;
  // set physical parameters
  for(it_type it=nlpars[0].begin();it!=nlpars[0].end();it++){
    if( it->second->fix == 1 ){
      upar.Add(it->second->nam,it->second->val);
      upar.Fix(it->second->nam);
    } else {
      //    upar.Add(it->second->nam,it->second->val,it->second->err,it->second->min,it->second->max);
      upar.Add(it->second->nam,it->second->min+(it->second->max-it->second->min)/3.,10*(it->second->max-it->second->min),it->second->min,it->second->max);
    }
  }
  // set other parameters
  for(it_type it=nlpars[1].begin();it!=nlpars[1].end();it++){
    if( it->second->fix == 1 ){
      upar.Add(it->second->nam,it->second->val);
      upar.Fix(it->second->nam);
    } else {
      //    upar.Add(it->second->nam,it->second->val,it->second->err,it->second->min,it->second->max);
      upar.Add(it->second->nam,it->second->min+(it->second->max-it->second->min)/3.,10*(it->second->max-it->second->min),it->second->min,it->second->max);
    }
  }
  // set parameters of all the lenses
  for(int k=0;k<collection->models.size();k++){
    for(it_type it=nlpars[2+k].begin();it!=nlpars[2+k].end();it++){
      if( it->second->fix == 1 ){
	upar.Add(it->second->nam,it->second->val);
	upar.Fix(it->second->nam);
      } else {
	//      upar.Add(it->second->nam,it->second->val,it->second->err,it->second->min,it->second->max);
	upar.Add(it->second->nam,it->second->min+(it->second->max-it->second->min)/3.,10*(it->second->max-it->second->min),it->second->min,it->second->max);
      }
    }
  }


  // get all parameters into one vector
  int ntot = 0;
  for(int i=0;i<nlpars.size();i++){
    ntot += nlpars[i].size();
  }
  std::vector<double> pars;
  for(int i=0;i<ntot;i++){
    pars.push_back(upar.Value(i));
  }


  // initialize user function and related structs
  myMinuitFCN myfcn(pars);
  myfcn.image      = image;
  myfcn.source     = source;
  myfcn.collection = collection;
  myfcn.nlpars     = nlpars;
  myfcn.matrices   = mat;
  myfcn.pcomp      = pcomp;


  myfcn(pars);


  
  ROOT::Minuit2::MnMigrad migrad(myfcn,upar);
  ROOT::Minuit2::FunctionMinimum min = migrad();
  //ROOT::Minuit2::MnSimplex minimizer(myfcn,upar);
  //ROOT::Minuit2::FunctionMinimum min = minimizer();

  std::cout << min << std::endl;

  char status;
  if( min.IsValid() ){
    status = 's';
  } else {
    status = 'f';
  }

  for(int i=0;i<pars.size();i++){
    std::cout << min.UserState().Value(i) << std::endl;
  }
}

