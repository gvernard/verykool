#ifndef SUPPORT_FUNCS_HPP
#define SUPPORT_FUNCS_HPP

#include <string>
#include <vector>
#include <map>

#include "json/json.h"

#include "nonLinearPars.hpp"


class ImagePlane;
class BaseSourcePlane;
class CollectionMassModels;
class mymatrices;
class precomp;

class Initialization {
public:
  std::string imgpath;        //path to image data
  std::string psfpath;        //path to PSF data
  std::string maskpath;       //path to mask
  std::string covpath;        //path to covariance matrix of the data
  std::string output;         //path to output directory

  std::string noise_flag;     //flag describing the kind of noise in the cov.dat file (uniform, noise map, correlated noise)
  std::string reg;            //the source regularization scheme
  std::string interp;         //interpolation scheme (bilinear,bicubic)
  int postmcmc;               //continue with MCMC sampling after finding a minimum/MAP value

  std::map<std::string,std::string> source;    //parameters of source plane and type (fixed,floating,adaptive)
  std::map<std::string,std::string> image;     //parameters of the image plane (x and y size and resolution)
  std::map<std::string,std::string> psf;       //parameters of the psf (x and y, total and cropped sizes)
  std::map<std::string,std::string> minimizer; //parameters for minimizer (MultiNest,Minuit2,etc)

  std::vector<std::string> mmodel;                           //mass models of the lenses
  std::vector<std::map<std::string,BaseNlpar*> > nlpars;     //non-linear parameters of the physical system [0], other [1], and of the mass models [2-*]


  Initialization(){};
  ~Initialization(){
    typedef std::map<std::string,BaseNlpar*>::iterator it_type;
    for(int i=0;i<this->nlpars.size();i++){
      for(it_type iterator=this->nlpars[i].begin();iterator!=this->nlpars[i].end();iterator++){
	delete this->nlpars[i][iterator->first];
      }
      this->nlpars[i].clear();
    }
    this->nlpars.clear();
  };

  void parseInputJSON(const char* path,const char* filename);


private:
  std::map<std::string,BaseNlpar*> nlparsFromJson(const Json::Value myjson);

};



void initialize_program(std::string path,std::string run,Initialization*& init,ImagePlane*& mydata,CollectionMassModels*& mycollection,BaseSourcePlane*& mysource,std::vector<myactive>& active,mymatrices*& matrices,precomp*& pcomp);
void finalize_program(Initialization* init,ImagePlane* mydata,CollectionMassModels* mycollection,BaseSourcePlane* mysource,mymatrices* matrices,precomp* pcomp);

void outputInitial(std::vector<std::map<std::string,BaseNlpar*> > nlpars,std::string output);
void outputGeneric(ImagePlane* image,BaseSourcePlane* source,std::vector<std::map<std::string,BaseNlpar*> > nlpars,precomp* pcomp,std::string output);

#endif /* SUPPORT_FUNCS_HPP */
