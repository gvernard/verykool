#ifndef SUPPORT_FUNCS_HPP
#define SUPPORT_FUNCS_HPP

#include <string>
#include <vector>
#include <map>

#include "json/json.h"

#include "nonLinearPars.hpp"


class BaseLikelihoodModel;
class ImagePlane;
class BaseSourcePlane;
class CollectionMassModels;

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

  std::map<std::string,std::string> source;    //parameters of source plane and type (fixed,floating,adaptive)
  std::map<std::string,std::string> image;     //parameters of the image plane (x and y size and resolution)
  std::map<std::string,std::string> psf;       //parameters of the psf (x and y, total and cropped sizes)
  std::map<std::string,std::string> minimizer; //parameters for minimizer (MultiNest,Minuit2,etc)

  std::vector<std::string> mmodel;                           //mass models of the lenses


  Initialization(){};
  ~Initialization(){};
  


  static void initialize_program(std::string path,std::string run,Initialization*& init,BaseLikelihoodModel*& mypars,ImagePlane*& mydata,CollectionMassModels*& mycollection,BaseSourcePlane*& mysource);
  static void finalize_program(Initialization* init,BaseLikelihoodModel* mypars,ImagePlane* mydata,CollectionMassModels* mycollection,BaseSourcePlane* mysource);
  static void outputGeneric(BaseLikelihoodModel* mypars,ImagePlane* image,BaseSourcePlane* source,std::string output);

  void parseInputJSON(std::string path,std::string filename);
  
private:
  void outputInitial(BaseLikelihoodModel* mypars);
};

#endif /* SUPPORT_FUNCS_HPP */
