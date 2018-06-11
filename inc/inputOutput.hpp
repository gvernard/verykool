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
class Pert;

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
  std::string smooth_like;    //likelihood model for the smooth mass model
  std::string pert_like;      //likelihood model for the perturbations

  std::map<std::string,std::string> source;           //parameters of source plane and type (fixed,floating,adaptive)
  std::map<std::string,std::string> image;            //parameters of the image plane (x and y size and resolution)
  std::map<std::string,std::string> psf;              //parameters of the psf (x and y, total and cropped sizes)
  std::map<std::string,std::string> perturbations;    //parameters for the perturbations
  std::map<std::string,std::string> smooth_minimizer; //parameters for the smooth model minimizer (MultiNest,Minuit2,etc)
  std::map<std::string,std::string> pert_minimizer;   //parameters for the perturbations minimizer (MultiNest,iterator,etc)

  std::vector<std::string> mmodel;                    //mass models of the lenses


  Initialization(){};
  ~Initialization(){};
  
  static void initialize_program(std::string path,std::string run,Initialization*& init,BaseLikelihoodModel*& smooth_like,ImagePlane*& mydata,CollectionMassModels*& mycollection,BaseSourcePlane*& mysource,BaseLikelihoodModel*& pert_like,Pert*& pert_mass_model);
  static void finalize_smooth(Initialization* init,BaseLikelihoodModel* smooth_like,ImagePlane* mydata,CollectionMassModels* mycollection,BaseSourcePlane* mysource);
  static void finalize_pert();
  static void outputGeneric(BaseLikelihoodModel* smooth_like,ImagePlane* image,BaseSourcePlane* source,std::string output);

private:
  void parseInputJSON(std::string path,std::string filename);
  void outputInitial(BaseLikelihoodModel* smooth_like);
};

#endif /* SUPPORT_FUNCS_HPP */
