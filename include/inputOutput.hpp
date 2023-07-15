#ifndef SUPPORT_FUNCS_HPP
#define SUPPORT_FUNCS_HPP

#include <string>
#include <vector>
#include <map>

#include "json/json.h"

#include "sourceProfile.hpp"

class Nlpar;
class BaseLikelihoodModel;
class ImagePlane;
class BaseSourcePlane;
class CollectionMassModels;
class Pert;
class BaseMinimizer;

class Initialization {
public:
  std::string imgpath;        //path to image data
  std::string psfpath;        //path to PSF data
  std::string maskpath;       //path to mask
  std::string noisepath;      //path to covariance matrix of the data
  std::string output;         //path to output directory

  std::string noise_flag;     //flag describing the kind of noise in the cov.dat file (uniform, noise map, correlated noise)
  std::string reg;            //the source regularization scheme
  std::string interp;         //interpolation scheme (bilinear,bicubic)
  std::string likeModel;      //likelihood model

  std::map<std::string,std::string> source;           //parameters of source plane and type (fixed,floating,adaptive)
  std::map<std::string,std::string> image;            //parameters of the image plane (x and y size and resolution)
  std::map<std::string,std::string> psf;              //parameters of the psf (x and y, total and cropped sizes)
  std::map<std::string,std::string> perturbations;    //parameters for the perturbations
  std::map<std::string,std::string> minimizer;        //parameters for the model minimizer (MultiNest,Minuit2,etc)

  std::vector<std::string> mmodel;                    //mass models of the lenses
  std::vector<std::string> lens_names;                //the name of each lens
  std::vector<std::string> src_names;                 //the name of each source
  std::vector< std::vector<Nlpar*> > nlpars_lenses;   //non-linear parameters for each lens model
  std::vector<Nlpar*> nlpars_physical;                //non-linear parameters physical
  std::vector<Nlpar*> nlpars_reg_s;                   //non-linear parameters for source regularization
  std::vector<Nlpar*> nlpars_reg_dpsi;                //non-linear parameters for dpsi regularization

  BaseProfile* prof_source0 = 0;

  Initialization(){};
  ~Initialization(){
    delete(prof_source0);
  };
  
  static void initialize_program(std::string path,std::string run,Initialization*& init,ImagePlane*& mydata,CollectionMassModels*& mycollection,BaseSourcePlane*& mysource,BaseSourcePlane*& source0,Pert*& pert_mass_model,BaseMinimizer*& minimizer);
  static void outputGeneric(BaseLikelihoodModel* smooth_like,ImagePlane* image,BaseSourcePlane* source,std::string output);

private:
  void parseInputJSON(std::string path,std::string filename);
};

#endif /* SUPPORT_FUNCS_HPP */
