#include "inputOutput.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "json/json.h"

#include "nonLinearPars.hpp"
#include "likelihoodModels.hpp"
#include "imagePlane.hpp"
#include "massModels.hpp"
#include "sourcePlane.hpp"
#include "sourceProfile.hpp"
#include "minimizers.hpp"


void Initialization::initialize_program(std::string path,std::string run,Initialization*& init,ImagePlane*& mydata,CollectionMassModels*& mycollection,BaseSourcePlane*& mysource,BaseSourcePlane*& source0,Pert*& pert_mass_model,BaseMinimizer*& minimizer){
  // Instantiate init object ---------------------------------------------------------------------------------------------------------------------------------------
  init = new Initialization();
  init->parseInputJSON(path,run);


  // Initializing components: image plane --------------------------------------------------------------------------------------------------------------------------
  mydata = new ImagePlane(init->imgpath,stoi(init->image["pix_x"]),stoi(init->image["pix_y"]),stof(init->image["siz_x"]),stof(init->image["siz_y"]));
  mydata->readB(init->psfpath,stoi(init->psf["pix_x"]),stoi(init->psf["pix_y"]),stoi(init->psf["crop_x"]),stoi(init->psf["crop_y"]));
  mydata->readC(init->noise_flag,init->noisepath);
  mydata->readS(init->maskpath);


  // Initializing components: mass model collection ----------------------------------------------------------------------------------------------------------------
  mycollection = new CollectionMassModels();


  // Update/initialize mass model parameters -----------------------------------------------------------------------------------------------------------------------
  mycollection->setPhysicalPars(init->nlpars_physical);
  mycollection->models.resize(init->mmodel.size());
  for(int k=0;k<mycollection->models.size();k++){
    std::vector<Nlpar*> lens = init->nlpars_lenses[k];
    mycollection->models[k] = FactoryMassModel::getInstance()->createMassModel(init->mmodel[k],lens);
  }
  mycollection->all_defl(mydata); // I need to deflect the image plane once I know the mass model parameters. Needed to create the adaptive grid.


  // Initializing components: source plane -------------------------------------------------------------------------------------------------------------------------
  mysource = FactorySourcePlane::getInstance()->createSourcePlane(init->source);


  // Update/initialize source --------------------------------------------------------------------------------------------------------------------------------------
  //initialize here, but also in each iteration for an adaptive grid                                                  [<---iteration dependent for adaptive source]
  mysource->inMask(mydata);
  if( mysource->type == "adaptive" ){
    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(mysource);
    ada->createAdaGrid(mydata,mycollection);
    ada->createDelaunay();
  }
  if( mysource->reg == "covariance_kernel" || mysource->reg == "covariance_kernel_in_identity_out" ){
    mysource->sample_reg = Nlpar::getSampleReg(init->nlpars_reg_s);
    std::string suffix = Nlpar::removeSuffix(init->nlpars_reg_s); // need to remove _dpsi suffix for the dpsi parameters to update the kernel, then I add it bac
    mysource->kernel = FactoryCovKernel::getInstance()->createCovKernel(init->source["kernel"],init->nlpars_reg_s);
    Nlpar::addSuffix(init->nlpars_reg_s,suffix);
  }
  mysource->constructH();


  // Initialize perturbations --------------------------------------------------------------------------------------------------------------------------------------
  if( init->likeModel == "perturbations_standard" || init->likeModel == "both" ){
    pert_mass_model = new Pert(std::stoi(init->perturbations["pix_x"]),std::stoi(init->perturbations["pix_y"]),mydata,init->perturbations["reg_dpsi"]);
    if( init->perturbations["reg_dpsi"] == "covariance_kernel" ){
      pert_mass_model->dpsi->sample_reg = Nlpar::getSampleReg(init->nlpars_reg_dpsi);
      std::string suffix = Nlpar::removeSuffix(init->nlpars_reg_dpsi); // need to remove _dpsi suffix for the dpsi parameters to update the kernel, then I add it bac
      pert_mass_model->dpsi->kernel = FactoryCovKernel::getInstance()->createCovKernel(init->perturbations["kernel_dpsi"],init->nlpars_reg_dpsi);
      Nlpar::addSuffix(init->nlpars_reg_dpsi,suffix);
    }
  }

  if( init->likeModel == "perturbations_standard" ){
    // HERE I CONSTRUCT source0 either as BaseSourcePlane or as BaseProfile.
    source0 = mysource->clone(); // x,y, and src are not copied (src is also empty at this stage)

    /*
    for(int i=0;i<source0->Sm;i++){
      source0->x[i]   = mysource->x[i];
      source0->y[i]   = mysource->y[i];
      source0->src[i] = init->prof_source0->value(source0->x[i],source0->y[i]);
    }
    */
    int i=0;
    std::string line;
    double xx,yy,vv;
    std::ifstream file("/net/callippus/data/users/gvernard/RESULTS/VKL_tests/test_standard/output/smooth_source_irregular.dat");
    while( std::getline(file,line) ){
      std::istringstream ss(line);
      ss >> vv >> xx >> yy;
      source0->x[i] = xx;
      source0->y[i] = yy;
      source0->src[i] = vv;
      i++;
    }

    AdaptiveSource* ada = dynamic_cast<AdaptiveSource*>(source0);
    ada->createDelaunay(); // required to properly set the Delaunay grid related CGAL variables

    source0->constructDerivatives();
    //    for(int i=0;i<source0->Sm;i++){
    //      printf("%20.5f %20.5f\n",source0->s_dx[i],source0->s_dy[i]);
    //    }
  }


  // Initialize minimizer ------------------------------------------------------------------------------------------------------------------------------------------
  minimizer = FactoryMinimizer::getInstance()->createMinimizer("dum",init->minimizer,init->output);
}




void Initialization::parseInputJSON(std::string path,std::string run){
  Json::Value root;
  Json::Value::Members jmembers;

  std::string runname = path + run;
  std::ifstream fin((runname+"vkl_input.json").c_str());
  fin >> root;


  //General parameters
  this->output = runname+root["output"].asString();
  if( root["maskpath"].asString() == "0" ){
    this->maskpath = "0";
  } else {
    this->maskpath = path+root["maskpath"].asString();
  }
  this->likeModel = root["parameter_model"].asString();


  //Parameters for the image
  this->imgpath = path+root["imgpath"].asString();
  jmembers = root["iplane"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->image[jmembers[i]] = root["iplane"][jmembers[i]].asString();
  }


  //Parameters of the source
  this->interp = root["interp"].asString();
  jmembers = root["splane"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->source[jmembers[i]] = root["splane"][jmembers[i]].asString();
  }


  //Calculate the number of adaptive source pixels based on the image data
  if( this->source["type"] == "adaptive" && ( this->source["mode"] == "image" || this->source["mode"] == "grid" ) ){
    int spacing = std::stoi(this->source["spacing"]);
    int Ni      = std::stoi(this->image["pix_y"]);
    int Nj      = std::stoi(this->image["pix_x"]);

    int count = 0;
    int i0    = (int) floor( (spacing-1)/2. );
    int j0    = (int) floor( (spacing-1)/2. );
    for(int i=i0;i<Ni;i=i+spacing){
      for(int j=j0;j<Nj;j=j+spacing){
	count++;
      }
    }

    this->source["sm"] = std::to_string( count );
    //    std::cout << count << std::endl;
    //    std::cout << (int) floor(Ni/spacing)*floor(Nj/spacing) << std::endl;
  }


  //Parameters for the psf
  this->psfpath = path+root["psfpath"].asString();
  jmembers = root["psf"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->psf[jmembers[i]] = root["psf"][jmembers[i]].asString();
  }


  //Parameters for the noise
  this->noise_flag = root["noise_flag"].asString();
  this->noisepath  = path+root["noisepath"].asString();


  //Parameters for the minimizer/sampler
  jmembers = root["minimizer"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->minimizer[jmembers[i]] = root["minimizer"][jmembers[i]].asString();
  }


  //NLPARS: Physical parameters (shear magnitude and direction)
  this->nlpars_physical = Nlpar::nlparsFromJsonVector(root["physical"]["nlpars"]);


  //NLPARS: Mass model for the lenses
  jmembers = root["lenses"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->mmodel.push_back( root["lenses"][jmembers[i]]["subtype"].asString() );
    this->lens_names.push_back( jmembers[i] );
    this->nlpars_lenses.push_back( Nlpar::nlparsFromJsonVector(root["lenses"][jmembers[i]]["nlpars"]) );
  }


  //NLPARS: Source regularization
  if( root["parameter_model"] == "standard" ){
    this->source["reg_s"] = root["reg_s"]["type"].asString();
    if( this->source["reg_s"] == "covariance_kernel" || this->source["reg_s"] == "covariance_kernel_in_identity_out" ){
      this->source["kernel"] = root["reg_s"]["subtype"].asString();
    }
    this->nlpars_reg_s = Nlpar::nlparsFromJsonVector(root["reg_s"]["nlpars"]);
  } else {
    this->source["reg_s"] = root["perturbations"]["reg_s"]["type"].asString();
    if( this->source["reg_s"] == "covariance_kernel" || this->source["reg_s"] == "covariance_kernel_in_identity_out" ){
      this->source["kernel"] = root["perturbations"]["reg_s"]["subtype"].asString();
    }
    this->nlpars_reg_s = Nlpar::nlparsFromJsonVector(root["perturbations"]["reg_s"]["nlpars"]);
  }
  for(int i=0;i<this->nlpars_reg_s.size();i++){
    this->nlpars_reg_s[i]->nam += "_s";
  }

  //Parameters for the perturbations
  if( root["parameter_model"] == "perturbations_standard" || root["parameter_model"] == "both" ){
    this->perturbations["pix_x"] = root["perturbations"]["pix_x"].asString();
    this->perturbations["pix_y"] = root["perturbations"]["pix_y"].asString();

    //NLPARS: Dpsi regularization
    this->perturbations["reg_dpsi"] = root["perturbations"]["reg_dpsi"]["type"].asString();
    if( this->perturbations["reg_dpsi"] == "covariance_kernel" ){
      this->perturbations["kernel_dpsi"] = root["perturbations"]["reg_dpsi"]["subtype"].asString();
    }
    this->nlpars_reg_dpsi = Nlpar::nlparsFromJsonVector(root["perturbations"]["reg_dpsi"]["nlpars"]);
    for(int i=0;i<this->nlpars_reg_dpsi.size();i++){
      this->nlpars_reg_dpsi[i]->nam += "_dpsi";
    }
  }

  if( root["parameter_model"] == "perturbations_standard" ){
    //source0
    Json::Value jsource0 = root["perturbations"]["source0"];
    if( jsource0["type"].asString() == "analytic" ){
      
      std::vector<std::string> names;
      std::vector<std::map<std::string,double> > all_pars;
      for(int i=0;i<jsource0["pars"].size();i++){
	std::string name = jsource0["pars"][i]["type"].asString();
	names.push_back(name);
	
	std::map<std::string,double> pars;
	const Json::Value::Members jpars = jsource0["pars"][i].getMemberNames();
	for(int j=0;j<jpars.size();j++){
	  if( jpars[j] != "type" ){
	    pars[jpars[j]] = jsource0["pars"][i][jpars[j]].asDouble();
	  }
	}
	all_pars.push_back(pars);
      }
      this->prof_source0 = new Analytic(names,all_pars);
      
    } else if( jsource0["type"].asString() == "delaunay" ){

      //      std::string filename = jsource0["filename"].asString();
      //      this->prof_source0 = new myDelaunay(filename);
      
    } else if( jsource0["type"].asString() == "fromfits" ){
      
      std::string filename = jsource0["filename"].asString();
      int Ni               = jsource0["Ni"].asInt();
      int Nj               = jsource0["Nj"].asInt();
      double height        = jsource0["height"].asDouble();
      double width         = jsource0["width"].asDouble();
      double x0            = jsource0["x0"].asDouble();
      double y0            = jsource0["y0"].asDouble();
      this->prof_source0 = new fromFITS(filename,Ni,Nj,height,width,x0,y0);
      
    }
  }
}
