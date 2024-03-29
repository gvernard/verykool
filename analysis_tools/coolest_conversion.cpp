#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <dirent.h>

#include "json/json.h"
#include <CCfits/CCfits>

#include "nonLinearPars.hpp"


// This code should run after a run of VKL has completed.
// Input:
// - The path to the case. We need the data, noise, mask, etc to produce the final .coolest files.
// - The name of the run.
// - The index of the output, e.g. could be the latest output.



void get_limits(Json::Value& mean,Json::Value& s1_lower,Json::Value& s1_upper,Json::Value limits,std::string pname,double pvalue){
  if( limits.isNull() ){
    mean     = Json::Value::null;;
    s1_lower = Json::Value::null;
    s1_upper = Json::Value::null;
  } else {
    mean     = limits[pname]["mean"].asDouble();
    s1_lower = limits[pname]["s1_lower"].asDouble();
    s1_upper = limits[pname]["s1_upper"].asDouble();
  }
}



int main(int argc,char* argv[]){
  std::string case_path = argv[1];
  std::string run_name = argv[2];
  
  std::string step;
  std::string fixed_file;

  if( argc == 4 ){
    // No time step given, use final output files
    step = "";
    fixed_file = argv[3];
  } else {
    step = argv[3];
    step = step + "_";
    fixed_file = argv[4];
  }

  // setting some nuissance variables
  Json::Value::Members jmembers;
  std::string runname = case_path + run_name;
  std::ifstream fin;
  std::map<std::string,std::string> names_lookup;
  Json::Value mass_model_1;
  std::string fname;
    
  // Read VKL input json file
  Json::Value input_json;
  fin.open((runname+"vkl_input.json").c_str(),std::ifstream::in);
  fin >> input_json;
  fin.close();

  // Getting/setting the initial values of some parameters
  std::string lname;
  if( input_json["parameter_model"] = "standard" ){
    lname = "smooth";
  } else {
    std::cout << "Unknown parameter model! Must be: 'standard', 'pert', 'both'" << std::endl;
    return 0;
  }

  // Read COOLEST fixed parameter file
  Json::Value coolest;
  fin.open(fixed_file.c_str(),std::ifstream::in);
  fin >> coolest;
  fin.close();
  double lens_redshift = 0.5;
  double source_redshift = 1.0;


  
  // Read VKL output parameter file
  Json::Value output_json;
  //fin.open((runname+"output/"+step+lname+"_output.json").c_str(),std::ifstream::in);
  fin.open((runname+"output/"+step+lname+"_minimizer_output.json").c_str(),std::ifstream::in);
  fin >> output_json;
  fin.close();
  //std::cout << output_json << std::endl;

  // Read limits file
  Json::Value limits;
  fin.open((runname+"analysis/"+step+"limits.json").c_str(),std::ifstream::in);
  if( fin.good() ){
    fin.open(fixed_file.c_str(),std::ifstream::in);
    fin >> limits;
    fin.close();
  } else {
    std::cout << "File: '" << runname << "analysis/" << step << "limits.json' not found, proceeding without limits" << std::endl;
    fin.close();
  }


  








  
 
  // Create output directory
  int dum;
  std::string outdir_name = (runname+"VKL_coolest");
  dum = system( ("rm -rf "+outdir_name).c_str() );
  mkdir(outdir_name.c_str(),0777);
 
  
  // Create lensing entities.
  // #############################################################################################################
  Json::Value lensing_entities = Json::Value(Json::arrayValue);

  Json::Value point_estimate;
  point_estimate["value"] = Json::Value::null;
  
  Json::Value posterior_stats;
  posterior_stats["mean"] = Json::Value::null;
  posterior_stats["median"] = Json::Value::null;
  posterior_stats["percentile_16th"] = Json::Value::null;
  posterior_stats["percentile_84th"] = Json::Value::null;
  
  Json::Value prior;
  prior["type"] = Json::Value::null;

  Json::Value mean,s1_lower,s1_upper;
  Json::Value param,this_point_estimate,this_posterior_stats,this_prior;


  
  
  // External shear is its own lensing entity
  Json::Value gamma;
  gamma["type"] = "external_shear";
  gamma["name"] = "An external shear";
  gamma["redshift"] = lens_redshift;
  gamma["mass_model"] = Json::Value(Json::arrayValue);
  mass_model_1["type"] = "ExternalShear";
  mass_model_1["parameters"] = Json::Value();
  names_lookup = {{"g","gamma_ext"},{"phi","phi_ext"}};
  
  std::vector<Nlpar*> nlpars_physical;
  nlpars_physical = Nlpar::nlparsFromJsonVector(input_json["physical"]["nlpars"]);
  for(int i=0;i<nlpars_physical.size();i++){
    this_point_estimate = point_estimate;
    if( nlpars_physical[i]->fix == 0 ){
      this_point_estimate["value"] = output_json["parameters"][nlpars_physical[i]->nam]["maxlike"].asDouble();
    } else {
      this_point_estimate["value"] = nlpars_physical[i]->val;
    }
    
    this_posterior_stats = posterior_stats;
    get_limits(mean,s1_lower,s1_upper,limits,nlpars_physical[i]->nam,nlpars_physical[i]->val);
    this_posterior_stats["mean"] = mean;
    this_posterior_stats["percentile_16th"] = s1_lower;
    this_posterior_stats["percentile_84th"] = s1_upper;

    this_prior = prior;

    param["point_estimate"] = this_point_estimate;
    param["posterior_stats"] = this_posterior_stats;
    param["prior"] = this_prior;
    mass_model_1["parameters"][names_lookup[nlpars_physical[i]->nam]] = param;
  }

  // Convert position angle from cartesian to East-of-North  
  mass_model_1["parameters"]["phi_ext"]["point_estimate"]["value"] = mass_model_1["parameters"]["phi_ext"]["point_estimate"]["value"].asDouble() - 90;
  if( !mass_model_1["parameters"]["phi_ext"]["posterior_stats"]["mean"].isNull() ){
    mass_model_1["parameters"]["phi_ext"]["posterior_stats"]["mean"] = mass_model_1["parameters"]["phi_ext"]["posterior_stats"]["mean"].asDouble() - 90;
  }
  if( !mass_model_1["parameters"]["phi_ext"]["posterior_stats"]["percentile_16th"].isNull() ){
    mass_model_1["parameters"]["phi_ext"]["posterior_stats"]["percentile_16th"] = mass_model_1["parameters"]["phi_ext"]["posterior_stats"]["percentile_16th"].asDouble() - 90;
  }
  if( !mass_model_1["parameters"]["phi_ext"]["posterior_stats"]["percentile_84th"].isNull() ){      
    mass_model_1["parameters"]["phi_ext"]["posterior_stats"]["percentile_84th"] = mass_model_1["parameters"]["phi_ext"]["posterior_stats"]["percentile_84th"].asDouble() - 90;
  }

  gamma["mass_model"].append( mass_model_1 );
  lensing_entities.append( gamma );



  
  // Loop over lenses and add them to the list of lensing entities
  jmembers = input_json["lenses"].getMemberNames();
  for(int m=0;m<jmembers.size();m++){
    std::vector<Nlpar*> nlpars_lens;
    nlpars_lens = Nlpar::nlparsFromJsonVector(input_json["lenses"][jmembers[m]]["nlpars"]);

    Json::Value lens;
    lens["type"] = "galaxy";
    lens["name"] = jmembers[m];
    lens["redshift"] = lens_redshift;
    lens["mass_model"] = Json::Value(Json::arrayValue);
    lens["light_model"] = Json::Value(Json::arrayValue);


    if( input_json["lenses"][jmembers[m]]["subtype"].asString() == "spemd" || input_json["lenses"][jmembers[m]]["subtype"].asString() == "pemd"){      

      if( input_json["lenses"][jmembers[m]]["subtype"].asString() == "spemd" ){
	mass_model_1["type"] = "SPEMD";
	names_lookup = {{"x0","center_x"},{"y0","center_y"},{"gam","gamma"},{"pa","phi"},{"q","q"},{"b","theta_E"},{"s","r_core"}};
      } else {
	mass_model_1["type"] = "PEMD";
	names_lookup = {{"x0","center_x"},{"y0","center_y"},{"gam","gamma"},{"pa","phi"},{"q","q"},{"b","theta_E"}};
      }
      mass_model_1["parameters"] = Json::Value();

      for(int i=0;i<nlpars_lens.size();i++){
	this_point_estimate = point_estimate;
	if( nlpars_lens[i]->fix == 0 ){
	  this_point_estimate["value"] = output_json["parameters"][jmembers[m]+"_"+nlpars_lens[i]->nam]["maxlike"].asDouble();
	} else {
	  this_point_estimate["value"] = nlpars_lens[i]->val;
	}	

	this_posterior_stats = posterior_stats;
	get_limits(mean,s1_lower,s1_upper,limits,jmembers[m]+"_"+nlpars_lens[i]->nam,nlpars_lens[i]->val);
	this_posterior_stats["mean"] = mean;
	this_posterior_stats["percentile_16th"] = s1_lower;
	this_posterior_stats["percentile_84th"] = s1_upper;

	this_prior = prior;

	param["point_estimate"] = this_point_estimate;
	param["posterior_stats"] = this_posterior_stats;
	param["prior"] = this_prior;
	mass_model_1["parameters"][names_lookup[nlpars_lens[i]->nam]] = param;
      }

      // Convert Einstein radius from VKL to COOLEST
      double vkl_b = mass_model_1["parameters"]["theta_E"]["point_estimate"]["value"].asDouble(); // It is still 'b'!
      double vkl_q = mass_model_1["parameters"]["q"]["point_estimate"]["value"].asDouble();
      mass_model_1["parameters"]["theta_E"]["point_estimate"]["value"] = vkl_b/sqrt(vkl_q);
      if( !mass_model_1["parameters"]["theta_E"]["posterior_stats"]["mean"].isNull() ){
	mass_model_1["parameters"]["theta_E"]["posterior_stats"]["mean"] = 1.0;
      }
      if( !mass_model_1["parameters"]["theta_E"]["posterior_stats"]["percentile_16th"].isNull() ){
	mass_model_1["parameters"]["theta_E"]["posterior_stats"]["percentile_16th"] = 1.0;
      }
      if( !mass_model_1["parameters"]["theta_E"]["posterior_stats"]["percentile_84th"].isNull() ){
	mass_model_1["parameters"]["theta_E"]["posterior_stats"]["percentile_84th"] = 1.0;
      }
      
      // Convert position angle from cartesian to East-of-North
      mass_model_1["parameters"]["phi"]["point_estimate"]["value"] = mass_model_1["parameters"]["phi"]["point_estimate"]["value"].asDouble() - 90;
      if( !mass_model_1["parameters"]["phi"]["posterior_stats"]["mean"].isNull() ){
	mass_model_1["parameters"]["phi"]["posterior_stats"]["mean"] = mass_model_1["parameters"]["phi"]["posterior_stats"]["mean"].asDouble() - 90;
      }
      if( !mass_model_1["parameters"]["phi"]["posterior_stats"]["percentile_16th"].isNull() ){
	mass_model_1["parameters"]["phi"]["posterior_stats"]["percentile_16th"] = mass_model_1["parameters"]["phi"]["posterior_stats"]["percentile_16th"].asDouble() - 90;
      }
      if( !mass_model_1["parameters"]["phi"]["posterior_stats"]["percentile_84th"].isNull() ){      
	mass_model_1["parameters"]["phi"]["posterior_stats"]["percentile_84th"] = mass_model_1["parameters"]["phi"]["posterior_stats"]["percentile_84th"].asDouble() - 90;
      }

      lens["mass_model"].append( mass_model_1 );

    } else if( input_json["lenses"][jmembers[m]]["subtype"].asString() == "sie" ){
      
      mass_model_1["type"] = "SIE";
      mass_model_1["parameters"] = Json::Value();
      names_lookup = {{"x0","center_x"},{"y0","center_y"},{"pa","phi"},{"q","q"},{"b","theta_E"}};
      
      for(int i=0;i<nlpars_lens.size();i++){
	this_point_estimate = point_estimate;
	if( nlpars_lens[i]->fix == 0 ){
	  this_point_estimate["value"] = output_json["parameters"][jmembers[m]+"_"+nlpars_lens[i]->nam]["maxlike"].asDouble();
	} else {
	  this_point_estimate["value"] = nlpars_lens[i]->val;
	}
	
	this_posterior_stats = posterior_stats;
	get_limits(mean,s1_lower,s1_upper,limits,jmembers[m]+"_"+nlpars_lens[i]->nam,nlpars_lens[i]->val);
	this_posterior_stats["mean"] = mean;
	this_posterior_stats["percentile_16th"] = s1_lower;
	this_posterior_stats["percentile_84th"] = s1_upper;

	this_prior = prior;

	param["point_estimate"] = this_point_estimate;
	param["posterior_stats"] = this_posterior_stats;
	param["prior"] = this_prior;	
	mass_model_1["parameters"][names_lookup[nlpars_lens[i]->nam]] = param;
      }      

      // Convert Einstein radius from VKL to COOLEST
      double vkl_b = mass_model_1["parameters"]["theta_E"]["point_estimate"]["value"].asDouble(); // It is still 'b'!
      double vkl_q = mass_model_1["parameters"]["q"]["point_estimate"]["value"].asDouble();
      mass_model_1["parameters"]["theta_E"]["point_estimate"]["value"] = vkl_b/sqrt(vkl_q);
      if( !mass_model_1["parameters"]["theta_E"]["posterior_stats"]["mean"].isNull() ){
	mass_model_1["parameters"]["theta_E"]["posterior_stats"]["mean"] = 1.0;
      }
      if( !mass_model_1["parameters"]["theta_E"]["posterior_stats"]["percentile_16th"].isNull() ){
	mass_model_1["parameters"]["theta_E"]["posterior_stats"]["percentile_16th"] = 1.0;
      }
      if( !mass_model_1["parameters"]["theta_E"]["posterior_stats"]["percentile_84th"].isNull() ){
	mass_model_1["parameters"]["theta_E"]["posterior_stats"]["percentile_84th"] = 1.0;
      }
      
      // Convert position angle from cartesian to East-of-North
      mass_model_1["parameters"]["phi"]["point_estimate"]["value"] = mass_model_1["parameters"]["phi"]["point_estimate"]["value"].asDouble() - 90;
      if( !mass_model_1["parameters"]["phi"]["posterior_stats"]["mean"].isNull() ){
	mass_model_1["parameters"]["phi"]["posterior_stats"]["mean"] = mass_model_1["parameters"]["phi"]["posterior_stats"]["mean"].asDouble() - 90;
      }
      if( !mass_model_1["parameters"]["phi"]["posterior_stats"]["percentile_16th"].isNull() ){
	mass_model_1["parameters"]["phi"]["posterior_stats"]["percentile_16th"] = mass_model_1["parameters"]["phi"]["posterior_stats"]["percentile_16th"].asDouble() - 90;
      }
      if( !mass_model_1["parameters"]["phi"]["posterior_stats"]["percentile_84th"].isNull() ){      
	mass_model_1["parameters"]["phi"]["posterior_stats"]["percentile_84th"] = mass_model_1["parameters"]["phi"]["posterior_stats"]["percentile_84th"].asDouble() - 90;
      }

      lens["mass_model"].append( mass_model_1 );

    } else {

      std::cout << "Unknown mass model! Must be: 'spemd', 'pemd', or 'sie'" << std::endl;
      return 0;

    }
    
    lensing_entities.append( lens );
  }
  

  // Add the source as a separate lensing entity
  
  // First create the .fits table for the irregular grid
  FILE* fh = fopen((runname+"output/"+step+lname+"_source_irregular.dat").c_str(),"r");
  std::vector<float> src_x;
  std::vector<float> src_y;
  std::vector<float> src_z;
  float dumx,dumy,dumz;
  while( fscanf(fh,"%20f %20f %20f\n",&dumz,&dumx,&dumy) != EOF ){
    src_x.push_back(dumx);
    src_y.push_back(dumy);
    src_z.push_back(dumz);
  }
  fclose(fh);
  int N_source = (int) src_x.size();

  fname = "my_VKL_source.fits";
  std::unique_ptr<CCfits::FITS> pFits(nullptr);
  pFits.reset( new CCfits::FITS("!"+outdir_name+"/"+fname,CCfits::Write) );
  
  std::string newName("NEW-EXTENSION");
  std::vector<std::string> ColFormats = {"E","E","E"};
  std::vector<std::string> ColNames = {"x","y","z"};
  std::vector<std::string> ColUnits = {"dum","dum","dum"};
  CCfits::Table* newTable = pFits->addTable(newName,N_source,ColNames,ColFormats,ColUnits);
  newTable->column("x").write(src_x,1);  
  newTable->column("y").write(src_y,1);
  newTable->column("z").write(src_z,1);

  // Then create the remaining json fields
  Json::Value source;
  source["type"] = "galaxy";
  source["name"] = "A VKL source";
  source["redshift"] = source_redshift;
  source["mass_model"] = Json::Value(Json::arrayValue);
  source["light_model"] = Json::Value(Json::arrayValue);
  
  Json::Value light_model_1;
  Json::Value pixels_irr;
  pixels_irr["field_of_view_x"] = Json::Value(Json::arrayValue);
  pixels_irr["field_of_view_x"].append(0);
  pixels_irr["field_of_view_x"].append(0);
  pixels_irr["field_of_view_y"] = Json::Value(Json::arrayValue);
  pixels_irr["field_of_view_y"].append(0);
  pixels_irr["field_of_view_y"].append(0);
  pixels_irr["num_pix"] = N_source;
  pixels_irr["fits_file"] = Json::Value();
  pixels_irr["fits_file"]["path"] = fname;  
  light_model_1["parameters"] = Json::Value();
  light_model_1["parameters"]["pixels"] = pixels_irr;
  light_model_1["type"] = "IrregularGrid";
  source["light_model"].append( light_model_1 );

  lensing_entities.append( source );


  coolest["lensing_entities"] = lensing_entities;




  // Create instrument
  // #############################################################################################################
  double pixel_size = input_json["iplane"]["siz_x"].asDouble()/input_json["iplane"]["pix_x"].asDouble();
  coolest["instrument"]["pixel_size"] = pixel_size;
  fname = "my_PSF.fits";
  
  Json::Value pixels_psf;
  pixels_psf["field_of_view_x"] = Json::Value(Json::arrayValue);
  pixels_psf["field_of_view_x"].append(0);
  pixels_psf["field_of_view_x"].append( pixel_size*input_json["psf"]["pix_x"].asDouble() );
  pixels_psf["field_of_view_y"] = Json::Value(Json::arrayValue);
  pixels_psf["field_of_view_y"].append(0);
  pixels_psf["field_of_view_y"].append( pixel_size*input_json["psf"]["pix_y"].asDouble() );
  pixels_psf["num_pix_x"] = input_json["psf"]["pix_x"].asInt();
  pixels_psf["num_pix_y"] = input_json["psf"]["pix_y"].asInt();
  pixels_psf["fits_file"] = Json::Value();
  pixels_psf["fits_file"]["path"] = fname;
  coolest["instrument"]["psf"]["pixels"] = pixels_psf;
  
  dum = system( ("cp "+case_path+"/"+input_json["psfpath"].asString()+" "+outdir_name+"/"+fname).c_str() );


  
  // Create observation
  // #############################################################################################################
  fname = "my_data.fits";
  
  Json::Value pixels_obs;
  pixels_obs["field_of_view_x"] = Json::Value(Json::arrayValue);
  pixels_obs["field_of_view_x"].append(0);
  pixels_obs["field_of_view_x"].append( input_json["iplane"]["siz_x"].asDouble() );
  pixels_obs["field_of_view_y"] = Json::Value(Json::arrayValue);
  pixels_obs["field_of_view_y"].append(0);
  pixels_obs["field_of_view_y"].append( input_json["iplane"]["siz_y"].asDouble() );
  pixels_obs["num_pix_x"] = input_json["iplane"]["pix_x"].asInt();
  pixels_obs["num_pix_y"] = input_json["iplane"]["pix_y"].asInt();
  pixels_obs["fits_file"] = Json::Value();
  pixels_obs["fits_file"]["path"] = fname;
  coolest["observation"]["pixels"] = pixels_obs;

  dum = system( ("cp "+case_path+"/"+input_json["imgpath"].asString()+" "+outdir_name+"/"+fname).c_str() );
  
  // Observation noise
  Json::Value noise;
  if( input_json["noise_flag"].asString() == "map" ){
    fname = "my_noise_map.fits";

    int Nx = input_json["iplane"]["pix_x"].asInt();
    int Ny = input_json["iplane"]["pix_y"].asInt();
    
    Json::Value pixels_noise;
    pixels_noise["field_of_view_x"] = Json::Value(Json::arrayValue);
    pixels_noise["field_of_view_x"].append(0);
    pixels_noise["field_of_view_x"].append( input_json["iplane"]["siz_x"].asDouble() );
    pixels_noise["field_of_view_y"] = Json::Value(Json::arrayValue);
    pixels_noise["field_of_view_y"].append(0);
    pixels_noise["field_of_view_y"].append( input_json["iplane"]["siz_y"].asDouble() );
    pixels_noise["num_pix_x"] = Nx;
    pixels_noise["num_pix_y"] = Ny;
    pixels_noise["fits_file"] = Json::Value();
    pixels_noise["fits_file"]["path"] = fname;
    noise["type"] = "NoiseMap";
    noise["noise_map"] = pixels_noise;
    
    std::valarray<float> contents(Nx*Ny);
    long naxis    = 2;
    long naxes[2] = {(long) Ny,(long) Nx};
    long Ntot = (long) Nx*Ny;
    std::unique_ptr<CCfits::FITS> pInfile( new CCfits::FITS(case_path+"/"+input_json["noisepath"].asString(),CCfits::Read,true) );
    std::unique_ptr<CCfits::FITS> pOutfile( new CCfits::FITS("!"+outdir_name+"/"+fname,FLOAT_IMG,naxis,naxes) );
    std::vector<long> extAx(2,(long) Ny);
    CCfits::ExtHDU* imageExt = pOutfile->addImage("NEW-EXTENSION",FLOAT_IMG,extAx);
    
    CCfits::PHDU& image = pInfile->pHDU();
    std::valarray<float> tmp;
    image.readAllKeys();
    image.read(tmp);
    for(int i=0;i<Ny;i++){
      for(int j=0;j<Nx;j++){
	contents[i*Nx+j] = sqrt(tmp[i*Nx+j]);
      }
    }
    
    long fpixel(1);
    imageExt->write(fpixel,Ntot,contents);
    pOutfile->pHDU().write(fpixel,Ntot,contents);
  } else if( input_json["noise_flag"].asString() == "uniform" ){
    double var;
    fin.open( (case_path+"/"+input_json["noisepath"].asString()).c_str() ,std::ifstream::in);
    fin >> var;
    fin.close();

    noise["type"] = "UniformGaussianNoise";
    noise["std_dev"] = sqrt(var);
  } else {
    std::cout << "Unknown noise type! Must be: 'map', 'uniform'" << std::endl;
    return 0;
  }
  coolest["observation"]["noise"] = noise;



  
  std::ofstream jsonfile(outdir_name+"/coolest_vkl.json");
  jsonfile << coolest;
  jsonfile.close();
  
  return 0;
}
