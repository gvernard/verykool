#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

#include "json/json.h"
#include <CCfits/CCfits>

#include "nonLinearPars.hpp"


// This code should run after a run of VKL has completed.
// Input:
// - The path to the case. We need the data, noise, mask, etc to produce the final .coolest files.
// - The name of the run.
// - The index of the output, e.g. could be the latest output.


int main(int argc,char* argv[]){
  std::string case_path = argv[1];
  std::string run_name = argv[2];
  std::string step = argv[3];

  // setting some nuissance variables
  Json::Value::Members jmembers;
  std::string runname = case_path + run_name;
  std::ifstream fin;
  std::map<std::string,std::string> names_lookup;
  Json::Value mass_model_1;
  
  // Read json input file
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
  double lens_redshift = 0.5;
  double source_redshift = 1.0;


  
  // Read output parameter file
  Json::Value output_json;
  fin.open((runname+"output/"+step+"_"+lname+"_output.json").c_str(),std::ifstream::in);
  fin >> output_json;
  fin.close();
  

  Json::Value lensing_entities = Json::Value(Json::arrayValue);

  
  // Create lensing entities.
  // Each one must have a light (empty for the lens) and mass (empty for the source) profile.


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
    Json::Value point_estimate;

    if( nlpars_physical[i]->fix == 0 ){
      // get value from output
      point_estimate["value"] = output_json["collapsed_active"][nlpars_physical[i]->nam]["map"];
    } else {
      // use fixed input value
      point_estimate["value"] = nlpars_physical[i]->val;
    }
    
    Json::Value param;
    param["point_estimate"] = point_estimate;

    mass_model_1["parameters"][names_lookup[nlpars_physical[i]->nam]] = param;
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



    if( input_json["lenses"][jmembers[m]]["subtype"].asString() == "spemd" ){
      mass_model_1["type"] = "SPEMD";
      mass_model_1["parameters"] = Json::Value();
      names_lookup = {{"x0","center_x"},{"y0","center_y"},{"gam","gamma"},{"pa","phi"},{"q","q"},{"b","theta_E"},{"s","r_core"}};

      for(int i=0;i<nlpars_lens.size();i++){
	Json::Value point_estimate;
	
	if( nlpars_lens[i]->fix == 0 ){
	  // get value from output
	  point_estimate["value"] = output_json["collapsed_active"][jmembers[m]+"_"+nlpars_lens[i]->nam]["map"];
	} else {
	  // use fixed input value
	  point_estimate["value"] = nlpars_lens[i]->val;
	}
	
	Json::Value param;
	param["point_estimate"] = point_estimate;
	
	mass_model_1["parameters"][names_lookup[nlpars_lens[i]->nam]] = param;
      }

      // NO conversion needed between the VKL and COOLEST parameters for the SPEMD profile.
      // This is especially true for the strength of the deflections, i.e. b_VKL.
      
      lens["mass_model"].append( mass_model_1 );

    } else if( input_json["lenses"][jmembers[m]]["subtype"].asString() == "sie" ){


    } else {
      std::cout << "Unknown mass model! Must be: 'spemd', 'sie'" << std::endl;
      return 0;
    }


    
    lensing_entities.append( lens );
  }
  



  // Add the source as a separate lensing entity
  
  // First create the .fits table for the irregular grid
  FILE* fh = fopen((runname+"output/"+step+"_"+lname+"_source_irregular.dat").c_str(),"r");
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

  std::string fname = "my_VKL_source.fits";
  std::unique_ptr<CCfits::FITS> pFits(nullptr);
  pFits.reset( new CCfits::FITS("!"+fname,CCfits::Write) );
  
  std::string newName("NEW-EXTENSION");
  std::vector<std::string> ColFormats = {"E","E","E"};
  std::vector<std::string> ColNames = {"x","y","z"};
  std::vector<std::string> ColUnits = {"dum","dum","dum"};
  CCfits::Table* newTable = pFits->addTable(newName,N_source,ColNames,ColFormats,ColUnits);
  newTable->column("x").write(src_x,1);  
  newTable->column("y").write(src_y,1);
  newTable->column("z").write(src_z,1);

 
  Json::Value source;
  source["type"] = "galaxy";
  source["name"] = "A VKL source";
  source["redshift"] = source_redshift;
  source["mass_model"] = Json::Value(Json::arrayValue);
  source["light_model"] = Json::Value(Json::arrayValue);
  
  Json::Value light_model_1;
  Json::Value pixels;
  pixels["field_of_view_x"] = Json::Value(Json::arrayValue);
  pixels["field_of_view_x"].append(0);
  pixels["field_of_view_x"].append(0);
  pixels["field_of_view_y"] = Json::Value(Json::arrayValue);
  pixels["field_of_view_y"].append(0);
  pixels["field_of_view_y"].append(0);
  pixels["num_pix"] = N_source;
  pixels["fits_file"] = fname;
  light_model_1["pixels"] = pixels;
  light_model_1["type"] = "IrregularGrid";
  source["light_model"].append( light_model_1 );

  lensing_entities.append( source );












  
  std::ofstream jsonfile("coolest_vkl.json");
  jsonfile << lensing_entities;
  jsonfile.close();

  
  return 0;
}
