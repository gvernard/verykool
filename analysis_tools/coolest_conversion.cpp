#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <dirent.h>
#include <type_traits>

#include "json/json.h"
#include <CCfits/CCfits>

#include "nonLinearPars.hpp"


// This code should run after a run of VKL has completed.
// Input:
// - The path to the case. We need the data, noise, mask, etc to produce the final .coolest files.
// - The name of the run.
// - The index of the output, e.g. could be the latest output.


class Param {
public:
  std::string name;
  double value;
  double mean;
  double p16;
  double p84;
  bool has_stats = false;

  Param(){};
  
  Param(std::string a,double b){
    this->name = a;
    this->value = b;
    this->has_stats = false;
  }
  
  Param(std::string a,double b,double c,double d,double e){
    this->name  = a;
    this->value = b;
    this->mean  = c;
    this->p16   = d;
    this->p84   = e;
    this->has_stats = true;
  }
};


class BaseConvModel {
public:
  std::string vkl_name;
  std::string coolest_name;
  std::map<std::string,std::string> names_lookup;

  BaseConvModel(std::string model_vkl_name){
    this->vkl_name = model_vkl_name;
  };
  
  std::map<std::string,std::string> get_names_lookup(){
    return names_lookup;
  }

  virtual std::vector<Param> convert_values(std::vector<Param> vkl_params) = 0;
};


class Shear: public BaseConvModel {
public:
  Shear(std::string model_vkl_name) : BaseConvModel(model_vkl_name) {
    this->coolest_name = "ExternalShear";
    this->names_lookup = {{"g","gamma_ext"},{"phi","phi_ext"}};
  };

  std::vector<Param> convert_values(std::vector<Param> vkl_params){
    std::vector<Param> coolest_params(vkl_params.size());

    for(int i=0;i<vkl_params.size();i++){
      std::string name = vkl_params[i].name;
      Param dum_param;
      dum_param.name  = this->names_lookup[name];
      if( name == "phi" ){
	// Convert position angle
	dum_param.value = vkl_params[i].value - 90.0;
	if( vkl_params[i].has_stats ){
	  dum_param.mean  = vkl_params[i].mean - 90.0;
	  dum_param.p16   = vkl_params[i].p16 - 90.0;
	  dum_param.p84   = vkl_params[i].p84 - 90.0;
	}
      } else {
	// No conversion (only the name)
	dum_param.value = vkl_params[i].value; 
	if( vkl_params[i].has_stats ){
	  dum_param.mean  = vkl_params[i].mean;
	  dum_param.p16   = vkl_params[i].p16;
	  dum_param.p84   = vkl_params[i].p84;
	}
      }
      coolest_params[i] = dum_param;
    }
    
    return coolest_params;
  }; 
};


class Sie: public BaseConvModel {
public:
  Sie(std::string model_vkl_name) : BaseConvModel(model_vkl_name) {
    this->coolest_name = "SIE";
    this->names_lookup = {{"x0","center_x"},{"y0","center_y"},{"pa","phi"},{"q","q"},{"b","theta_E"}};
  };

  std::vector<Param> convert_values(std::vector<Param> vkl_params){
    std::vector<Param> coolest_params(vkl_params.size());

    for(int i=0;i<vkl_params.size();i++){
      std::string name = vkl_params[i].name;
      Param dum_param;
      dum_param.name  = this->names_lookup[name];
      if( name == "pa" ){
	// Convert position angle
	dum_param.value = vkl_params[i].value - 90.0;
	if( vkl_params[i].has_stats ){
	  dum_param.mean  = vkl_params[i].mean - 90.0;
	  dum_param.p16   = vkl_params[i].p16 - 90.0;
	  dum_param.p84   = vkl_params[i].p84 - 90.0;
	}
      } else if( name == "b" ){
	// Convert 'b' to 'theta_E', for which we need the value of "q" first
	double fac;
	for(int k=0;k<vkl_params.size();k++){
	  if( vkl_params[k].name == "q" ){
	    fac = sqrt(vkl_params[k].value);
	    break;
	  }
	}
	dum_param.value = vkl_params[i].value/fac;
	if( vkl_params[i].has_stats ){
	  dum_param.mean  = vkl_params[i].mean/fac;
	  dum_param.p16   = vkl_params[i].p16/fac;
	  dum_param.p84   = vkl_params[i].p84/fac;
	}
      } else {
	// No conversion (only the name)
	dum_param.value = vkl_params[i].value; 
	if( vkl_params[i].has_stats ){
	  dum_param.mean  = vkl_params[i].mean;
	  dum_param.p16   = vkl_params[i].p16;
	  dum_param.p84   = vkl_params[i].p84;
	}
      }
      coolest_params[i] = dum_param;
    }
    
    return coolest_params;
  };
};


class Spemd: public Sie {
public:
  Spemd(std::string model_vkl_name) : Sie(model_vkl_name) {
    this->coolest_name = "SPEMD";
    this->names_lookup = {{"x0","center_x"},{"y0","center_y"},{"gam","gamma"},{"pa","phi"},{"q","q"},{"b","theta_E"},{"s","r_core"}};
  };
};


class Pemd: public Sie {
public:
  Pemd(std::string model_vkl_name) : Sie(model_vkl_name) {
    this->coolest_name = "PEMD";
    this->names_lookup = {{"x0","center_x"},{"y0","center_y"},{"gam","gamma"},{"pa","phi"},{"q","q"},{"b","theta_E"}};
  };
};


class Source: public BaseConvModel {
public:
  Source(std::string model_vkl_name) : BaseConvModel(model_vkl_name) {
    this->coolest_name = "IrregularGrid";
    this->names_lookup = {{"lambda","lambda"},{"sdev","l_corr"},{"eta","eta"}};
  };

  std::vector<Param> convert_values(std::vector<Param> vkl_params){
    std::vector<Param> coolest_params(vkl_params.size());

    for(int i=0;i<vkl_params.size();i++){
      // No conversion (only the name)
      std::string name = vkl_params[i].name;
      Param dum_param;
      dum_param.name  = this->names_lookup[name];
      dum_param.value = vkl_params[i].value; 
      if( vkl_params[i].has_stats ){
	dum_param.mean  = vkl_params[i].mean;
	dum_param.p16   = vkl_params[i].p16;
	dum_param.p84   = vkl_params[i].p84;
      }
      coolest_params[i] = dum_param;
    }
    
    return coolest_params;
  };
};


class Dpsi: public Source {
public:
  Dpsi(std::string model_vkl_name) : Source(model_vkl_name) {
    this->coolest_name = "PixelatedRegularGridPotential";
  };
};


class FactoryConvModel {//This is a singleton class.
public:
  FactoryConvModel(FactoryConvModel const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryConvModel const&) = delete;
  
  static FactoryConvModel* getInstance(){
    static FactoryConvModel dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }
  
  BaseConvModel* createConvModel_from_COOLEST_name(const std::string model_coolest_name){
    std::string model_vkl_name = this->models_lookup[model_coolest_name];
    return this->createConvModel_from_VKL_name(model_vkl_name);
  }
  
  BaseConvModel* createConvModel_from_VKL_name(const std::string model_vkl_name){
    if( model_vkl_name == "shear" ){
      return new Shear(model_vkl_name);
    } else if( model_vkl_name == "sie" ){
      return new Sie(model_vkl_name);
    } else if( model_vkl_name == "spemd" ){
      return new Spemd(model_vkl_name);
    } else if( model_vkl_name == "pemd" ){
      return new Pemd(model_vkl_name);
    } else if( model_vkl_name == "source" ){
      return new Source(model_vkl_name);  
    } else if( model_vkl_name == "dpsi" ){
      return new Dpsi(model_vkl_name);  
    } else {
      std::cout << "DANGER! Not recognized VKL type of lensing entity: " << model_vkl_name <<std::endl;
      return NULL;
    }
  }
  
private:
  std::map<std::string,std::string> models_lookup = {
    {"ExternalShear","shear"},
    {"SIE","sie"},
    {"SPEMD","spemd"},
    {"PEMD","pemd"},
    {"IrregularGrid","source"},
    {"PixelatedRegularGridPotential","dpsi"}
  };
  FactoryConvModel(){};
};




bool get_limits(Json::Value& mean,Json::Value& s1_lower,Json::Value& s1_upper,Json::Value limits,std::string pname,double pvalue){
  if( limits.isNull() ){
    mean     = Json::Value::null;;
    s1_lower = Json::Value::null;
    s1_upper = Json::Value::null;
    return false;
  } else {
    mean     = limits[pname]["mean"].asDouble();
    s1_lower = limits[pname]["s1_lower"].asDouble();
    s1_upper = limits[pname]["s1_upper"].asDouble();
    return true;
  }
}

template<typename mytype>
std::string myjoin(std::vector<mytype> input,std::string delim){
    if( input.empty() ){
      return std::string();
    }
    
    std::vector<std::string> dum(input.size());
    for(int i=0;i<input.size();i++){
      dum[i] = std::to_string(input[i]);
    }
        
    return std::accumulate(dum.begin() + 1, dum.end(), dum[0],
			   [&delim](std::string x, std::string y) {
			     return x + delim + y;
			   }
			   );
}
template <> std::string myjoin<std::string>(std::vector<std::string> input,std::string delim){
    if( input.empty() ){
      return std::string();
    }

    return std::accumulate(input.begin() + 1, input.end(), input[0],
			 [&delim](std::string x, std::string y) {
			   return x + delim + y;
			 }
			 );
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
  std::string fname;
    
  // Read VKL input json file
  Json::Value input_json;
  fin.open((runname+"vkl_input.json").c_str(),std::ifstream::in);
  fin >> input_json;
  fin.close();

  
  // Getting/setting the initial values of some parameters
  std::string lname;
  if( input_json["parameter_model"] == "standard" ){
    lname = "smooth";
  } else if( input_json["parameter_model"] == "both" ){
    lname = "both";    
  } else {
    std::cout << "Unknown parameter model! Must be: 'standard', 'pert', 'both'" << std::endl;
    return 0;
  }
  std::string minimizer = input_json["minimizer"]["type"].asString();
  
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
  if( minimizer == "multinest" ){
    fin.open((runname+"analysis/"+step+"limits.json").c_str(),std::ifstream::in);
    if( fin.good() ){
      fin.open(fixed_file.c_str(),std::ifstream::in);
      fin >> limits;
    } else {
      std::cout << "File: '" << runname << "analysis/" << step << "limits.json' not found, proceeding without limits" << std::endl;
    }
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

  Json::Value a_mass_model;
  Json::Value mean,s1_lower,s1_upper;
  Json::Value param,this_point_estimate,this_posterior_stats,this_prior;
  bool has_limits;
  double val;
  BaseConvModel* conv_model;

  
  
  
  // External shear is its own lensing entity
  // #############################################################################################################
  std::vector<Nlpar*> nlpars_physical;
  nlpars_physical = Nlpar::nlparsFromJsonVector(input_json["physical"]["nlpars"]);
  

  // Construct vector of Param instances
  std::vector<Param> vkl_params(nlpars_physical.size());
  for(int i=0;i<nlpars_physical.size();i++){
    if( nlpars_physical[i]->fix == 0 ){
      val = output_json["parameters"]["shear_"+nlpars_physical[i]->nam]["maxlike"].asDouble();
    } else {
      val = nlpars_physical[i]->val;
    }
    has_limits = get_limits(mean,s1_lower,s1_upper,limits,nlpars_physical[i]->nam,nlpars_physical[i]->val);
    if( has_limits ){
      vkl_params[i] = Param(nlpars_physical[i]->nam,val,mean.asDouble(),s1_lower.asDouble(),s1_upper.asDouble());
    } else {
      vkl_params[i] = Param(nlpars_physical[i]->nam,val);
    }
  }
    
  // Apply conversions
  conv_model = FactoryConvModel::getInstance()->createConvModel_from_VKL_name("shear");
  std::vector<Param> coolest_params = conv_model->convert_values(vkl_params);
  // for(int k=0;k<vkl_params.size();k++){
  //   std::cout << vkl_params[k].name << " " << vkl_params[k].value << "       " << coolest_params[k].name << " " << coolest_params[k].value << std::endl;
  // }  
  
  // Create Lensing Entity json object
  a_mass_model["type"] = "ExternalShear";
  a_mass_model["parameters"] = Json::Value();
  for(int i=0;i<coolest_params.size();i++){
    this_point_estimate = point_estimate;
    this_point_estimate["value"] = coolest_params[i].value;
    
    this_posterior_stats = posterior_stats;
    if( coolest_params[i].has_stats ){
      this_posterior_stats["mean"] = coolest_params[i].mean;
      this_posterior_stats["percentile_16th"] = coolest_params[i].p16;
      this_posterior_stats["percentile_84th"] = coolest_params[i].p84;
    }
    
    this_prior = prior;

    param["point_estimate"]  = this_point_estimate;
    param["posterior_stats"] = this_posterior_stats;
    param["prior"]           = this_prior;
    a_mass_model["parameters"][coolest_params[i].name] = param;
  }

  Json::Value gamma;
  gamma["type"] = "MassField";
  gamma["name"] = "An external shear";
  gamma["redshift"] = lens_redshift;
  gamma["mass_model"] = Json::Value(Json::arrayValue);
  gamma["mass_model"].append( a_mass_model );
  lensing_entities.append( gamma );


  
  // Loop over lenses and add them to the list of lensing entities
  // #############################################################################################################
  jmembers = input_json["lenses"].getMemberNames();
  for(int m=0;m<jmembers.size();m++){
    std::vector<Nlpar*> nlpars_lens;
    nlpars_lens = Nlpar::nlparsFromJsonVector(input_json["lenses"][jmembers[m]]["nlpars"]);

    // Construct vector of Param instances
    std::vector<Param> vkl_params(nlpars_lens.size());
    for(int i=0;i<nlpars_lens.size();i++){
      if( nlpars_lens[i]->fix == 0 ){
	val = output_json["parameters"][jmembers[m]+"_"+nlpars_lens[i]->nam]["maxlike"].asDouble();
      } else {
	val = nlpars_physical[i]->val;
      }
      has_limits = get_limits(mean,s1_lower,s1_upper,limits,nlpars_lens[i]->nam,nlpars_lens[i]->val);
      if( has_limits ){
	vkl_params[i] = Param(nlpars_lens[i]->nam,val,mean.asDouble(),s1_lower.asDouble(),s1_upper.asDouble());
      } else {
	vkl_params[i] = Param(nlpars_lens[i]->nam,val);
      }
    }

    // Apply conversions
    std::string model_vkl_name = input_json["lenses"][jmembers[m]]["subtype"].asString();
    conv_model = FactoryConvModel::getInstance()->createConvModel_from_VKL_name(model_vkl_name);
    std::vector<Param> coolest_params = conv_model->convert_values(vkl_params);
    // for(int k=0;k<vkl_params.size();k++){
    //   std::cout << vkl_params[k].name << " " << vkl_params[k].value << "       " << coolest_params[k].name << " " << coolest_params[k].value << std::endl;
    // }
    
    // Create Lensing Entity json object
    a_mass_model["type"] = conv_model->coolest_name;
    a_mass_model["parameters"] = Json::Value();
    for(int i=0;i<coolest_params.size();i++){
      this_point_estimate = point_estimate;
      this_point_estimate["value"] = coolest_params[i].value;
    
      this_posterior_stats = posterior_stats;
      if( coolest_params[i].has_stats ){
	this_posterior_stats["mean"] = coolest_params[i].mean;
	this_posterior_stats["percentile_16th"] = coolest_params[i].p16;
	this_posterior_stats["percentile_84th"] = coolest_params[i].p84;
      }
      
      this_prior = prior;
      
      param["point_estimate"]  = this_point_estimate;
      param["posterior_stats"] = this_posterior_stats;
      param["prior"]           = this_prior;
      a_mass_model["parameters"][coolest_params[i].name] = param;
    }    
    Json::Value lens;
    lens["type"] = "Galaxy";
    lens["name"] = jmembers[m];
    lens["redshift"] = lens_redshift;
    lens["mass_model"] = Json::Value(Json::arrayValue);
    lens["mass_model"].append( a_mass_model );
    lens["light_model"] = Json::Value(Json::arrayValue);

    
    lensing_entities.append( lens );
  }




  // Add the Dpsi as a separate lensing entity
  // #############################################################################################################  
  if( lname == "both" ){

    int N_dpsi_x = input_json["dpsi"]["pix_x"].asInt();
    int N_dpsi_y = input_json["dpsi"]["pix_y"].asInt();
    std::string fname = "vkl_dpsi.fits";

    // Copy perturbations .fits file
    std::ifstream src((runname+"output/" + step + lname + "_dpsi.fits").c_str(),std::ios::binary);
    std::ofstream dst((outdir_name + "/" + fname).c_str(), std::ios::binary);
    dst << src.rdbuf();

    // Create the remaining json fields
    Json::Value dpsi;
    dpsi["type"] = "Galaxy";
    dpsi["name"] = "A field of VKL pixellated potential corrections";
    dpsi["redshift"] = lens_redshift;
    dpsi["mass_model"] = Json::Value(Json::arrayValue);
    dpsi["light_model"] = Json::Value(Json::arrayValue);
    dpsi["name"] = "dpsi";
  
    Json::Value a_mass_model;
    Json::Value pixels_reg;
    pixels_reg["field_of_view_x"] = Json::Value(Json::arrayValue);
    pixels_reg["field_of_view_x"].append(0);
    pixels_reg["field_of_view_x"].append(0);
    pixels_reg["field_of_view_y"] = Json::Value(Json::arrayValue);
    pixels_reg["field_of_view_y"].append(0);
    pixels_reg["field_of_view_y"].append(0);
    pixels_reg["num_pix_x"] = N_dpsi_x;
    pixels_reg["num_pix_y"] = N_dpsi_y;
    pixels_reg["fits_file"] = Json::Value();
    pixels_reg["fits_file"]["path"] = fname;  
    a_mass_model["parameters"] = Json::Value();
    a_mass_model["parameters"]["pixels"] = pixels_reg;
    a_mass_model["type"] = "PixelatedRegularGridPotential";
    dpsi["mass_model"].append( a_mass_model );
    
    lensing_entities.append( dpsi );
  }

  
  

  
  // Add the source as a separate lensing entity
  // #############################################################################################################  
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
  source["type"] = "Galaxy";
  source["name"] = "A VKL source";
  source["redshift"] = source_redshift;
  source["mass_model"] = Json::Value(Json::arrayValue);
  source["light_model"] = Json::Value(Json::arrayValue);

  jmembers = input_json["sources"].getMemberNames();
  source["name"] = jmembers[0];

  
  Json::Value a_light_model;
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
  a_light_model["parameters"] = Json::Value();
  a_light_model["parameters"]["pixels"] = pixels_irr;
  a_light_model["type"] = "IrregularGrid";
  source["light_model"].append( a_light_model );

  lensing_entities.append( source );
  int ind_source = lensing_entities.size()-1;

  coolest["lensing_entities"] = lensing_entities;


  


  
  // This maps the VKL full parameter names, e.g. 'dum_pa', to the COOLEST convention <entity index>-<massfield/galaxy>-<mass/light>-<profile index>-<profile type>-<parameter name> 
  // #############################################################################################################
  std::map<std::string,std::string> mydict;

  // Mass models
  for(int i=0;i<lensing_entities.size();i++){
    
    for(int j=0;j<lensing_entities[i]["mass_model"].size();j++){
      std::string model_coolest_name = lensing_entities[i]["mass_model"][j]["type"].asString();
      conv_model = FactoryConvModel::getInstance()->createConvModel_from_COOLEST_name(model_coolest_name);
      std::map<std::string,std::string> names_lookup = conv_model->get_names_lookup();
      if( model_coolest_name == "ExternalShear" ){
	// External shear
	for(std::map<std::string,std::string>::iterator it = names_lookup.begin();it != names_lookup.end();it++){
	  mydict.insert( {"shear_"+it->first,std::to_string(i)+"-massfield-mass-"+std::to_string(j)+"-ExternalShear-"+it->second} );
	}
      } else {
	// Any other mass model	
	for(std::map<std::string,std::string>::iterator it = names_lookup.begin();it != names_lookup.end();it++){
	  mydict.insert(
			{
			  lensing_entities[i]["name"].asString()+"_"+it->first,
			  std::to_string(i)+"-galaxy-mass-"+std::to_string(j)+"-"+lensing_entities[i]["mass_model"][j]["type"].asString()+"-"+it->second
			}
			);
	  
	}
      }
      
    }
  }

  // Light models, currently only for the source
  int i = ind_source;
  conv_model = FactoryConvModel::getInstance()->createConvModel_from_VKL_name("source");
  std::map<std::string,std::string> names_lookup = conv_model->get_names_lookup();
  for(int j=0;j<lensing_entities[i]["light_model"].size();j++){
    for(std::map<std::string,std::string>::iterator it = names_lookup.begin();it != names_lookup.end();it++){
      mydict.insert(
		    {
		      lensing_entities[i]["name"].asString()+"_"+it->first,
		      std::to_string(i)+"-galaxy-light-"+std::to_string(j)+"-"+lensing_entities[i]["light_model"][j]["type"].asString()+"-"+it->second
		    }
		    );
    }
  }
    
  // for(std::map<std::string,std::string>::iterator it = mydict.begin();it != mydict.end();it++){
  //   std::cout << it->first << " " << it->second << std::endl;
  // }



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
  if( input_json["noise_flag"].asString() == "sigma_map" ){
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
  } else if( input_json["noise_flag"].asString() == "sigma_uniform" ){
    double var;
    fin.open( (case_path+"/"+input_json["noisepath"].asString()).c_str() ,std::ifstream::in);
    fin >> var;
    fin.close();

    noise["type"] = "UniformGaussianNoise";
    noise["std_dev"] = sqrt(var);
  } else {
    std::cout << "Unknown noise type! Must be: 'sigma_map', 'sigma_uniform'" << std::endl;
    return 0;
  }
  coolest["observation"]["noise"] = noise;



  



  // Creating chain file
  if( minimizer == "multinest" ){

    // Open output file
    //FILE* fout = fopen((outdir_name+"/coolest_chain.csv").c_str(),"w");
    std::string chain_file_name = "coolest_chain.csv";
    std::ofstream fout((outdir_name+"/"+chain_file_name).c_str(),std::ifstream::out);
    std::string out_line;

    
    // Read in sampled parameter names
    std::vector<std::string> vkl_names;
    char vkl_name[30];
    char tag[30];
    FILE* fp = fopen((runname+"output/"+lname+"_postdist.paramnames").c_str(),"r");
    while( fscanf(fp,"%s %s",vkl_name,tag) != EOF ){
      vkl_names.push_back(vkl_name);
    }
    fclose(fp);    
    int npar = vkl_names.size();


    // Match VKL to COOLEST parameter names
    std::vector<std::string> coolest_names(npar);
    for(int i=0;i<npar;i++){
      coolest_names[i] = mydict[vkl_names[i]];
      //std::cout << vkl_names[i] << " " << coolest_names[i] << std::endl;
    }
      

    // Find unique COOLEST model names
    std::vector<std::string> trunc(npar);
    for(int i=0;i<npar;i++){
      int last = coolest_names[i].find_last_of("-");
      trunc[i] = coolest_names[i].substr(0,last);
    }
    std::vector<std::string> uniq = trunc; 
    sort(uniq.begin(),uniq.end());
    std::vector<std::string>::iterator it = std::unique(uniq.begin(),uniq.end());
    uniq.resize( std::distance(uniq.begin(),it) );
    int Nuniq = uniq.size();
    
    // Create vector of Conversion Models (ConvModels) that matches the unique COOLEST model names
    std::vector<BaseConvModel*> conv_models(Nuniq);
    for(int i=0;i<Nuniq;i++){
      int last = uniq[i].find_last_of("-");
      std::string model_coolest_name = uniq[i].substr(last+1);
      conv_models[i] = FactoryConvModel::getInstance()->createConvModel_from_COOLEST_name(model_coolest_name);
    }

    // Match the unique models to the final parameters, i.e. create a vector of index vectors.
    std::string substr;
    std::vector< std::vector<int> > indices(Nuniq);
    for(int i=0;i<Nuniq;i++){
      for(int j=0;j<npar;j++){
	int last = coolest_names[j].find_last_of("-");
	substr = coolest_names[j].substr(0,last);
	if( substr == uniq[i] ){
	  indices[i].push_back(j);
	}
      }
    }


    // Now we are ready to process the postdist file
    
    // Write parameter names
    out_line = myjoin(coolest_names,",");
    fout << out_line << ",probability_weights" << std::endl;
    //out_line = myjoin(vkl_names,"\t");
    //fout << out_line << "," << "probabil" <<  std::endl;

    
    // Read VKL multinest postdist file
    std::string line;
    std::ifstream post((runname+"output/"+lname+"_postdist.txt").c_str(),std::ifstream::in);
    if( post.is_open() ){
      while( std::getline(post,line) ){


	// Read all the parameters into the par_vals vector
	std::stringstream ss(line);
	std::string buf;
	double prob;
	std::vector<double> par_vals(npar);
	ss >> buf;
	ss >> buf;
	prob = std::stold(buf);
	for(int i=0;i<npar;i++){
	  ss >> buf;
	  par_vals[i] = std::stold(buf);
	}	

	//out_line = myjoin(par_vals,"\t");
	//fout << out_line << "," << prob <<  std::endl;
	
	// Convert each par_vals entry (valus from the postdist file) through its ConvModel
	for(int i=0;i<Nuniq;i++){
	  int N_model_pars = indices[i].size();

	  // Create vector<Param> holding the VKL parameter values
	  std::vector<Param> vkl_dum(N_model_pars);
	  for(int j=0;j<N_model_pars;j++){
	    int index = indices[i][j];
	    // keep only the VKL parameter name, without the prefix ('shear', <lens_name>, <source_name>, etc)
	    int last = vkl_names[index].find("_");
	    std::string par_name = vkl_names[index].substr(last+1);
	    vkl_dum[j] = Param(par_name,par_vals[index]);
	  }

	  // Apply conversion
	  std::vector<Param> conv_pars = conv_models[i]->convert_values(vkl_dum);

	  // Update par_vals
	  for(int j=0;j<N_model_pars;j++){
	    int index = indices[i][j];
	    par_vals[index] = conv_pars[j].value;
	  }

	}

	// Finally, output the line from postdist into the new file
	out_line = myjoin(par_vals,",");
	fout << out_line << "," << prob <<  std::endl;	
      }      
    }
    post.close();

    fout.close();


    // Add chain file path to meta
    if( coolest["meta"].isNull() ){
      coolest["meta"] = Json::Value::null;
    }
    coolest["meta"]["chain_file_name"] = chain_file_name;
  }
  




  std::ofstream jsonfile(outdir_name+"/coolest_vkl.json");
  jsonfile << coolest;
  jsonfile.close();
  
  return 0;
}
