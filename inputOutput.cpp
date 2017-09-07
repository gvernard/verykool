#include "inputOutput.hpp"

void outputInitial(std::vector<std::map<std::string,BaseNlpar*> > nlpars,std::string output){

  //Set active parameters from vector of sets of non-linear parameters (nlpars; order matters!!!)
  std::vector<myactive> active;
  myactive tmp;
  std::vector<std::string> names;

  //set active parameters (and periodicity) in order: physical, other, lens parameters for each lens
  for(int j=0;j<nlpars.size();j++){
    names = BaseNlpar::getActive(nlpars[j]);
    for(int i=0;i<names.size();i++){
      tmp = {j,names[i],nlpars[j][names[i]]->per};
      active.push_back(tmp);
    }
    names.clear();
  }


  // File with the parameter names
  std::ofstream f_names(output+"corner.paramnames",std::ofstream::out);
  for(int i=0;i<active.size();i++){
    f_names << nlpars[active[i].index][active[i].nam]->nam;
    f_names << " ";
    f_names << nlpars[active[i].index][active[i].nam]->nam;
    f_names << std::endl;
  }
  f_names.close();


  // File with the parameter ranges
  std::ofstream f_ranges(output+"corner.ranges",std::ofstream::out);
  for(int i=0;i<active.size();i++){
    f_ranges << nlpars[active[i].index][active[i].nam]->nam;
    f_ranges << " ";
    f_ranges << nlpars[active[i].index][active[i].nam]->min;
    f_ranges << " ";
    f_ranges << nlpars[active[i].index][active[i].nam]->max;
    f_ranges << std::endl;
  }
  f_ranges.close();

}




void outputGeneric(ImagePlane* image,BaseSourcePlane* source,std::vector<std::map<std::string,BaseNlpar*> > nlpars,precomp* pcomp,std::string output){
  
  // Print MAP parameters from init.nlpars
  typedef std::map<std::string,BaseNlpar*>::iterator it_nlpar;
  for(int i=0;i<nlpars.size();i++){
    for(it_nlpar myiterator=nlpars[i].begin();myiterator!=nlpars[i].end();myiterator++){
      std::cout << myiterator->first << " " << myiterator->second->val << std::endl;
      //    printf("%8s: %7.3f\n",iterator->second->nam,iterator->second->val);
    }
  }
  
  // Output reconstructed source
  source->outputSource(output + "vkl_voronoi.dat");
  
  // Output errors of reconstructed source
  double* errors = (double*) calloc(source->Sm,sizeof(double));
  getSourceErrors(source->Sm,errors,pcomp);
  source->outputSourceErrors(errors,output + "vkl_voronoi_errors.dat");
  
  // Create mock data (lensed MAP source)
  ImagePlane mock_data(image->Nj,image->Ni,image->width,image->height);
  getMockData(&mock_data,source,pcomp);
  
  // Output image of the model (model)
  mock_data.writeImage(output + "vkl_image.fits");
  
  // Output residual image (diference between my_data and mock_data)
  ImagePlane residual(image->Nj,image->Ni,image->width,image->height);
  for(int i=0;i<image->Ni;i++){
    for(int j=0;j<image->Nj;j++){
      residual.img[i*image->Nj+j] = image->img[i*image->Nj+j] - mock_data.img[i*image->Nj+j];
    }
  }
  residual.writeImage(output + "vkl_residual.fits");



  // Output general parameters, and final non-linear parameter values in json
  Json::Value json_output;
  

  Json::Value pars;
  typedef std::map<std::string,BaseNlpar*>::iterator it_nlpar;
  for(int i=0;i<nlpars.size();i++){
    for(it_nlpar myiterator=nlpars[i].begin();myiterator!=nlpars[i].end();myiterator++){
      pars[myiterator->first] = myiterator->second->val;
    }
  }
  json_output["pars"] = pars;


  Json::Value other;
  other["Nsource"] = source->Sm;
  other["Ndata"]   = image->Nm;
  other["Nmask"]   = static_cast<int>( image->lookup.size() );
  other["Psize"]   = image->width/image->Ni;
  json_output["other"]  = other;


  Json::Value terms;
  double pi          =  3.14159265358979323846;
  terms["ED"]        =  pcomp->chi2/2.;
  terms["ES"]        =  pcomp->reg/2.;
  terms["logdetA"]   = -pcomp->detA/2.;
  terms["logdetHtH"] =  pcomp->detHtH/2.;
  terms["logdetC"]   =  pcomp->detC/2.;
  terms["logl"]      =  source->Sm*log10(nlpars[1]["lambda"]->val)/2.;
  terms["log2pi"]    = -image->Nm*log10(2*pi)/2.;
  json_output["terms"] = terms;


  std::ofstream jsonfile(output+"vkl_output.json");
  jsonfile << json_output;
  jsonfile.close();
}




int Initialization::parseInputJSON(const char* c_path,const char* c_filename){
  Json::Value root;
  Json::Value::Members jmembers;

  std::string path(c_path);
  std::string filename(c_path);
  filename.append(c_filename);
  std::ifstream fin(filename);
  fin >> root;

  //General parameters
  this->output = path+root["output"].asString();
  if( root["maskpath"].asString() == "0" ){
    this->maskpath = "0";
  } else {
    this->maskpath = path+root["maskpath"].asString();
  }
  

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
  this->source["reg"] = root["reg"].asString();
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
  this->covpath    = path+root["covpath"].asString();

  //Parameters for the minimizer/sampler
  this->postmcmc = root["postmcmc"].asInt();
  this->minimizer["type"] = root["minimizer_type"].asString();
  jmembers = root[ this->minimizer["type"] ].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    this->minimizer[jmembers[i]] = root[ this->minimizer["type"] ][jmembers[i]].asString();
  }

  //Non-linear parameters: read the physical parameters in the class variables, always at position [0] nlpars.
  this->nlpars.push_back( nlparsFromJson(root["physical"]) );

  //Non-linear parameters: Read the other parameters (regularization lambda, etc) in the class variables, always at position [1] of nlpars.
  this->nlpars.push_back( nlparsFromJson(root["other"]) );

  //Non-linear parameters: read the mass model type and its parameters in the class variables, always at positions [2-*] of nlpars.
  const Json::Value lenses = root["lenses"];
  for(int i=0;i<lenses.size();i++){
    this->mmodel.push_back( lenses[i]["subtype"].asString() );
    this->nlpars.push_back( nlparsFromJson(lenses[i]["nlpars"]) );
  }



  return 0;
}


std::map<std::string,BaseNlpar*> Initialization::nlparsFromJson(const Json::Value myjson){
  std::map<std::string,BaseNlpar*> mymap;
  BaseNlpar* ptr;
  for(int i=0;i<myjson.size();i++){
    std::string nam     = myjson[i]["nam"].asString();
    int fix             = myjson[i]["fix"].asInt();
    int per             = myjson[i]["per"].asInt();
    double val          = myjson[i]["val"].asDouble();
    double err          = myjson[i]["err"].asDouble();
    double min          = myjson[i]["min"].asDouble();
    double max          = myjson[i]["max"].asDouble();
    std::string pri_nam = myjson[i]["pri_nam"].asString();

    std::map<std::string,double> pri_par;
    Json::Value::Members jmembers = myjson[i]["pri_par"].getMemberNames();
    for(int i=0;i<jmembers.size();i++){
      pri_par[jmembers[i]] = myjson[i][jmembers[i]].asDouble();
    }

    ptr = FactoryNlpar::getInstance()->createNlpar(nam,fix,per,val,err,min,max,pri_nam,pri_par);
    mymap[ptr->nam] = ptr;
  }
  
  return mymap;
}
