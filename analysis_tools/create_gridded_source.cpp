#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "massModels.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"

#include "sourceProfile.hpp"
#include "fitsHeader.hpp" // from FProject


int main(int argc,char* argv[]){
  std::string path = argv[1];
  std::string run  = argv[2];
  int res = std::stoi(argv[3]);
  std::string step;
  std::string lmodel;
  std::string out_path;
  if( argc == 5 ){
    lmodel = argv[4];
    step = "";
    out_path = path + run + "output/" + lmodel;
  } else if( argc == 6 ){
    lmodel = argv[4];
    step = argv[5];
    out_path = path + run + "output/" + step + "_" + lmodel;
  } else {
    std::cout <<  "Either 4 or 5 command line arguments required: path, run, resolution, lmodel, <step>" << std::endl;
    std::cout << argc << " provided, exiting!!!" << std::endl;
  }



  myDelaunay* source = new myDelaunay(out_path + "_source_irregular.dat");
  double size = source->sourceExtent();






  FixedSource* splane = new FixedSource(res,res,size,"identity");
  
  //set a fixed source profile
  for(int i=0;i<splane->Si;i++){
    for(int j=0;j<splane->Sj;j++){
      splane->src[i*splane->Sj+j] = source->value(splane->x[j],splane->y[i]);
    }
  }
  
  // output source profile
  splane->outputSource("source_model.fits");
  
  // add size (in physical units) of the source plane fits file to its header
  std::map<std::string,std::string> header;
  char buffer [50];
  sprintf(buffer,"%f",size);
  header["size"] = buffer;
  header["WIDTH"] = buffer;
  header["HEIGHT"] = buffer;
  addFitsHeader("source_model.fits",header);
  
  
  delete(source);
  delete(splane);

  return 0;
}
