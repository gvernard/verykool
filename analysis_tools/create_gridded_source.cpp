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
  std::string fname = argv[1];
  int res = std::stoi(argv[2]);


  myDelaunay* source = new myDelaunay(fname);
  double size = source->sourceExtent();


  ImagePlane* splane = new ImagePlane(res,res,size,size);
  //set a fixed source profile
  for(int i=0;i<splane->Nm;i++){
    splane->img[i] = source->value(splane->x[i],splane->y[i]);
  }
  
  // output source profile
  splane->writeImage("source_model.fits");

  
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
