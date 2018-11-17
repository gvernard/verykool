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
  std::string out_path = argv[1];
  int res = std::stoi(argv[2]);

  myDelaunay* source = new myDelaunay(out_path + "vkl_source_irregular.dat");



  // define size of the source plane here
  double sum = 0.0;
  for(int i=0;i<source->N;i++){
    sum += source->src[i];
  }
  double limit = 0.99*sum;
  
  double psum,xmin,xmax,ymin,ymax;
  double hsize = 0.0;
  double dsize = 0.25;
  int jmax = 20;
  for(int j=1;j<jmax;j++){
    hsize = j*dsize;

    xmin = -hsize;
    xmax =  hsize;
    ymin = -hsize;
    ymax =  hsize;

    psum = 0.0;
    for(int i=0;i<source->N;i++){
      if( xmin < source->x[i] && source->x[i] < xmax && ymin < source->y[i] && source->y[i] < ymax ){
	psum += source->src[i];
      }
    }
    //    std::cout << hsize << " " << psum << std::endl;
    if( psum > limit ){
      break;
    }
  }

  if( hsize == jmax*dsize ){
    std::cout << "Source larger than " << jmax*dsize << " arcsec!!!" << std::endl;
    return 0;
  }
  double size = 2*hsize; // arcsec






  FixedSource* splane = new FixedSource(res,res,size,"identity");
  
  //set a fixed source profile
  source->profile(splane->Sj,splane->Si,splane->x,splane->y,splane->src);
  
  // output source profile
  splane->outputSource(""); // writes to "vkl_source.fits" at the current directory
  
  // add size (in physical units) of the source plane fits file to its header
  std::map<std::string,std::string> header;
  char buffer [50];
  sprintf(buffer,"%f",size);
  header["size"] = buffer;
  addFitsHeader("vkl_source.fits",header);
  
  
  delete(source);
  delete(splane);

  return 0;
}
