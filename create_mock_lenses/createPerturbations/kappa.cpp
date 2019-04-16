#include "imagePlane.hpp"
#include <iostream>
#include <fstream>
#include "json/json.h"

int main(int argc,char* argv[]){
  //Read in the JSON input
  Json::Value root;
  std::ifstream fin(argv[1]);
  fin >> root;

  std::string path = root["output"].asString();
  int Nj = root["Nwidth"].asInt();      // pixels in x
  int Ni = root["Nheight"].asInt();     // pixels in y
  double w = root["width"].asDouble();  // size in x
  double h = root["height"].asDouble(); // size in y
  double ddx,ddy;

  std::string inname  = path + "dpsi.fits";
  std::string outname = path + "kappa.fits";

  ImagePlane dpsi(inname,Ni,Nj,w,h);
  ImagePlane kappa(dpsi);


  double dx2 = pow( fabs(dpsi.x[1] - dpsi.x[0]), 2);
  double dy2 = pow( fabs(dpsi.y[Nj] - dpsi.y[0]), 2);
  //  std::cout << dx2 << " " << dy2 << std::endl;


  // Calculate second derivatives:

  // first row
  ddx = 1.0*dpsi.img[0] - 2.0*dpsi.img[1] + 1.0*dpsi.img[2];
  ddy = 1.0*dpsi.img[0] - 2.0*dpsi.img[Nj] + 1.0*dpsi.img[2*Nj];
  kappa.img[0] = 0.5*(ddx/dx2 + ddy/dy2);
  for(int j=1;j<Nj-1;j++){
    ddx = 1.0*dpsi.img[j-1] - 2.0*dpsi.img[j] + 1.0*dpsi.img[j+1];
    ddy = 1.0*dpsi.img[j] - 2.0*dpsi.img[Nj+j] + 1.0*dpsi.img[2*Nj+j];
    kappa.img[j] = 0.5*(ddx/dx2 + ddy/dy2);
  }
  ddx = 1.0*dpsi.img[Nj-3] - 2.0*dpsi.img[Nj-2] + 1.0*dpsi.img[Nj-1];
  ddy = 1.0*dpsi.img[Nj-1] - 2.0*dpsi.img[2*Nj-1] + 1.0*dpsi.img[3*Nj-1];
  kappa.img[Nj-1] = 0.5*(ddx/dx2 + ddy/dy2);


  // in-between rows
  for(int i=1;i<Ni-1;i++){
    ddx = 1.0*dpsi.img[i*Nj] - 2.0*dpsi.img[i*Nj+1] + 1.0*dpsi.img[i*Nj+2];
    ddy = 1.0*dpsi.img[(i-1)*Nj] - 2.0*dpsi.img[i*Nj] + 1.0*dpsi.img[(i+1)*Nj];
    kappa.img[i*Nj] = 0.5*(ddx/dx2 + ddy/dy2);
    for(int j=1;j<Nj-1;j++){
      ddx = 1.0*dpsi.img[i*Nj+j-1] - 2.0*dpsi.img[i*Nj+j] + 1.0*dpsi.img[i*Nj+j+1];
      ddy = 1.0*dpsi.img[(i-1)*Nj+j] - 2.0*dpsi.img[i*Nj+j] + 1.0*dpsi.img[(i+1)*Nj+j];
      kappa.img[i*Nj+j] = 0.5*(ddx/dx2 + ddy/dy2);
    }
    ddx = 1.0*dpsi.img[i*Nj+Nj-3] - 2.0*dpsi.img[i*Nj+Nj-2] + 1.0*dpsi.img[i*Nj+Nj-1];
    ddy = 1.0*dpsi.img[(i-1)*Nj+Nj-1] - 2.0*dpsi.img[i*Nj+Nj-1] + 1.0*dpsi.img[(i+1)*Nj+Nj-1];
    kappa.img[i*Nj+Nj-1] = 0.5*(ddx/dx2 + ddy/dy2);
  }


  // last row
  ddx = 1.0*dpsi.img[(Ni-1)*Nj] - 2.0*dpsi.img[(Ni-1)*Nj+1] + 1.0*dpsi.img[(Ni-1)*Nj+2];
  ddy = 1.0*dpsi.img[(Ni-1)*Nj] - 2.0*dpsi.img[(Ni-2)*Nj] + 1.0*dpsi.img[(Ni-3)*Nj];
  kappa.img[(Ni-1)*Nj] = 0.5*(ddx/dx2 + ddy/dy2);
  for(int j=1;j<Nj-1;j++){
    ddx = 1.0*dpsi.img[(Ni-1)*Nj+j-1] - 2.0*dpsi.img[(Ni-1)*Nj+j] + 1.0*dpsi.img[(Ni-1)*Nj+j+1];
    ddy = 1.0*dpsi.img[(Ni-1)*Nj+j] - 2.0*dpsi.img[(Ni-2)*Nj+j] + 1.0*dpsi.img[(Ni-3)*Nj+j];
    kappa.img[(Ni-1)*Nj+j] = 0.5*(ddx/dx2 + ddy/dy2);
  }
  ddx = 1.0*dpsi.img[(Ni-1)*Nj+Nj-3] - 2.0*dpsi.img[(Ni-1)*Nj+Nj-2] + 1.0*dpsi.img[(Ni-1)*Nj+Nj-1];
  ddy = 1.0*dpsi.img[(Ni-1)*Nj+Nj-1] - 2.0*dpsi.img[(Ni-2)*Nj+Nj-1] + 1.0*dpsi.img[(Ni-3)*Nj+Nj-1];
  kappa.img[(Ni-1)*Nj+Nj-1] = 0.5*(ddx/dx2 + ddy/dy2);



  kappa.writeImage(outname);


  return 0;
}
