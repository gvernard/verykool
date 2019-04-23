#include <algorithm>
#include <string>
#include <cmath>
#include <fstream>
#include <memory>
#include <iostream>
#include <vector>
#include <map>
#include <valarray>

#include <fftw3.h>

#include "json/json.h"

#include "imagePlane.hpp"



void addUniNoiz(int seed,double sigma,int Ni,int Nj,double* data){
  const double two_pi = 2.0*3.14159265358979323846;
  double z1,z2,u1,u2;
  srand48(seed);
  //  double* noiz = (double*) calloc(Ni*Nj,sizeof(double));
  
  //Applying the Box-Muller transformation
  for(int i=0;i<Ni*Nj;i++){
    u1 = drand48();
    u2 = drand48();
    z1 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    //    z2 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    
    //    noiz[i] = z1*sigma;
    data[i] += z1*sigma;
  }
}

void addMapNoiz(int seed,double* sigmas,int Ni,int Nj,double* data){
  const double two_pi = 2.0*3.14159265358979323846;
  double z1,z2,u1,u2;
  srand48(seed);
  //  double* noiz = (double*) calloc(Ni*Nj,sizeof(double));
  
  //Applying the Box-Muller transformation
  for(int i=0;i<Ni*Nj;i++){
    u1 = drand48();
    u2 = drand48();
    z1 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    //    z2 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    
    //    noiz[i] = z1*sigma;
    data[i] += z1*sigmas[i];
  }
}

void addMapRealization(double factor,double* noise,int Ni,int Nj,double* data){
  for(int i=0;i<Ni*Nj;i++){
    data[i] += factor*noise[i];
  }
}

int main(int argc,char* argv[]){

  //=============== BEGIN:PARSE INPUT =======================
  // Read in the JSON input
  Json::Value root;
  if( argc > 1 ){
    std::ifstream fin(argv[1]);
    fin >> root;
  }

  std::string output     = root["output"].asString();
  std::string imgpath    = root["imgpath"].asString();
  std::string psfpath    = root["psfpath"].asString();

  // Read the noise properties as a map of strings
  std::map<std::string,std::string> noise;
  const Json::Value::Members jmembers = root["noise"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    noise[jmembers[i]] = root["noise"][jmembers[i]].asString();
  }
  std::string noise_flag = noise["noise_flag"];

  // Read the image properties as a map of strings
  std::map<std::string,std::string> image;
  const Json::Value iplane = root["iplane"];
  image["pix_x"] = iplane["pix_x"].asString();
  image["pix_y"] = iplane["pix_y"].asString();
  image["width"] = iplane["width"].asString();
  image["height"] = iplane["height"].asString();
  image["inf_x"] = iplane["inf_x"].asString();
  image["inf_y"] = iplane["inf_y"].asString();

  // Read the psf properties as a map of strings
  std::map<std::string,std::string> psf;
  const Json::Value jpsf = root["psf"];
  psf["pix_x"]  = jpsf["pix_x"].asString();
  psf["pix_y"]  = jpsf["pix_y"].asString();
  psf["crop_x"] = jpsf["crop_x"].asString();
  psf["crop_y"] = jpsf["crop_y"].asString();
  //================= END:PARSE INPUT =======================




  //=============== BEGIN:INITIALIZATION =======================
  // Read clean image data
  ImagePlane mydata(imgpath,stoi(image["inf_x"]),stoi(image["inf_y"]),stof(image["width"]),stof(image["height"]));
  //================= END:INITIALIZATION =======================




  //=============== BEGIN:PROCESS IMAGE =======================
  // Convolve image with psf
  if( psfpath != "0" ){
    int Ni = mydata.Ni;
    int Nj = mydata.Nj;
    ImagePlane mypsf(psfpath,stoi(psf["pix_x"]),stoi(psf["pix_y"]),1.0,1.0); // last two arguments are dummy

    // Create psf kernel
    int Pi     = stoi(psf["pix_x"]);
    int Pj     = stoi(psf["pix_y"]);
    int Ncropx = stoi(psf["crop_x"]);
    int Ncropy = stoi(psf["crop_y"]);
    int loffx,roffx,toffy,boffy;

    loffx = floor(Ncropx/2.0);
    if( Ncropx % 2 == 0 ){
      roffx = floor(Ncropx/2.0);
    } else {
      roffx = floor(Ncropx/2.0) + 1;
    }
    toffy = floor(Ncropy/2.0);
    if( Ncropy % 2 == 0 ){
      boffy = floor(Ncropy/2.0);
    } else {
      boffy = floor(Ncropy/2.0) + 1;
    }

    double* blur = (double*) calloc(Ncropx*Ncropy,sizeof(double));
    int offset = (floor(Pi/2.0)-toffy)*Pi + (floor(Pj/2.0)-loffx);
    for(int i=0;i<Ncropy;i++){
      for (int j=0;j<Ncropx;j++){
	blur[i*Ncropx+j] = mypsf.img[offset+i*Pi+j];
      }
    }
        
    int bNx = Ncropx/2.0;
    int bNy = Ncropy/2.0;
    double* kernel = (double*) calloc(Ni*Nj,sizeof(double));
    for(int j=0;j<bNy;j++){
      for(int i=0;i<bNx;i++){
	kernel[j*Ni+i]                    = blur[bNy*2*bNx+bNx+j*2*bNx+i];
	kernel[Ni-bNx+j*Ni+i]             = blur[bNy*2*bNx+j*2*bNx+i];
	kernel[Ni*(Nj-bNy)+j*Ni+i]        = blur[bNx+2*bNx*j+i];
	kernel[Ni*(Nj-bNy)+Ni-bNx+j*Ni+i] = blur[2*bNx*j+i];
      }
    }
    

    // Convolve with psf kernel
    fftw_complex* f_image  = (fftw_complex*) fftw_malloc(Ni*Nj*sizeof(fftw_complex));
    fftw_complex* f_kernel = (fftw_complex*) fftw_malloc(Ni*Nj*sizeof(fftw_complex));
    
    fftw_plan p1;
    p1 = fftw_plan_dft_r2c_2d(Ni,Nj,kernel,f_kernel,FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_destroy_plan(p1);
    
    p1 = fftw_plan_dft_r2c_2d(Ni,Nj,mydata.img,f_image,FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_destroy_plan(p1);
        
    double dum1,dum2;
    for(int i=0;i<Ni;i++) {
      for(int j=0;j<Nj;j++) {
	dum1 = f_image[i*Nj+j][0]*f_kernel[i*Nj+j][0] - f_image[i*Nj+j][1]*f_kernel[i*Nj+j][1];
	dum2 = f_image[i*Nj+j][0]*f_kernel[i*Nj+j][1] + f_image[i*Nj+j][1]*f_kernel[i*Nj+j][0];
	f_image[i*Nj+j][0] = dum1;
	f_image[i*Nj+j][1] = dum2;
      }
    }
        
    p1 = fftw_plan_dft_c2r_2d(Ni,Nj,f_image,mydata.img,FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_destroy_plan(p1);
    
    fftw_free(f_image);
    fftw_free(f_kernel);
    
    // Normalize output
    for(int i=0;i<Ni*Nj;i++){
      mydata.img[i] /= (Ni*Nj);
    }
  }



  // Bin image from 'infinite' to observed resolution
  ImagePlane obs_img(stoi(image["pix_x"]),stoi(image["pix_y"]),stof(image["width"]),stof(image["height"]));
  double inf_dx = stof(image["width"])/stoi(image["inf_x"]);
  double inf_dy = stof(image["height"])/stoi(image["inf_y"]);
  double obs_dx = stof(image["width"])/stoi(image["pix_x"]);
  double obs_dy = stof(image["height"])/stoi(image["pix_y"]);
  for(int i=0;i<mydata.Ni;i++){
    int ii = (int) floor(i*inf_dy/obs_dy);
    for(int j=0;j<mydata.Nj;j++){
      int jj = (int) floor(j*inf_dx/obs_dx);
      obs_img.img[ii*obs_img.Nj + jj] += mydata.img[i*mydata.Nj + j];
    }
  }

    


  // Add different kinds of noise to the image
  if( noise_flag == "uniform" ){
    double maxdata = *std::max_element(obs_img.img,obs_img.img+obs_img.Ni*obs_img.Nj);
    double sigma = maxdata/stof(noise["sn"]);
    addUniNoiz(stoi(noise["seed"]),sigma,obs_img.Ni,obs_img.Nj,obs_img.img);
    std::ofstream myfile(output+"noise.dat",std::ios::out);
    myfile << sigma << std::endl;
    myfile.close();
  } else if( noise_flag == "map" ){
    ImagePlane sigmas(noise["file"],stoi(image["pix_x"]),stoi(image["pix_y"]),stof(image["width"]),stof(image["height"]));
    addMapNoiz(stoi(noise["seed"]),sigmas.img,obs_img.Ni,obs_img.Nj,obs_img.img);
    sigmas.writeImage(output+"noise.fits");
  } else if( noise_flag == "realization" ){
    ImagePlane noise_realization(noise["file"],stoi(image["pix_x"]),stoi(image["pix_y"]),stof(image["width"]),stof(image["height"]));
    addMapRealization(stof(noise["factor"]),noise_realization.img,obs_img.Ni,obs_img.Nj,obs_img.img);
    noise_realization.writeImage(output+"noise.fits");
  } else if( noise_flag == "correlated" ){

  } else {

  }
  //================= END:PROCESS IMAGE =======================




  //=============== BEGIN:OUTPUT =======================
  // Write image with added noise
  obs_img.writeImage(output+"image.fits");
  //================= END:OUTPUT =======================




  return 0;
}
