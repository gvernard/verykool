#include <algorithm>
#include <string>
#include <cmath>
#include <fstream>
#include <memory>
#include <iostream>
#include <vector>
#include <map>

#include <fftw3.h>

#include "json/json.h"

#include "imagePlane.hpp"



int main(int argc,char* argv[]){

  //=============== BEGIN:PARSE INPUT =======================
  //Read in the JSON input
  Json::Value root;
  if( argc > 1 ){//read json object from file
    std::ifstream fin(argv[1]);
    fin >> root;
  }

  std::string outpath    = root["output"].asString();
  double smear           = root["smear"].asDouble();      // The sdev of the Gaussian in arcsec.
  double threshold       = root["threshold"].asDouble();  // A multiplying factor of the maximum noise.
  std::string noise_flag;
  if( root.isMember("noise_flag") ){
    noise_flag = root["noise_flag"].asString();
  } else {
    noise_flag = "none";
  }

  //Read the image properties as a map of strings
  std::map<std::string,std::string> image;
  const Json::Value iplane = root["iplane"];
  image["pix_x"] = iplane["pix_x"].asString();
  image["pix_y"] = iplane["pix_y"].asString();
  image["width"] = iplane["width"].asString();
  image["height"] = iplane["height"].asString();
  //================= END:PARSE INPUT =======================




  //=============== BEGIN:INITIALIZATION =======================
  //Read clean image data
  ImagePlane mydata(outpath+"image.fits",stoi(image["pix_x"]),stoi(image["pix_y"]),stof(image["width"]),stof(image["height"]));

  double img_max = 0.0;
  for(int i=0;i<mydata->Nm;i++){
    if( mydata.img[i] > img_max ){
      img_max = mydata.img[i];
    }
  }
  double threshold_noise = threshold*img_max;
  

  /*
  double threshold_noise;
  if( noise_flag == "none" ){
    threshold_noise = threshold;
  } else if( noise_flag == "uniform" ){
    std::fstream noisefile(outpath+"noise.dat",std::ios_base::in);
    double sigma;
    noisefile >> sigma;
    threshold_noise = threshold*sigma;
  } else {
    ImagePlane noise(outpath+"noise.fits",stoi(image["pix_x"]),stoi(image["pix_y"]),stof(image["width"]),stof(image["height"]));
    double sigma_max = 0.0;
    for(int k=0;k<noise.Nm;k++){
      if( noise.img[k] > sigma_max ){
	sigma_max = noise.img[k];
      }
    }
    threshold_noise = threshold*sigma_max;    
  }
  */
  //  std::cout << threshold_noise << std::endl;
  //================= END:INITIALIZATION =======================
  



  //=============== BEGIN:CREATE KERNEL =======================
  ImagePlane blur = mydata;
  int Ni = mydata.Ni;
  int Nj = mydata.Nj;
  
  double dx = mydata.width/mydata.Nj;
  double dy = mydata.height/mydata.Ni;
  double factor1 = 1.0/(2.0*M_PI*pow(smear,2));
  double factor2 = 1.0/(2.0*pow(smear,2));
  for(int i=0;i<Ni;i++){
    for(int j=0;j<Nj;j++){
      double x = (j-Nj/2)*dx;
      double y = (i-Ni/2)*dy;
      double e = -factor2*(pow(x,2) + pow(y,2));
      blur.img[i*Nj+j] = factor1*exp(e);
    }
  }

  double blur_sum = 0.0;
  for(int i=0;i<Ni*Nj;i++){
    blur_sum += blur.img[i];
  }
  for(int i=0;i<Ni*Nj;i++){
    blur.img[i] /= blur_sum;
  }

  int bNx = Nj/2.0;
  int bNy = Ni/2.0;
  double* kernel = (double*) calloc(mydata.Ni*mydata.Nj,sizeof(double));
  for(int j=0;j<bNy;j++){
    for(int i=0;i<bNx;i++){
      kernel[j*Ni+i]                    = blur.img[bNy*2*bNx+bNx+j*2*bNx+i];
      kernel[Ni-bNx+j*Ni+i]             = blur.img[bNy*2*bNx+j*2*bNx+i];
      kernel[Ni*(Nj-bNy)+j*Ni+i]        = blur.img[bNx+2*bNx*j+i];
      kernel[Ni*(Nj-bNy)+Ni-bNx+j*Ni+i] = blur.img[2*bNx*j+i];
    }
  }
  //=============== END:CREATE KERNEL =========================

    


  //=============== BEGIN:CONVOLUTION =========================
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
  free(kernel);
  
  //Normalize output
  for(int i=0;i<Ni*Nj;i++){
    mydata.img[i] /= (Ni*Nj);
  }
  //=============== END:CONVOLUTION ===========================

    

  //=============== BEGIN:OUTPUT =======================
  ImagePlane mask = mydata;
  for(int i=0;i<Ni*Nj;i++){
    if( mydata.img[i] > threshold_noise ){
      mask.img[i] = 1;
    } else {
      mask.img[i] = 0;
    }
  }

  mask.writeImage(outpath+"mask.fits");
  //================= END:OUTPUT =======================




  return 0;
}
