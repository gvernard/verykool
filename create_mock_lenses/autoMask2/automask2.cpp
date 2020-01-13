#define _USE_MATH_DEFINES

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
  double threshold       = root["threshold"].asDouble();  // % of the maximum brightness below which set flux to zero
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
  //Read clean data
  ImagePlane mydata(outpath+"image.fits",stoi(image["pix_x"]),stoi(image["pix_y"]),stof(image["width"]),stof(image["height"]));
  ImagePlane mynoise(outpath+"noise_realization.fits",stoi(image["pix_x"]),stoi(image["pix_y"]),stof(image["width"]),stof(image["height"]));

  double img_max = 0.0;
  for(int i=0;i<mydata.Nm;i++){
    mydata.img[i] -= mynoise.img[i];
    if( mydata.img[i] > img_max ){
      img_max = mydata.img[i];
    }
  }

  double threshold_brightness = img_max*threshold;
  for(int i=0;i<mydata.Nm;i++){  
    if( mydata.img[i] > threshold_brightness ){
      mydata.img[i] = 1;
    } else {
      mydata.img[i] = 0;
    }
  }
  
  mydata.writeImage("test.fits");
  //================= END:INITIALIZATION =======================
  



  //=============== BEGIN:FIND RING RADII =======================
  int i = 2; // starting with a 2 pixel radius
  int imax = (int) floor(mydata.Ni/2.0);
  bool clause = true;
  while( clause && i<imax ){
    double R = i*mydata.height/mydata.Ni;
    int Ntheta = (int) ceil( M_PI*R/(2.0*mydata.height/mydata.Ni) );
    double dtheta = M_PI/(2.0*Ntheta);
    for(int j=0;j<Ntheta;j++){
      double theta = j*dtheta;
      int iy = (int) floor(sin(theta)*R/(mydata.height/mydata.Ni));
      int ix = (int) floor(cos(theta)*R/(mydata.height/mydata.Ni));
      
      int index1 = (imax - iy)*mydata.Nj + mydata.Nj/2 + ix;
      int index2 = (imax - iy)*mydata.Nj + mydata.Nj/2 - ix;
      int index3 = (imax + iy)*mydata.Nj + mydata.Nj/2 + ix;
      int index4 = (imax + iy)*mydata.Nj + mydata.Nj/2 - ix;
      
      if( mydata.img[index1] == 1 || mydata.img[index2] == 1 || mydata.img[index3] == 1 || mydata.img[index4] == 1 ){
	clause = false;
	imax = i;
	break;
      }
    }
    
    i++;
  }
  double inner_radius = imax*mydata.height/mydata.Ni;
  //std::cout << "Inner radius is: " << imax << std::endl;



  
  //  FILE* fhh = fopen("indices.dat","w");

  i =  imax+2;  // starting with the inner radius
  imax = (int) floor(mydata.Ni/2.0);
  clause = true;
  while( clause && i<imax ){
    clause = false;
    double R = i*mydata.height/mydata.Ni;
    int Ntheta = (int) ceil( M_PI*R/(2.0*mydata.height/mydata.Ni) );
    double dtheta = M_PI/(2.0*Ntheta);
    for(int j=0;j<Ntheta;j++){
      double theta = j*dtheta;
      int iy = (int) floor(sin(theta)*R/(mydata.height/mydata.Ni));
      int ix = (int) floor(cos(theta)*R/(mydata.height/mydata.Ni));
      
      int index1 = (imax - iy)*mydata.Nj + mydata.Nj/2 + ix;
      int index2 = (imax - iy)*mydata.Nj + mydata.Nj/2 - ix;
      int index3 = (imax + iy)*mydata.Nj + mydata.Nj/2 + ix;
      int index4 = (imax + iy)*mydata.Nj + mydata.Nj/2 - ix;

      //fprintf(fhh,"%d %d %d %d %f %f %f %f\n",index1,index2,index3,index4,mydata.img[index1],mydata.img[index2],mydata.img[index3],mydata.img[index4]);
      //fprintf(fhh,"%d %d\n",ix,iy);
      
      if( mydata.img[index1] == 1 || mydata.img[index2] == 1 || mydata.img[index3] == 1 || mydata.img[index4] == 1 ){
	clause = true;
	break;
      }
    }

    i++;
  }
  double outer_radius = (i-1)*mydata.height/mydata.Ni;
  //  std::cout << "Outer radius is: " << i*mydata.height/mydata.Ni << std::endl;
  //fclose(fhh);
  
  
  //  FILE* fh = fopen("radii.dat","w");
  //  fprintf(fh,"%f\n",inner_radius);
  //  fprintf(fh,"%f\n",outer_radius);
  //  fclose(fh);
  //=============== END:FIND RING RADII =======================


  


  //=============== BEGIN:CREATE KERNEL =======================
  ImagePlane blur = mydata;
  int Ni = mydata.Ni;
  int Nj = mydata.Nj;

  smear = smear*(outer_radius - inner_radius)/3.0; // set the 3 sigma of the Guassian
  
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
  img_max = 0.0;
  for(int i=0;i<mydata.Nm;i++){
    if( mydata.img[i] > img_max ){
      img_max = mydata.img[i];
    }
  }

  ImagePlane mask = mydata;
  threshold_brightness = img_max*threshold;
  for(int i=0;i<mydata.Nm;i++){  
    //    if( mydata.img[i] > threshold_brightness ){
    if( mydata.img[i] > 0.00001 ){
      mask.img[i] = 1;
    } else {
      mask.img[i] = 0;
    }
  }

  mask.writeImage(outpath+"mask.fits");
  //================= END:OUTPUT =======================




  return 0;
}
