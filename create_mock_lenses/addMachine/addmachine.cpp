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



void addUniNoiz(int seed,double sigma,int Ni,int Nj,double* noise){
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
    noise[i] = z1*sigma;
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



class PSF {
public:
  ImagePlane* original_psf = NULL;
  ImagePlane* scaled_psf = NULL;

  PSF(std::string fname,int pix_x,int pix_y,double width,double height,ImagePlane* mydata){
    this->original_psf = new ImagePlane(fname,pix_x,pix_y,width,height);
    interpolatePSF(mydata);
  }

  ~PSF(){
    delete(original_psf);
    delete(scaled_psf);
  }

  void interpolatePSF(ImagePlane* mydata){
    //    double newPixSize  = (mydata->xmax - mydata->xmin)/mydata->Nj;
    double newPixSize  = (mydata->width)/(mydata->Nj);
    double origPixSize = this->original_psf->width/this->original_psf->Nj;

    // Decide on the profile width and height in pixels based on the input profile
    int newNj,newNi;
    if( mydata->width < this->original_psf->width ){
      newNj = mydata->Nj;
      newNi = mydata->Ni;
    } else {
      newNj = floor( mydata->Nj*(0.901*this->original_psf->width)/mydata->width );
      newNi = floor( mydata->Ni*(0.901*this->original_psf->height)/mydata->height );
    }

    double neww    = newNj*newPixSize;
    double xoffset = (this->original_psf->width - neww)/2.0;
    double newh    = newNi*newPixSize;
    double yoffset = (this->original_psf->height - newh)/2.0;

    this->scaled_psf = new ImagePlane(newNi,newNj,newh,neww);

    
    //    double sum = 0.0;
    double x,y,xp,yp,dx,dy,ddx,ddy,w00,w10,w01,w11,f00,f10,f01,f11;
    int ii,jj;

    for(int i=0;i<this->scaled_psf->Ni;i++){
      y  = yoffset+i*newPixSize;
      ii = floor( y/origPixSize );
      yp = ii*origPixSize;
      dy = (y - yp)/origPixSize;
      ddy = (1.0 - dy);

      for(int j=0;j<this->scaled_psf->Nj;j++){
	x  = xoffset+j*newPixSize;
	jj = floor( x/origPixSize );
	xp = jj*origPixSize;
	dx = (x - xp)/origPixSize;
	ddx = (1.0 - dx);

	// first index: i (y direction) second index: j (x direction)
	w00 = ddx*ddy;
	w01 = dx*ddy;
	w10 = dy*ddx;
	w11 = dx*dy;

	f00 = this->original_psf->img[ii*this->original_psf->Nj+jj];
	f01 = this->original_psf->img[ii*this->original_psf->Nj+jj+1];
	f10 = this->original_psf->img[(ii+1)*this->original_psf->Nj+jj];
	f11 = this->original_psf->img[(ii+1)*this->original_psf->Nj+jj+1];

	this->scaled_psf->img[i*this->scaled_psf->Nj+j] = f00*w00 + f10*w10 + f01*w01 + f11*w11;
	//	sum += this->scaled_psf->img[i*this->scaled_psf->Nj+j];
      }
    }

    //    for(int i=0;i<this->scaled_psf->Nm;i++){
    //      this->scaled_psf->img[i] /= sum;
    //    }

  }

};










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
  image["pix_x"]  = iplane["pix_x"].asString();
  image["pix_y"]  = iplane["pix_y"].asString();
  image["width"]  = iplane["width"].asString();
  image["height"] = iplane["height"].asString();
  image["inf_x"]  = iplane["inf_x"].asString();
  image["inf_y"]  = iplane["inf_y"].asString();

  // Read the psf properties as a map of strings
  std::map<std::string,std::string> psf;
  const Json::Value jpsf = root["psf"];
  psf["pix_x"]  = jpsf["pix_x"].asString();
  psf["pix_y"]  = jpsf["pix_y"].asString();
  psf["width"]  = jpsf["width"].asString();
  psf["height"] = jpsf["height"].asString();
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
    PSF mypsf(psfpath,stoi(psf["pix_x"]),stoi(psf["pix_y"]),stof(psf["width"]),stof(psf["height"]),&mydata);
    //    mypsf.scaled_psf->writeImage(output+"psf_scaled.fits");

    // Create psf kernel. Cropped PSF always even in dimensions.
    int Ncropx = floor( stof(psf["crop_x"])/(mypsf.scaled_psf->width/mypsf.scaled_psf->Nj) ); // crop x in arcsec divided by the data pixel size
    int Ncropy = floor( stof(psf["crop_y"])/(mypsf.scaled_psf->height/mypsf.scaled_psf->Ni) );
    if( Ncropx % 2 != 0 ){
      Ncropx += 1;
    }
    if( Ncropy % 2 != 0 ){
      Ncropy += 1;
    }
    int loffx,roffx,toffy,boffy;
    loffx = floor(Ncropx/2.0);
    roffx = floor(Ncropx/2.0);
    toffy = floor(Ncropy/2.0);
    boffy = floor(Ncropy/2.0);

    double sum = 0.0;
    double* blur = (double*) calloc(Ncropx*Ncropy,sizeof(double));
    int offset = (floor(mypsf.scaled_psf->Ni/2.0)-toffy)*mypsf.scaled_psf->Ni + (floor(mypsf.scaled_psf->Nj/2.0)-loffx);
    for(int i=0;i<Ncropy;i++){
      for(int j=0;j<Ncropx;j++){
	blur[i*Ncropx+j] = mypsf.scaled_psf->img[offset+i*mypsf.scaled_psf->Ni+j];
	sum += blur[i*Ncropx+j];
      }
    }
    for(int i=0;i<Ncropx*Ncropy;i++){
      blur[i] /= sum;
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
  //  mydata.writeImage(output+"psf_image.fits");


  // Bin image from 'infinite' to observed resolution
  ImagePlane obs_img(stoi(image["pix_x"]),stoi(image["pix_y"]),stof(image["width"]),stof(image["height"]));
  int* counts = (int*) calloc(obs_img.Nm,sizeof(int));
  double inf_dx = stof(image["width"])/stoi(image["inf_x"]);
  double inf_dy = stof(image["height"])/stoi(image["inf_y"]);
  double obs_dx = stof(image["width"])/stoi(image["pix_x"]);
  double obs_dy = stof(image["height"])/stoi(image["pix_y"]);
  for(int i=0;i<mydata.Ni;i++){
    int ii = (int) floor(i*inf_dy/obs_dy);
    for(int j=0;j<mydata.Nj;j++){
      int jj = (int) floor(j*inf_dx/obs_dx);
      obs_img.img[ii*obs_img.Nj + jj] += mydata.img[i*mydata.Nj + j];
      counts[ii*obs_img.Nj + jj] += 1;
    }
  }
  for(int i=0;i<obs_img.Nm;i++){
    //    obs_img.img[i] = obs_img.img[i]*obs_dx*obs_dy/counts[i];
    obs_img.img[i] = obs_img.img[i]/counts[i];
  }
  //  obs_img.writeImage(output+"binned_image.fits");


  double sum = 0.0;
  for(int i=0;i<obs_img.Nm;i++){
    sum += obs_img.img[i]*obs_dx*obs_dy;
  }
  //  std::cout << sum << std::endl;



  // Add different kinds of noise to the image
  if( noise_flag == "uniform" ){
    double maxdata = *std::max_element(obs_img.img,obs_img.img+obs_img.Ni*obs_img.Nj);
    double sigma = maxdata/stof(noise["sn"]);
    ImagePlane noise_realization(stoi(image["pix_x"]),stoi(image["pix_y"]),stof(image["width"]),stof(image["height"]));
    addUniNoiz(stoi(noise["seed"]),sigma,obs_img.Ni,obs_img.Nj,noise_realization.img);
    for(int i=0;i<obs_img.Nm;i++){
      obs_img.img[i] += noise_realization.img[i];
    }
    noise_realization.writeImage(output+"noise_realization.fits");
    std::ofstream myfile(output+"noise.dat",std::ios::out);
    myfile << sigma*sigma << std::endl;
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
