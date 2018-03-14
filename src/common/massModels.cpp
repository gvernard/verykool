#include "massModels.hpp"

#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <CCfits/CCfits>

//Abstract class: BaseMassModel
//===============================================================================================================
void BaseMassModel::setMassPars(std::vector<Nlpar*> nlpars){
  for(int i=0;i<nlpars.size();i++){
    this->mpars[nlpars[i]->nam] = nlpars[i]->val;
  }
}

void BaseMassModel::printMassPars(){
  typedef std::map<std::string,double>::iterator some_type;
  for(some_type iterator=this->mpars.begin();iterator!=this->mpars.end();iterator++){
    std::cout << iterator->first << " " << iterator->second << std::endl;
    //    printf("%8s: %7.3f\n",iterator->second->nam,iterator->second->val);
  }
}


//Derived class from BaseMassModel: Sie (Singular Isothermal Ellipsoid)
//===============================================================================================================
Sie::Sie(std::vector<Nlpar*> nlpars){
  this->n = 5;
  setMassPars(nlpars);
}

void Sie::defl(double xin,double yin,double& xout,double& yout){
  double b  = this->mpars["b"];
  double q  = this->mpars["q"];
  double pa = this->mpars["pa"] * 0.01745329251;//in rad
  double x0 = this->mpars["x0"];
  double y0 = this->mpars["y0"];

  //rotate the coordinate system according to position angle and translate to the lens center
  double x_t =  (xin-x0)*cos(pa) + (yin-y0)*sin(pa);
  double y_t = -(xin-x0)*sin(pa) + (yin-y0)*cos(pa);

  if( fabs(x_t) < 0.0001 && fabs(y_t) < 0.0001 ){
    if( std::signbit(x_t) ){
      x_t = -0.0001;
    } else {
      x_t =  0.0001;
    }
    if( std::signbit(y_t) ){
      y_t = -0.0001;
    } else {
      y_t =  0.0001;
    }
  }

  double fac   = 1.0-q*q;
  double omega = q*q*x_t*x_t + y_t*y_t; // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.
  double fac2  = sqrt(fac/omega);

  double ax_t = (b/sqrt(fac))*atan(x_t*fac2);
  double ay_t = (b/sqrt(fac))*atanh(y_t*fac2);
  
  //rotate back according to position angle, no need to translate (this is equivalent to rotating by -pa using the same equations as above)
  double ax =  ax_t*cos(pa) - ay_t*sin(pa);
  double ay =  ax_t*sin(pa) + ay_t*cos(pa);
  
  xout = ax;
  yout = ay;
}

//Derived class from BaseMassModel: Spemd (Softened Power-law Elliptical Mass Density)
//===============================================================================================================
Spemd::Spemd(std::vector<Nlpar*> nlpars){
  this->n = 7;
  setMassPars(nlpars);
}

void Spemd::defl(double xin,double yin,double& xout,double& yout){
  double b  = this->mpars["b"];
  double q  = this->mpars["q"];
  double pa = this->mpars["pa"] * 0.01745329251;//in rad
  double x0 = this->mpars["x0"];
  double y0 = this->mpars["y0"];
  double e  = this->mpars["e"];
  double s2 = this->mpars["s"] * this->mpars["s"];
  double defl[2] = {0.,0.};

  //rotate according to position angle and translate to the lens center
  double x_t =  (xin-x0)*cos(pa) + (yin-y0)*sin(pa);
  double y_t = -(xin-x0)*sin(pa) + (yin-y0)*cos(pa);

  fastelldefl_(&x_t,&y_t,&b,&e,&q,&s2,defl);

  double ax_t = defl[0];
  double ay_t = defl[1];

  //rotate back according to position angle, no need to translate
  double ax =  ax_t*cos(pa) - ay_t*sin(pa);
  double ay =  ax_t*sin(pa) + ay_t*cos(pa);
  
  xout = ax;
  yout = ay;
}

//Derived class from BaseMassModel: Pert (perturbations on a grid)
//===============================================================================================================
Pert::Pert(std::string filepath,int a,int b,double c,double d){
  Ni     = a;
  Nj     = b;
  width  = c;
  height = d;
  int Nm = Ni*Nj;


  // Read the .fits file with the potential corrections
  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  image.readAllKeys();

  std::valarray<float> psi(image.axis(0)*image.axis(1));
  std::valarray<float> tmp;
  image.read(tmp);


  //convert FITS standard (bottom to top) to the one used in this code (top to bottom)
  for(int j=0;j<Nj;j++){
    for(int i=0;i<Ni;i++){
      psi[j*Ni+i] = tmp[(Nj-j-1)*Ni+i];
    }
  }


  // Set the pixel coordinates
  i0 = Ni/2.0;
  j0 = Nj/2.0;
  di = width/(Ni);
  dj = height/(Nj);


  // Calculate gradients
  dpsidx = (double*) calloc(Nm,sizeof(double));
  dpsidy = (double*) calloc(Nm,sizeof(double));

  // first row
  dpsidx[0] = (psi[1]-psi[0])/dj;
  dpsidy[0] = (psi[Nj]-psi[0])/di;
  for(int j=1;j<Nj-1;j++){
    dpsidx[j] = (psi[j+1]-psi[j-1])/(2*dj);
    dpsidy[j] = (psi[Nj+j]-psi[j])/di;
  }
  dpsidx[Nj-1] = (psi[Nj-1]-psi[Nj-2])/dj;
  dpsidy[Nj-1] = (psi[Nj+Nj-1]-psi[Nj-1])/di;

  // in-between rows
  for(int i=1;i<Ni-1;i++){
    dpsidx[i*Nj+0] = (psi[i*Nj+1]-psi[i*Nj+0])/dj;
    dpsidy[i*Nj+0] = (psi[(i+1)*Nj+0]-psi[(i-1)*Nj+0])/(2*di);
    for(int j=1;j<Nj-1;j++){
      dpsidx[i*Nj+j] = (psi[i*Nj+j+1]-psi[i*Nj+j-1])/(2*dj);
      dpsidy[i*Nj+j] = (psi[(i+1)*Nj+j]-psi[(i-1)*Nj+j])/(2*di);
    }
    dpsidx[i*Nj+Nj-1] = (psi[i*Nj+Nj-1]-psi[i*Nj+Nj-2])/dj;
    dpsidy[i*Nj+Nj-1] = (psi[(i+1)*Nj+Nj-1]-psi[(i-1)*Nj+Nj-1])/(2*di);
  }

  // last row
  dpsidx[(Ni-1)*Nj+0] = (psi[(Ni-1)*Nj+1]-psi[(Ni-1)*Nj+0])/dj;
  dpsidy[(Ni-1)*Nj+0] = (psi[(Ni-1)*Nj+0]-psi[(Ni-2)*Nj+0])/di;
  for(int j=1;j<Nj-1;j++){
    dpsidx[(Ni-1)*Nj+j] = (psi[(Ni-1)*Nj+j+1]-psi[(Ni-1)*Nj+j-1])/(2*dj);
    dpsidy[(Ni-1)*Nj+j] = (psi[(Ni-1)*Nj+j]-psi[(Ni-2)*Nj+j])/di;
  }
  dpsidx[(Ni-1)*Nj+Nj-1] = (psi[(Ni-1)*Nj+Nj-1]-psi[(Ni-1)*Nj+Nj-2])/dj;
  dpsidy[(Ni-1)*Nj+Nj-1] = (psi[(Ni-1)*Nj+Nj-1]-psi[(Ni-2)*Nj+Nj-1])/di;

}

void Pert::defl(double xin,double yin,double& xout,double& yout){
  int j = floor(xin/dj+j0);
  int i = floor(-yin/di+i0);

  double ax = dpsidx[i*Nj+j];
  double ay = dpsidy[i*Nj+j];
  //double ax = 0;
  //double ay = 0;

  xout = ax;
  yout = ay;
}
