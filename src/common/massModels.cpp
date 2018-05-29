#include "massModels.hpp"

#include <string>
#include <vector>
#include <map>
#include <cstdlib>

#include "imagePlane.hpp"
#include "tableDefinition.hpp"


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
  double q  = this->mpars["q"];
  double e  = this->mpars["e"];
  double b  = pow(this->mpars["b"],2.0*e)*(2.0-2.0*e)/(q*2.0);
  double pa = this->mpars["pa"] * 0.01745329251;//in rad
  double x0 = this->mpars["x0"];
  double y0 = this->mpars["y0"];
  double s2 = this->mpars["s"] * this->mpars["s"];
  double defl[2] = {0.0,0.0};

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
  this->dpsi = new ImagePlane(filepath,a,b,c,d);

  this->dpsi_dx = (double*) calloc(this->dpsi->Nm,sizeof(double));
  this->dpsi_dy = (double*) calloc(this->dpsi->Nm,sizeof(double));

  this->Ddpsi.Ti = 2*this->dpsi->Nm;
  this->Ddpsi.Tj = this->dpsi->Nm;

  this->i0 = this->dpsi->Ni/2.0;
  this->j0 = this->dpsi->Nj/2.0;
  this->di = this->dpsi->width/(this->dpsi->Ni);
  this->dj = this->dpsi->height/(this->dpsi->Nj);

  updateDpsi(this->dpsi->img);
}

Pert::Pert(ImagePlane* new_dpsi){
  // make a deep copy of the ImagePlane object

  this->dpsi_dx = (double*) calloc(this->dpsi->Nm,sizeof(double));
  this->dpsi_dy = (double*) calloc(this->dpsi->Nm,sizeof(double));

  this->Ddpsi.Ti = 2*this->dpsi->Nm;
  this->Ddpsi.Tj = this->dpsi->Nm;

  this->i0 = this->dpsi->Ni/2.0;
  this->j0 = this->dpsi->Nj/2.0;
  this->di = this->dpsi->width/(this->dpsi->Ni);
  this->dj = this->dpsi->height/(this->dpsi->Nj);

  updateDpsi(this->dpsi->img);
}

void Pert::updateDpsi(double* new_dpsi){
  int Ni = this->dpsi->Ni;
  int Nj = this->dpsi->Nj;

  for(int i=0;i<this->dpsi->Nm;i++){
    this->dpsi->img[i] = new_dpsi[i];
  }  

  // Calculate derivatives:
  // first row
  this->dpsi_dx[0] = (this->dpsi->img[1]-this->dpsi->img[0])/dj;
  this->dpsi_dy[0] = (this->dpsi->img[Nj]-this->dpsi->img[0])/di;
  for(int j=1;j<Nj-1;j++){
    this->dpsi_dx[j] = (this->dpsi->img[j+1]-this->dpsi->img[j-1])/(2*dj);
    this->dpsi_dy[j] = (this->dpsi->img[Nj+j]-this->dpsi->img[j])/di;
  }
  this->dpsi_dx[Nj-1] = (this->dpsi->img[Nj-1]-this->dpsi->img[Nj-2])/dj;
  this->dpsi_dy[Nj-1] = (this->dpsi->img[Nj+Nj-1]-this->dpsi->img[Nj-1])/di;

  // in-between rows
  for(int i=1;i<Ni-1;i++){
    this->dpsi_dx[i*Nj+0] = (this->dpsi->img[i*Nj+1]-this->dpsi->img[i*Nj+0])/dj;
    this->dpsi_dy[i*Nj+0] = (this->dpsi->img[(i+1)*Nj+0]-this->dpsi->img[(i-1)*Nj+0])/(2*di);
    for(int j=1;j<Nj-1;j++){
      this->dpsi_dx[i*Nj+j] = (this->dpsi->img[i*Nj+j+1]-this->dpsi->img[i*Nj+j-1])/(2*dj);
      this->dpsi_dy[i*Nj+j] = (this->dpsi->img[(i+1)*Nj+j]-this->dpsi->img[(i-1)*Nj+j])/(2*di);
    }
    this->dpsi_dx[i*Nj+Nj-1] = (this->dpsi->img[i*Nj+Nj-1]-this->dpsi->img[i*Nj+Nj-2])/dj;
    this->dpsi_dy[i*Nj+Nj-1] = (this->dpsi->img[(i+1)*Nj+Nj-1]-this->dpsi->img[(i-1)*Nj+Nj-1])/(2*di);
  }

  // last row
  this->dpsi_dx[(Ni-1)*Nj+0] = (this->dpsi->img[(Ni-1)*Nj+1]-this->dpsi->img[(Ni-1)*Nj+0])/dj;
  this->dpsi_dy[(Ni-1)*Nj+0] = (this->dpsi->img[(Ni-1)*Nj+0]-this->dpsi->img[(Ni-2)*Nj+0])/di;
  for(int j=1;j<Nj-1;j++){
    this->dpsi_dx[(Ni-1)*Nj+j] = (this->dpsi->img[(Ni-1)*Nj+j+1]-this->dpsi->img[(Ni-1)*Nj+j-1])/(2*dj);
    this->dpsi_dy[(Ni-1)*Nj+j] = (this->dpsi->img[(Ni-1)*Nj+j]-this->dpsi->img[(Ni-2)*Nj+j])/di;
  }
  this->dpsi_dx[(Ni-1)*Nj+Nj-1] = (this->dpsi->img[(Ni-1)*Nj+Nj-1]-this->dpsi->img[(Ni-1)*Nj+Nj-2])/dj;
  this->dpsi_dy[(Ni-1)*Nj+Nj-1] = (this->dpsi->img[(Ni-1)*Nj+Nj-1]-this->dpsi->img[(Ni-2)*Nj+Nj-1])/di;


  // Create Ddpsi table
  std::vector<mytriplet> tmp;//need to make sure that the Ddpsi triplet vector is a new one

  for(int i=0;i<this->dpsi->Nm;i++){
    tmp.push_back({2*i,i,this->dpsi_dx[i]});
    tmp.push_back({2*i+1,i,this->dpsi_dy[i]});
  }

  this->Ddpsi.tri.swap(tmp);
}

void Pert::defl(double xin,double yin,double& xout,double& yout){
  int j = floor(xin/dj+j0);
  int i = floor(-yin/di+i0);

  double ax = this->dpsi_dx[i*this->dpsi->Nj+j];
  double ay = this->dpsi_dy[i*this->dpsi->Nj+j];

  xout = ax;
  yout = ay;
}
