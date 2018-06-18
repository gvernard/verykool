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
Pert::Pert(int a,int b,double c,double d){
  this->dpsi = new ImagePlane(a,b,c,d);

  this->dpsi_dx = (double*) calloc(this->dpsi->Nm,sizeof(double));
  this->dpsi_dy = (double*) calloc(this->dpsi->Nm,sizeof(double));

  this->Bdev.Ti = 2*this->dpsi->Nm;
  this->Bdev.Tj = this->dpsi->Nm;

  this->di = this->dpsi->height/(this->dpsi->Ni);
  this->dj = this->dpsi->width/(this->dpsi->Nj);

  createBdev();  
}

Pert::Pert(std::string filepath,int a,int b,double c,double d){
  this->dpsi = new ImagePlane(filepath,a,b,c,d);

  this->dpsi_dx = (double*) calloc(this->dpsi->Nm,sizeof(double));
  this->dpsi_dy = (double*) calloc(this->dpsi->Nm,sizeof(double));

  this->Bdev.Ti = 2*this->dpsi->Nm;
  this->Bdev.Tj = this->dpsi->Nm;

  this->di = this->dpsi->height/(this->dpsi->Ni);
  this->dj = this->dpsi->width/(this->dpsi->Nj);

  updateDpsi(this->dpsi->img);
  createBdev();
}

Pert::Pert(ImagePlane* new_dpsi){
  // make a deep copy of the ImagePlane object

  this->dpsi_dx = (double*) calloc(this->dpsi->Nm,sizeof(double));
  this->dpsi_dy = (double*) calloc(this->dpsi->Nm,sizeof(double));

  this->Bdev.Ti = 2*this->dpsi->Nm;
  this->Bdev.Tj = this->dpsi->Nm;

  this->di = this->dpsi->height/(this->dpsi->Ni);
  this->dj = this->dpsi->width/(this->dpsi->Nj);

  updateDpsi(this->dpsi->img);
  createBdev();
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

}

void Pert::createBdev(){
  // This function creates the table of finite difference coefficients for the x and y derivatives of dpsi.
  // The result is a Ndpsi x Ndpsi sparse matrix.
  // The coefficients are the same as in function Pert::updateDpsi
  int Nj = this->dpsi->Nj;
  int Ni = this->dpsi->Ni;

  std::vector<mytriplet> tmp;//need to make sure that the Ddpsi triplet vector is a new one

  // Calculate derivatives:
  // first row
  tmp.push_back({  0,  0,         -1.0/dj });
  tmp.push_back({  0,  1,          1.0/dj });
  tmp.push_back({  1,  0,         -1.0/di });
  tmp.push_back({  1, Nj,          1.0/di });
  for(int j=1;j<Nj-1;j++){
    tmp.push_back({   2*j,  j-1, -1.0/(2.0*dj) });
    tmp.push_back({   2*j,  j+1,  1.0/(2.0*dj) });
    tmp.push_back({ 2*j+1,    j,       -1.0/di });
    tmp.push_back({ 2*j+1, j+Nj,        1.0/di });
  }
  tmp.push_back({   2*(Nj-1),    Nj-2,      -1.0/dj });
  tmp.push_back({   2*(Nj-1),    Nj-1,       1.0/dj });
  tmp.push_back({ 2*(Nj-1)+1,    Nj-1,      -1.0/di });
  tmp.push_back({ 2*(Nj-1)+1, Nj+Nj-1,       1.0/di });

  // in-between rows
  for(int i=1;i<Ni-1;i++){
    tmp.push_back({   2*(i*Nj),     i*Nj,       -1.0/dj });
    tmp.push_back({   2*(i*Nj),   i*Nj+1,        1.0/dj });
    tmp.push_back({ 2*(i*Nj)+1, (i-1)*Nj, -1.0/(2.0*di) });
    tmp.push_back({ 2*(i*Nj)+1, (i+1)*Nj,  1.0/(2.0*di) });
    for(int j=1;j<Nj-1;j++){
      tmp.push_back({   2*(i*Nj+j),   i*Nj+j-1, -1.0/(2.0*dj) });
      tmp.push_back({   2*(i*Nj+j),   i*Nj+j+1,  1.0/(2.0*dj) });
      tmp.push_back({ 2*(i*Nj+j)+1, (i-1)*Nj+j, -1.0/(2.0*di) });
      tmp.push_back({ 2*(i*Nj+j)+1, (i+1)*Nj+j,  1.0/(2.0*di) });
    }
    tmp.push_back({   2*(i*Nj+Nj-1),     i*Nj+Nj-2,       -1.0/dj });
    tmp.push_back({   2*(i*Nj+Nj-1),     i*Nj+Nj-1,        1.0/dj });
    tmp.push_back({ 2*(i*Nj+Nj-1)+1, (i-1)*Nj+Nj-1, -1.0/(2.0*di) });
    tmp.push_back({ 2*(i*Nj+Nj-1)+1, (i+1)*Nj+Nj-1,  1.0/(2.0*di) });
  }

  // last row
  tmp.push_back({   2*(Ni-1)*Nj,   (Ni-1)*Nj,       -1.0/dj });
  tmp.push_back({   2*(Ni-1)*Nj, (Ni-1)*Nj+1,        1.0/dj });
  tmp.push_back({ 2*(Ni-1)*Nj+1,   (Ni-2)*Nj,       -1.0/di });
  tmp.push_back({ 2*(Ni-1)*Nj+1,   (Ni-1)*Nj,        1.0/di });
  for(int j=1;j<Nj-1;j++){
    tmp.push_back({   2*((Ni-1)*Nj+j), (Ni-1)*Nj+j-1, -1.0/(2.0*dj) });
    tmp.push_back({   2*((Ni-1)*Nj+j), (Ni-1)*Nj+j+1,  1.0/(2.0*dj) });
    tmp.push_back({ 2*((Ni-1)*Nj+j)+1,   (Ni-2)*Nj+j,       -1.0/di });
    tmp.push_back({ 2*((Ni-1)*Nj+j)+1,   (Ni-1)*Nj+j,        1.0/di });
  }
  tmp.push_back({   2*((Ni-1)*Nj+Nj-1), (Ni-1)*Nj+Nj-2,   -1.0/dj });
  tmp.push_back({   2*((Ni-1)*Nj+Nj-1), (Ni-1)*Nj+Nj-1,    1.0/dj });
  tmp.push_back({ 2*((Ni-1)*Nj+Nj-1)+1, (Ni-2)*Nj+Nj-1,   -1.0/di });
  tmp.push_back({ 2*((Ni-1)*Nj+Nj-1)+1, (Ni-1)*Nj+Nj-1,    1.0/di });

  this->Bdev.tri.swap(tmp);
}



void Pert::createAint(ImagePlane* data){
  this->Aint.Ti = 2*data->Nm;
  this->Aint.Tj = 2*this->dpsi->Nm;

  std::vector<mytriplet> tmp;

  int i,j;
  double xa,ya,xb,yb,w00,w10,w01,w11,f00,f10,f01,f11;
  double den = this->di*this->dj;

  int Nj = this->dpsi->Nj;

  for(int k=0;k<data->Nm;k++){
    int j = floor( (data->x[k]+this->dpsi->width/2.0)/dj );
    int i = floor( (data->y[k]+this->dpsi->height/2.0)/di );  

    if( j == this->dpsi->Nj-1 ){
      j = j-2;
    }
    if( i == this->dpsi->Ni-1 ){
      i = i-2;
    }

    ya  = data->y[k] - this->dpsi->y[i+1];
    yb  = this->dpsi->y[i] - data->y[k];
    xa  = data->x[k] - this->dpsi->x[j];
    xb  = this->dpsi->x[j+1] - data->x[k];

    w00 = xb*ya/den;
    w10 = xb*yb/den;
    w01 = xa*ya/den;
    w11 = xa*yb/den;

    tmp.push_back({ 2*k,       2*(i*Nj+j),  w00 });
    tmp.push_back({ 2*k,   2*((i+1)*Nj+j),  w10 });
    tmp.push_back({ 2*k,     2*(i*Nj+j+1),  w01 });
    tmp.push_back({ 2*k, 2*((i+1)*Nj+j+1),  w11 });

    tmp.push_back({ 2*k+1,       2*(i*Nj+j)+1,  w00 });
    tmp.push_back({ 2*k+1,   2*((i+1)*Nj+j)+1,  w10 });
    tmp.push_back({ 2*k+1,     2*(i*Nj+j+1)+1,  w01 });
    tmp.push_back({ 2*k+1, 2*((i+1)*Nj+j+1)+1,  w11 });
  }

  this->Aint.tri.swap(tmp);
}




void Pert::defl(double xin,double yin,double& xout,double& yout){
  int j = floor( (xin+this->dpsi->width/2.0)/dj );
  int i = floor( (yin+this->dpsi->height/2.0)/di );  

  if( j == this->dpsi->Nj-1 ){
    j = j-2;
  }
  if( i == this->dpsi->Ni-1 ){
    i = i-2;
  }
  
  double den,xa,ya,xb,yb,w00,w10,w01,w11,f00,f10,f01,f11;

  // Be careful: the indices count from top left in the interpolation scheme below
  ya  = yin - this->dpsi->y[i+1];
  yb  = this->dpsi->y[i] - yin;
  xa  = xin - this->dpsi->x[j];
  xb  = this->dpsi->x[j+1] - xin;
  den = this->di*this->dj;

  w00 = xb*ya;
  w10 = xb*yb;
  w01 = xa*ya;
  w11 = xa*yb;

  // Derivative on x
  f00 = this->dpsi_dx[i*this->dpsi->Nj+j];
  f10 = this->dpsi_dx[(i+1)*this->dpsi->Nj+j];
  f01 = this->dpsi_dx[i*this->dpsi->Nj+j+1];
  f11 = this->dpsi_dx[(i+1)*this->dpsi->Nj+j+1];
  double ax = (f00*w00 + f10*w10 + f01*w01 + f11*w11)/den;

  // Derivative on y
  f00 = this->dpsi_dy[i*this->dpsi->Nj+j];
  f10 = this->dpsi_dy[(i+1)*this->dpsi->Nj+j];
  f01 = this->dpsi_dy[i*this->dpsi->Nj+j+1];
  f11 = this->dpsi_dy[(i+1)*this->dpsi->Nj+j+1];
  double ay = (f00*w00 + f10*w10 + f01*w01 + f11*w11)/den;

  xout = ax;
  yout = ay;

}

//Class: CollectionMassModels
//===============================================================================================================
CollectionMassModels::CollectionMassModels(){
  this->mpars["g1"] = 0.0;
  this->mpars["g2"] = 0.0;
}
CollectionMassModels::CollectionMassModels(std::vector<Nlpar*> nlpars){
  this->setPhysicalPars(nlpars);
};
CollectionMassModels::~CollectionMassModels(){
  for(int i=0;i<this->models.size();i++){
    delete(this->models[i]);
  }
  mpars.clear();
};

void CollectionMassModels::setPhysicalPars(std::vector<Nlpar*> nlpars){
  for(int i=0;i<nlpars.size();i++){
    this->mpars[nlpars[i]->nam] = nlpars[i]->val;
  }
  this->mpars["phi"] *= 0.01745329251;
  this->mpars["g1"]  = this->mpars["g"]*cos(2*this->mpars["phi"]);
  this->mpars["g2"]  = this->mpars["g"]*sin(2*this->mpars["phi"]);
}

void CollectionMassModels::all_defl(double xin,double yin,double& xout,double& yout){
  double ax   = 0.0;
  double ay   = 0.0;
  double dumx = 0.0;
  double dumy = 0.0;
  for(int i=0;i<this->models.size();i++){
    this->models[i]->defl(xin,yin,dumx,dumy);
    ax += dumx;
    ay += dumy;
  }
  xout = (1.0-this->mpars["g1"])*xin - this->mpars["g2"]*yin - ax;
  yout = (1.0+this->mpars["g1"])*yin - this->mpars["g2"]*xin - ay;
}

void CollectionMassModels::all_defl(ImagePlane* image){
  double dumx = 0.0;
  double dumy = 0.0;
  double xin,yin;
  
  for(int j=0;j<image->Nm;j++){
    xin = image->x[j];
    yin = image->y[j];
    double ax   = 0.0;
    double ay   = 0.0;
    for(int i=0;i<this->models.size();i++){
      this->models[i]->defl(xin,yin,dumx,dumy);
      ax += dumx;
      ay += dumy;
    }
    image->defl_x[j] = (1.0-this->mpars["g1"])*xin - this->mpars["g2"]*yin - ax;
    image->defl_y[j] = (1.0+this->mpars["g1"])*yin - this->mpars["g2"]*xin - ay;
  }
}
