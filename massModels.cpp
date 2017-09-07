#include "massModels.hpp"

//Abstract class: BaseMassModel
//===============================================================================================================
void BaseMassModel::setMassPars(std::map<std::string,BaseNlpar*> nlpars){
  typedef std::map<std::string,BaseNlpar*>::iterator it_type;
  for(it_type iterator=nlpars.begin();iterator!=nlpars.end();iterator++){
    this->mpars[iterator->first] = iterator->second->val;
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
Sie::Sie(std::map<std::string,BaseNlpar*> nlpars){
  this->n = 5;
  setMassPars(nlpars);
}

void Sie::defl(double xin,double yin,double& xout,double& yout){
  double b  = this->mpars["b"];
  double q  = this->mpars["q"];
  double pa = this->mpars["pa"] * 0.01745329251;//in rad
  double x0 = this->mpars["x0"];
  double y0 = this->mpars["y0"];;

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

  double fac   = 1.-q*q;
  double omega = q*q*x_t*x_t + y_t*y_t; // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.
  double fac2  = sqrt(fac/omega);
  
  double ax_t = (b/sqrt(fac))*atan(x_t*fac2);
  double ay_t = (b/sqrt(fac))*atanh(y_t*fac2);
  
  //rotate back according to position angle, no need to translate
  double ax =  ax_t*cos(-pa) + ay_t*sin(-pa);
  double ay = -ax_t*sin(-pa) + ay_t*cos(-pa);
  
  //  xout = xin - ax;
  //  yout = yin - ay;
  xout = ax;
  yout = ay;
}

/*
//Derived class from BaseMassModel: Sie_g (Singular Isothermal Ellipsoid + external shear)
//===============================================================================================================
Sie_g::Sie_g(std::map<std::string,BaseNlpar*> nlpars){
  this->n = 7;
  this->mpars["b"]   = nlpars["b"]->val;
  this->mpars["q"]   = nlpars["q"]->val;
  this->mpars["pa"]  = nlpars["pa"]->val;
  this->mpars["x0"]  = nlpars["x0"]->val;
  this->mpars["y0"]  = nlpars["y0"]->val;
  this->mpars["g"]   = nlpars["g"]->val;
  this->mpars["phi"] = nlpars["phi"]->val;
}

void Sie_g::defl(double xin,double yin,double& xout,double& yout){
  double b   = this->mpars["b"];
  double q   = this->mpars["q"];
  double pa  = this->mpars["pa"] * 0.01745329251;//in rad
  double x0  = this->mpars["x0"];
  double y0  = this->mpars["y0"];
  double g   = this->mpars["g"];
  double phi = this->mpars["phi"] * 0.01745329251;//in rad

  //rotate according to position angle and translate to the lens center
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

  double fac   = sqrt(1.-q*q);
  double omega = sqrt(q*q*x_t*x_t + y_t*y_t);
  double fac2  = fac/omega;
  
  double ax_t = (b/fac)*atan(x_t*fac2);
  double ay_t = (b/fac)*atanh(y_t*fac2);
  
  //rotate back according to position angle, no need to translate
  double ax =  ax_t*cos(-pa) + ay_t*sin(-pa);
  double ay = -ax_t*sin(-pa) + ay_t*cos(-pa);
  
  //add deflection from external shear
  double g1 = g*cos(2*phi);
  double g2 = g*sin(2*phi);
  
  xout = (1.-g1)*xin + g2*yin - ax;
  yout = g2*xin + (1.+g1)*yin - ay;
}
*/


//Derived class from BaseMassModel: Spemd (Softened Power-law Elliptical Mass Density)
//===============================================================================================================
Spemd::Spemd(std::map<std::string,BaseNlpar*> nlpars){
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
  double ax =  ax_t*cos(-pa) + ay_t*sin(-pa);
  double ay = -ax_t*sin(-pa) + ay_t*cos(-pa);
  
  xout = xin - ax;
  yout = yin - ay;
}

/*
//Derived class from BaseMassModel: Spemd_g (Softened Power-law Elliptical Mass Density + external shear)
//===============================================================================================================
Spemd_g::Spemd_g(std::map<std::string,BaseNlpar*> nlpars){
  this->n = 9;
  this->mpars["b"]   = nlpars["b"]->val;
  this->mpars["q"]   = nlpars["q"]->val;
  this->mpars["pa"]  = nlpars["pa"]->val;
  this->mpars["x0"]  = nlpars["x0"]->val;
  this->mpars["y0"]  = nlpars["y0"]->val;
  this->mpars["e"]   = nlpars["e"]->val;
  this->mpars["s"]   = nlpars["s"]->val;
  this->mpars["g"]   = nlpars["g"]->val;
  this->mpars["phi"] = nlpars["phi"]->val;
}

void Spemd_g::defl(double xin,double yin,double& xout,double& yout){
  double b   = this->mpars["b"];
  double q   = this->mpars["q"];
  double pa  = this->mpars["pa"] * 0.01745329251;//in rad
  double x0  = this->mpars["x0"];
  double y0  = this->mpars["y0"];
  double e   = this->mpars["e"];
  double s2  = this->mpars["s"] * this->mpars["s"];
  double g   = this->mpars["g"];
  double phi = this->mpars["phi"] * 0.01745329251;//in rad
  double defl[2] = {0.,0.};

  //rotate according to position angle and translate to the lens center
  double x_t =  (xin-x0)*cos(pa) + (yin-y0)*sin(pa);
  double y_t = -(xin-x0)*sin(pa) + (yin-y0)*cos(pa);

  fastelldefl_(&x_t,&y_t,&b,&e,&q,&s2,defl);

  double ax_t = defl[0];
  double ay_t = defl[1];

  //rotate back according to position angle, no need to translate
  double ax =  ax_t*cos(-pa) + ay_t*sin(-pa);
  double ay = -ax_t*sin(-pa) + ay_t*cos(-pa);
  
  //add deflection from external shear
  double g1 = g*cos(2*phi);
  double g2 = g*sin(2*phi);

  xout = (1.-g1)*xin + g2*yin - ax;
  yout = g2*xin + (1.+g1)*yin - ay;
}
*/
