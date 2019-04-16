#ifndef SOURCE_PROFILE_HPP
#define SOURCE_PROFILE_HPP

#include <cmath>
#include <string>
#include <map>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2_algorithms.h>

#include "nonLinearPars.hpp"
#include "imagePlane.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

class BaseProfile {
public:
  std::map<std::string,double> pars;
  int npoly = 100;
  double* x;
  double* y;
  const double pi  = 3.14159265358979323846;
  const double fac = 0.01745329251;
  
  BaseProfile(){};
  ~BaseProfile(){
    free(x);
    free(y);
    pars.clear();
  }

  virtual double value(double x,double y) = 0;

  void profile(int Sj,int Si,double* sx,double* sy,double* s);
  void ellipticalContour();
};

  
class Sersic: public BaseProfile {
public:
  double index;

  Sersic(std::map<std::string,Nlpar*> nlpars);
  double value(double x,double y);
};


class ProGauss: public BaseProfile {
public:
  ProGauss(std::map<std::string,Nlpar*> nlpars);
  double value(double x,double y);
};


class myVoronoi: public BaseProfile {
public:
  int Ncells;
  std::vector<int> sizes;
  std::vector<Point*> cells;
  std::vector<double> values;
  
  myVoronoi(std::string filename);
  ~myVoronoi(){
    for(int i=0;i<Ncells;i++){
      delete(cells[i]);
    }
  }
  double value(double x,double y);
};


class myDelaunay: public BaseProfile {
public:
  int N;
  double* src;
  
  myDelaunay(std::string filename);
  ~myDelaunay(){
    free(src);
    free(convex_hull);
  }
  double value(double x,double y);

private:
  struct atriangle {
    int a;
    int b;
    int c;
  };
  std::vector<atriangle> triangles;
  Point* convex_hull;
  int ch_size;
};


class fromFITS: public BaseProfile {
public:
  fromFITS(std::string filename,int Ni,int Nj,double height,double width,double x0,double y0);
  ~fromFITS(){
    delete(mySource);
  }
  double value(double x,double y);

private:
  ImagePlane* mySource;
};



class FactoryProfile{//This is a singleton class.
public:
  FactoryProfile(FactoryProfile const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryProfile const&) = delete;

  static FactoryProfile* getInstance(){
    static FactoryProfile dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseProfile* createProfile(const std::string &modelname,std::map<std::string,Nlpar*> nlpars){
    if( modelname == "sersic" ){
      return new Sersic(nlpars);
    } else if ( modelname == "gauss" ){
      return new ProGauss(nlpars);
    } else {
      return NULL;
    }
  }

  BaseProfile* createProfile(const std::string &modelname,std::map<std::string,std::string> pars){
    if( modelname == "voronoi" ){
      return new myVoronoi(pars["filename"]);
    } else if( modelname == "delaunay" ){
      return new myDelaunay(pars["filename"]);
    } else if( modelname == "fromfits" ){
      return new fromFITS(pars["filename"],std::stoi(pars["Ni"]),std::stoi(pars["Nj"]),std::stof(pars["height"]),std::stof(pars["width"]),std::stof(pars["x0"]),std::stof(pars["y0"]));
    } else {
      return NULL;
    }
  }

private:
  FactoryProfile(){};
};

#endif /* SOURCE_PROFILE_HPP */
