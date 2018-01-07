#ifndef SOURCE_PLANE_HPP
#define SOURCE_PLANE_HPP

#include <string>
#include <vector>
#include <cstdlib>
#include <map>

/*
#include <boost/geometry.hpp> 
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/algorithms/append.hpp> 
BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)
*/

#include "covKernels.hpp"


typedef std::map<int,double>::iterator it_int_double;

struct mytable;
class Nlpar;
class ImagePlane;
class CollectionMassModels;


class BaseSourcePlane {
public:
  std::string type;               // adaptive or regular
  int Si;                         // pixels in x, if adaptive then Si is set to 0
  int Sj;                         // pixels in y, if adaptive then Sj is set to 0
  int Sm;                         // total number of pixels
  double* src;                    // values of the pixels
  double* x;                      // pixel x-coordinates
  double* y;                      // pixel y-coordinates
  std::vector<double> bound_x;    // embedding polygon vertex x-coordinates
  std::vector<double> bound_y;    // embedding polygon vertex y-coordinates
  std::string reg;                // name of the regularization scheme
  int eigenSparseMemoryAllocForH; // estimate of the non-zero elements per row of the regularization matrix H
  bool sample_reg = false;        // sampling regularization matrix related parameters
  BaseCovKernel* kernel;          // pointer to kernel class
  
  BaseSourcePlane(){};

  ~BaseSourcePlane(){
    free(src);
    free(x);
    free(y);
    if( this->reg == "covariance_kernel" ){
      delete this->kernel;
    }
  };


  //virtual members
  virtual void constructL(ImagePlane* image,CollectionMassModels* mycollection,mytable* L) = 0;
  virtual void constructH(mytable* H) = 0;
  virtual void outputSource(const std::string filename) = 0;
  virtual void outputSourceErrors(double* errors,const std::string filename) = 0;


  //non-virtual members
  void normalize();
  void setSourceCovariance(std::vector<Nlpar> reg_pars){};
};


class FixedSource: public BaseSourcePlane {
public:

  //virtual members
  FixedSource(int source_i,int source_j,double size,std::string reg_scheme);
  ~FixedSource();
  void constructL(ImagePlane* image,CollectionMassModels* mycollection,mytable* L);
  void constructH(mytable* H);
  void outputSource(const std::string filename);
  void outputSourceErrors(double* errors,const std::string filename);

  //non-virtual members
  void setGrid(std::map<std::string,std::string> pars);
  void boundPolygon();
  bool pointInPolygon(double x,double y);

private:
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};

class FloatingSource: public BaseSourcePlane {
public:

  //virtual members
  FloatingSource(int source_i,int source_j,double size,double x0,double y0,std::string reg_scheme);
  ~FloatingSource();
  void constructL(ImagePlane* image,CollectionMassModels* mycollection,mytable* L);
  void constructH(mytable* H);
  void outputSource(const std::string filename);
  void outputSourceErrors(double* errors,const std::string filename);

  //non-virtual members
  void setGrid(std::map<std::string,std::string> pars);
  void boundPolygon();
  bool pointInPolygon(double x,double y);

private:
  double x0;
  double y0;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};

class AdaptiveSource: public BaseSourcePlane {
public:

  //virtual members
  AdaptiveSource(int Sm,std::string reg_scheme);
  AdaptiveSource(std::string mode,int Sm,int spacing,std::string reg_scheme);
  ~AdaptiveSource();
  void constructL(ImagePlane* image,CollectionMassModels* mycollection,mytable* L);
  void constructH(mytable* H);
  void outputSource(const std::string filename);
  void outputSourceErrors(double* errors,const std::string filename);

  //non-virtual members
  void createAdaGrid(ImagePlane* image,CollectionMassModels* mycollection);
  void createDelaunay();
  void writeTriangles();
  void setGrid(std::map<std::string,std::string> pars){};
  bool pointInPolygon(double x,double y){};
  //  void boundPolygon();

private:

  struct xypoint {
    double x;
    double y;
  };
  
  struct a_triangle {
    int a;
    int b;
    int c;
  };

  xypoint intersection_point_x(xypoint p0,xypoint p1,xypoint p2);
  xypoint intersection_point_y(xypoint p0,xypoint p1,xypoint p2);

  std::vector<int> used;
  std::string mode;
  int spacing;
  int n_triangles;
  std::vector<a_triangle> triangles;
  std::vector< std::vector<int> > opposite_edges_per_vertex;
};



class FactorySourcePlane{//This is a singleton class.
public:
  FactorySourcePlane(FactorySourcePlane const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactorySourcePlane const&) = delete;

  static FactorySourcePlane* getInstance(){
    static FactorySourcePlane dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseSourcePlane* createSourcePlane(std::map<std::string,std::string> source){
    if( source["type"] == "fixed" ){
      return new FixedSource(stoi(source["sx"]),stoi(source["sy"]),stof(source["size"]),source["reg"]);
    } else if( source["type"] == "floating" ){
      return new FloatingSource(stoi(source["sx"]),stoi(source["sy"]),stof(source["size"]),stof(source["x0"]),stof(source["y0"]),source["reg"]);
    } else if( source["type"] == "adaptive" ){
      if( source["mode"] == "random" ){
	return new AdaptiveSource(stoi(source["sm"]),source["reg"]);
      } else if( source["mode"] == "image" || source["mode"] == "grid" ){
	return new AdaptiveSource(source["mode"],stoi(source["sm"]),stoi(source["spacing"]),source["reg"]);
      }
    } else {
      return NULL;
    }
  }

private:
  FactorySourcePlane(){};
};


#endif /* SOURCE_PLANE_HPP */
