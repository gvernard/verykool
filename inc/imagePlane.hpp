#ifndef IMAGE_PLANE_HPP
#define IMAGE_PLANE_HPP

#include <valarray>
#include <string>
#include <vector>
#include <map>

#include "tableDefinition.hpp"


class ImagePlane {
public:
  int Ni;                    //pixels in x direction
  int Nj;                    //pixels in y direction
  int Nm;                    //total pixels in the image data
  std::map<int,int> lookup;  //matching the indices of the data pixels to the indices of the mask pixels
  double width;              //in arcsec
  double height;             //in arcsec
  double* img;               //values of the pixels
  double* x;                 //pixel x-coordinates in arcsec
  double* y;                 //pixel y-coordinates in arcsec
  int* active;               //active image pixels used in the construction of the adaptive grid
  mytable B;
  mytable C;
  mytable S;

  ImagePlane(){};
  ImagePlane(const std::string filepath,int i,int j,double w,double h);  //used to read images
  ImagePlane(int i,int j,double w,double h);                             //used to create images
  ImagePlane(int i,double w,double h);                                   //used to create masked images
  ImagePlane(const ImagePlane& image);
  ~ImagePlane();

  void writeImage(const std::string filename);
  void writeBin(const std::string filename);
  void readFits(const std::string filename,std::valarray<float>& contents);
  void readB(const std::string filepath,int i,int j,int ci,int cj);
  void readC(const std::string flag,const std::string filepath);
  void readS(const std::string filepath);
  //  void setMaskedC(mytable* Cout,mytable* S,mytable* C);
  void maskData(std::map<int,int> lookup,ImagePlane* masked);
};

#endif /* IMAGE_PLANE_HPP */
