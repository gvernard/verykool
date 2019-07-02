#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "json/json.h"

#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "sourceProfile.hpp"
#include "nonLinearPars.hpp"
#include "covKernels.hpp"
#include "likelihoodModels.hpp"
#include "eigenAlgebra.hpp"
#include "tableDefinition.hpp"


int main(int argc,char* argv[]){
  std::string src_path  = argv[1];
  std::string json_file = argv[2];
  //  std::string src_true  = argv[3];


  
  Json::Value json_input;
  std::ifstream fin((json_file).c_str());
  fin >> json_input;



  // Read in a Delaunay source
  myDelaunay* delaunay = new myDelaunay(src_path);



  // Define the different adaptive sources and their regularization parameters/schemes
  AdaptiveSource* s_ada_curv = new AdaptiveSource(delaunay->N,"curvature");


  AdaptiveSource* s_ada_modgauss  = new AdaptiveSource(delaunay->N,"covariance_kernel");
  std::vector<Nlpar*> modgauss_nlpars = FactoryLikelihoodModel::getInstance()->nlparsFromJsonVector(json_input["modgauss"]);
  BaseCovKernel* modgauss_kernel = FactoryCovKernel::getInstance()->createCovKernel("modgauss",modgauss_nlpars);
  s_ada_modgauss->kernel = modgauss_kernel;
  for(int i=0;i<modgauss_nlpars.size();i++){
    delete(modgauss_nlpars[i]);
  }
  
  AdaptiveSource* s_ada_gauss  = new AdaptiveSource(delaunay->N,"covariance_kernel");
  std::vector<Nlpar*> gauss_nlpars = FactoryLikelihoodModel::getInstance()->nlparsFromJsonVector(json_input["gauss"]);
  BaseCovKernel* gauss_kernel = FactoryCovKernel::getInstance()->createCovKernel("gauss",gauss_nlpars);
  s_ada_gauss->kernel = gauss_kernel;
  for(int i=0;i<gauss_nlpars.size();i++){
    delete(gauss_nlpars[i]);
  }
  


  // Set the parameters of the various sources to the input one
  for(int i=0;i<delaunay->N;i++){
    s_ada_curv->x[i]     = delaunay->x[i];
    s_ada_curv->y[i]     = delaunay->y[i];
    s_ada_curv->src[i]   = delaunay->src[i];

    s_ada_modgauss->x[i]   = delaunay->x[i];
    s_ada_modgauss->y[i]   = delaunay->y[i];
    s_ada_modgauss->src[i] = delaunay->src[i];

    s_ada_gauss->x[i]   = delaunay->x[i];
    s_ada_gauss->y[i]   = delaunay->y[i];
    s_ada_gauss->src[i] = delaunay->src[i];
  }


  // Getting the curvature covariance matrix, i.e. the inverse of HtH
  s_ada_curv->createDelaunay();
  s_ada_curv->constructH();
  Eigen::SparseMatrix<double> H(s_ada_curv->H.Ti,s_ada_curv->H.Tj);
  H.reserve(Eigen::VectorXi::Constant(s_ada_curv->Sm,8));
  for(int i=0;i<s_ada_curv->H.tri.size();i++){
    H.insert(s_ada_curv->H.tri[i].i,s_ada_curv->H.tri[i].j) = s_ada_curv->H.tri[i].v;
  }
  Eigen::SparseMatrix<double> HtH(s_ada_curv->Sm,s_ada_curv->Sm);
  HtH = (H.transpose() * H);
  Eigen::SparseMatrix<double> HtHinv(s_ada_curv->Sm,s_ada_curv->Sm);
  double detHtH = 0.0;
  BaseAlgebra* algebra = new BaseAlgebra();
  algebra->getInverseAndDet(HtH,HtHinv,detHtH);
  H.resize(0,0);
  HtH.resize(0,0);
  std::vector<mytriplet> HtHinv_tri;
  for(int k=0;k<HtHinv.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(HtHinv,k);it;++it){
      HtHinv_tri.push_back({it.row(),it.col(),it.value()});
      //HtHinv_tri.push_back({it.col(),it.row(),it.value()});
    }
  }
  HtHinv.resize(0,0);

  s_ada_modgauss->createDelaunay();
  s_ada_modgauss->constructH();
  
  s_ada_gauss->createDelaunay();
  s_ada_gauss->constructH();


  
  // Write the regularization matrix of the source as an image
  /*
  ImagePlane* matrix_curv = new ImagePlane(s_ada_curv->Sm,s_ada_curv->Sm,1.0,1.0);
  for(int i=0;i<s_ada_curv->H.tri.size();i++){
    int nx = s_ada_curv->H.tri[i].i;
    int ny = s_ada_curv->H.tri[i].j;
    matrix_curv->img[nx*s_ada_curv->Sm + ny] = s_ada_curv->H.tri[i].v;
  }
  matrix_curv->writeImage("matrix_H.fits");
  delete(matrix_curv);
  */

  ImagePlane* matrix_HtHinv = new ImagePlane(s_ada_curv->Sm,s_ada_curv->Sm,1.0,1.0);
  for(int i=0;i<HtHinv_tri.size();i++){
    int nx = HtHinv_tri[i].i;
    int ny = HtHinv_tri[i].j;
    matrix_HtHinv->img[nx*s_ada_curv->Sm + ny] = HtHinv_tri[i].v;
  }
  matrix_HtHinv->writeImage("matrix_HtHinv.fits");
  delete(matrix_HtHinv);

  ImagePlane* matrix_modgauss = new ImagePlane(s_ada_modgauss->Sm,s_ada_modgauss->Sm,1.0,1.0);
  for(int i=0;i<s_ada_modgauss->H.tri.size();i++){
    int nx = s_ada_modgauss->H.tri[i].i;
    int ny = s_ada_modgauss->H.tri[i].j;
    matrix_modgauss->img[nx*s_ada_modgauss->Sm + ny] = s_ada_modgauss->H.tri[i].v;
  }
  matrix_modgauss->writeImage("matrix_modgauss.fits");
  delete(matrix_modgauss);

  ImagePlane* matrix_gauss = new ImagePlane(s_ada_gauss->Sm,s_ada_gauss->Sm,1.0,1.0);
  for(int i=0;i<s_ada_gauss->H.tri.size();i++){
    int nx = s_ada_gauss->H.tri[i].i;
    int ny = s_ada_gauss->H.tri[i].j;
    matrix_gauss->img[nx*s_ada_gauss->Sm + ny] = s_ada_gauss->H.tri[i].v;
  }
  matrix_gauss->writeImage("matrix_gauss.fits");
  delete(matrix_gauss);




  double Nbins = 50;
  double rmax = 1.0; // arcsec
  double rmin = 0.0; // arcsec
  double dr = (rmax-rmin)/Nbins;
  double* bins = (double*) calloc(Nbins,sizeof(double));
  for(int i=0;i<Nbins;i++){
    bins[i] = rmin + dr*i;
  }

  // Binning by r the covariance matrix for curvature
  double* vals_curv = (double*) calloc(Nbins,sizeof(double));
  int* counts_curv  = (int*) calloc(Nbins,sizeof(int));
  for(int k=0;k<HtHinv_tri.size();k++){
    int ia = HtHinv_tri[k].i;
    int ib = HtHinv_tri[k].j;
    double r = hypot(s_ada_curv->x[ia]-s_ada_curv->x[ib],s_ada_curv->y[ia]-s_ada_curv->y[ib]);
    if( r < rmax ){
      int index = (int) floor(r/dr);
      vals_curv[index] += HtHinv_tri[k].v;
      counts_curv[index]++;
    }
  }
  for(int i=0;i<Nbins;i++){
    if( counts_curv[i] > 0 ){
      vals_curv[i] /= counts_curv[i];
    }
  }
  free(counts_curv);

  // Binning by r the covariance matrix for a modgauss kernel
  double* vals_modgauss = (double*) calloc(Nbins,sizeof(double));
  int* counts_modgauss  = (int*) calloc(Nbins,sizeof(int));
  for(int k=0;k<s_ada_modgauss->H.tri.size();k++){
    int ia = s_ada_modgauss->H.tri[k].i;
    int ib = s_ada_modgauss->H.tri[k].j;
    double r = hypot(s_ada_modgauss->x[ia]-s_ada_modgauss->x[ib],s_ada_modgauss->y[ia]-s_ada_modgauss->y[ib]);
    if( r < rmax ){
      int index = (int) floor(r/dr);
      vals_modgauss[index] += s_ada_modgauss->H.tri[k].v;
      counts_modgauss[index]++;
    }
  }
  for(int i=0;i<Nbins;i++){
    if( counts_modgauss[i] > 0 ){
      vals_modgauss[i] /= counts_modgauss[i];
    }
  }
  free(counts_modgauss);


  // Binning by r the covariance matrix for a gauss kernel
  double* vals_gauss = (double*) calloc(Nbins,sizeof(double));
  int* counts_gauss  = (int*) calloc(Nbins,sizeof(int));
  for(int k=0;k<s_ada_gauss->H.tri.size();k++){
    int ia = s_ada_gauss->H.tri[k].i;
    int ib = s_ada_gauss->H.tri[k].j;
    double r = hypot(s_ada_gauss->x[ia]-s_ada_gauss->x[ib],s_ada_gauss->y[ia]-s_ada_gauss->y[ib]);
    if( r < rmax ){
      int index = (int) floor(r/dr);
      vals_gauss[index] += s_ada_gauss->H.tri[k].v;
      counts_gauss[index]++;
    }
  }
  for(int i=0;i<Nbins;i++){
    if( counts_gauss[i] > 0 ){
      vals_gauss[i] /= counts_gauss[i];
    }
  }
  free(counts_gauss);

  /*
  // Binning by r the covariance matrix of the reconstructed source
  double* vals_s = (double*) calloc(Nbins,sizeof(double));
  int* counts_s  = (int*) calloc(Nbins,sizeof(int));
  for(int k=0;k<s_ada_curv->Sm;k++){
    for(int q=k+1;q<s_ada_curv->Sm;q++){
      double r = hypot(s_ada_curv->x[q]-s_ada_curv->x[k],s_ada_curv->y[q]-s_ada_curv->y[k]);
      if( r < rmax ){
	int index = (int) floor(r/dr);
	vals_s[index] += s_ada_curv->src[q]*s_ada_curv->src[k];
	counts_s[index]++; 
      }  
    }
  }
  for(int i=0;i<Nbins;i++){
    if( counts_s[i] > 0 ){
      vals_s[i] /= counts_s[i];
    }
  }
  free(counts_s);
  */

  
  // Binning by r the covariance matrix of the true source
  //  ImagePlane* s_true = new ImagePlane(src_true,700,700,1.0,1.0);
  //  double* vals_s_true = (double*) calloc(Nbins,sizeof(double));
  //  int* counts_s_true  = (int*) calloc(Nbins,sizeof(int));
  //  for(int k=0;k<s_true->Nm;k++){
  //    for(int q=k+1;q<s_true->Nm;q++){
  //      double r = hypot(s_true->x[q]-s_true->x[k],s_true->y[q]-s_true->y[k]);
  //      if( r < rmax ){
  //	int index = (int) floor(r/dr);
  //	vals_s_true[index] += s_true->img[q]*s_true->img[k];
  //	counts_s_true[index]++; 
  //      }  
  //    }
  //  }
  //  for(int i=0;i<Nbins;i++){
  //    if( counts_s_true[i] > 0 ){
  //      vals_s_true[i] /= counts_s_true[i];
  //    }
  //  }
  



  FILE* fh = fopen("corr_theo.dat","w");
  for(int i=0;i<Nbins;i++){
    fprintf(fh,"%5.3f %8.3f %8.3f %8.3f\n",bins[i],vals_curv[i],vals_modgauss[i],vals_gauss[i]);
  }
  fclose(fh);



  free(bins);
  free(vals_curv);
  free(vals_modgauss);
  free(vals_gauss);
  //  free(vals_s);
  //  free(vals_s_true);
  //  free(counts_s_true);

  
  delete(delaunay);
  delete(s_ada_curv);
  delete(s_ada_modgauss);
  delete(s_ada_gauss);
  //  delete(s_true);

  return 0;
}
