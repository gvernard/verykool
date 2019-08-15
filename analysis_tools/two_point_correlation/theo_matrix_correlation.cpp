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
  std::string src_path  = argv[1]; // The reconstructed source is needed only to read in the adaptive grid.
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


  // Getting the covariance and inverse covariance matrices
  BaseAlgebra* algebra = new BaseAlgebra();


  s_ada_curv->createDelaunay();
  s_ada_curv->constructH();
  Eigen::SparseMatrix<double> H(s_ada_curv->H.Ti,s_ada_curv->H.Tj);
  H.reserve(Eigen::VectorXi::Constant(s_ada_curv->Sm,8));
  for(int i=0;i<s_ada_curv->H.tri.size();i++){
    H.insert(s_ada_curv->H.tri[i].i,s_ada_curv->H.tri[i].j) = s_ada_curv->H.tri[i].v;
  }
  Eigen::SparseMatrix<double> C_curv_inv(s_ada_curv->Sm,s_ada_curv->Sm); // HtH
  C_curv_inv = (H.transpose() * H);
  H.resize(0,0);
  Eigen::SparseMatrix<double> C_curv(s_ada_curv->Sm,s_ada_curv->Sm); // HtH_inv
  double dum_detHtH = 0.0;
  algebra->getInverseAndDet(C_curv_inv,C_curv,dum_detHtH);
  Eigen::SparseMatrix<double> curv_ide(s_ada_modgauss->Sm,s_ada_modgauss->Sm);
  curv_ide = C_curv_inv*C_curv;


  /*
  // Replace the source pixel values with their variance for curvature
  for(int i=0;i<s_ada_curv->Sm;i++){
    s_ada_curv->mask_vertices[i] = 1;
  }
  s_ada_curv->outputSource("");

  for(int k=0;k<HtH.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(HtH,k);it;++it){
      if( it.row() == it.col() ){
	s_ada_curv->src[it.row()] = it.value();
      }
    }
  }
  s_ada_curv->outputSource("curv_variance");
  */





  // Modgauss
  s_ada_modgauss->createDelaunay();
  s_ada_modgauss->constructH();

  Eigen::SparseMatrix<double> C_modgauss(s_ada_modgauss->H.Ti,s_ada_modgauss->H.Tj);
  C_modgauss.reserve(Eigen::VectorXi::Constant(s_ada_modgauss->Sm,100));
  for(int i=0;i<s_ada_modgauss->H.tri.size();i++){
    C_modgauss.insert(s_ada_modgauss->H.tri[i].i,s_ada_modgauss->H.tri[i].j) = s_ada_modgauss->H.tri[i].v;
  }
  Eigen::SparseMatrix<double> C_modgauss_inv(s_ada_modgauss->Sm,s_ada_modgauss->Sm);
  double dum_detCmodgauss = 0.0;
  algebra->getInverseAndDet(C_modgauss,C_modgauss_inv,dum_detCmodgauss);
  Eigen::SparseMatrix<double> modgauss_ide(s_ada_modgauss->Sm,s_ada_modgauss->Sm);
  modgauss_ide = C_modgauss_inv*C_modgauss;




  // Gauss
  s_ada_gauss->createDelaunay();
  s_ada_gauss->constructH();
  
  Eigen::SparseMatrix<double> C_gauss(s_ada_gauss->H.Ti,s_ada_gauss->H.Tj);
  C_gauss.reserve(Eigen::VectorXi::Constant(s_ada_gauss->Sm,100));
  for(int i=0;i<s_ada_gauss->H.tri.size();i++){
    C_gauss.insert(s_ada_gauss->H.tri[i].i,s_ada_gauss->H.tri[i].j) = s_ada_gauss->H.tri[i].v;
  }
  Eigen::SparseMatrix<double> C_gauss_inv(s_ada_gauss->Sm,s_ada_gauss->Sm);
  double dum_detCgauss = 0.0;
  algebra->getInverseAndDet(C_gauss,C_gauss_inv,dum_detCgauss);
  Eigen::SparseMatrix<double> gauss_ide(s_ada_modgauss->Sm,s_ada_modgauss->Sm);
  gauss_ide = C_gauss_inv*C_gauss;
  


  

  // Curvature
  ImagePlane* matrix_curv_inv = new ImagePlane(s_ada_curv->Sm,s_ada_curv->Sm,1.0,1.0);
  for(int k=0;k<C_curv_inv.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(C_curv_inv,k);it;++it){
      int nx = it.row();
      int ny = it.col();
      matrix_curv_inv->img[nx*s_ada_curv->Sm + ny] = it.value();
    }
  }
  matrix_curv_inv->writeImage("inv_matrix_curv.fits");
  delete(matrix_curv_inv);

  ImagePlane* matrix_curv = new ImagePlane(s_ada_curv->Sm,s_ada_curv->Sm,1.0,1.0);
  for(int k=0;k<C_curv.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(C_curv,k);it;++it){
      int nx = it.row();
      int ny = it.col();
      matrix_curv->img[nx*s_ada_curv->Sm + ny] = it.value();
    }
  }
  matrix_curv->writeImage("matrix_curv.fits");
  delete(matrix_curv);

  ImagePlane* matrix_curv_ide = new ImagePlane(s_ada_curv->Sm,s_ada_curv->Sm,1.0,1.0);
  for(int k=0;k<curv_ide.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(curv_ide,k);it;++it){
      int nx = it.row();
      int ny = it.col();
      matrix_curv_ide->img[nx*s_ada_curv->Sm + ny] = it.value();
    }
  }
  matrix_curv_ide->writeImage("matrix_curv_ide.fits");
  delete(matrix_curv_ide);




  // Modgauss
  ImagePlane* matrix_modgauss = new ImagePlane(s_ada_modgauss->Sm,s_ada_modgauss->Sm,1.0,1.0);
  for(int k=0;k<C_modgauss.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(C_modgauss,k);it;++it){
      int nx = it.row();
      int ny = it.col();
      matrix_modgauss->img[nx*s_ada_modgauss->Sm + ny] = it.value();
    }
  }
  matrix_modgauss->writeImage("matrix_modgauss.fits");
  delete(matrix_modgauss);

  ImagePlane* matrix_modgauss_inv = new ImagePlane(s_ada_modgauss->Sm,s_ada_modgauss->Sm,1.0,1.0);
  for(int k=0;k<C_modgauss_inv.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(C_modgauss_inv,k);it;++it){
      int nx = it.row();
      int ny = it.col();
      matrix_modgauss_inv->img[nx*s_ada_modgauss->Sm + ny] = it.value();
    }
  }
  matrix_modgauss_inv->writeImage("inv_matrix_modgauss.fits");
  delete(matrix_modgauss_inv);

  ImagePlane* matrix_modgauss_ide = new ImagePlane(s_ada_curv->Sm,s_ada_curv->Sm,1.0,1.0);
  for(int k=0;k<modgauss_ide.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(modgauss_ide,k);it;++it){
      int nx = it.row();
      int ny = it.col();
      matrix_modgauss_ide->img[nx*s_ada_curv->Sm + ny] = it.value();
    }
  }
  matrix_modgauss_ide->writeImage("matrix_modgauss_ide.fits");
  delete(matrix_modgauss_ide);




  // Modgauss
  ImagePlane* matrix_gauss = new ImagePlane(s_ada_gauss->Sm,s_ada_gauss->Sm,1.0,1.0);
  for(int k=0;k<C_gauss.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(C_gauss,k);it;++it){
      int nx = it.row();
      int ny = it.col();
      matrix_gauss->img[nx*s_ada_gauss->Sm + ny] = it.value();
    }
  }
  matrix_gauss->writeImage("matrix_gauss.fits");
  delete(matrix_gauss);

  ImagePlane* matrix_gauss_inv = new ImagePlane(s_ada_gauss->Sm,s_ada_gauss->Sm,1.0,1.0);
  for(int k=0;k<C_gauss_inv.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(C_gauss_inv,k);it;++it){
      int nx = it.row();
      int ny = it.col();
      matrix_gauss_inv->img[nx*s_ada_gauss->Sm + ny] = it.value();
    }
  }
  matrix_gauss_inv->writeImage("inv_matrix_gauss.fits");
  delete(matrix_gauss_inv);

  ImagePlane* matrix_gauss_ide = new ImagePlane(s_ada_curv->Sm,s_ada_curv->Sm,1.0,1.0);
  for(int k=0;k<gauss_ide.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(gauss_ide,k);it;++it){
      int nx = it.row();
      int ny = it.col();
      matrix_gauss_ide->img[nx*s_ada_curv->Sm + ny] = it.value();
    }
  }
  matrix_gauss_ide->writeImage("matrix_gauss_ide.fits");
  delete(matrix_gauss_ide);









  double Nbins = 100;
  double rmax = 2.0; // arcsec
  double rmin = 0.0; // arcsec
  double dr = (rmax-rmin)/Nbins;
  double* bins = (double*) calloc(Nbins,sizeof(double));
  for(int i=0;i<Nbins;i++){
    bins[i] = rmin + dr*i;
  }
  
  // Binning by r the covariance matrix for curvature
  double* vals_curv = (double*) calloc(Nbins,sizeof(double));
  int* counts_curv  = (int*) calloc(Nbins,sizeof(int));
  for(int k=0;k<C_curv.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(C_curv,k);it;++it){
      int ia = it.row();
      int ib = it.col();
      double r = hypot(s_ada_curv->x[ia]-s_ada_curv->x[ib],s_ada_curv->y[ia]-s_ada_curv->y[ib]);
      if( r < rmax && ia != ib ){
	int index = (int) floor(r/dr);
	vals_curv[index] += it.value();
	counts_curv[index]++;
      }
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
  for(int k=0;k<C_modgauss.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(C_modgauss,k);it;++it){
      int ia = it.row();
      int ib = it.col();
      double r = hypot(s_ada_modgauss->x[ia]-s_ada_modgauss->x[ib],s_ada_modgauss->y[ia]-s_ada_modgauss->y[ib]);
      if( r < rmax && ia != ib ){
	int index = (int) floor(r/dr);
	vals_modgauss[index] += it.value();
	counts_modgauss[index]++;
      }
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
  for(int k=0;k<C_gauss.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(C_gauss,k);it;++it){
      int ia = it.row();
      int ib = it.col();
      double r = hypot(s_ada_gauss->x[ia]-s_ada_gauss->x[ib],s_ada_gauss->y[ia]-s_ada_gauss->y[ib]);
      if( r < rmax && ia != ib ){
	int index = (int) floor(r/dr);
	vals_gauss[index] += it.value();
	counts_gauss[index]++;
      }
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

  /*  
  // Binning by r the covariance matrix of the true source
  ImagePlane* s_true = new ImagePlane(src_true,700,700,1.0,1.0);
  double* vals_s_true = (double*) calloc(Nbins,sizeof(double));
  int* counts_s_true  = (int*) calloc(Nbins,sizeof(int));
  for(int k=0;k<s_true->Nm;k++){
    for(int q=k+1;q<s_true->Nm;q++){
      double r = hypot(s_true->x[q]-s_true->x[k],s_true->y[q]-s_true->y[k]);
      if( r < rmax ){
  	int index = (int) floor(r/dr);
  	vals_s_true[index] += s_true->img[q]*s_true->img[k];
  	counts_s_true[index]++; 
      }  
    }
  }
  for(int i=0;i<Nbins;i++){
    if( counts_s_true[i] > 0 ){
      vals_s_true[i] /= counts_s_true[i];
    }
  }
  */  



  
  FILE* fh = fopen("corr_theo.dat","w");
  for(int i=0;i<Nbins;i++){
    fprintf(fh,"%5.3f %8.3f %8.3f %8.3f\n",bins[i],vals_curv[i],vals_modgauss[i],vals_gauss[i]);
    //fprintf(fh,"%5.3f %8.3f\n",bins[i],vals_curv[i]);
  }
  fclose(fh);
  


  C_curv.resize(0,0);
  C_curv_inv.resize(0,0);
  C_modgauss.resize(0,0);
  C_modgauss_inv.resize(0,0);
  C_gauss.resize(0,0);
  C_gauss_inv.resize(0,0);



  
  free(bins);
  free(vals_curv);
  free(vals_modgauss);
  free(vals_gauss);


  //  free(vals_s_true);
  //  free(counts_s_true);

  
  delete(delaunay);
  delete(s_ada_curv);
  delete(s_ada_modgauss);
  delete(s_ada_gauss);
  //  delete(s_true);

  return 0;
}
