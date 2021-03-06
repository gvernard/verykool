#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "json/json.h"

#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "massModels.hpp"
#include "nonLinearPars.hpp"
#include "covKernels.hpp"
#include "eigenAlgebra.hpp"
#include "tableDefinition.hpp"



void writeMatrix(Eigen::SparseMatrix<double> matrix,BaseSourcePlane* src,std::string fname){
  ImagePlane* img_matrix = new ImagePlane(src->Sm,src->Sm,1.0,1.0);
  for(int k=0;k<matrix.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(matrix,k);it;++it){
      int nx = it.row();
      int ny = it.col();
      img_matrix->img[nx*src->Sm + ny] = it.value();
    }
  }
  img_matrix->writeImage(fname);
  delete(img_matrix);
}


void maskMatrix(ImagePlane* mask,Eigen::SparseMatrix<double>& matrix){
  for(int k=0;k<matrix.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(matrix,k);it;++it){
      int nx = it.row();
      int ny = it.col();
      if( mask->img[nx*mask->Nm + ny] == 0 ){
	matrix.coeffRef(nx,ny) = 0.0;
      }
    }
  }
}




int main(int argc,char* argv[]){
  std::string json_file = argv[1];
  //  std::string src_true  = argv[3];


  
  Json::Value json_input;
  std::ifstream fin((json_file).c_str());
  fin >> json_input;

  int Ni = json_input["iplane"]["pix_x"].asInt();
  int Nj = json_input["iplane"]["pix_y"].asInt();
  double w = json_input["iplane"]["siz_x"].asDouble();
  double h = json_input["iplane"]["siz_y"].asDouble();


  

  // Define the different adaptive sources and their regularization parameters/schemes
  std::cout << "Initialization" << std::endl;
  Pert* pert_curv = new Pert(Ni,Nj,w,h,"curvature");


  Pert* pert_modgauss = new Pert(Ni,Nj,w,h,"covariance_kernel");
  std::vector<Nlpar*> modgauss_nlpars = Nlpar::nlparsFromJsonVector(json_input["modgauss"]);
  BaseCovKernel* modgauss_kernel = FactoryCovKernel::getInstance()->createCovKernel("modgauss",modgauss_nlpars);
  pert_modgauss->dpsi->kernel = modgauss_kernel;
  for(int i=0;i<modgauss_nlpars.size();i++){
    delete(modgauss_nlpars[i]);
  }


  Pert* pert_gauss = new Pert(Ni,Nj,w,h,"covariance_kernel");
  std::vector<Nlpar*> gauss_nlpars = Nlpar::nlparsFromJsonVector(json_input["gauss"]);
  BaseCovKernel* gauss_kernel = FactoryCovKernel::getInstance()->createCovKernel("gauss",gauss_nlpars);
  pert_gauss->dpsi->kernel = gauss_kernel;
  for(int i=0;i<gauss_nlpars.size();i++){
    delete(gauss_nlpars[i]);
  }




  

  // Getting the covariance and inverse covariance matrices
  std::cout << "Getting the covariance and inverse covariance matrices" << std::endl;

  BaseAlgebra* algebra = new BaseAlgebra();

  // Curvature
  std::cout << "Curvature...";
  pert_curv->dpsi->constructH();
  Eigen::SparseMatrix<double> H(pert_curv->dpsi->H.Ti,pert_curv->dpsi->H.Tj);
  H.reserve(Eigen::VectorXi::Constant(pert_curv->dpsi->Sm,8));
  for(int i=0;i<pert_curv->dpsi->H.tri.size();i++){
    H.insert(pert_curv->dpsi->H.tri[i].i,pert_curv->dpsi->H.tri[i].j) = pert_curv->dpsi->H.tri[i].v;
  }
  Eigen::SparseMatrix<double> C_curv_inv(pert_curv->dpsi->Sm,pert_curv->dpsi->Sm); // HtH
  C_curv_inv = (H.transpose() * H);
  H.resize(0,0);
  Eigen::SparseMatrix<double> C_curv(pert_curv->dpsi->Sm,pert_curv->dpsi->Sm); // HtH_inv
  double dum_detHtH = 0.0;
  algebra->getInverseAndDet(C_curv_inv,C_curv,dum_detHtH);
  //  Eigen::SparseMatrix<double> curv_ide(pert_curv->dpsi->Sm,pert_curv->dpsi->Sm);
  //  curv_ide = C_curv_inv*C_curv;
  std::cout << "done" << std::endl;


  // Modgauss
  std::cout << "Modgauss...";
  pert_modgauss->dpsi->constructH();
  Eigen::SparseMatrix<double> C_modgauss(pert_modgauss->dpsi->H.Ti,pert_modgauss->dpsi->H.Tj);
  C_modgauss.reserve(Eigen::VectorXi::Constant(pert_modgauss->dpsi->Sm,1000));
  for(int i=0;i<pert_modgauss->dpsi->H.tri.size();i++){
    C_modgauss.insert(pert_modgauss->dpsi->H.tri[i].i,pert_modgauss->dpsi->H.tri[i].j) = pert_modgauss->dpsi->H.tri[i].v;
  }
  //  Eigen::SparseMatrix<double> C_modgauss_inv(pert_modgauss->dpsi->Sm,pert_modgauss->dpsi->Sm);
  //  double dum_detCmodgauss = 0.0;
  //  algebra->getInverseAndDet(C_modgauss,C_modgauss_inv,dum_detCmodgauss);
  //  Eigen::SparseMatrix<double> modgauss_ide(pert_modgauss->dpsi->Sm,pert_modgauss->dpsi->Sm);
  //  modgauss_ide = C_modgauss_inv*C_modgauss;
  std::cout << "done" << std::endl;


  // Gauss
  std::cout << "Gauss...";
  pert_gauss->dpsi->constructH();
  Eigen::SparseMatrix<double> C_gauss(pert_gauss->dpsi->H.Ti,pert_gauss->dpsi->H.Tj);
  C_gauss.reserve(Eigen::VectorXi::Constant(pert_gauss->dpsi->Sm,1000));
  for(int i=0;i<pert_gauss->dpsi->H.tri.size();i++){
    C_gauss.insert(pert_gauss->dpsi->H.tri[i].i,pert_gauss->dpsi->H.tri[i].j) = pert_gauss->dpsi->H.tri[i].v;
  }
  //  Eigen::SparseMatrix<double> C_gauss_inv(pert_gauss->dpsi->Sm,pert_gauss->dpsi->Sm);
  //  double dum_detCgauss = 0.0;
  //  algebra->getInverseAndDet(C_gauss,C_gauss_inv,dum_detCgauss);
  //  Eigen::SparseMatrix<double> gauss_ide(pert_modgauss->dpsi->Sm,pert_modgauss->dpsi->Sm);
  //  gauss_ide = C_gauss_inv*C_gauss;
  std::cout << "done" << std::endl;




  // Mask covariance matrices
  if( json_input.isMember("maskPath") ){
    std::cout << "Masking covariance matrices" << std::endl;
    ImagePlane* mask = new ImagePlane(json_input["maskPath"].asString(),30,30,3.5,3.5);
    maskMatrix(mask,C_curv);
    maskMatrix(mask,C_gauss);
    maskMatrix(mask,C_modgauss);
    std::cout << "done" << std::endl;
  }




  
  // Write the covariance, inverse covariance, and their product (identity) matrices
  std::cout << "Write the covariance, inverse covariance, and their product (identity) matrices" << std::endl;

  // Curvature
  std::cout << "Curvature...";
  writeMatrix(C_curv,pert_curv->dpsi,"matrix_curv.fits");
  //  writeMatrix(C_curv_inv,pert_curv->dpsi,"matrix_curv_inv.fits");
  //  writeMatrix(curv_ide,pert_curv->dpsi,"matrix_curv_ide.fits");
  std::cout << "done" << std::endl;
  
  // Modgauss
  std::cout << "Modgauss...";
  writeMatrix(C_modgauss,pert_curv->dpsi,"matrix_modgauss.fits");
  //  writeMatrix(C_modgauss_inv,pert_curv->dpsi,"matrix_modgauss_inv.fits");
  //  writeMatrix(modgauss_ide,pert_curv->dpsi,"matrix_modgauss_ide.fits");
  std::cout << "done" << std::endl;
  
  // Gauss
  std::cout << "Gauss...";
  writeMatrix(C_gauss,pert_gauss->dpsi,"matrix_gauss.fits");
  //  writeMatrix(C_gauss_inv,pert_gauss->dpsi,"matrix_gauss_inv.fits");
  //  writeMatrix(gauss_ide,pert_gauss->dpsi,"matrix_gauss_ide.fits");
  std::cout << "done" << std::endl;







  // Extract the two-point correlation function
  std::cout << "Extract the two-point correlation function" << std::endl;
  
  double Nbins = 50;
  double rmax = 4.0; // arcsec
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
      double r = hypot(pert_curv->dpsi->x[ia]-pert_curv->dpsi->x[ib],pert_curv->dpsi->y[ia]-pert_curv->dpsi->y[ib]);
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
      double r = hypot(pert_modgauss->dpsi->x[ia]-pert_modgauss->dpsi->x[ib],pert_modgauss->dpsi->y[ia]-pert_modgauss->dpsi->y[ib]);
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
      double r = hypot(pert_gauss->dpsi->x[ia]-pert_gauss->dpsi->x[ib],pert_gauss->dpsi->y[ia]-pert_gauss->dpsi->y[ib]);
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

  
  FILE* fh = fopen("corr_theo.dat","w");
  for(int i=0;i<Nbins;i++){
    fprintf(fh,"%5.3f %8.3f %8.3f %8.3f\n",bins[i],vals_curv[i],vals_modgauss[i],vals_gauss[i]);
    //fprintf(fh,"%5.3f %8.3f\n",bins[i],vals_curv[i]);
  }
  fclose(fh);



  C_curv.resize(0,0);
  //  C_curv_inv.resize(0,0);
  C_modgauss.resize(0,0);
  //  C_modgauss_inv.resize(0,0);
  C_gauss.resize(0,0);
  //  C_gauss_inv.resize(0,0);



  
  free(bins);
  free(vals_curv);
  free(vals_modgauss);
  free(vals_gauss);


  delete(pert_curv);
  delete(pert_modgauss);
  delete(pert_gauss);

  return 0;
}
