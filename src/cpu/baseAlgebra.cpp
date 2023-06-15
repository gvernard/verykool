#include "baseAlgebra.hpp"

#include "sourcePlane.hpp"



// Class: BaseAlgebra
//===============================================================================================================
void BaseAlgebra::setAlgebraField(BaseSourcePlane* source,Eigen::SparseMatrix<double> mat_in,Eigen::SparseMatrix<double>& mat_out,double& det_out){
  Eigen::SparseMatrix<double> out(source->Sm,source->Sm);

  if( source->reg == "covariance_kernel" || source->reg == "identity" ){
    this->getInverseAndDet(mat_in,out,det_out);
  } else {
    // calculate the product HtH of a derivative based H matrix
    Eigen::SparseMatrix<double> mat_in_t(mat_in.rows(),mat_in.cols());
    mat_in_t = mat_in.transpose();
    out = (mat_in_t * mat_in);    
    mat_in_t.resize(0,0);
    det_out = -this->getDeterminant(out); // the minus sign is because Cs is actually (HtH)^-1
  }
  
  mat_out = out;
  out.resize(0,0);
}

void BaseAlgebra::getInverseAndDet(Eigen::SparseMatrix<double> mat_in,Eigen::SparseMatrix<double>& mat_out,double& det_out){
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> inv(mat_in.rows(),mat_in.cols());
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  Eigen::VectorXd idc(mat_in.cols()),c(mat_in.cols());
  double det = 1.0;

  
  // get inverse
  solver.compute(mat_in);
  //std::cout << "Factorization (" << mat_in.cols() << ") done" << std::endl;

  
  for(int i=0;i<mat_in.rows();i++){
    for(int j=0;j<idc.size();j++){
      idc[j] = 0;
    }
    idc[i] = 1;
    //this needs to be done in two steps, otherwise it does not work when using LU
    c = solver.solve(idc);
    inv.col(i) = c;
  }
  
  //Eigen::MatrixXd b = Eigen::MatrixXd::Identity(mat_in.cols(),mat_in.cols());
  //inv = solver.solve(b);

  //std::cout << "Inversion done" << std::endl;

  
  // get determinant of mat_in (NOT its inverse)
  Eigen::VectorXd diag = solver.vectorD();
  diag = solver.vectorD();
  for(int i=0;i<diag.size();i++){
    det += log( *(diag.data()+i) );
  }
  //std::cout << "Determinant done" << std::endl;

  
  det_out = det;
  mat_out = inv.sparseView();
  inv.resize(0,0);
}

double BaseAlgebra::getDeterminant(Eigen::SparseMatrix<double> mat){
  Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
  Eigen::VectorXd diag;
  double det = 1.0;

  solver.analyzePattern(mat);
  solver.factorize(mat);
  diag = solver.vectorD();
  for(int i=0;i<diag.size();i++){
    det += log( *(diag.data()+i) );
  }
  
  return det;
}

