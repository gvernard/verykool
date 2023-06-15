#include "baseAlgebra.hpp"

#include "sourcePlane.hpp"

#include "cusparse_cholesky_solver.h"

using Ordering = Eigen::AMDOrdering<Eigen::SparseMatrix<double,Eigen::StorageOptions::RowMajor>::StorageIndex>;
using PermutationMatrix = Ordering::PermutationType;

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
  Eigen::SparseMatrix<double,Eigen::StorageOptions::RowMajor> Acsr = mat_in; // solver supports CSR format
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> inv(mat_in.rows(),mat_in.cols());
  Eigen::VectorXd idc(mat_in.cols()),diag(mat_in.cols()),c(mat_in.cols());
  int nnz = mat_in.nonZeros();
  double det = 1.0;


  auto solver = CuSparseCholeskySolver<double>::create(mat_in.cols());
  
  // compute permutation
  PermutationMatrix P;
  Ordering ordering;
  ordering(Acsr.selfadjointView<Eigen::Upper>(), P);
  // set permutation to solver
  solver->setPermutaion(mat_in.cols(), P.indices().data());

  solver->analyze(nnz, Acsr.outerIndexPtr(), Acsr.innerIndexPtr());
  solver->factorize(Acsr.valuePtr());
  std::cout << "Factorization (" << mat_in.cols() << ") done" << std::endl;

  for(int i=0;i<mat_in.rows();i++){
    for(int j=0;j<idc.size();j++){
      idc[j] = 0;
    }
    idc[i] = 1;
    //this needs to be done in two steps, otherwise it does not work when using LU
    solver->solve(idc.data(),c.data());
    inv.col(i) = c;
  }
  std::cout << "Inversion done" << std::endl;

  /*
  solver->diagonal(diag.data());
  for(int i=0;i<diag.size();i++){
    det += log( *(diag.data()+i) );
  }
  */
  std::cout << "Determinant done" << std::endl;
  
  det_out = det;
  mat_out = inv.sparseView();
  inv.resize(0,0);
}

double BaseAlgebra::getDeterminant(Eigen::SparseMatrix<double> mat){
  Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
  Eigen::VectorXd diag;
  double det = 1.0;

  /*
  auto solver = CuSparseCholeskySolver<double>::create(mat.cols());
  
  // compute permutation
  PermutationMatrix P;
  Ordering ordering;
  ordering(Acsr.selfadjointView<Eigen::Upper>(), P);
  // set permutation to solver
  solver->setPermutaion(mat.cols(), P.indices().data());

  solver->analyze(mat.nonZeros(), Acsr.outerIndexPtr(), Acsr.innerIndexPtr());
  solver->factorize(Acsr.valuePtr());
  std::cout << "Factorization done" << std::endl;
  
  solver->diagonal(diag.data());
  for(int i=0;i<diag.size();i++){
    det += log( *(diag.data()+i) );
  }
  */
  
  return det;
}

