#ifndef BASE_ALGEBRA_HPP
#define BASE_ALGEBRA_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

class BaseSourcePlane;

class BaseAlgebra {
public:
  Eigen::VectorXd d;
  Eigen::SparseMatrix<double> StCS;
  Eigen::SparseMatrix<double> C_inv;
  Eigen::SparseMatrix<double> B;

  BaseAlgebra(){};
  ~BaseAlgebra(){};

  void setAlgebraField(BaseSourcePlane* source,Eigen::SparseMatrix<double> mat_in,Eigen::SparseMatrix<double>& mat_out,double& det_out);
  void getInverseAndDet(Eigen::SparseMatrix<double> mat_in,Eigen::SparseMatrix<double>& mat_out,double& det_out);
  double getDeterminant(Eigen::SparseMatrix<double> mat);
};

#endif /* BASE_ALGEBRA_HPP */
