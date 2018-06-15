#ifndef ALGEBRA_HPP
#define ALGEBRA_HPP

// Eigen implementation of all the linear algebra operations and related variables

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

class ImagePlane;
class BaseSourcePlane;
class BaseLikelihoodModel;
class Pert;

class StandardAlgebra {
public:
  StandardAlgebra(BaseLikelihoodModel* a);
  StandardAlgebra(){};
  ~StandardAlgebra(){};

  void setAlgebraInit(ImagePlane* mydata,BaseSourcePlane* mysource);
  void setAlgebraRuntime(ImagePlane* image,BaseSourcePlane* source,double lambda);
  void solveLinearSparseS(ImagePlane* image,BaseSourcePlane* source);
  void getSourceErrors(int Sm,double* errors);
  void getMockData(ImagePlane* mockdata,BaseSourcePlane* source);

  // this could be private, but I still want it inherited by any child class (like PerturbationsAlgebra)
  BaseLikelihoodModel* likeModel;
  Eigen::VectorXd d;
  Eigen::SparseMatrix<double> StCS;
  Eigen::SparseMatrix<double> C;
  Eigen::SparseMatrix<double> B;
  Eigen::SparseMatrix<double> M;
  Eigen::SparseMatrix<double> Mt;
  Eigen::SparseMatrix<double> A;
  Eigen::SparseMatrix<double> HtH;
};


class PerturbationsAlgebra : public StandardAlgebra {
public:
  PerturbationsAlgebra(BaseLikelihoodModel* a);
  ~PerturbationsAlgebra(){};

  void setAlgebraInitPert(ImagePlane* mydata,Pert* mypert,double* res);
  void setAlgebraRuntimePert(ImagePlane* image,BaseSourcePlane* source,Pert* mypert,BaseLikelihoodModel* smooth_like,double lambda);
  void solvePert(ImagePlane* image,BaseSourcePlane* dpsi,BaseLikelihoodModel* smooth_like);

  Eigen::VectorXd dd;
  Eigen::SparseMatrix<double> Dpsi;
  Eigen::SparseMatrix<double> M_pert;
  Eigen::SparseMatrix<double> Mt_pert;
  Eigen::SparseMatrix<double> HtH_pert;
  Eigen::SparseMatrix<double> A_pert;
};


#endif /* ALGEBRA_HPP */
