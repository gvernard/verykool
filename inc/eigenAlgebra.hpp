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
class StandardLikelihood;
class PerturbationsLikelihood;
class Pert;

class StandardAlgebra {
public:
  StandardAlgebra(StandardLikelihood* a);
  ~StandardAlgebra(){};

  void setAlgebraInit(ImagePlane* mydata,BaseSourcePlane* mysource);
  void setAlgebraRuntime(ImagePlane* image,BaseSourcePlane* source,double lambda);
  void solveLinearSparseS(ImagePlane* image,BaseSourcePlane* source);
  void getSourceErrors(int Sm,double* errors);
  void getMockData(ImagePlane* mockdata,BaseSourcePlane* source);


  StandardLikelihood* likeModel;
  Eigen::VectorXd d;
  Eigen::SparseMatrix<double> StCS;
  Eigen::SparseMatrix<double> C;
  Eigen::SparseMatrix<double> B;
  Eigen::SparseMatrix<double> M;
  Eigen::SparseMatrix<double> Mt;
  Eigen::SparseMatrix<double> A;
  Eigen::SparseMatrix<double> HtH;
};


class PerturbationsAlgebra {
public:
  PerturbationsAlgebra(PerturbationsLikelihood* a);
  ~PerturbationsAlgebra(){};

  void setAlgebraInit(ImagePlane* mydata,Pert* mypert,double* res);
  void setAlgebraRuntime(ImagePlane* image,BaseSourcePlane* source,Pert* pert_mass_model,BaseLikelihoodModel* smooth_like,double lambda);
  void solvePert(ImagePlane* image,Pert* pert_mass_model,BaseLikelihoodModel* smooth_like);


  PerturbationsLikelihood* likeModel;
  Eigen::VectorXd dd;
  Eigen::SparseMatrix<double> Dpsi;
  Eigen::SparseMatrix<double> HtH_pert;
  Eigen::SparseMatrix<double> M_pert;
  Eigen::SparseMatrix<double> Mt_pert;
  Eigen::SparseMatrix<double> A_pert;
};


#endif /* ALGEBRA_HPP */
