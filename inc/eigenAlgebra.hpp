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
class SmoothLikelihood;
class PertLikelihood;
class PertIterationLikelihood;
class Pert;



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



class SmoothAlgebra : public BaseAlgebra {
public:
  SmoothAlgebra(SmoothLikelihood* a);
  ~SmoothAlgebra(){};

  void setAlgebraInit(ImagePlane* image,BaseSourcePlane* source);
  void setAlgebraRuntime(ImagePlane* image,BaseSourcePlane* source);
  void solveSource(BaseSourcePlane* source);
  void getSourceErrors(int Sm,double* errors);
  void getMockData(ImagePlane* mockdata,BaseSourcePlane* source);
  void solveDpsi(Pert* pert_mass_model,ImagePlane* image,BaseSourcePlane* source);

  SmoothLikelihood* likeModel;
  Eigen::SparseMatrix<double> M;
  Eigen::SparseMatrix<double> Mt;
  Eigen::SparseMatrix<double> A;
  Eigen::SparseMatrix<double> Cs_inv;
  Eigen::SparseMatrix<double> block;
};



class PertAlgebra : public BaseAlgebra {
public:
  PertAlgebra(PertLikelihood* a);
  ~PertAlgebra(){};

  void setAlgebraInit(ImagePlane* image,BaseSourcePlane* source,BaseSourcePlane* source0,Pert* pert_mass_model);
  void setAlgebraRuntime(BaseSourcePlane* source,BaseSourcePlane* source0,Pert* pert_mass_model);
  void solveSourcePert(BaseSourcePlane* source,Pert* pert_mass_model);
  void constructDsDpsi(ImagePlane* image,BaseSourcePlane* source,Pert* pert_mass_model);
  void constructNormalizingJmatrix(BaseSourcePlane* source,Pert* pert_mass_model,Eigen::SparseMatrix<double>& J_out,double lambda_s,double lambda_dpsi);
  void solvePerturbationsResiduals(std::string output,Pert* pert_mass_model);
  void getMockData(ImagePlane* mockdata,BaseSourcePlane* source,Eigen::SparseMatrix<double> B);

  PertLikelihood* likeModel;
  Eigen::SparseMatrix<double> DsDpsi;
  Eigen::SparseMatrix<double> RtR;
  Eigen::SparseMatrix<double> M_r;
  Eigen::SparseMatrix<double> Mt_r;
  Eigen::SparseMatrix<double> A_r;
  Eigen::SparseMatrix<double> Cs_inv;
  Eigen::SparseMatrix<double> Cp_inv;
  Eigen::SparseMatrix<double> J;
  Eigen::SparseMatrix<double> block_s;
  Eigen::SparseMatrix<double> block_p;
};


#endif /* ALGEBRA_HPP */
