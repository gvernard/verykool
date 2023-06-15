#ifndef ALGEBRA_HPP
#define ALGEBRA_HPP

// Eigen implementation of all the linear algebra operations and related variables

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

#include "baseAlgebra.hpp"

class ImagePlane;
class BaseSourcePlane;
class BaseLikelihoodModel;
class SmoothLikelihood;
class PertLikelihood;
class BothLikelihood;
class Pert;


class SmoothAlgebra : public BaseAlgebra {
public:
  SmoothAlgebra(BaseLikelihoodModel* a);
  ~SmoothAlgebra(){};

  void setAlgebraInit(ImagePlane* image,BaseSourcePlane* source,double lambda);
  void setAlgebraRuntime(ImagePlane* image,BaseSourcePlane* source,double lambda);
  void solveSource(BaseSourcePlane* source,double lambda);
  void getSourceErrors(int Sm,double* errors);
  void getMockData(ImagePlane* mockdata,BaseSourcePlane* source);
  void solveDpsi(Pert* pert_mass_model,ImagePlane* image,BaseSourcePlane* source);

  BaseLikelihoodModel* likeModel;
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


class BothAlgebra : BaseAlgebra {
  // The public member variables of the base class are not used here, only the functions
public:
  BothAlgebra(BothLikelihood*);
  ~BothAlgebra(){};
  void constructDsDpsi(ImagePlane* image,BaseSourcePlane* source,Pert* pert_mass_model);
  void setAlgebraInit(BaseSourcePlane* source,Pert* pert_mass_model);
  void setAlgebraRuntime(BaseSourcePlane* source,Pert* pert_mass_model);
  void solveSourcePert(BaseSourcePlane* source,Pert* pert_mass_model);
  void getMockData(ImagePlane* mockdata,BaseSourcePlane* source,Eigen::SparseMatrix<double> B);
  void getMockData(ImagePlane* mockdata,BaseSourcePlane* source,Pert* pert_mass_model);

  
  BothLikelihood* likeModel;
  Eigen::SparseMatrix<double> DsDpsi;
  Eigen::SparseMatrix<double> RtR;
  Eigen::SparseMatrix<double> M_r;
  Eigen::SparseMatrix<double> Mt_r;
  Eigen::SparseMatrix<double> A_r;
  Eigen::SparseMatrix<double> Cp_inv;
};


#endif /* ALGEBRA_HPP */
