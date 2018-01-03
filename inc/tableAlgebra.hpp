#ifndef TABLE_ALGEBRA_HPP
#define TABLE_ALGEBRA_HPP

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"



class BaseSourcePlane;
class ImagePlane;
class BaseNlpar;

struct mytriplet {
  int i;
  int j;
  double v;
};

struct mytable {
  int Ti;
  int Tj;
  std::vector<mytriplet> tri;
};


class mymatrices {
public:
  mytable S;
  mytable C;
  mytable B;
  mytable H;
  mytable L;

  mymatrices(){};
  ~mymatrices(){};

  void initMatrices(ImagePlane* image,BaseSourcePlane* source,std::string maskpath,std::string noise_flag,std::string covpath,std::string psfpath,std::map<std::string,std::string> psf);
};

class precomp {
public:
  Eigen::VectorXd d;
  Eigen::SparseMatrix<double> StCS;
  Eigen::SparseMatrix<double> C;
  Eigen::SparseMatrix<double> B;
  Eigen::SparseMatrix<double> M;
  Eigen::SparseMatrix<double> Mt;
  Eigen::SparseMatrix<double> A;
  Eigen::SparseMatrix<double> HtH;
  double detC;
  double detHtH;
  double detA;
  double chi2;
  double reg;

  precomp(){};
  precomp(ImagePlane* image,BaseSourcePlane* source);
  ~precomp();
};

void setAlgebraInit(ImagePlane* image,BaseSourcePlane* source,mymatrices* mat,precomp* pcomp);
void setAlgebraRuntime(ImagePlane* image,BaseSourcePlane* source,std::map<std::string,BaseNlpar*> reg_pars,mymatrices* mat,precomp* pcomp);
void solveLinearSparseS(ImagePlane* image,BaseSourcePlane* source,precomp* pcomp);
double getLogLike(ImagePlane* image,BaseSourcePlane* source,precomp* pcomp,std::vector<std::map<std::string,BaseNlpar*> > nlpars);

void getSourceErrors(int Sm,double* errors,precomp* pcomp);
void getMockData(ImagePlane* mockdata,BaseSourcePlane* source,precomp* pcomp);
void getMin(ImagePlane* image,BaseSourcePlane* source,precomp* pcomp);



class gaussKernel : public HODLR_Matrix {
public:
  //  gaussKernel(int Sm,double* x,double* y,std::map<std::string,BaseNlpar*> pars) : pSm(Sm), px(x), py(y), ppars(pars) {};
  gaussKernel(int dum){};
  
  double get_Matrix_Entry(const unsigned i,const unsigned j){
    //double d = sqrt( pow(px[i]-px[j],2) + pow(py[i]-py[j],2) );
    //return ppars["spec_ampl"]->val*exp( -0.5*d*d/ppars["spec_sig"]->val );
    return 1.0;
  }

private:
  int pSm;
  double* px;
  double* py;
  std::map<std::string,BaseNlpar*> ppars;  
};

//void getInverseCovarianceMatrix(BaseSourcePlane* source,std::map<std::string,BaseNlpar*> pars,precomp* pcomp);
void getInverseCovarianceMatrixHODLR(mymatrices* mat,std::map<std::string,BaseNlpar*> pars,precomp* pcomp);


#endif /* TABLE_ALGEBRA_HPP */
