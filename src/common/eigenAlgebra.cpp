#include "eigenAlgebra.hpp"

#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "likelihoodModels.hpp"

#include <iostream>

StandardAlgebra::StandardAlgebra(StandardLikelihood* a){
  this->likeModel = a;
}



void StandardAlgebra::setAlgebraInit(ImagePlane* image,BaseSourcePlane* source){
  int Nm = image->Nm;
  int Sm = source->Sm;

  Eigen::Map<Eigen::VectorXd> d(image->img,Nm);
  
  Eigen::SparseMatrix<double> SS(image->S.Ti,image->S.Tj);
  SS.reserve(Eigen::VectorXi::Constant(image->S.Ti,1));//overestimating the mask matrix number of entries per row
  for(int i=0;i<image->S.tri.size();i++){  SS.insert(image->S.tri[i].i,image->S.tri[i].j) = image->S.tri[i].v;  }

  Eigen::SparseMatrix<double> CC(image->C.Ti,image->C.Tj);
  CC.reserve(Eigen::VectorXi::Constant(image->C.Ti,8));//overestimating the covariance matrix number of entries per row
  for(int i=0;i<image->C.tri.size();i++){  CC.insert(image->C.tri[i].i,image->C.tri[i].j) = image->C.tri[i].v;  }

  Eigen::SparseMatrix<double> BB(image->B.Ti,image->B.Tj);
  BB.reserve(Eigen::VectorXi::Constant(image->B.Ti,900));//overestimating the number of entries per row of the blurring matrix BB (900 for a 30x30 psf)
  //  for(int i=0;i<imageB.tri.size();i++){  BB.insert(mat->B.tri[i].i,mat->B.tri[i].j) = mat->B.tri[i].v;  }//for row-major
  for(int i=0;i<image->B.tri.size();i++){  BB.insert(image->B.tri[i].j,image->B.tri[i].i) = image->B.tri[i].v;  }//for column-major (swap the i and j indices of tri struct)

  Eigen::SparseMatrix<double> HH(Sm,Sm);
  HH.reserve(Eigen::VectorXi::Constant(Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH row (different number for 1st,2nd order derivative etc)
  for(int i=0;i<source->H.tri.size();i++){  HH.insert(source->H.tri[i].i,source->H.tri[i].j) = source->H.tri[i].v;  }


  /*
  //calculate determinant of CCmasked
  Eigen::SparseMatrix<double> CCtmp(mat->S.Ti,mat->S.Tj);
  CCtmp = (CC * SS);
  //now fix the indices to be in the reduced ranged of CCmasked
  Eigen::SparseMatrix<double> CCmasked(image->lookup.size(),image->lookup.size());
  for(int k=0;k<CCtmp.outerSize();k++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(CCtmp,k);it;++it){
	CCmasked.insert( image->lookup[it.row()], image->lookup[it.col()] ) = it.value();
    }
  }
  //  CCtmp.resize(0,0);
  */

  Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
  Eigen::VectorXd diag;
  double dum;


  //calculate determinant of CC
  double detC = 0.;
  solver.analyzePattern(CC);
  solver.factorize(CC);
  diag = solver.vectorD();
  for(int i=0;i<diag.size();i++){
    dum = *(diag.data()+i);
    if( dum > 1.e-20 ){
      detC += log10( dum );
    }
  }


  //Calculate StCS
  Eigen::SparseMatrix<double> SSt(Nm,Nm);
  Eigen::SparseMatrix<double> StCS(Nm,Nm);
  SSt   = SS.transpose();
  StCS  = (SSt * CC * SS);


  // Calculate HtH
  Eigen::SparseMatrix<double> HtH(Sm,Sm);
  if( source->reg == "covariance_kernel" ){
    // calculate the inverse of the source covariance matrix stored in H, and store it in HtH
  } else {
    // calculate the product HtH of a derivative based H matrix
    Eigen::SparseMatrix<double> HHt(Sm,Sm);
    HHt = HH.transpose();
    HtH = (HHt * HH);    
    HHt.resize(0,0);
  }

  //calculate determinant of HtH
  double detHtH = 0.;
  solver.analyzePattern(HtH);
  solver.factorize(HtH);
  diag = solver.vectorD();
  for(int i=0;i<diag.size();i++){
    dum = *(diag.data()+i);
    if( dum > 1.e-20 ){
      detHtH += log10( dum );
    }
  }





  this->d       = d;
  this->C       = CC;
  this->StCS    = StCS;
  this->B       = BB;
  this->HtH     = HtH;

  this->likeModel->terms["detC"]    = detC/2.0;
  this->likeModel->terms["detHtH"]  = detHtH/2.0;

  SS.resize(0,0);
  SSt.resize(0,0);
  CC.resize(0,0);
  BB.resize(0,0);
  HH.resize(0,0);
  HtH.resize(0,0);
}



// Calculate tables and related quantities at every iteration
void StandardAlgebra::setAlgebraRuntime(ImagePlane* image,BaseSourcePlane* source,double lambda){
  int Nm = image->Nm;
  int Sm = source->Sm;


  Eigen::SparseMatrix<double> LL(source->L.Ti,source->L.Tj);
  LL.reserve(Eigen::VectorXi::Constant(source->L.Ti,4));//setting the non-zero coefficients per row of the lensing matrix (4 for bi-linear interpolation scheme, etc)
  for(int i=0;i<source->L.tri.size();i++){  LL.insert(source->L.tri[i].i,source->L.tri[i].j) = source->L.tri[i].v;  }
  
  Eigen::SparseMatrix<double> M(Nm,Sm);
  Eigen::SparseMatrix<double> Mt(Sm,Nm);

  M         = this->B*LL;
  Mt        = M.transpose();
  this->M  = M;
  this->Mt = Mt;


  //I need also to calculate HtH, and detHtH for an adaptive source
  if( source->type == "adaptive" || source->sample_reg ){

    Eigen::SparseMatrix<double> HH(Sm,Sm);
    HH.reserve(Eigen::VectorXi::Constant(Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH row (different number for 1st,2nd order derivative etc)
    for(int i=0;i<source->H.tri.size();i++){  HH.insert(source->H.tri[i].i,source->H.tri[i].j) = source->H.tri[i].v;  }
    Eigen::SparseMatrix<double> HtH(Sm,Sm);

    // Calculate HtH
    if( source->sample_reg ){
      // calculate the inverse of the source covariance matrix stored in H, and store it in HtH

      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> inv(Sm,Sm);
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
      solver.compute(HH);
      Eigen::VectorXd idc(HH.cols()),c(HH.cols());
      for(int i=0;i<Sm;i++){
	for(int j=0;j<idc.size();j++){
	  idc[j] = 0;
	}
	idc[i] = 1;
	//this needs to be done in two steps, otherwise it does not work when using LU
	c = solver.solve(idc);
	inv.col(i) = c;
      }
      HtH = inv.sparseView();

      inv.resize(0,0);
      
    } else {
      // calculate the product HtH of a derivative based H matrix
      Eigen::SparseMatrix<double> HHt(Sm,Sm);
      HHt = HH.transpose();
      HtH = (HHt * HH);    
      HHt.resize(0,0);
    }
    HH.resize(0,0);   // no longer needed, all information has been passed to HtH
    this->HtH = HtH;
    
    
    // Calculate the determinant of HtH
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
    Eigen::VectorXd diag;
    double dum;
    double detHtH = 0.;
    solver.analyzePattern(HtH);
    solver.factorize(HtH);
    diag = solver.vectorD();
    for(int i=0;i<diag.size();i++){
      dum = *(diag.data()+i);
      if( dum > 1.e-20 ){
	detHtH += log10( dum );
      }
    }
    this->likeModel->terms["detHtH"] = detHtH/2.0;


    HtH.resize(0,0);
  }
  
  Eigen::SparseMatrix<double> A(Sm,Sm);
  if( source->reg == "covariance_kernel" ){
    A = Mt*this->C*M + this->HtH;
  } else {
    A = Mt*this->C*M + lambda*this->HtH;
  }
  this->A = A;



  LL.resize(0,0);
  M.resize(0,0);
  Mt.resize(0,0);
  A.resize(0,0);
}


// Solve for the source
void StandardAlgebra::solveLinearSparseS(ImagePlane* image,BaseSourcePlane* source){
  int Nm = image->Nm;
  int Sm = source->Sm;

  Eigen::SparseMatrix<double> A(Sm,Sm);
  A = this->A;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> inv(Sm,Sm);
  
  //Get the inverse
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  solver.compute(A);
  Eigen::VectorXd idc(A.cols()),c(A.cols());
  for(int i=0;i<Sm;i++){
    for(int j=0;j<idc.size();j++){
      idc[j] = 0;
    }
    idc[i] = 1;
    //this needs to be done in two steps, otherwise it does not work when using LU
    c = solver.solve(idc);
    inv.col(i) = c;
  }

  //Calculate determinant of A from the LDLT(Cholesky) decomposition
  double detA = 0.;
  Eigen::VectorXd diag = solver.vectorD();
  for(int i=0;i<diag.size();i++){
    detA += log10( *(diag.data()+i) );
  }
  this->likeModel->terms["detA"] = -detA/2.0;

  //Get the Most Probable solution for the source
  Eigen::VectorXd s(Sm);
  s = inv*(this->Mt*this->C*this->d);
  Eigen::Map<Eigen::VectorXd>(source->src,s.size()) = s;

  //Get the chi-squared term
  Eigen::VectorXd st = s.transpose();
  Eigen::VectorXd y  = this->d - (this->M*s);
  Eigen::VectorXd yt = y.transpose();
  double chi2        = yt.dot(this->StCS*y);
  double reg         = st.dot(this->HtH*s);

  this->likeModel->terms["chi2"] = -chi2/2.0;
  this->likeModel->terms["reg"]  = -reg/2.0;


  inv.resize(0,0);
  A.resize(0,0);
}


void StandardAlgebra::getMockData(ImagePlane* mockdata,BaseSourcePlane* source){
  Eigen::Map<Eigen::VectorXd> s(source->src,source->Sm);
  Eigen::VectorXd tmp = this->M*s;
  Eigen::Map<Eigen::VectorXd>(mockdata->img,tmp.size()) = tmp;
}


void StandardAlgebra::getSourceErrors(int Sm,double* errors){
  Eigen::SparseMatrix<double> A(Sm,Sm);
  A = this->A;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> inv(Sm,Sm);
  
  //Get the inverse
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  solver.compute(A);
  Eigen::VectorXd idc(A.cols()),c(A.cols());
  for(int i=0;i<Sm;i++){
    for(int j=0;j<idc.size();j++){
      idc[j] = 0;
    }
    idc[i] = 1;
    //this needs to be done in two steps, otherwise it does not work when using LU
    c = solver.solve(idc);
    inv.col(i) = c;
  }


  //  get the diagonal elements of matrix inv as source values
  Eigen::VectorXd diag = inv.diagonal();
  for(int i=0;i<diag.size();i++){
    errors[i] = diag[i];
  }

}

