#include "tableAlgebra.hpp"

#include <iostream>

#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "nonLinearPars.hpp"

precomp::precomp(ImagePlane* image,BaseSourcePlane* source){
  int Nm = image->Nm;  
  int Sm = source->Sm;
  
  d.resize(Nm);
  C.resize(Nm,Nm);
  B.resize(Nm,Nm);
  M.resize(Nm,Sm);
  Mt.resize(Sm,Nm);
  A.resize(Sm,Sm);
  HtH.resize(Sm,Sm);
  detC   = 1.;
  detHtH = 1.;
  detA   = 1.;
}

precomp::~precomp(){
  d.resize(0);
  C.resize(0,0);
  B.resize(0,0);
  M.resize(0,0);
  Mt.resize(0,0);
  A.resize(0,0);
  HtH.resize(0,0);
}


//This is independent of the specific linear algebra package used.
// S: mask,
// C: covariance matrix,
// B: blurring matrix,
// L: lensing matrix,
// H: source regularization matrix
void mymatrices::initMatrices(ImagePlane* image,BaseSourcePlane* source,std::string maskpath,std::string noise_flag,std::string covpath,std::string psfpath,std::map<std::string,std::string> psf){
  //Read the mask (maximum allowed size Nm x Nm). The indices of the pixels in the data are matched to the indices of the pixels in the mask by the 'lookup' map.
  this->S.Ti = image->Nm;
  this->S.Tj = image->Nm;
  image->readS(&this->S,maskpath);

  //Read data covariance matrix (full size = Ni*Nj x Ni*Nj or Nm x Nm, but actual non-zero elements depend on the noise_flag)
  //ATTENTION: in the following C is always the inverse of the covariance matrix.
  this->C.Ti = image->Nm;
  this->C.Tj = image->Nm;
  image->readC(&this->C,noise_flag,covpath);

  //Construct blurring matrix (Ni*Nj x Ni*Nj) or (Nm x Nm)
  this->B.Ti = image->Nm;
  this->B.Tj = image->Nm;
  image->readB(&this->B,psfpath,stoi(psf["pix_x"]),stoi(psf["pix_y"]),stoi(psf["crop_x"]),stoi(psf["crop_y"]));

  //Lensing matrix (Ni*Nj x Si*Sj) or (Nm x Sm) initialize here, but also in each iteration                                               [<---iteration dependent]
  this->L.Ti = image->Nm;
  this->L.Tj = source->Sm;
  //mysource->constructL(&my_data,mycollection,&matrices.L);//No need to construct the lensing matrix here

  //Regularization matrix (Si*Sj x Si*Sj) or (Sm x Sm)                                                                [<---iteration dependent for adaptive source]
  this->H.Ti = source->Sm;
  this->H.Tj = source->Sm;
  source->constructH(&this->H);
}


// This function generates the model image from the given source and nlpars using tables
void getMockData(ImagePlane* mockdata,BaseSourcePlane* source,precomp* pcomp){
  Eigen::Map<Eigen::VectorXd> s(source->src,source->Sm);
  Eigen::VectorXd tmp = pcomp->M*s;
  Eigen::Map<Eigen::VectorXd>(mockdata->img,tmp.size()) = tmp;
}

// Initialize tables and related quantities
void setAlgebraInit(ImagePlane* image,BaseSourcePlane* source,mymatrices* mat,precomp* pcomp){
  int Nm = image->Nm;
  int Sm = source->Sm;

  Eigen::Map<Eigen::VectorXd> d(image->img,Nm);
  
  Eigen::SparseMatrix<double> SS(mat->S.Ti,mat->S.Tj);
  SS.reserve(Eigen::VectorXi::Constant(mat->S.Ti,1));//overestimating the covariance matrix number of entries per row
  for(int i=0;i<mat->S.tri.size();i++){  SS.insert(mat->S.tri[i].i,mat->S.tri[i].j) = mat->S.tri[i].v;  }

  Eigen::SparseMatrix<double> CC(mat->C.Ti,mat->C.Tj);
  CC.reserve(Eigen::VectorXi::Constant(mat->C.Ti,8));//overestimating the covariance matrix number of entries per row
  for(int i=0;i<mat->C.tri.size();i++){  CC.insert(mat->C.tri[i].i,mat->C.tri[i].j) = mat->C.tri[i].v;  }

  Eigen::SparseMatrix<double> BB(mat->B.Ti,mat->B.Tj);
  BB.reserve(Eigen::VectorXi::Constant(mat->B.Ti,900));//overestimating the number of entries per row of the blurring matrix BB (900 for a 30x30 psf)
  //  for(int i=0;i<mat->B.tri.size();i++){  BB.insert(mat->B.tri[i].i,mat->B.tri[i].j) = mat->B.tri[i].v;  }//for row-major
  for(int i=0;i<mat->B.tri.size();i++){  BB.insert(mat->B.tri[i].j,mat->B.tri[i].i) = mat->B.tri[i].v;  }//for column-major (swap the i and j indices of tri struct)

  Eigen::SparseMatrix<double> HH(Sm,Sm);
  HH.reserve(Eigen::VectorXi::Constant(Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH row (different number for 1st,2nd order derivative etc)
  for(int i=0;i<mat->H.tri.size();i++){  HH.insert(mat->H.tri[i].i,mat->H.tri[i].j) = mat->H.tri[i].v;  }


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





  pcomp->d       = d;
  pcomp->C       = CC;
  pcomp->StCS    = StCS;
  pcomp->detC    = detC;
  pcomp->B       = BB;
  pcomp->HtH     = HtH;
  pcomp->detHtH  = detHtH;

  SS.resize(0,0);
  SSt.resize(0,0);
  CC.resize(0,0);
  BB.resize(0,0);
  HH.resize(0,0);
  HtH.resize(0,0);
}


// Calculate tables and related quantities at every iteration
void setAlgebraRuntime(ImagePlane* image,BaseSourcePlane* source,std::map<std::string,BaseNlpar*> reg_pars,mymatrices* mat,precomp* pcomp){
  int Nm = image->Nm;
  int Sm = source->Sm;


  Eigen::SparseMatrix<double> LL(mat->L.Ti,mat->L.Tj);
  LL.reserve(Eigen::VectorXi::Constant(mat->L.Ti,4));//setting the non-zero coefficients per row of the lensing matrix (4 for bi-linear interpolation scheme, etc)
  for(int i=0;i<mat->L.tri.size();i++){  LL.insert(mat->L.tri[i].i,mat->L.tri[i].j) = mat->L.tri[i].v;  }
  
  Eigen::SparseMatrix<double> M(Nm,Sm);
  Eigen::SparseMatrix<double> Mt(Sm,Nm);

  M         = pcomp->B*LL;
  Mt        = M.transpose();
  pcomp->M  = M;
  pcomp->Mt = Mt;


  //I need also to calculate HtH, and detHtH for an adaptive source
  if( source->type == "adaptive" || source->sample_reg ){

    Eigen::SparseMatrix<double> HH(Sm,Sm);
    HH.reserve(Eigen::VectorXi::Constant(Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH row (different number for 1st,2nd order derivative etc)
    for(int i=0;i<mat->H.tri.size();i++){  HH.insert(mat->H.tri[i].i,mat->H.tri[i].j) = mat->H.tri[i].v;  }
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
    pcomp->HtH = HtH;
    
    
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
    pcomp->detHtH  = detHtH;


    HtH.resize(0,0);
  }
  
  Eigen::SparseMatrix<double> A(Sm,Sm);
  if( source->reg == "covariance_kernel" ){
    A = Mt*pcomp->C*M + pcomp->HtH;
  } else {
    A = Mt*pcomp->C*M + reg_pars["lambda"]->val*pcomp->HtH;
  }
  pcomp->A = A;


  LL.resize(0,0);
  M.resize(0,0);
  Mt.resize(0,0);
  A.resize(0,0);
}

// Solve for the source
void solveLinearSparseS(ImagePlane* image,BaseSourcePlane* source,precomp* pcomp){
  int Nm = image->Nm;
  int Sm = source->Sm;

  Eigen::SparseMatrix<double> A(Sm,Sm);
  A = pcomp->A;
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
  pcomp->detA = detA;

  //Get the Most Probable solution for the source
  Eigen::VectorXd s(Sm);
  s = inv*(pcomp->Mt*pcomp->C*pcomp->d);
  Eigen::Map<Eigen::VectorXd>(source->src,s.size()) = s;

  inv.resize(0,0);
  A.resize(0,0);
}




void getSourceErrors(int Sm,double* errors,precomp* pcomp){
  Eigen::SparseMatrix<double> A(Sm,Sm);
  A = pcomp->A;
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




// Calculate the chi squared and regularization terms using table algebra
void getMin(ImagePlane* image,BaseSourcePlane* source,precomp* pcomp){
  Eigen::Map<Eigen::VectorXd> s(source->src,source->Sm);
  Eigen::VectorXd st = s.transpose();

  Eigen::VectorXd y  = pcomp->d - (pcomp->M*s);
  Eigen::VectorXd yt = y.transpose();
  double chi2        = yt.dot(pcomp->StCS*y);
  double reg         = st.dot(pcomp->HtH*s);

  pcomp->chi2 = chi2;
  pcomp->reg  = reg;
}


double getLogLike(ImagePlane* image,BaseSourcePlane* source,precomp* pcomp,std::vector<std::map<std::string,BaseNlpar*> > nlpars){
  double pi  = 3.14159265358979323846;

  if( source->reg == "covariance_kernel" ){
    getMin(image,source,pcomp);
    double g   = pcomp->chi2/2. + pcomp->reg/2.;
    double f1  = image->lookup.size()*log10(2*pi)/2.;
    //double f2  = source->Sm*log10(nlpars[1]["lambda"]->val)/2.;
    double val = -g -f1 +(pcomp->detC)/2. +(pcomp->detHtH)/2. -(pcomp->detA)/2.;
    
    for(int i=0;i<nlpars.size();i++){
      for(auto it = nlpars[i].cbegin(); it != nlpars[i].cend(); ++it){
	printf(" %12.8f",it->second->val);
      }
    }
    printf("\n");
    printf("%16.7f  %16.7f  %16.7f  %16.7f  %16.7f  %16.7f  %16.7f\n",pcomp->chi2,pcomp->reg,-g,-f1,pcomp->detC/2.,pcomp->detHtH/2.,-pcomp->detA/2.);
    printf("%16.7f\n\n",val);

    return val;
  } else {
    getMin(image,source,pcomp);
    double g   = pcomp->chi2/2.0 + nlpars[1]["lambda"]->val*pcomp->reg/2.0;
    double f1  = image->lookup.size()*log10(2*pi)/2.0;
    double f2  = source->Sm*log10(nlpars[1]["lambda"]->val)/2.0;
    double val = -g -f1 +f2 +(pcomp->detC)/2.0 +(pcomp->detHtH)/2.0 -(pcomp->detA)/2.0;
    
    for(int i=0;i<nlpars.size();i++){
      for(auto it = nlpars[i].cbegin(); it != nlpars[i].cend(); ++it){
	printf(" %12.8f",it->second->val);
      }
    }
    printf("\n");
    printf("%16.7f  %16.7f  %16.7f  %16.7f  %16.7f  %16.7f  %16.7f  %16.7f  %16.7f\n",pcomp->chi2,pcomp->reg,nlpars[1]["lambda"]->val*pcomp->reg,-g,-f1,f2,pcomp->detC/2.0,pcomp->detHtH/2.0,-pcomp->detA/2.0);
    printf("%16.7f\n\n",val);

    return val;
  }

}



//void getInverseCovarianceMatrix(BaseSourcePlane* source,std::map<std::string,BaseNlpar*> pars,precomp* pcomp){
void getInverseCovarianceMatrixHODLR(mymatrices* mat,std::map<std::string,BaseNlpar*> pars,precomp* pcomp){

  // Set up the kernel.
  //  gaussKernel kernel(source->Sm,source->x,source->y,pars);
  gaussKernel kernel(3);

  // Settings for the solver.
  //  unsigned N       = source->Sm;
  unsigned N       = mat->C.Ti;
  unsigned nLeaf   = 50;
  double tolerance = 1e-14;

  // Build the RHS matrix.
  Eigen::MatrixXd b(N,1),x;
  Eigen::VectorXd yvar(N);
  for(int i=0;i<N;i++){
    b(i,0) = 1.0;
    //    b(i,0) = source->src[i];
    //    yvar(i) = _ferr[i] * _ferr[i];
    //    yvar(i) = 1.0;
    yvar(i) = sqrt(1.0/mat->C.tri[i].v);
  }

  // Setting things up
  HODLR_Tree<gaussKernel> * A = new HODLR_Tree<gaussKernel>(&kernel,N,nLeaf);

  // Assembling the matrix in HODLR form
  A->assemble_Matrix(yvar,tolerance,'s');

  // Factoring the matrix
  A->compute_Factor();

  // Solving the system
  A->solve(b,x);
  
  // Computing the determinant
  double determinant;
  A->compute_Determinant(determinant);

  //Get the elements of the inverse
  A->computeInverse();


  pcomp->detC = determinant;
  std::cout << "det: " << determinant << std::endl;
  //  pcomp->C    = C;
}

