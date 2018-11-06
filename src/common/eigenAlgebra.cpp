#include "eigenAlgebra.hpp"

#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "massModels.hpp"
#include "likelihoodModels.hpp"

#include <iostream>
#include <cmath>


// Class: BaseAlgebra
//===============================================================================================================
void BaseAlgebra::setAlgebraField(BaseSourcePlane* source,Eigen::SparseMatrix<double> mat_in,Eigen::SparseMatrix<double>& mat_out,double& det_out){
  Eigen::SparseMatrix<double> out(source->Sm,source->Sm);

  if( source->reg == "covariance_kernel" ){
    this->getInverseAndDet(mat_in,out,det_out);
  } else {
    // calculate the product HtH of a derivative based H matrix
    Eigen::SparseMatrix<double> mat_in_t(mat_in.rows(),mat_in.cols());
    mat_in_t = mat_in.transpose();
    out = (mat_in_t * mat_in);    
    mat_in_t.resize(0,0);
    det_out = this->getDeterminant(out);
  }
  
  mat_out = out;
  out.resize(0,0);
}

void BaseAlgebra::getInverseAndDet(Eigen::SparseMatrix<double> mat_in,Eigen::SparseMatrix<double>& mat_out,double& det_out){
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> inv(mat_in.rows(),mat_in.cols());
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  Eigen::VectorXd idc(mat_in.cols()),c(mat_in.cols());
  double det = 0.0;

  // get inverse
  solver.compute(mat_in);
  for(int i=0;i<mat_in.rows();i++){
    for(int j=0;j<idc.size();j++){
      idc[j] = 0;
    }
    idc[i] = 1;
    //this needs to be done in two steps, otherwise it does not work when using LU
    c = solver.solve(idc);
    inv.col(i) = c;
  }

  // get determinant
  Eigen::VectorXd diag = solver.vectorD();
  diag = solver.vectorD();
  for(int i=0;i<diag.size();i++){
    det += log10( *(diag.data()+i) );
    //    dum = *(diag.data()+i);
    //    if( dum > 1.e-20 ){
    //      det += log10( dum );
    //    }
  }

  det_out = det;
  mat_out = inv.sparseView();
  inv.resize(0,0);
}

double BaseAlgebra::getDeterminant(Eigen::SparseMatrix<double> mat){
  Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
  Eigen::VectorXd diag;
  double dum;
  
  double det = 0.0;
  solver.analyzePattern(mat);
  solver.factorize(mat);
  diag = solver.vectorD();
  for(int i=0;i<diag.size();i++){
    //    std::cout << " " << *(diag.data()+i);
    det += log10( *(diag.data()+i) );
    //    dum = *(diag.data()+i);
    //    if( dum > 1.e-20 ){
    //      det += log10( dum );
    //    }
  }
  //  std::cout << std::endl << std::endl;

  return det;
}


// Class: SmoothAlgebra
//===============================================================================================================
SmoothAlgebra::SmoothAlgebra(SmoothLikelihood* a){
  this->likeModel = a;
}


// Calculate tables and related quantities at initialization
void SmoothAlgebra::setAlgebraInit(ImagePlane* image,BaseSourcePlane* source){
  // Read the data
  Eigen::Map<Eigen::VectorXd> d(image->img,image->Nm);
  this->d = d;
  
  // Read CC and calculate detCs
  Eigen::SparseMatrix<double> CC(image->C.Ti,image->C.Tj);
  CC.reserve(Eigen::VectorXi::Constant(image->C.Ti,8));//overestimating the covariance matrix number of entries per row
  for(int i=0;i<image->C.tri.size();i++){  CC.insert(image->C.tri[i].i,image->C.tri[i].j) = image->C.tri[i].v;  }
  this->C = CC;
  double detC = this->getDeterminant(this->C);
  this->likeModel->terms["detC"] = detC/2.0;
  CC.resize(0,0);

  // Read mask and calculate StCS
  Eigen::SparseMatrix<double> SS(image->S.Ti,image->S.Tj);
  SS.reserve(Eigen::VectorXi::Constant(image->S.Ti,1));//overestimating the mask matrix number of entries per row
  for(int i=0;i<image->S.tri.size();i++){  SS.insert(image->S.tri[i].i,image->S.tri[i].j) = image->S.tri[i].v;  }
  Eigen::SparseMatrix<double> SSt(image->Nm,image->Nm);
  Eigen::SparseMatrix<double> StCS(image->Nm,image->Nm);
  SSt  = SS.transpose();
  StCS = (SSt * this->C * SS);
  this->StCS = StCS;
  SS.resize(0,0);
  SSt.resize(0,0);

  // Read B
  Eigen::SparseMatrix<double> BB(image->B.Ti,image->B.Tj);
  BB.reserve(Eigen::VectorXi::Constant(image->B.Ti,900));//overestimating the number of entries per row of the blurring matrix BB (900 for a 30x30 psf)
  //  for(int i=0;i<imageB.tri.size();i++){  BB.insert(mat->B.tri[i].i,mat->B.tri[i].j) = mat->B.tri[i].v;  }//for row-major
  for(int i=0;i<image->B.tri.size();i++){  BB.insert(image->B.tri[i].j,image->B.tri[i].i) = image->B.tri[i].v;  }//for column-major (swap the i and j indices of tri struct)
  this->B = BB;
  BB.resize(0,0);

  // Read HH and calculate Cs and detCs
  Eigen::SparseMatrix<double> HH(source->Sm,source->Sm);
  HH.reserve(Eigen::VectorXi::Constant(source->Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH row (different number for 1st,2nd order derivative etc)
  for(int i=0;i<source->H.tri.size();i++){  HH.insert(source->H.tri[i].i,source->H.tri[i].j) = source->H.tri[i].v;  }
  double detCs = 0.0;
  this->setAlgebraField(source,HH,this->Cs,detCs);
  this->likeModel->terms["detCs"] = detCs/2.0;
  HH.resize(0,0);

  // Calculate two likelihood terms
  this->likeModel->terms["Nilog2p"] = -(image->lookup.size()*log10(2*M_PI)/2.0); // for some reason i need the outer parenthsis here, otherwise there is a memory problem
  this->likeModel->terms["Nslogl"]  = source->Sm*log10(Nlpar::getValueByName("lambda",this->likeModel->reg))/2.0;
}


// Calculate tables and related quantities at every iteration
void SmoothAlgebra::setAlgebraRuntime(ImagePlane* image,BaseSourcePlane* source){
  // Update Mt and M matrices based on the new lensing matrix L
  Eigen::SparseMatrix<double> LL(source->L.Ti,source->L.Tj);
  LL.reserve(Eigen::VectorXi::Constant(source->L.Ti,4));//setting the non-zero coefficients per row of the lensing matrix (4 for bi-linear interpolation scheme, etc)
  for(int i=0;i<source->L.tri.size();i++){  LL.insert(source->L.tri[i].i,source->L.tri[i].j) = source->L.tri[i].v;  }
  Eigen::SparseMatrix<double> M(image->Nm,source->Sm);
  Eigen::SparseMatrix<double> Mt(source->Sm,image->Nm);
  M        = this->B*LL;
  Mt       = M.transpose();
  this->M  = M;
  this->Mt = Mt;
  LL.resize(0,0);
  M.resize(0,0);
  Mt.resize(0,0);

  // If the source is adaptive, or if the source covariance matrix has changed, calculate Cs and detCs
  if( source->type == "adaptive" || source->sample_reg ){
    Eigen::SparseMatrix<double> HH(source->Sm,source->Sm);
    HH.reserve(Eigen::VectorXi::Constant(source->Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH row (different number for 1st,2nd order derivative etc)
    for(int i=0;i<source->H.tri.size();i++){  HH.insert(source->H.tri[i].i,source->H.tri[i].j) = source->H.tri[i].v;  }

    double detCs = 0.0;
    this->setAlgebraField(source,HH,this->Cs,detCs);
    this->likeModel->terms["detCs"] = detCs/2.0;
    HH.resize(0,0);
  }

  // Calculate the new A matrix
  Eigen::SparseMatrix<double> A(source->Sm,source->Sm);
  A = this->Mt*this->C*this->M + Nlpar::getValueByName("lambda",this->likeModel->reg)*this->Cs;
  this->A = A;
  A.resize(0,0);
}


// Solve for the source and calculate likelihood terms at every iteration
void SmoothAlgebra::solveSource(BaseSourcePlane* source){
  // Get the inverse and the det of A
  Eigen::SparseMatrix<double> inv(source->Sm,source->Sm);
  double detA = 0.0;
  this->getInverseAndDet(this->A,inv,detA);
  this->likeModel->terms["detA"] = -detA/2.0;

  //Get the Most Probable solution for the source
  Eigen::VectorXd s(source->Sm);
  s = inv*(this->Mt*this->C*this->d);
  Eigen::Map<Eigen::VectorXd>(source->src,s.size()) = s;
  inv.resize(0,0);

  //Get the chi-squared and reg terms
  Eigen::VectorXd st = s.transpose();
  Eigen::VectorXd y  = this->d - (this->M*s);
  Eigen::VectorXd yt = y.transpose();
  double chi2        = yt.dot(this->StCS*y);
  double reg         = st.dot(this->Cs*s);
  this->likeModel->terms["chi2"] = -chi2/2.0;
  this->likeModel->terms["reg"]  = -Nlpar::getValueByName("lambda",this->likeModel->reg)*reg/2.0;
}



void SmoothAlgebra::getMockData(ImagePlane* mockdata,BaseSourcePlane* source){
  Eigen::Map<Eigen::VectorXd> s(source->src,source->Sm);
  Eigen::VectorXd tmp = this->M*s;
  Eigen::Map<Eigen::VectorXd>(mockdata->img,tmp.size()) = tmp;
}


void SmoothAlgebra::getSourceErrors(int Sm,double* errors){
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

void SmoothAlgebra::solveDpsi(Pert* pert_mass_model,ImagePlane* image,BaseSourcePlane* source){
  // ATTENTION!!!  >>>  It seems that solving for the inverse od Dpsi does not work  <<<  ATTENTION!!!

  // Get the inverse of the Dpsi matrix (a gradient regularization matrix)
  Eigen::SparseMatrix<double> Dpsi(pert_mass_model->dpsi->Sm,pert_mass_model->dpsi->Sm);
  Dpsi.reserve(Eigen::VectorXi::Constant(pert_mass_model->dpsi->Sm,pert_mass_model->dpsi->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH row (different number for 1st,2nd order derivative etc)
  for(int i=0;i<pert_mass_model->dpsi->H.tri.size();i++){  Dpsi.insert(pert_mass_model->dpsi->H.tri[i].i,pert_mass_model->dpsi->H.tri[i].j) = pert_mass_model->dpsi->H.tri[i].v;  }
  

  Eigen::SparseMatrix<double> Dpsi_inv(source->Sm,source->Sm);
  double det = 0.0;
  this->getInverseAndDet(Dpsi,Dpsi_inv,det);
  Dpsi.resize(0,0);
  
  // Get the vector of residuals multiplied by the inverse of Ds, which is a diagonal matrix
  Eigen::VectorXd dd(image->Nm);
  for(int i=0;i<image->Nm;i++){
    dd[i] = this->likeModel->res->img[i]/(source->Ds.tri[2*i].v + source->Ds.tri[2*i+1].v);
  }

  Eigen::VectorXd dpsi(pert_mass_model->dpsi->Sm);
  dpsi = Dpsi_inv*dd;
  Eigen::Map<Eigen::VectorXd>(pert_mass_model->dpsi->src,dpsi.size()) = dpsi;


  Dpsi_inv.resize(0,0);
}







// Class: PertAlgebra
//===============================================================================================================
PertAlgebra::PertAlgebra(PertLikelihood* a){
  this->likeModel = a;
}

void PertAlgebra::setAlgebraInit(BaseSourcePlane* source,Pert* pert_mass_model){
  // Read source regularization matrix
  Eigen::SparseMatrix<double> HH_s(source->Sm,source->Sm);
  HH_s.reserve(Eigen::VectorXi::Constant(source->Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH_s row (different number for 1st,2nd order derivative etc)
  for(int i=0;i<source->H.tri.size();i++){  HH_s.insert(source->H.tri[i].i,source->H.tri[i].j) = source->H.tri[i].v;  }

  Eigen::SparseMatrix<double> Cs_dum(source->Sm,source->Sm);
  double detCs = 0.0;
  this->setAlgebraField(source,HH_s,Cs_dum,detCs);
  this->likeModel->terms["detCs"] = detCs/2.0;
  this->Cs = Cs_dum;
  Cs_dum.resize(0,0);
  HH_s.resize(0,0);


  // Read perturbations regularization matrix
  Eigen::SparseMatrix<double> HH_dpsi(pert_mass_model->dpsi->Sm,pert_mass_model->dpsi->Sm);
  HH_dpsi.reserve(Eigen::VectorXi::Constant(pert_mass_model->dpsi->Sm,pert_mass_model->dpsi->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH_dpsi row (different number for 1st,2nd order derivative etc)
  for(int i=0;i<pert_mass_model->dpsi->H.tri.size();i++){  HH_dpsi.insert(pert_mass_model->dpsi->H.tri[i].i,pert_mass_model->dpsi->H.tri[i].j) = pert_mass_model->dpsi->H.tri[i].v;  }

  Eigen::SparseMatrix<double> Cp_dum(pert_mass_model->dpsi->Sm,pert_mass_model->dpsi->Sm);
  double detCp = 0.0;
  this->setAlgebraField(source,HH_dpsi,Cp_dum,detCp);
  this->likeModel->terms["detCp"] = detCp/2.0;
  this->Cp = Cp_dum;
  Cp_dum.resize(0,0);
  HH_dpsi.resize(0,0);

  /*
  // Create J "normalizing" matrix and its inverse
  Eigen::SparseMatrix<double> J_dum(source->Sm+pert_mass_model->dpsi->Sm,source->Sm+pert_mass_model->dpsi->Sm);
  constructNormalizingJmatrix(source,pert_mass_model,J_dum,Nlpar::getValueByName("lambda_s",this->likeModel->reg_s),Nlpar::getValueByName("lambda_dpsi",this->likeModel->reg_dpsi));
  this->J = J_dum;
  J_dum.resize(0,0);
  */

  // Create DsDpsi sparse matrix directly based on source->Ds (varying at each call) and image->crosses (fixed at initialization)
  this->constructDsDpsi(source,pert_mass_model);

  this->likeModel->terms["detCd"]   = this->likeModel->smooth_like->terms["detC"];
  this->likeModel->terms["Nslogls"] = source->Sm*log10(Nlpar::getValueByName("lambda_s",this->likeModel->reg_s))/2.0;
  this->likeModel->terms["Nploglp"] = pert_mass_model->dpsi->Sm*log10(Nlpar::getValueByName("lambda_dpsi",this->likeModel->reg_dpsi))/2.0;
  this->likeModel->terms["Nilog2p"] = -(this->likeModel->smooth_like->image->lookup.size()*log10(2*M_PI)/2.0); // for some reason i need the outer parenthsis here, otherwise there is a memory problem
}

void PertAlgebra::constructNormalizingJmatrix(BaseSourcePlane* source,Pert* pert_mass_model,Eigen::SparseMatrix<double>& J_out,double lambda_s,double lambda_dpsi){
  Eigen::SparseMatrix<double> J_dum(source->Sm+pert_mass_model->dpsi->Sm,source->Sm+pert_mass_model->dpsi->Sm);

  for(int i=0;i<source->Sm;i++){
    J_dum.insert(i,i) = 1.0;
  }
  //  double factor = lambda_s/lambda_dpsi;
  double factor = 1.0/lambda_dpsi;

  for(int i=0;i<pert_mass_model->dpsi->Sm;i++){
    J_dum.insert(source->Sm+i,source->Sm+i) = factor;
    //J_dum.insert(source->Sm+i,source->Sm+i) = 1.0;
  }

  /*
  int Si = pert_mass_model->dpsi->Si;
  int Sj = pert_mass_model->dpsi->Sj;
  for(int j=0;j<Sj;j++){
    J_dum.insert( source->Sm+j, source->Sm+j ) = 1.0;
  }
  for(int i=1;i<Si-1;i++){
    //First pixel of each image row: unity
    J_dum.insert( source->Sm+i*Sj, source->Sm+i*Sj ) = 1.0;
   
    //central 2nd derivative in both X and Y directions
    for(int j=1;j<Sj-1;j++){
      J_dum.insert( source->Sm+i*Sj+j, source->Sm+i*Sj+j ) = factor;
    }
    
    //Last pixel of each image row: unity
    J_dum.insert( source->Sm+i*Sj+Sj-1, source->Sm+i*Sj+Sj-1 ) = 1.0;
  }
  for(int j=0;j<Sj;j++){
    J_dum.insert( source->Sm+(Si-1)*Sj+j, source->Sm+(Si-1)*Sj+j ) = 1.0;
  }
  */  

  J_out = J_dum;
  J_dum.resize(0,0);
}

void PertAlgebra::constructDsDpsi(BaseSourcePlane* source,Pert* pert_mass_model){
  Eigen::SparseMatrix<double> DsDpsi_dum(likeModel->smooth_like->image->Nm,pert_mass_model->dpsi->Sm);
  DsDpsi_dum.reserve(Eigen::VectorXi::Constant(pert_mass_model->dpsi->Sm,12));
  int i0,j0,j_index;
  double dsx,dsy;
  int rel_i[12] = {-1,-1,0,0,0,0,1,1,1,1,2,2};
  int rel_j[12] = {0,1,-1,0,1,2,-1,0,1,2,0,1};
  double coeffs[12] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  for(int h=0;h<likeModel->smooth_like->image->Nm;h++){
    dsx = source->Ds.tri[2*h].v;
    dsy = source->Ds.tri[2*h+1].v;

    i0 = likeModel->smooth_like->image->crosses[h]->i0;
    j0 = likeModel->smooth_like->image->crosses[h]->j0;

    coeffs[0]  = likeModel->smooth_like->image->crosses[h]->coeff_y[0]*dsy;
    coeffs[1]  = likeModel->smooth_like->image->crosses[h]->coeff_y[4]*dsy;
    coeffs[2]  = likeModel->smooth_like->image->crosses[h]->coeff_x[0]*dsx;
    coeffs[3]  = likeModel->smooth_like->image->crosses[h]->coeff_y[1]*dsy + likeModel->smooth_like->image->crosses[h]->coeff_x[1]*dsx;
    coeffs[4]  = likeModel->smooth_like->image->crosses[h]->coeff_y[5]*dsy + likeModel->smooth_like->image->crosses[h]->coeff_x[2]*dsx;
    coeffs[5]  = likeModel->smooth_like->image->crosses[h]->coeff_x[3]*dsx;
    coeffs[6]  = likeModel->smooth_like->image->crosses[h]->coeff_x[4]*dsx;
    coeffs[7]  = likeModel->smooth_like->image->crosses[h]->coeff_y[2]*dsy + likeModel->smooth_like->image->crosses[h]->coeff_x[5]*dsx;
    coeffs[8]  = likeModel->smooth_like->image->crosses[h]->coeff_y[6]*dsy + likeModel->smooth_like->image->crosses[h]->coeff_x[6]*dsx;
    coeffs[9]  = likeModel->smooth_like->image->crosses[h]->coeff_x[7]*dsx;
    coeffs[10] = likeModel->smooth_like->image->crosses[h]->coeff_y[3]*dsy;
    coeffs[11] = likeModel->smooth_like->image->crosses[h]->coeff_y[7]*dsy;

    for(int q=0;q<12;q++){
      if( coeffs[q] != 0.0 ){
	j_index = (i0 + rel_i[q])*pert_mass_model->dpsi->Sj + (j0 + rel_j[q]);
	DsDpsi_dum.insert( h, j_index) = coeffs[q];
      }
      coeffs[q] = 0.0; // need to reset the vector "coeffs"
    }
  }
  this->DsDpsi = DsDpsi_dum;
  DsDpsi_dum.resize(0,0);
}


void PertAlgebra::setAlgebraRuntime(BaseSourcePlane* source,Pert* pert_mass_model,SmoothLikelihood* smooth_like){
  // required before constructing M_r
  if( likeModel->name == "pert_iter" ){
    this->constructDsDpsi(source,pert_mass_model);
  }

  Eigen::SparseMatrix<double> M_r_dum(smooth_like->image->Nm,pert_mass_model->dpsi->Sm+source->Sm);
  Eigen::SparseMatrix<double> Mt_r_dum(smooth_like->image->Nm,pert_mass_model->dpsi->Sm+source->Sm);

  Eigen::SparseMatrix<double> block(smooth_like->image->Nm,pert_mass_model->dpsi->Sm+source->Sm);
  block.reserve(Eigen::VectorXi::Constant(smooth_like->image->Nm,30));//overestimating the number of non-zero coefficients per HH row (different number for 1st,2nd order derivative etc)
  for(int i=0;i<source->L.tri.size();i++){  block.insert(source->L.tri[i].i,source->L.tri[i].j) = source->L.tri[i].v;  }
  for(int k=0;k<this->DsDpsi.outerSize();++k){
    for(Eigen::SparseMatrix<double>::InnerIterator it(this->DsDpsi,k);it;++it){
      block.insert(it.row(),source->Sm+it.col()) = -it.value();
    }
  }
  M_r_dum  = (smooth_like->algebra->B) * block;
  Mt_r_dum = M_r_dum.transpose();
  this->M_r = M_r_dum;
  this->Mt_r = Mt_r_dum;
  block.resize(0,0);
  M_r_dum.resize(0,0);
  Mt_r_dum.resize(0,0);


  // The following steps are needed to calculate RtR for covariance based regularizations
  // Calculate source regularization
  if( source->sample_reg && this->likeModel->reg_s_type == "covariance_kernel" ){
    // Read updated source regularization matrix
    Eigen::SparseMatrix<double> HH_s(source->Sm,source->Sm);
    HH_s.reserve(Eigen::VectorXi::Constant(source->Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH_s row (different number for 1st,2nd order derivative etc)
    for(int i=0;i<source->H.tri.size();i++){  HH_s.insert(source->H.tri[i].i,source->H.tri[i].j) = source->H.tri[i].v;  }
    
    Eigen::SparseMatrix<double> Cs(source->Sm,source->Sm);
    double detCs = 0.0;
    this->setAlgebraField(source,HH_s,Cs,detCs);
    this->likeModel->terms["detCs"] = detCs/2.0;
    this->Cs = Cs;
    Cs.resize(0,0);
    HH_s.resize(0,0);
  }
  
  
  // Calculate perturbation regularization
  if( pert_mass_model->dpsi->sample_reg && this->likeModel->reg_dpsi_type == "covariance_kernel" ){
    // Read updated perturbations regularization matrix
    Eigen::SparseMatrix<double> HH_dpsi(pert_mass_model->dpsi->Sm,pert_mass_model->dpsi->Sm);
    HH_dpsi.reserve(Eigen::VectorXi::Constant(pert_mass_model->dpsi->Sm,pert_mass_model->dpsi->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH_dpsi row (different number for 1st,2nd order derivative etc)
    for(int i=0;i<pert_mass_model->dpsi->H.tri.size();i++){  HH_dpsi.insert(pert_mass_model->dpsi->H.tri[i].i,pert_mass_model->dpsi->H.tri[i].j) = pert_mass_model->dpsi->H.tri[i].v;  }

    Eigen::SparseMatrix<double> Cp(pert_mass_model->dpsi->Sm,pert_mass_model->dpsi->Sm);
    double detCp = 0.0;
    this->setAlgebraField(source,HH_dpsi,Cp,detCp);
    this->likeModel->terms["detCp"] = detCp/2.0;
    this->Cp = Cp;
    Cp.resize(0,0);
    HH_dpsi.resize(0,0);
  }
  

  // Calculate RtR
  Eigen::SparseMatrix<double> RtR_dum(source->Sm+pert_mass_model->dpsi->Sm,source->Sm+pert_mass_model->dpsi->Sm);
  double lambda_s = Nlpar::getValueByName("lambda_s",this->likeModel->reg_s);
  for(int k=0;k<this->Cs.outerSize();++k){
    for(Eigen::SparseMatrix<double>::InnerIterator it(this->Cs,k);it;++it){
      RtR_dum.insert(it.row(),it.col()) = lambda_s * it.value();
    }
  }
  double lambda_dpsi = Nlpar::getValueByName("lambda_dpsi",this->likeModel->reg_dpsi);
  for(int k=0;k<this->Cp.outerSize();++k){
    for(Eigen::SparseMatrix<double>::InnerIterator it(this->Cp,k);it;++it){
      RtR_dum.insert(source->Sm+it.row(),source->Sm+it.col()) = lambda_dpsi * it.value();
    }
  }
  this->RtR = RtR_dum;
  RtR_dum.resize(0,0);

  /*
  // Calculate J matrix
  Eigen::SparseMatrix<double> J_dum(source->Sm+pert_mass_model->dpsi->Sm,source->Sm+pert_mass_model->dpsi->Sm);
  constructNormalizingJmatrix(source,pert_mass_model,J_dum,lambda_s,lambda_dpsi);
  this->J = J_dum;
  J_dum.resize(0,0);
  */

  //std::cout << this->RtR << std::endl;
  //std::cout << this->J*this->RtR << std::endl;


  Eigen::SparseMatrix<double> A_r_dum(smooth_like->image->Nm,pert_mass_model->dpsi->Sm+source->Sm);
  //A_r_dum = this->J*(this->Mt_r*smooth_like->algebra->C*this->M_r + this->RtR);
  A_r_dum = this->Mt_r*smooth_like->algebra->C*this->M_r + this->RtR;
  this->A_r = A_r_dum;

  A_r_dum.resize(0,0);
}


void PertAlgebra::solveSourcePert(BaseSourcePlane* source,Pert* pert_mass_model,SmoothLikelihood* smooth_like){
  // Get the inverse and det of A_r
  Eigen::SparseMatrix<double> inv(source->Sm+pert_mass_model->dpsi->Sm,source->Sm+pert_mass_model->dpsi->Sm);
  double detA = 0.0;
  this->getInverseAndDet(this->A_r,inv,detA);
  double lambda_s    = Nlpar::getValueByName("lambda_s",this->likeModel->reg_s);
  double lambda_dpsi = Nlpar::getValueByName("lambda_dpsi",this->likeModel->reg_dpsi);
  //this->likeModel->terms["detA"] = -(detA - pert_mass_model->dpsi->Sm*log10(1.0/lambda_dpsi) )/2.0;
  this->likeModel->terms["detA"] = -detA/2.0;
  
  //Get the Most Probable solution for the source and potential
  Eigen::VectorXd r(source->Sm+pert_mass_model->dpsi->Sm);
  
  //r = (inv*this->J)*(this->Mt_r*this->likeModel->smooth_like->algebra->C*this->likeModel->smooth_like->algebra->d);
  r = inv*(this->Mt_r*this->likeModel->smooth_like->algebra->C*this->likeModel->smooth_like->algebra->d);
  //std::cout << r << std::endl;
  inv.resize(0,0);

  /*
  double* dum = (double*) calloc(source->Sm+pert_mass_model->dpsi->Sm,sizeof(double));
  Eigen::Map<Eigen::VectorXd>(dum,r.size()) = r;
  for(int i=0;i<source->Sm;i++){
    std::cout << " " << dum[i];
  }
  std::cout << std::endl;
  for(int i=0;i<pert_mass_model->dpsi->Sm;i++){
    std::cout << " " << dum[source->Sm+i];
  }
  std::cout << std::endl;
  free(dum);
  */

  //Get the chi-squared term
  Eigen::VectorXd y  = this->likeModel->smooth_like->algebra->d - (this->M_r*r);
  Eigen::VectorXd yt = y.transpose();
  double chi2        = yt.dot(this->likeModel->smooth_like->algebra->StCS*y);
  this->likeModel->terms["chi2"] = -chi2/2.0;

  // Get the regularization term
  Eigen::VectorXd rt = r.transpose();
  double reg         = rt.dot(this->RtR*r);
  this->likeModel->terms["reg"] = -reg/2.0;

  //  Eigen::Map<Eigen::VectorXd>(source->src,source->Sm) = r;
  //  Eigen::Map<Eigen::VectorXd>(pert_mass_model->dpsi->src,dpsi.size()) = dpsi;
  for(int i=0;i<source->Sm;i++){
    source->src[i] = r[i];
    //source->src[i] = -r[i];
  }
  for(int i=0;i<pert_mass_model->dpsi->Sm;i++){
    pert_mass_model->dpsi->src[i] = r[source->Sm+i];
  }
}

