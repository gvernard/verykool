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

  // get determinant of mat_in (NOT its inverse)
  Eigen::VectorXd diag = solver.vectorD();
  diag = solver.vectorD();
  for(int i=0;i<diag.size();i++){
    det += log( *(diag.data()+i) );
    //    dum = *(diag.data()+i);
    //    if( dum > 1.e-20 ){
    //      det += log( dum );
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
    det += log( *(diag.data()+i) );
    //    dum = *(diag.data()+i);
    //    if( dum > 1.e-20 ){
    //      det += log( dum );
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
  
  // Read CC and calculate detC
  if( image->noise_flag == "correlated" ){
    Eigen::SparseMatrix<double> CC(image->C.Ti,image->C.Tj);
    CC.reserve(Eigen::VectorXi::Constant(image->C.Ti,8));//overestimating the covariance matrix number of entries per row
    for(int i=0;i<image->C.tri.size();i++){  CC.insert(image->C.tri[i].i,image->C.tri[i].j) = image->C.tri[i].v;  }
    double detC = 0.0;
    this->getInverseAndDet(CC,this->C_inv,detC);
    this->likeModel->terms["detCd"] = -detC/2.0;
    CC.resize(0,0);
  } else {
    // Data covariance matrix is diagonal
    Eigen::SparseMatrix<double> Cdum(image->C.Ti,image->C.Tj);
    Cdum.reserve(Eigen::VectorXi::Constant(image->C.Ti,8));//overestimating the covariance matrix number of entries per row
    double det_Cd = 0.0;
    for(int i=0;i<image->C.tri.size();i++){
      Cdum.insert(image->C.tri[i].i,image->C.tri[i].j) = 1.0/image->C.tri[i].v;
      det_Cd += log(image->C.tri[i].v);
    }
    this->C_inv = Cdum;
    this->likeModel->terms["detCd"] = -det_Cd/2.0;
    Cdum.resize(0,0);
  }

  // Read mask and calculate StCS
  Eigen::SparseMatrix<double> SS(image->S.Ti,image->S.Tj);
  SS.reserve(Eigen::VectorXi::Constant(image->S.Ti,1));//overestimating the mask matrix number of entries per row
  for(int i=0;i<image->S.tri.size();i++){  SS.insert(image->S.tri[i].i,image->S.tri[i].j) = image->S.tri[i].v;  }
  Eigen::SparseMatrix<double> SSt(image->Nm,image->Nm);
  Eigen::SparseMatrix<double> StCS(image->Nm,image->Nm);
  SSt  = SS.transpose();
  StCS = (SSt * this->C_inv * SS);
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

  // Read HH and calculate Cs_inv and detCs_inv
  Eigen::SparseMatrix<double> HH(source->H.Ti,source->H.Tj);
  HH.reserve(Eigen::VectorXi::Constant(source->Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH row (different number for 1st,2nd order derivative etc)
  for(int i=0;i<source->H.tri.size();i++){  HH.insert(source->H.tri[i].i,source->H.tri[i].j) = source->H.tri[i].v;  }
  double detCs = 0.0;
  this->setAlgebraField(source,HH,this->Cs_inv,detCs);
  this->likeModel->terms["detCs"] = -detCs/2.0;
  HH.resize(0,0);

  // Calculate the data related constant likelihood term
  this->likeModel->terms["Nilog2pi"] = -(image->Nmask*log(2*M_PI)/2.0); // for some reason i need the outer parenthsis here, otherwise there is a memory problem

  // Calculate the source related constant likelihood term (dependent on lambda)
  if( this->likeModel->source->reg == "curvature_in_identity_out" || this->likeModel->source->reg == "covariance_kernel_in_identity_out" ){
    Eigen::SparseMatrix<double> block_dum(source->Sm,source->Sm);
    block_dum.reserve(Eigen::VectorXi::Constant(source->Sm,10));
    for(int i=0;i<source->Sm;i++){
      if( source->mask_vertices[i] == 0 ){
	block_dum.insert(i,i) = source->lambda_out[i];
      }
    }
    this->block = block_dum;
    block_dum.resize(0,0);

    this->likeModel->terms["Noutlog2pi"] = -(source->Sm-source->Smask)*log(2*M_PI)/2.0;
    this->likeModel->terms["detLout"] = -source->lambda_out_sum/2.0;
  }
  this->likeModel->terms["Nslogl"]    =  source->Sm*log(Nlpar::getValueByName("lambda",this->likeModel->reg))/2.0;
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

  // If the source is adaptive, or if the source covariance matrix has changed, calculate Cs_inv and detCs
  if( source->type == "adaptive" || source->sample_reg ){
    Eigen::SparseMatrix<double> HH(source->H.Ti,source->H.Tj);
    HH.reserve(Eigen::VectorXi::Constant(source->Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH row (different number for 1st,2nd order derivative etc)
    for(int i=0;i<source->H.tri.size();i++){  HH.insert(source->H.tri[i].i,source->H.tri[i].j) = source->H.tri[i].v;  }
    double detCs = 0.0;
    this->setAlgebraField(source,HH,this->Cs_inv,detCs);
    this->likeModel->terms["detCs"] = -detCs/2.0;
    HH.resize(0,0);
  }

  Eigen::SparseMatrix<double> A_dum(source->Sm,source->Sm);
  if( this->likeModel->source->reg == "curvature_in_identity_out" || this->likeModel->source->reg == "covariance_kernel_in_identity_out" ){
    /*
    // construct block regularization matrix
    Eigen::SparseMatrix<double> block_dum(source->Sm,source->Sm);
    block_dum.reserve(Eigen::VectorXi::Constant(source->Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH row (different number for 1st,2nd order derivative etc)
    for(int i=0;i<source->Sm;i++){
      if( source->mask_vertices[i] == 0 ){
	block_dum.insert(i,i) = source->lambda_out[i];
      }
    }
    double lambda_in = Nlpar::getValueByName("lambda",this->likeModel->reg);
    for(int i=0;i<source->H.tri.size();i++){
      block_dum.insert( source->in_total[source->H.tri[i].i], source->in_total[source->H.tri[i].j]) = lambda_in * source->H.tri[i].v;
    }
    this->block = block_dum;
    block_dum.resize(0,0);
    A_dum = this->Mt*this->StCS*this->M + this->block;
    */
    A_dum = this->Mt*this->StCS*this->M + Nlpar::getValueByName("lambda",this->likeModel->reg)*this->Cs_inv + this->block;

  } else {
    //  A = this->Mt*this->C*this->M + Nlpar::getValueByName("lambda",this->likeModel->reg)*this->Cs;
    A_dum = this->Mt*this->StCS*this->M + Nlpar::getValueByName("lambda",this->likeModel->reg)*this->Cs_inv;
  }


  this->A = A_dum;
  A_dum.resize(0,0);
}


// Solve for the source and calculate likelihood terms at every iteration
void SmoothAlgebra::solveSource(BaseSourcePlane* source){
  // Get the inverse and the det of A
  Eigen::SparseMatrix<double> A_inv(source->Sm,source->Sm);
  double detA = 0.0;
  this->getInverseAndDet(this->A,A_inv,detA);
  this->likeModel->terms["detA"] = -detA/2.0;

  //Get the Most Probable solution for the source
  Eigen::VectorXd s(source->Sm);
  //  s = inv*(this->Mt*this->C*this->d);
  s = A_inv*(this->Mt*this->StCS*this->d);
  Eigen::Map<Eigen::VectorXd>(source->src,s.size()) = s;
  A_inv.resize(0,0);

  /*
  // Read in a matching source
  FILE* fh2 = fopen("/net/argo/data/users/gvernard/RESULTS/VKL_MODELS/test_modelling_perturbations/test_6_nopert/base_run/output/vkl_source_irregular.dat","r");
  float dumaa,dumbb,dumcc;
  for(int i=0;i<source->Sm;i++){
    fscanf(fh2,"%f%f%f\n",&dumaa,&dumbb,&dumcc);
    source->src[i] = dumaa;
  }  
  fclose(fh2);

  Eigen::Map<Eigen::VectorXd> dum1(source->src,source->Sm);
  Eigen::VectorXd dum = this->Cs*dum1;
  //  Eigen::VectorXd dum = dum1;
  dum[0] = 10.0;
  s = dum;
  Eigen::Map<Eigen::VectorXd>(source->src,s.size()) = s;
  */


  /*
  // Calculating the Wiener filter
  std::cout << "Starting Wiener..." << std::endl;
  Eigen::SparseMatrix<double> Cs(source->Sm,source->Sm);
  double dum = 0;
  this->getInverseAndDet(this->Cs_inv,Cs,dum);
  std::cout << "Cs done" << std::endl;

  Eigen::SparseMatrix<double> Cd(this->likeModel->image->C.Ti,this->likeModel->image->C.Tj);
  Cd.reserve(Eigen::VectorXi::Constant(this->likeModel->image->C.Ti,8));//overestimating the covariance matrix number of entries per row
  for(int i=0;i<this->likeModel->image->C.tri.size();i++){  Cd.insert(this->likeModel->image->C.tri[i].i,this->likeModel->image->C.tri[i].j) = this->likeModel->image->C.tri[i].v;  }
  std::cout << "Cd done" << std::endl;

  Eigen::SparseMatrix<double> Q(this->likeModel->image->Nm,this->likeModel->image->Nm);
  Eigen::SparseMatrix<double> mat(this->likeModel->image->Nm,source->Sm);
  mat = this->M*Cs;
  Q = mat*this->Mt;
  Q += Cd;
  //  Q = this->M*Cs*this->Mt + Cd;
  std::cout << "Q done" << std::endl;

  Eigen::SparseMatrix<double> Q_inv(this->likeModel->image->Nm,this->likeModel->image->Nm);
  this->getInverseAndDet(Q,Q_inv,dum);
  std::cout << "Q_inv done" << std::endl;

  Eigen::VectorXd sw(source->Sm);
  sw = Cs*this->Mt*Q_inv*this->d;
  std::cout << "sw done" << std::endl;

  std::cout << "Wiener OK!!!" << std::endl;
  */



  // Get the chi-squared term
  Eigen::VectorXd y  = this->d - (this->M*s);
  Eigen::VectorXd yt = y.transpose();
  double chi2        = yt.dot(this->StCS*y);
  this->likeModel->terms["chi2"] = -chi2/2.0;

  // Get the regularization term
  Eigen::VectorXd st = s.transpose();
  if( this->likeModel->source->reg == "curvature_in_identity_out" || this->likeModel->source->reg == "covariance_kernel_in_identity_out" ){
    //    double reg                    =  st.dot(this->block*s);
    double reg1                   =  st.dot(this->Cs_inv*s);
    double reg2                   =  st.dot(this->block*s);
    this->likeModel->terms["reg"] = -Nlpar::getValueByName("lambda",this->likeModel->reg)*reg1/2.0;
    this->likeModel->terms["reg_out"] = -reg2/2.0;
  } else {
    double reg                     = st.dot(this->Cs_inv*s);
    this->likeModel->terms["reg"]  = -Nlpar::getValueByName("lambda",this->likeModel->reg)*reg/2.0;
  }
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
  Eigen::SparseMatrix<double> HH_s(source->H.Ti,source->H.Tj);
  HH_s.reserve(Eigen::VectorXi::Constant(source->Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH_s row (different number for 1st,2nd order derivative etc)
  for(int i=0;i<source->H.tri.size();i++){  HH_s.insert(source->H.tri[i].i,source->H.tri[i].j) = source->H.tri[i].v;  }
  double detCs = 0.0;
  this->setAlgebraField(source,HH_s,this->Cs_inv,detCs);
  this->likeModel->terms["detCs"] = -detCs/2.0;
  HH_s.resize(0,0);

  // Read perturbations regularization matrix
  Eigen::SparseMatrix<double> HH_dpsi(pert_mass_model->dpsi->Sm,pert_mass_model->dpsi->Sm);
  HH_dpsi.reserve(Eigen::VectorXi::Constant(pert_mass_model->dpsi->Sm,pert_mass_model->dpsi->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH_dpsi row (different number for 1st,2nd order derivative etc)
  for(int i=0;i<pert_mass_model->dpsi->H.tri.size();i++){  HH_dpsi.insert(pert_mass_model->dpsi->H.tri[i].i,pert_mass_model->dpsi->H.tri[i].j) = pert_mass_model->dpsi->H.tri[i].v;  }

  double detCp = 0.0;
  this->setAlgebraField(source,HH_dpsi,this->Cp_inv,detCp);
  this->likeModel->terms["detCp"] = -detCp/2.0;
  HH_dpsi.resize(0,0);

  /*
  // Create J "normalizing" matrix and its inverse
  Eigen::SparseMatrix<double> J_dum(source->Sm+pert_mass_model->dpsi->Sm,source->Sm+pert_mass_model->dpsi->Sm);
  constructNormalizingJmatrix(source,pert_mass_model,J_dum,Nlpar::getValueByName("lambda_s",this->likeModel->reg_s),Nlpar::getValueByName("lambda_dpsi",this->likeModel->reg_dpsi));
  this->J = J_dum;
  J_dum.resize(0,0);
  */

  // Create DsDpsi sparse matrix directly based on source->Ds (varying at each call) and image->crosses (fixed at initialization)
  this->constructDsDpsi(this->likeModel->smooth_like->image,source,pert_mass_model);

  this->likeModel->terms["detCd"]    = this->likeModel->smooth_like->terms["detCd"];
  this->likeModel->terms["Nilog2pi"] = -(this->likeModel->smooth_like->image->Nmask*log(2*M_PI)/2.0); // for some reason i need the outer parenthsis here, otherwise there is a memory problem

  if( source->reg == "curvature_in_identity_out" || source->reg == "covariance_kernel_in_identity_out" ){
    Eigen::SparseMatrix<double> block_dum(source->Sm+pert_mass_model->dpsi->Sm,source->Sm+pert_mass_model->dpsi->Sm);
    block_dum.reserve(Eigen::VectorXi::Constant(source->Sm+pert_mass_model->dpsi->Sm,10));
    for(int i=0;i<source->Sm;i++){
      if( source->mask_vertices[i] == 0 ){
	block_dum.insert(i,i) = source->lambda_out[i];
      }
    }
    this->block_s = block_dum;
    block_dum.resize(0,0);

    this->likeModel->terms["Noutlog2pi_s"] = -(source->Sm-source->Smask)*log(2*M_PI)/2.0;
    this->likeModel->terms["detLout_s"] = -source->lambda_out_sum/2.0;
  }
  this->likeModel->terms["Nslogl_s"] = source->Sm*log(Nlpar::getValueByName("lambda_s",this->likeModel->reg_s))/2.0;

  if( pert_mass_model->dpsi->reg == "curvature_in_identity_out" || pert_mass_model->dpsi->reg == "covariance_kernel_in_identity_out" ){
    Eigen::SparseMatrix<double> block_dum(source->Sm+pert_mass_model->dpsi->Sm,source->Sm+pert_mass_model->dpsi->Sm);
    block_dum.reserve(Eigen::VectorXi::Constant(source->Sm+pert_mass_model->dpsi->Sm,10));
    for(int i=0;i<pert_mass_model->dpsi->Sm;i++){
      if( pert_mass_model->dpsi->mask_vertices[i] == 0 ){
	block_dum.insert(source->Sm+i,source->Sm+i) = pert_mass_model->dpsi->lambda_out[i];
      }
    }
    this->block_p = block_dum;
    block_dum.resize(0,0);

    this->likeModel->terms["Noutlog2pi_p"] = -(pert_mass_model->dpsi->Sm-pert_mass_model->dpsi->Smask)*log(2*M_PI)/2.0;
    this->likeModel->terms["detLout_p"] = -pert_mass_model->dpsi->lambda_out_sum/2.0;
  }
  this->likeModel->terms["Nplogl_p"] = pert_mass_model->dpsi->Sm*log(Nlpar::getValueByName("lambda_dpsi",this->likeModel->reg_dpsi))/2.0;
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

void PertAlgebra::constructDsDpsi(ImagePlane* image,BaseSourcePlane* source,Pert* pert_mass_model){
  Eigen::SparseMatrix<double> DsDpsi_dum(image->Nm,pert_mass_model->dpsi->Sm);
  DsDpsi_dum.reserve(Eigen::VectorXi::Constant(pert_mass_model->dpsi->Sm,24));
  int i0,j0,j_index;
  double dsx,dsy;
  int rel_i[12] = {-1,-1,0,0,0,0,1,1,1,1,2,2};
  int rel_j[12] = {0,1,-1,0,1,2,-1,0,1,2,0,1};
  double coeffs[12] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  
  for(int h=0;h<image->Nm;h++){
    dsx = source->Ds.tri[2*h].v;
    dsy = source->Ds.tri[2*h+1].v;

    i0 = image->crosses[h]->i0;
    j0 = image->crosses[h]->j0;

    coeffs[0]  = image->crosses[h]->coeff_y[0]*dsy;
    coeffs[1]  = image->crosses[h]->coeff_y[4]*dsy;
    coeffs[2]  = image->crosses[h]->coeff_x[0]*dsx;
    coeffs[3]  = image->crosses[h]->coeff_y[1]*dsy + image->crosses[h]->coeff_x[1]*dsx;
    coeffs[4]  = image->crosses[h]->coeff_y[5]*dsy + image->crosses[h]->coeff_x[2]*dsx;
    coeffs[5]  = image->crosses[h]->coeff_x[3]*dsx;
    coeffs[6]  = image->crosses[h]->coeff_x[4]*dsx;
    coeffs[7]  = image->crosses[h]->coeff_y[2]*dsy + image->crosses[h]->coeff_x[5]*dsx;
    coeffs[8]  = image->crosses[h]->coeff_y[6]*dsy + image->crosses[h]->coeff_x[6]*dsx;
    coeffs[9]  = image->crosses[h]->coeff_x[7]*dsx;
    coeffs[10] = image->crosses[h]->coeff_y[3]*dsy;
    coeffs[11] = image->crosses[h]->coeff_y[7]*dsy;
    
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


void PertAlgebra::setAlgebraRuntime(BaseSourcePlane* source,Pert* pert_mass_model){
  // required before constructing M_r
  if( likeModel->name == "pert_iter" ){
    this->constructDsDpsi(this->likeModel->smooth_like->image,source,pert_mass_model);
  }

  Eigen::SparseMatrix<double> M_r_dum(this->likeModel->smooth_like->image->Nm,pert_mass_model->dpsi->Sm+source->Sm);
  Eigen::SparseMatrix<double> Mt_r_dum(this->likeModel->smooth_like->image->Nm,pert_mass_model->dpsi->Sm+source->Sm);

  Eigen::SparseMatrix<double> block(this->likeModel->smooth_like->image->Nm,pert_mass_model->dpsi->Sm+source->Sm);
  block.reserve(Eigen::VectorXi::Constant(this->likeModel->smooth_like->image->Nm,30));//overestimating the number of non-zero coefficients per HH row (different number for 1st,2nd order derivative etc)
  for(int i=0;i<source->L.tri.size();i++){  block.insert(source->L.tri[i].i,source->L.tri[i].j) = source->L.tri[i].v;  }
  for(int k=0;k<this->DsDpsi.outerSize();++k){
    for(Eigen::SparseMatrix<double>::InnerIterator it(this->DsDpsi,k);it;++it){
      block.insert(it.row(),source->Sm+it.col()) = -it.value();
    }
  }
  M_r_dum  = (this->likeModel->smooth_like->algebra->B) * block;
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
    Eigen::SparseMatrix<double> HH_s(source->H.Ti,source->H.Tj);
    HH_s.reserve(Eigen::VectorXi::Constant(source->Sm,source->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH_s row (different number for 1st,2nd order derivative etc)
    for(int i=0;i<source->H.tri.size();i++){  HH_s.insert(source->H.tri[i].i,source->H.tri[i].j) = source->H.tri[i].v;  }
    double detCs = 0.0;
    this->setAlgebraField(source,HH_s,this->Cs_inv,detCs);
    this->likeModel->terms["detCs"] = -detCs/2.0;
    HH_s.resize(0,0);
  }
  
  
  // Calculate perturbation regularization
  if( pert_mass_model->dpsi->sample_reg && this->likeModel->reg_dpsi_type == "covariance_kernel" ){
    // Read updated perturbations regularization matrix
    Eigen::SparseMatrix<double> HH_dpsi(pert_mass_model->dpsi->Sm,pert_mass_model->dpsi->Sm);
    HH_dpsi.reserve(Eigen::VectorXi::Constant(pert_mass_model->dpsi->Sm,pert_mass_model->dpsi->eigenSparseMemoryAllocForH));//overestimating the number of non-zero coefficients per HH_dpsi row (different number for 1st,2nd order derivative etc)
    for(int i=0;i<pert_mass_model->dpsi->H.tri.size();i++){  HH_dpsi.insert(pert_mass_model->dpsi->H.tri[i].i,pert_mass_model->dpsi->H.tri[i].j) = pert_mass_model->dpsi->H.tri[i].v;  }
    double detCp = 0.0;
    this->setAlgebraField(source,HH_dpsi,this->Cp_inv,detCp);
    this->likeModel->terms["detCp"] = -detCp/2.0;
    HH_dpsi.resize(0,0);
  }
  

  // Calculate RtR
  Eigen::SparseMatrix<double> RtR_dum(source->Sm+pert_mass_model->dpsi->Sm,source->Sm+pert_mass_model->dpsi->Sm);
  double lambda_s = Nlpar::getValueByName("lambda_s",this->likeModel->reg_s);
  for(int k=0;k<this->Cs_inv.outerSize();++k){
    for(Eigen::SparseMatrix<double>::InnerIterator it(this->Cs_inv,k);it;++it){
      RtR_dum.insert(it.row(),it.col()) = lambda_s * it.value();
    }
  }
  if( source->reg == "curvature_in_identity_out" || source->reg == "covariance_kernel_in_identity_out" ){
    RtR_dum += this->block_s;
  }
  double lambda_dpsi = Nlpar::getValueByName("lambda_dpsi",this->likeModel->reg_dpsi);
  for(int k=0;k<this->Cp_inv.outerSize();++k){
    for(Eigen::SparseMatrix<double>::InnerIterator it(this->Cp_inv,k);it;++it){
      RtR_dum.insert(source->Sm+it.row(),source->Sm+it.col()) = lambda_dpsi * it.value();
    }
  }
  if( pert_mass_model->dpsi->reg == "curvature_in_identity_out" || pert_mass_model->dpsi->reg == "covariance_kernel_in_identity_out" ){
    RtR_dum += this->block_p;
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


  Eigen::SparseMatrix<double> A_r_dum(this->likeModel->smooth_like->image->Nm,pert_mass_model->dpsi->Sm+source->Sm);
  //A_r_dum = this->J*(this->Mt_r*this->likeModel->smooth_like->algebra->C*this->M_r + this->RtR);
  //A_r_dum = this->Mt_r*this->likeModel->smooth_like->algebra->C*this->M_r + this->RtR;
  A_r_dum = this->Mt_r*this->likeModel->smooth_like->algebra->StCS*this->M_r + this->RtR;
  this->A_r = A_r_dum;

  A_r_dum.resize(0,0);
}


void PertAlgebra::solveSourcePert(BaseSourcePlane* source,Pert* pert_mass_model){
  // Get the inverse and det of A_r
  Eigen::SparseMatrix<double> A_inv(source->Sm+pert_mass_model->dpsi->Sm,source->Sm+pert_mass_model->dpsi->Sm);
  double detA = 0.0;
  this->getInverseAndDet(this->A_r,A_inv,detA);
  double lambda_s    = Nlpar::getValueByName("lambda_s",this->likeModel->reg_s);
  double lambda_dpsi = Nlpar::getValueByName("lambda_dpsi",this->likeModel->reg_dpsi);
  //this->likeModel->terms["detA"] = -(detA - pert_mass_model->dpsi->Sm*log(1.0/lambda_dpsi) )/2.0;
  this->likeModel->terms["detA"] = -detA/2.0;
  
  // Get the Most Probable solution for the source and potential
  Eigen::VectorXd r(source->Sm+pert_mass_model->dpsi->Sm);
  
  //r = (inv*this->J)*(this->Mt_r*this->likeModel->smooth_like->algebra->C*this->likeModel->smooth_like->algebra->d);
  //r = inv*(this->Mt_r*this->likeModel->smooth_like->algebra->C*this->likeModel->smooth_like->algebra->d);
  r = A_inv*(this->Mt_r*this->likeModel->smooth_like->algebra->StCS*this->likeModel->smooth_like->algebra->d);
  //std::cout << r << std::endl;
  A_inv.resize(0,0);


  // Set the edge values of dpsi to zero
  for(int j=0;j<pert_mass_model->dpsi->Sj;j++){
    r[source->Sm+j] = 0.0;
  }
  for(int i=1;i<pert_mass_model->dpsi->Si-1;i++){
    r[source->Sm+i*pert_mass_model->dpsi->Sj] = 0.0;
    r[source->Sm+i*pert_mass_model->dpsi->Sj+pert_mass_model->dpsi->Sj-1] = 0.0;
  }
  for(int j=0;j<pert_mass_model->dpsi->Sj;j++){
    r[source->Sm+(pert_mass_model->dpsi->Si-1)*pert_mass_model->dpsi->Sj+j] = 0.0;
  }


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

  // Get the chi-squared term
  Eigen::VectorXd y  = this->likeModel->smooth_like->algebra->d - (this->M_r*r);
  Eigen::VectorXd yt = y.transpose();
  double chi2        = yt.dot(this->likeModel->smooth_like->algebra->StCS*y);
  this->likeModel->terms["chi2"] = -chi2/2.0;

  // Get the regularization term (adding the block matrix is already incorporated in RtR)
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


void PertAlgebra::solvePerturbationsResiduals(std::string output,Pert* pert_mass_model){
  Eigen::Map<Eigen::VectorXd> dpsi(pert_mass_model->dpsi->src,pert_mass_model->dpsi->Sm);
  //  Eigen::VectorXd res = -this->likeModel->smooth_like->algebra->B * this->DsDpsi * dpsi;
  Eigen::VectorXd res = -this->likeModel->smooth_like->algebra->B * this->DsDpsi * dpsi;

  ImagePlane res_img(this->likeModel->smooth_like->image->Ni,this->likeModel->smooth_like->image->Nj,this->likeModel->smooth_like->image->width,this->likeModel->smooth_like->image->height);
  
  for(int i=0;i<res_img.Ni;i++){
    for(int j=0;j<res_img.Nj;j++){
      res_img.img[i*res_img.Nj+j] = res[i*res_img.Nj+j];
      this->likeModel->smooth_like->algebra->d[i*res_img.Nj+j] = res[i*res_img.Nj+j];      
      this->likeModel->smooth_like->image->img[i*res_img.Nj+j] = res[i*res_img.Nj+j];
    }
  }

  res_img.writeImage(output + "res_dpsi.fits");
}

void PertAlgebra::getMockData(ImagePlane* mockdata,BaseSourcePlane* source,Eigen::SparseMatrix<double> B){
  Eigen::Map<Eigen::VectorXd> s(source->src,source->Sm);

  Eigen::SparseMatrix<double> LL(source->L.Ti,source->L.Tj);
  LL.reserve(Eigen::VectorXi::Constant(source->L.Ti,4));//setting the non-zero coefficients per row of the lensing matrix (4 for bi-linear interpolation scheme, etc)
  for(int i=0;i<source->L.tri.size();i++){  LL.insert(source->L.tri[i].i,source->L.tri[i].j) = source->L.tri[i].v;  }

  Eigen::VectorXd tmp = B*LL*s;
  Eigen::Map<Eigen::VectorXd>(mockdata->img,tmp.size()) = tmp;
}
