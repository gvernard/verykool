#include "sourcePlane.hpp"

//Abstract class: BaseSourcePlane
//===============================================================================================================
void BaseSourcePlane::normalize(){
  double sum = 0.;
  for(int j=0;j<this->Sm;j++){
      sum += this->src[j];
  }
  
  for(int j=0;j<this->Sm;j++){
    this->src[j] /= sum;
  }
}
//===============================================================================================================


//Derived class from BaseSourcePlane: FixedSource
//===============================================================================================================

//virtual
FixedSource::FixedSource(int i,int j,double size,std::string reg_scheme){
  type = "fixed";
  reg  = reg_scheme;
  Si   = i;
  Sj   = j;
  Sm   = Si*Sj;

  src = (double*) calloc(Sm,sizeof(double));
  x   = (double*) calloc(Sm,sizeof(double));
  y   = (double*) calloc(Sm,sizeof(double));

  std::map<std::string,std::string> pars;
  pars["size"] = std::to_string(size);

  this->setGrid(pars);
  this->boundPolygon();
}

//virtual
FixedSource::~FixedSource(){}

//non-virtual
void FixedSource::setGrid(std::map<std::string,std::string> pars){
  double size = stof(pars["size"]);

  this->xmax =  size;
  this->xmin = -size;
  this->ymax =  size;
  this->ymin = -size;

  double dx = (this->xmax - this->xmin)/this->Sj;
  double dy = (this->ymax - this->ymin)/this->Si;

  this->xmax -= dx;
  this->ymax -= dy;

  for(int j=0;j<this->Sj;j++){
    for(int i=0;i<this->Si;i++){
      this->x[j*this->Si+i] = this->xmin + i*dx;
      this->y[j*this->Si+i] = this->ymin + (this->Sj-1-j)*dy;//reflect y-axis
    }
  }
}

//non-virtual
void FixedSource::boundPolygon(){
  std::vector<double> tx {this->xmin,this->xmin,this->xmax,this->xmax,this->xmin};
  std::vector<double> ty {this->ymin,this->ymax,this->ymax,this->ymin,this->ymin};
  this->bound_x  = tx;
  this->bound_y  = ty;
}

//non-virtual
bool FixedSource::pointInPolygon(double x,double y){
  if( x <= this->xmin || x >= this->xmax || y <= this->ymin || y >= this->ymax ){
    return false;
  } else {
    return true;
  }
}

//virtual
void FixedSource::constructL(ImagePlane* image,CollectionMassModels* mycollection,mytable* L){
  int Nm = image->Nm;
  int Si = this->Si;
  int Sj = this->Sj;

  long excluded = 0;
  long included = 0;
  double xp     = 0;
  double yp     = 0;
  double dx     = (this->xmax - this->xmin)/(Si - 1);
  double dy     = (this->ymax - this->ymin)/(Sj - 1);
  double norm   = 1./(dx*dy);
  int ic        = 0;
  int jc        = 0;
  double g1=0,g2=0,g3=0,g4=0;
  std::vector<mytriplet> tmp;//need to make sure that the L triplet vector is a new one


  for(int i=0;i<Nm;i++){

    mycollection->all_defl(image->x[i],image->y[i],xp,yp);
    
    if( this->pointInPolygon(xp,yp) ){
      included++;
      
      //	this->interpolateSource(i,j,xdefl,ydefl,&L);
      
      //Indices corresponding to the top left pixel
      jc = (int) floor((xp - this->xmin)/dx);
      ic = (int) floor((yp - this->ymin)/dy);
      
      //Now interpolate between neighbouring pixels and add entries to L matrix.
      //The interpolation function could return an array of column indices in the row of the L matrix, and the corresponding weights.
      //The following is for bi-linear interpolation:
      g1 = (ic+1)*dy - (yp - this->ymin);
      g2 = (jc+1)*dx - (xp - this->xmin);
      g3 = (yp - this->ymin) - ic*dy;
      g4 = (xp - this->xmin) - jc*dx;
      
      tmp.push_back({i, (Si-2-ic)*Sj   + (jc),       g1*g2*norm});
      tmp.push_back({i, (Si-2-ic)*Sj   + (jc)+1,     g1*g4*norm});
      tmp.push_back({i, (Si-2-ic+1)*Sj + (jc),       g3*g2*norm});
      tmp.push_back({i, (Si-2-ic+1)*Sj + (jc)+1,     g3*g4*norm});
    } else {
      excluded++;
    } 

  }

  L->tri.swap(tmp);
  
  //  std::cout << excluded << " image pixels lie outside the source grid." << std::endl;
  //  std::cout << included << " image pixels lie inside the source grid." << std::endl;
}


//virtual
void FixedSource::constructH(mytable* H){
  std::vector<mytriplet> tmp;//need to make sure that the H triplet vector is a new one

  if( this->reg == "identity" ){//---------------------------> zero order
    
    for(int i=0;i<H->Ti;i++){
      tmp.push_back({i,i,1});
    }

  } else if ( this->reg == "gradient" ){//-------------------> first order
    
    int Si = this->Si;
    int Sj = this->Sj;
    double ddx = 1./fabs(this->x[1]  - this->x[0]);
    double ddy = 1./fabs(this->y[Sj] - this->y[0]);

    for(int i=1;i<Si-1;i++){
      for(int j=1;j<Sj-1;j++){
	tmp.push_back({  j*Si+i,  (j-1)*Si+i,  ddy});
	tmp.push_back({  j*Si+i,    j*Si+i-1,  ddx});
	tmp.push_back({  j*Si+i,    j*Si+i+1,  ddx});
	tmp.push_back({  j*Si+i,  (j+1)*Si+i,  ddy});
      }
    }

  } else if ( this->reg == "curvature" ){//-------------------> second order

    int Si = this->Si;
    int Sj = this->Sj;
    double dx = fabs(this->x[1]  - this->x[0]);
    double dy = fabs(this->y[Sj] - this->y[0]);
    double ddx2 =  1.0/(dx*dx);
    double ddy2 =  1.0/(dy*dy);
    double fac  = -2.0*(ddx2+ddy2);
    double fac2 =  2.0*(ddy2-ddx2);
    double fac3 =  4.0*(ddx2+ddy2);
    //    std::cout << fac << " " << fac2 << " " << fac3 << std::endl;

    //First pixel of first image row:  forward 2nd derivative in X, forward 2nd derivative in Y, both of 2nd order accuracy


    for(int j=0;j<Sj;j++){
      tmp.push_back({  j,  j,           fac   });
    }

    /*
    H->tri.push_back({  0,  0,         fac3   });
    H->tri.push_back({  0,  1*Sj,  -5.*ddy2   });
    H->tri.push_back({  0,  2*Sj,   4.*ddy2   });
    H->tri.push_back({  0,  3*Sj,  -1.*ddy2   });
    H->tri.push_back({  0,  1,     -5.*ddx2   });
    H->tri.push_back({  0,  2,      4.*ddx2   });
    H->tri.push_back({  0,  3,     -1.*ddx2   });

    //First row of image pixels: central 2nd derivative in X, forward 2nd derivative in Y, both of 2nd order accuracy
    for(int j=1;j<Sj-1;j++){
      H->tri.push_back({  j,  j,           fac2   });
      H->tri.push_back({  j,  1*Sj+j,  -5.*ddy2   });
      H->tri.push_back({  j,  2*Sj+j,   4.*ddy2   });
      H->tri.push_back({  j,  3*Sj+j,  -1.*ddy2   });
      H->tri.push_back({  j,  j-1,      1.*ddx2   });
      H->tri.push_back({  j,  j+1,      1.*ddx2   });
    }

    //Last pixel of first image row:  backward 2nd derivative in X, forward 2nd derivative in Y, both of 2nd order accuracy
    H->tri.push_back({  Sj-1,  Sj-1,          fac3   });
    H->tri.push_back({  Sj-1,  2*Sj-1,    -5.*ddy2   });
    H->tri.push_back({  Sj-1,  3*Sj-1,     4.*ddy2   });
    H->tri.push_back({  Sj-1,  4*Sj-1,    -1.*ddy2   });
    H->tri.push_back({  Sj-1,  Sj-2,      -5.*ddx2   });
    H->tri.push_back({  Sj-1,  Sj-3,       4.*ddx2   });
    H->tri.push_back({  Sj-1,  Sj-4,      -1.*ddx2   });
    */

    for(int i=1;i<Si-1;i++){

      //First pixel of each image row: forward 2nd derivative in X-direction, central 2nd derivative in Y, both of 2nd order accuracy
      tmp.push_back({  i*Sj,  i*Sj,           fac   });

      /*
      H->tri.push_back({  i*Sj,  i*Sj,           fac2   });
      H->tri.push_back({  i*Sj,  i*Sj+1,     -5.*ddx2   });
      H->tri.push_back({  i*Sj,  i*Sj+2,      4.*ddx2   });
      H->tri.push_back({  i*Sj,  i*Sj+3,     -1.*ddx2   });
      H->tri.push_back({  i*Sj,  (i-1)*Sj,    1.*ddy2   });
      H->tri.push_back({  i*Sj,  (i+1)*Sj,    1.*ddy2   });
      */
      //central 2nd derivative of 2nd order accuracy in both X and Y directions
      for(int j=1;j<Sj-1;j++){
	tmp.push_back({  i*Sj+j,        i*Sj+j,        fac   });
	tmp.push_back({  i*Sj+j,      i*Sj+j-1,    1.*ddx2   });
	tmp.push_back({  i*Sj+j,      i*Sj+j+1,    1.*ddx2   });
	tmp.push_back({  i*Sj+j,    (i-1)*Sj+j,    1.*ddy2   });
	tmp.push_back({  i*Sj+j,    (i+1)*Sj+j,    1.*ddy2   });
      }

      //Last pixel of each image row: backward 2nd derivative in X-direction, central 2nd derivative in Y, both of 2nd order accuracy
      tmp.push_back({  (i+1)*Sj-1,  (i+1)*Sj-1,            fac   });
      /*
      H->tri.push_back({  (i+1)*Sj-1,  (i+1)*Sj-1,            fac2   });
      H->tri.push_back({  (i+1)*Sj-1,  (i+1)*Sj-2,        -5.*ddx2   });
      H->tri.push_back({  (i+1)*Sj-1,  (i+1)*Sj-3,         4.*ddx2   });
      H->tri.push_back({  (i+1)*Sj-1,  (i+1)*Sj-4,        -1.*ddx2   });
      H->tri.push_back({  (i+1)*Sj-1,  i*Sj-1,             1.*ddy2   });
      H->tri.push_back({  (i+1)*Sj-1,  (i+2)*Sj-1,         1.*ddy2   });
      */
    }

    for(int j=0;j<Sj;j++){
      tmp.push_back({  (Si-1)*Sj+j,  (Si-1)*Sj+j,       fac   });
    }

    //First pixel of last image row:  forward 2nd derivative in X, backward 2nd derivative in Y, both of 2nd order accuracy
    /*
    H->tri.push_back({  (Si-1)*Sj,  (Si-1)*Sj,            fac3   });
    H->tri.push_back({  (Si-1)*Sj,  (Si-2)*Sj,        -5.*ddy2   });
    H->tri.push_back({  (Si-1)*Sj,  (Si-3)*Sj,         4.*ddy2   });
    H->tri.push_back({  (Si-1)*Sj,  (Si-4)*Sj,        -1.*ddy2   });
    H->tri.push_back({  (Si-1)*Sj,  (Si-1)*Sj+1,      -5.*ddx2   });
    H->tri.push_back({  (Si-1)*Sj,  (Si-1)*Sj+2,       4.*ddx2   });
    H->tri.push_back({  (Si-1)*Sj,  (Si-1)*Sj+3,      -1.*ddx2   });

    //Last row of image pixels:  central 2nd derivative in X, backward 2nd derivative in Y, both of 2nd order accuracy
    for(int j=1;j<Sj-1;j++){
      H->tri.push_back({  (Si-1)*Sj+j,  (Si-1)*Sj+j,       fac2   });
      H->tri.push_back({  (Si-1)*Sj+j,  (Si-2)*Sj+j,   -5.*ddy2   });
      H->tri.push_back({  (Si-1)*Sj+j,  (Si-3)*Sj+j,    4.*ddy2   });
      H->tri.push_back({  (Si-1)*Sj+j,  (Si-4)*Sj+j,   -1.*ddy2   });
      H->tri.push_back({  (Si-1)*Sj+j,  (Si-1)*Sj+j-1,  1.*ddx2   });
      H->tri.push_back({  (Si-1)*Sj+j,  (Si-1)*Sj+j+1,  1.*ddx2   });
    }

    //Last pixel of last image row:  backward 2nd derivative in X, backward 2nd derivative in Y, both of 2nd order accuracy
    H->tri.push_back({  Si*Sj-1,  Si*Sj-1,            fac3   });
    H->tri.push_back({  Si*Sj-1,  (Si-1)*Sj-1,    -5.*ddy2   });
    H->tri.push_back({  Si*Sj-1,  (Si-2)*Sj-1,     4.*ddy2   });
    H->tri.push_back({  Si*Sj-1,  (Si-3)*Sj-1,    -1.*ddy2   });
    H->tri.push_back({  Si*Sj-1,  Si*Sj-2,        -5.*ddx2   });
    H->tri.push_back({  Si*Sj-1,  Si*Sj-3,         4.*ddx2   });
    H->tri.push_back({  Si*Sj-1,  Si*Sj-4,        -1.*ddx2   });
    */

  }

  H->tri.swap(tmp);
}


//virtual
void FixedSource::outputSource(const std::string filename){
  //Write PNG:
  //writeArrayPngPP(filename,this->Sj,this->Si,this->src);

  //Write FITS:
  long naxis    = 2;
  long naxes[2] = {(long) this->Si,(long) this->Sj};
  long Ntot = (long) this->Si*this->Sj;

  std::unique_ptr<CCfits::FITS> pFits(nullptr);
  pFits.reset( new CCfits::FITS("!"+filename,DOUBLE_IMG,naxis,naxes) );
  
  std::vector<long> extAx(2,(long) this->Si);
  std::string newName("NEW-EXTENSION");
  CCfits::ExtHDU* imageExt = pFits->addImage(newName,DOUBLE_IMG,extAx);

  std::valarray<double> array(Ntot);
  long count = 0;
  for(int j=0;j<this->Sj;j++){
    for(int i=0;i<this->Si;i++){
      array[(this->Sj-1-j)*this->Si+i] = this->src[count];
      //      array[j*this->Si+i] = this->src[count];
      count++;
    }
  }

  long fpixel(1);
  imageExt->write(fpixel,(long) Ntot,array);
  pFits->pHDU().addKey("EXPOSURE",13,"Total Exposure Time"); 
  pFits->pHDU().write(fpixel,Ntot,array); 
}

//virtual
void FixedSource::outputSourceErrors(double* errors,const std::string filename){}


//Derived class from BaseSourcePlane: FloatingSource
//===============================================================================================================

//virtual
FloatingSource::FloatingSource(int i,int j,double size,double x0,double y0,std::string reg_scheme){
  type = "floating";
  reg  = reg_scheme;
  Si   = i;
  Sj   = j;

  src = (double*) calloc(i*j,sizeof(double));
  x   = (double*) calloc(i*j,sizeof(double));
  y   = (double*) calloc(i*j,sizeof(double));

  std::map<std::string,std::string> pars;
  pars["size"] = std::to_string(size);
  pars["x0"]   = std::to_string(x0);
  pars["y0"]   = std::to_string(y0);

  this->setGrid(pars);
  this->boundPolygon();
}

//virtual
FloatingSource::~FloatingSource(){}

//non-virtual
void FloatingSource::setGrid(std::map<std::string,std::string> pars){
  double size = stof(pars["size"]);
  double x0   = stof(pars["x0"]);
  double y0   = stof(pars["y0"]);

  this->x0   = x0;
  this->y0   = y0;
  this->xmax = x0 + size;
  this->xmin = x0 - size;
  this->ymax = y0 + size;
  this->ymin = y0 - size;

  double dx = (this->xmax - this->xmin)/this->Sj;
  double dy = (this->ymax - this->ymin)/this->Si;

  this->xmax -= dx;
  this->ymax -= dy;

  for(int j=0;j<this->Sj;j++){
    for(int i=0;i<this->Si;i++){
      this->x[j*this->Si+i] = this->xmin + i*dx;
      this->y[j*this->Si+i] = this->ymin + (this->Sj-1-j)*dy;//reflect y-axis
    }
  }
}

//non-virtual
void FloatingSource::boundPolygon(){
  std::vector<double> tx {this->xmin,this->xmin,this->xmax,this->xmax,this->xmin};
  std::vector<double> ty {this->ymin,this->ymax,this->ymax,this->ymin,this->ymin};
  this->bound_x  = tx;
  this->bound_y  = ty;
}

//non-virtual
bool FloatingSource::pointInPolygon(double x,double y){
  if( x <= this->xmin || x >= this->xmax || y <= this->ymin || y >= this->ymax ){
    return false;
  } else {
    return true;
  }
}

//virtual
void FloatingSource::constructL(ImagePlane* image,CollectionMassModels* mycollection,mytable* L){
  int Ni = image->Ni;
  int Nj = image->Nj;
  int Si = this->Si;
  int Sj = this->Sj;

  long excluded = 0;
  long included = 0;
  double xp     = 0;
  double yp     = 0;
  double dx     = (this->xmax - this->xmin)/(Si - 1);
  double dy     = (this->ymax - this->ymin)/(Sj - 1);
  double norm   = 1./(dx*dy);
  int ic        = 0;
  int jc        = 0;
  double g1=0,g2=0,g3=0,g4=0;
  std::vector<mytriplet> tmp;//need to make sure that the L triplet vector is a new one

  for(int i=0;i<Ni;i++){
    for(int j=0;j<Nj;j++){

      mycollection->all_defl(image->x[i*Nj+j],image->y[i*Nj+j],xp,yp);

      if( this->pointInPolygon(xp,yp) ){
	included++;

	//	this->interpolateSource(i,j,xdefl,ydefl,&L);

	//Indices corresponding to the top left pixel
	jc = (int) floor((xp - this->xmin)/dx);
	ic = (int) floor((yp - this->ymin)/dy);
	
	//Now interpolate between neighbouring pixels and add entries to L matrix.
	//The interpolation function could return an array of column indices in the row of the L matrix, and the corresponding weights.
	//The following is for bi-linear interpolation:
	g1 = (ic+1)*dy - (yp - this->ymin);
	g2 = (jc+1)*dx - (xp - this->xmin);
	g3 = (yp - this->ymin) - ic*dy;
	g4 = (xp - this->xmin) - jc*dx;
	
	tmp.push_back({i*Nj+j, (Si-2-ic)*Sj   + (jc),       g1*g2*norm});
	tmp.push_back({i*Nj+j, (Si-2-ic)*Sj   + (jc)+1,     g1*g4*norm});
	tmp.push_back({i*Nj+j, (Si-2-ic+1)*Sj + (jc),       g3*g2*norm});
	tmp.push_back({i*Nj+j, (Si-2-ic+1)*Sj + (jc)+1,     g3*g4*norm});
      } else {
	excluded++;
      } 
    }
  }

  L->tri.swap(tmp);

  //  std::cout << excluded << " image pixels lie outside the source grid." << std::endl;
  //  std::cout << included << " image pixels lie inside the source grid." << std::endl;
}

//virtual
void FloatingSource::constructH(mytable* H){
  std::vector<mytriplet> tmp;//need to make sure that the L triplet vector is a new one

  if( this->reg == "identity" ){//---------------------------> zero order
    
    for(int i=0;i<H->Ti;i++){
      tmp.push_back({i,i,1});
    }

  } else if ( this->reg == "gradient" ){//-------------------> first order
    
    int Si = this->Si;
    int Sj = this->Sj;
    double ddx = 1./fabs(this->x[1] - this->x[0]);
    double ddy = 1./fabs(this->y[Sj] - this->y[0]);

    for(int i=1;i<Si-1;i++){
      for(int j=1;j<Sj-1;j++){
	tmp.push_back({  j*Si+i,  (j-1)*Si+i,  ddy});
	tmp.push_back({  j*Si+i,    j*Si+i-1,  ddx});
	tmp.push_back({  j*Si+i,    j*Si+i+1,  ddx});
	tmp.push_back({  j*Si+i,  (j+1)*Si+i,  ddy});
      }
    }

  } else if ( this->reg == "curvature" ){//-------------------> second order

    int Si = this->Si;
    int Sj = this->Sj;
    double dx = fabs(this->x[1] - this->x[0]);
    double dy = fabs(this->y[Sj] - this->y[0]);
    double ddx2 =  1.0/(dx*dx);
    double ddy2 =  1.0/(dy*dy);
    double sum  = -2.0*(ddx2+ddy2);
    int ii=0,jj=0;

    for(int i=0;i<Si;i++){
      tmp.push_back({i,i,1.});
    }
    for(int j=1;j<Sj-1;j++){
      tmp.push_back({  j*Si,  j*Si,  1.});
      for(int i=1;i<Si-1;i++){
	tmp.push_back({  j*Si+i,  (j-1)*Si+i,  ddy2});
	tmp.push_back({  j*Si+i,    j*Si+i-1,  ddx2});
	tmp.push_back({  j*Si+i,      j*Si+i,  sum});
	tmp.push_back({  j*Si+i,    j*Si+i+1,  ddx2});
	tmp.push_back({  j*Si+i,  (j+1)*Si+i,  ddy2});
      }
      tmp.push_back({  (j+1)*Si-1,  (j+1)*Si-1,  1.});
    }
    for(int i=Si*(Sj-1);i<Si*Sj;i++){
      tmp.push_back({i,i,1.});
    }

  }

  H->tri.swap(tmp);
}

//virtual
void FloatingSource::outputSource(const std::string filename){
  //Write PNG:
  //writeArrayPngPP(filename,this->Sj,this->Si,this->src);

  //Write FITS:
  long naxis    = 2;
  long naxes[2] = {(long) this->Si,(long) this->Sj};
  long Ntot = (long) this->Si*this->Sj;

  std::unique_ptr<CCfits::FITS> pFits(nullptr);
  pFits.reset( new CCfits::FITS("!"+filename,DOUBLE_IMG,naxis,naxes) );
  
  std::vector<long> extAx(2,(long) this->Si);
  std::string newName("NEW-EXTENSION");
  CCfits::ExtHDU* imageExt = pFits->addImage(newName,DOUBLE_IMG,extAx);

  std::valarray<double> array(Ntot);
  long count = 0;
  for(int j=0;j<this->Sj;j++){
    for(int i=0;i<this->Si;i++){
      array[(this->Sj-1-j)*this->Si+i] = this->src[count];
      //      array[j*this->Si+i] = this->src[count];
      count++;
    }
  }

  long fpixel(1);
  imageExt->write(fpixel,(long) Ntot,array);
  pFits->pHDU().addKey("EXPOSURE",13,"Total Exposure Time"); 
  pFits->pHDU().write(fpixel,Ntot,array); 
}

//virtual
void FloatingSource::outputSourceErrors(double* errors,const std::string filename){}


//Derived class from BaseSourcePlane: AdaptiveSource
//===============================================================================================================

//virtual
AdaptiveSource::AdaptiveSource(int a,std::string reg_scheme){
  type    = "adaptive";
  mode    = "random";
  spacing = 0;
  reg     = reg_scheme;
  Sm      = a;
  Si      = a;
  Sj      = 0;

  src = (double*) calloc(Sm,sizeof(double));
  x   = (double*) calloc(Sm,sizeof(double));
  y   = (double*) calloc(Sm,sizeof(double));
}

//virtual
AdaptiveSource::AdaptiveSource(std::string m,int a,int b,std::string reg_scheme){
  type    = "adaptive";
  mode    = m;
  spacing = b;
  reg     = reg_scheme;
  Sm      = a;// + 8;
  Si      = a;// + 8;
  Sj      = 0;

  src = (double*) calloc(Sm,sizeof(double));
  x   = (double*) calloc(Sm,sizeof(double));
  y   = (double*) calloc(Sm,sizeof(double));
}

//virtual
AdaptiveSource::~AdaptiveSource(){
  std::vector<a_triangle> ().swap(triangles);
}

//virtual
void AdaptiveSource::constructL(ImagePlane* image,CollectionMassModels* mycollection,mytable* L){
  double xp,yp;
  double wa,wb,wc;
  double ybc,xac,xcb,yac,xxc,yyc,den;
  int ja,jb,jc;
  a_triangle triangle;
  std::vector<mytriplet> tmp;//need to make sure that the L triplet vector is a new one

  int included = 0;
  int excluded = 0;

  //  FILE* fh = fopen("rays.dat","w");

  for(int i=0;i<image->Nm;i++){

    if( image->active[i] == 1 ){

      mycollection->all_defl(image->x[i],image->y[i],xp,yp);
      //    fprintf(fh,"%12.4f %12.4f\n",xp,yp);
      
      int flag = 0;
      for(int j=0;j<this->n_triangles;j++){
	triangle = this->triangles[j];
	
	ybc = this->y[triangle.b] - this->y[triangle.c];//(yb-yc)
	xac = this->x[triangle.a] - this->x[triangle.c];//(xa-xc)
	xcb = this->x[triangle.c] - this->x[triangle.b];//(xc-xb)
	yac = this->y[triangle.a] - this->y[triangle.c];//(ya-yc)
	xxc = xp                  - this->x[triangle.c];//(x -xc)
	yyc = yp                  - this->y[triangle.c];//(y -yc)
	den = ybc*xac + xcb*yac;
	
	wa = ( ybc*xxc+xcb*yyc)/den;
	wb = (-yac*xxc+xac*yyc)/den;
	wc = 1. - wa - wb;
	
	if( wa <= 1. && wb <= 1. && wc <= 1. ){
	  tmp.push_back({i,   triangle.a,   wa});
	  tmp.push_back({i,   triangle.b,   wb});
	  tmp.push_back({i,   triangle.c,   wc});

	  flag = 1;
	  break;
	}
      }
      
      if( flag == 0 ){
	excluded++;
      } else {
	included++;
      }
      
    }
  }

  for(int i=0;i<this->used.size();i++){
    tmp.push_back({this->used[i],   i,   1.});
  }


  //fclose(fh);

  //std::cout << included << " " << excluded << std::endl;
  //std::cout << tmp.size() << std::endl;
  
  L->tri.swap(tmp);
}

//non-virtual
void AdaptiveSource::createAdaGrid(ImagePlane* image,CollectionMassModels* mycollection){
  if( this->mode == "random" ){
    
    double xmax =  image->width/2.;
    double xmin = -image->width/2.;
    double ymax =  image->height/2.;
    double ymin = -image->height/2.;
    
    double xtmp,ytmp;
    srand48(time(NULL));
    for(int i=0;i<this->Sm;i++){
      xtmp = drand48()*(xmax - xmin) + xmin;
      ytmp = drand48()*(ymax - ymin) + ymin;
      mycollection->all_defl(xtmp,ytmp,this->x[i],this->y[i]);   
    }

    for(int i=0;i<image->Nm;i++){
      image->active[i] = 1;
    }

  } else if( this->mode == "image" ){

    for(int i=0;i<image->Nm;i++){
      image->active[i] = 1;
    }

    std::vector<int> used;
    double xtmp,ytmp;
    int i0    = (int) floor( (this->spacing-1)/2. );
    int j0    = (int) floor( (this->spacing-1)/2. );
    int count = 0;//must go up to Sm
    for(int i=i0;i<image->Ni;i=i+this->spacing){
      for(int j=j0;j<image->Nj;j=j+this->spacing){
	xtmp = image->x[i*image->Nj+j];
	ytmp = image->y[i*image->Nj+j];
	image->active[i*image->Nj+j] = 0;
	used.push_back(i*image->Nj+j);
	mycollection->all_defl(xtmp,ytmp,this->x[count],this->y[count]);
	count++;
      }
    }

    this->used.swap(used);
    //    std::cout << "The number of pixels in the adaptive source is " << count << " and it must be " << this->Sm << std::endl;

  } else if( this->mode == "grid" ){

    int Nj = image->Nj;
    int Ni = image->Ni;

    int i0    = floor(Ni/2);
    int j0    = floor(Nj/2);
    double di = image->width/(Ni-1);
    double dj = image->height/(Nj-1);
    
    int count = 0;
    double xtmp,ytmp;
    //    FILE* fh1 = fopen("ada_grid.dat","w");
    for(int ii=0;ii<Ni;ii=ii+this->spacing){
      for(int jj=0;jj<Nj;jj=jj+this->spacing){
	xtmp   =  (jj-j0)*di;
	ytmp   = -(ii-i0)*dj;//reflect y-axis
	mycollection->all_defl(xtmp,ytmp,this->x[count],this->y[count]);
	//	fprintf(fh1,"%12.4f %12.4f\n",this->x[count],this->y[count]);
	count++;
      }
    }
    //    fclose(fh1);
    //    std::cout << "The number of pixels in the adaptive source is " << count << " and it must be " << this->Sm << std::endl;


    std::vector<xypoint> frame(8);
    double d_out = 5.0;
    double d_in  = 0.000001;

    frame[0].x = image->x[0]         - d_out;
    frame[0].y = image->y[0]         + d_out;
    frame[1].x = image->x[Nj-1]      + d_out;
    frame[1].y = image->y[0]         + d_out;
    frame[2].x = image->x[0]         - d_out;
    frame[2].y = image->y[(Ni-1)*Nj] - d_out;
    frame[3].x = image->x[Nj-1]      + d_out;
    frame[3].y = image->y[(Ni-1)*Nj] - d_out;

    frame[4].x = mycollection->models[0]->mpars["x0"] - d_in;
    frame[4].y = mycollection->models[0]->mpars["y0"] + d_in;
    frame[5].x = mycollection->models[0]->mpars["x0"] + d_in;
    frame[5].y = mycollection->models[0]->mpars["y0"] + d_in;
    frame[6].x = mycollection->models[0]->mpars["x0"] - d_in;
    frame[6].y = mycollection->models[0]->mpars["y0"] - d_in;
    frame[7].x = mycollection->models[0]->mpars["x0"] + d_in;
    frame[7].y = mycollection->models[0]->mpars["y0"] - d_in;


    //    FILE* fh = fopen("outer_grid.dat","w");
    for(int i=0;i<frame.size();i++){
      mycollection->all_defl(frame[i].x,frame[i].y,this->x[count+i],this->y[count+i]);
      //      fprintf(fh,"%12.4f %12.4f\n",this->x[count+i],this->y[count+i]);
    }
    //    fclose(fh);

    for(int i=0;i<image->Nm;i++){
      image->active[i] = 1;
    }

  }

}

//non-virtual
void AdaptiveSource::createDelaunay(){
  typedef CGAL::Exact_predicates_inexact_constructions_kernel            K;
  typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int,K>    Vb;
  typedef CGAL::Triangulation_face_base_with_info_2<unsigned int,K>      Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                    Tds;
  typedef CGAL::Delaunay_triangulation_2<K,Tds>                          Delaunay;
  typedef K::Point_2                                                     Point;

  std::vector< std::pair<Point,int> > points;
  for(int i=0;i<this->Sm;i++){
    points.push_back( std::make_pair(Point(x[i],y[i]),i) );
  }

  Delaunay triangulation;
  triangulation.insert(points.begin(),points.end());


  //Get each Delaunay triangle in my own struct, and give them an ID
  //[constructing this->triangles]
  Delaunay::Finite_faces_iterator fit;
  Delaunay::Face_handle face;
  a_triangle triangle;
  int face_id = 0;
  this->n_triangles = triangulation.number_of_faces();
  this->triangles.resize( triangulation.number_of_faces() );

  for(fit=triangulation.finite_faces_begin();fit!=triangulation.finite_faces_end();fit++){
    face = fit;
    face->info() = face_id;

    triangle.a = (int) face->vertex(0)->info();
    triangle.b = (int) face->vertex(1)->info();
    triangle.c = (int) face->vertex(2)->info();
    
    this->triangles[face_id] = triangle;
    face_id++;
  }


  //Find which triangle IDs around each vertex. This is zero for the convex-hull vertices (having more than 1 infinite incident faces) 
  //[constructing this->face_ids_per_vertex]
  Delaunay::Face_circulator fc;
  Delaunay::Vertex_circulator vc;
  Delaunay::Finite_vertices_iterator fvi;
  this->opposite_edges_per_vertex.resize( this->Sm );
  for(fvi=triangulation.finite_vertices_begin();fvi!=triangulation.finite_vertices_end();fvi++){

    int inf = 0;
    std::vector<int> opposite_indices;
    fc = triangulation.incident_faces(fvi);
    do{
      //      vertex_face_ids.push_back( fc->info() );
      vc = triangulation.incident_vertices(fvi,fc);
      opposite_indices.push_back( vc->info() );
      vc++;
      opposite_indices.push_back( vc->info() );      
      if( triangulation.is_infinite(fc) ){
	inf++;
      }
    }while( ++fc != triangulation.incident_faces(fvi) );

    if( inf > 1 ){
      opposite_indices.resize(0);
    }

    this->opposite_edges_per_vertex[ fvi->info() ] = opposite_indices;
  }


}


//virtual
void AdaptiveSource::constructH(mytable* H){
  std::vector<mytriplet> tmp;//need to make sure that the L triplet vector is a new one

  if( this->reg == "identity" ){//---------------------------> zero order
    
    for(int i=0;i<H->Ti;i++){
      tmp.push_back({i,i,1});
    }

  } else if ( this->reg == "gradient" ){//-------------------> first order





    
    for(int i=0;i<this->Sm;i++){
      
      xypoint p0 = {this->x[i],this->y[i]};
      std::vector<int> indices = this->opposite_edges_per_vertex[i];
      
      if( indices.size() == 0 ){

	//we are on a convex-hull point that has no smoothing
	tmp.push_back({i,i,1});

      } else {
	
	std::map<int,double> weights;
	std::map<int,double> weights_x;
	std::map<int,double> weights_y;
	std::vector<double> points_x;
	std::vector<double> points_y;
	  double l0,l1,l2,y20,xm0,x02,ym0,y01,x10,den;

	for(int j=0;j<indices.size();j=j+2){
	  
	  xypoint p1 = {this->x[indices[j]],this->y[indices[j]]};
	  xypoint p2 = {this->x[indices[j+1]],this->y[indices[j+1]]};

	  //Smoothing in the x-direction
	  double ymin = 0.;
	  double ymax = 0.;
	  if( p1.y > p2.y ){
	    ymax = p1.y;
	    ymin = p2.y;
	  } else {
	    ymax = p2.y;
	    ymin = p1.y;
	  }
	  
	  //Find intersection point (if it exists), interpolate to get the weights, and add them to <map>weights
	  if( ymin <= p0.y && p0.y <= ymax ){
	    xypoint pint = intersection_point_x(p0,p1,p2);
	    points_x.push_back( pint.x );

	    xypoint pmid = {(p0.x+pint.x)/2.,p0.y};
	    y20 = p2.y    - p0.y;
	    xm0 = pmid.x  - p0.x;
	    x02 = p0.x    - p2.x;
	    ym0 = pmid.y  - p0.y;
	    y01 = p0.y    - p1.y;
	    x10 = p1.x    - p0.x;
	    den = y20*x10 - y01*x02;//to make y10 appear without calculating it again
	    
	    l1 = (y20*xm0 + x02*ym0)/den;
	    l2 = (y01*xm0 + x10*ym0)/den;
	    l0 = 1. - l1 - l2;

	    weights_x[ i            ] -= l0;
	    weights_x[ indices[j]   ] += l1;
	    weights_x[ indices[j+1] ] += l2;
	  }


	  //Smoothing in the y-direction
	  double xmin = 0.;
	  double xmax = 0.;

	  if( p1.x > p2.x ){
	    xmax = p1.x;
	    xmin = p2.x;
	  } else {
	    xmax = p2.x;
	    xmin = p1.x;
	  }

	  //Find intersection point (if it exists), interpolate to get the weights, and add them to <map>weights
	  if( xmin <= p0.x && p0.x <= xmax ){
	    xypoint pint = intersection_point_y(p0,p1,p2);
	    points_y.push_back( pint.y );

	    xypoint pmid = {p0.x,(p0.y+pint.y)/2.};
	    y20 = p2.y    - p0.y;
	    xm0 = pmid.x  - p0.x;
	    x02 = p0.x    - p2.x;
	    ym0 = pmid.y  - p0.y;
	    y01 = p0.y    - p1.y;
	    x10 = p1.x    - p0.x;
	    den = y20*x10 - y01*x02;//to make y10 appear without calculating it again
	    
	    l1 = (y20*xm0 + x02*ym0)/den;
	    l2 = (y01*xm0 + x10*ym0)/den;
	    l0 = 1. - l1 - l2;

	    weights_y[ i            ] -= l0;
	    weights_y[ indices[j]   ] += l1;
	    weights_y[ indices[j+1] ] += l2;
	  }

	}


	double dx = (points_x[0] + p0.x)/2. - (points_x[1] + p0.x)/2.;
	for(it_int_double iterator=weights_x.begin();iterator!=weights_x.end();iterator++){
	  weights[iterator->first] += iterator->second/dx;
	}
	weights_x.clear();

	double dy = (points_y[0] + p0.y)/2. - (points_y[1] + p0.y)/2.;
	for(it_int_double iterator=weights_y.begin();iterator!=weights_y.end();iterator++){
	  weights[iterator->first] += iterator->second/dy;
	}
	weights_y.clear();
	
	for(it_int_double iterator=weights.begin();iterator!=weights.end();iterator++){
	  tmp.push_back({i,iterator->first,iterator->second});
	}
	weights.clear();

      }

    }







  } else if ( this->reg == "curvature" ){//-------------------> second order

    for(int i=0;i<this->Sm;i++){
      
      xypoint p0 = {this->x[i],this->y[i]};
      std::vector<int> indices = this->opposite_edges_per_vertex[i];
      
      if( indices.size() == 0 ){

	//we are on a convex-hull point that has no smoothing
	tmp.push_back({i,i,1});

      } else {
	
	std::map<int,double> weights;

	for(int j=0;j<indices.size();j=j+2){
	  double ymin = 0.;
	  double ymax = 0.;
	  double xmin = 0.;
	  double xmax = 0.;
	  
	  xypoint p1 = {this->x[indices[j]],this->y[indices[j]]};
	  xypoint p2 = {this->x[indices[j+1]],this->y[indices[j+1]]};
	  
	  //Smoothing in the x-direction
	  if( p1.y > p2.y ){
	    ymax = p1.y;
	    ymin = p2.y;
	  } else {
	    ymax = p2.y;
	    ymin = p1.y;
	  }
	  
	  //Find intersection point (if it exists), interpolate to get the weights, and add them to <map>weights
	  if( ymin <= p0.y && p0.y <= ymax ){
	    xypoint pint = intersection_point_x(p0,p1,p2);
	    double l2 = ((p1.y-p0.y)*(pint.x-p1.x)+(p0.x-p1.x)*(pint.y-p1.y))/((p2.y-p1.y)*(p0.x-p1.x)+(p1.x-p2.x)*(p0.y-p1.y));
	    double l1 = 1.- l2;
	    //	    double d12 = pow(p1.x-p2.x,2) + pow(p1.y-p2.y,2);
	    //	    double di2 = pow(pint.x-p2.x,2) + pow(pint.y-p2.y,2);
	    //	    double l1 = sqrt(di2/d12);
	    //	    double l2 = 1. - l1;
	    double d  = fabs(p0.x-pint.x);
	    weights[ i            ] -= 1./d;
	    weights[ indices[j]   ] += l1/d;
	    weights[ indices[j+1] ] += l2/d;
	  }

	  //Smoothing in the y-direction
	  if( p1.x > p2.x ){
	    xmax = p1.x;
	    xmin = p2.x;
	  } else {
	    xmax = p2.x;
	    xmin = p1.x;
	  }

	  //Find intersection point (if it exists), interpolate to get the weights, and add them to <map>weights
	  if( xmin <= p0.x && p0.x <= xmax ){
	    xypoint pint = intersection_point_y(p0,p1,p2);
	    double l2 = ((p1.y-p0.y)*(pint.x-p1.x)+(p0.x-p1.x)*(pint.y-p1.y))/((p2.y-p1.y)*(p0.x-p1.x)+(p1.x-p2.x)*(p0.y-p1.y));
	    double l1 = 1.- l2;
	    //	    double d12 = pow(p1.x-p2.x,2) + pow(p1.y-p2.y,2);
	    //	    double di2 = pow(pint.x-p2.x,2) + pow(pint.y-p2.y,2);
	    //	    double l1 = sqrt(di2/d12);
	    //	    double l2 = 1. - l1;
	    double d  = fabs(p0.y-pint.y);
	    weights[ i            ] -= 1./d;
	    weights[ indices[j]   ] += l1/d;
	    weights[ indices[j+1] ] += l2/d;
	  }

	}
	
	//Loop through weights and add entries to H
	for(it_int_double iterator=weights.begin();iterator!=weights.end();iterator++){
	  tmp.push_back({i,iterator->first,iterator->second});
	}
	weights.clear();

      }

    }

    
  }

  H->tri.swap(tmp);
}


//private
AdaptiveSource::xypoint AdaptiveSource::intersection_point_x(xypoint p0,xypoint p1,xypoint p2){
  xypoint pint;

  if( p1.x == p2.x ){
    pint.x = p1.x;
    pint.y = p0.y;
  } else {
    pint.x = p1.x + (p0.y-p1.y)*(p2.x-p1.x)/(p2.y-p1.y);
    pint.y = p0.y;
  }

  return pint;
}

AdaptiveSource::xypoint AdaptiveSource::intersection_point_y(xypoint p0,xypoint p1,xypoint p2){
  xypoint pint;

  if( p1.y == p2.y ){
    pint.x = p0.x;
    pint.y = p1.y;
  } else {
    pint.x = p0.x;
    pint.y = (p2.y-p1.y)*(p0.x-p1.x)/(p2.x-p1.x) + p1.y;
  }

  return pint;
}


//non-virtual
void AdaptiveSource::writeTriangles(){
  a_triangle t;
  FILE* fh = fopen("triangles.dat","w");
  for(int i=0;i<this->n_triangles;i++){
    t = this->triangles[i];
    fprintf(fh,"%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",this->x[t.a],this->y[t.a],this->x[t.b],this->y[t.b],this->x[t.c],this->y[t.c]);
  }
  fclose(fh);
}

//virtual
void AdaptiveSource::outputSource(const std::string filename){
  typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
  typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
  typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
  typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
  typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    Voronoi;
  typedef Voronoi::Locate_result                                               Locate_result;
  typedef Voronoi::Face_handle                                                 Face_handle;
  typedef Voronoi::Ccb_halfedge_circulator                                     Ccb_halfedge_circulator;
  typedef K::Point_2                                                           Point;

  //Scale the source solution between 0 and 1
  double* ss = (double*) calloc(this->Sm,sizeof(double));

  double max = 0.;
  for(int i=0;i<this->Sm;i++){
    if( this->src[i] > max ){
      max = this->src[i];
    }
  }
  double min = 100.;
  for(int i=0;i<this->Sm;i++){
    if( this->src[i] < min ){
      min = this->src[i];
    }
  }
  for(int i=0;i<this->Sm;i++){
    ss[i] = (src[i]-min)/(max-min);
  }

  //Calculate the Voronoi graph
  Voronoi voronoi;
  std::vector<Point> points;
  for(int i=0;i<this->Sm;i++){
    points.push_back( Point(this->x[i],this->y[i]) );
  }
  voronoi.insert(points.begin(),points.end());

  //Locate the centres of the Voronoi cells and then get the vertices of the surrounding face
  //This code is based on an online example from the CGAL voronoi package
  FILE* fh = fopen(filename.c_str(),"w");
  for(int i=0;i<this->Sm;i++){

    if( this->opposite_edges_per_vertex[i].size() != 0 ){
      fprintf(fh,"%12.5f",ss[i]);

      Locate_result f   = voronoi.locate(Point(this->x[i],this->y[i]));
      Face_handle* face = boost::get<Face_handle>(&f);
    
      Ccb_halfedge_circulator ec_start = (*face)->ccb();
      Ccb_halfedge_circulator ec = ec_start;
      do {
	Point p = ec->source()->point();
	fprintf(fh,"%12.5f %12.5f",p.x(),p.y());
      } while( ++ec != ec_start );
      fprintf(fh,"\n");
    }

  }
  fclose(fh);

  /*
  //Set up an image grid, match the value of its pixels to the voronoi cells, and output it as fits
  std::string fitsname = path + "vkl_source.fits";

  double width  = 1.2;
  double height = 1.2;

  ImagePlane myfits(50,50,width,height);
  double xmin = - width/2.;
  double xmax =   width/2.;
  double ymin = - height/2.;
  double ymax =   height/2.;

  for(int i=0;i<this->Sm;i++){
    if( xmin < this->x[i] && this->x[i] < xmax && ymin < this->y[i] && this->y[i] < ymax ){
      std::cout << this->x[i] << " " << this->y[i] << std::endl;
      Locate_result ri = voronoi.locate(Point(this->x[i],this->y[i]));

      for(int m=0;m<myfits.Nm;m++){
	Locate_result rm = voronoi.locate(Point(myfits.x[m],myfits.y[m]));
	if( rm == ri ){
	  myfits.img[m] = ss[i];
	  break;
	}
      }

    }
  }
  myfits.writeImage(fitsname);
  */


  free(ss);
}


//virtual
void AdaptiveSource::outputSourceErrors(double* errors,const std::string filename){
  typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
  typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
  typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
  typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
  typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    Voronoi;
  typedef Voronoi::Locate_result                                               Locate_result;
  typedef Voronoi::Face_handle                                                 Face_handle;
  typedef Voronoi::Ccb_halfedge_circulator                                     Ccb_halfedge_circulator;
  typedef K::Point_2                                                           Point;

  //Scale the errors between 0 and 1 (same is done for the source)
  double* ee = (double*) calloc(this->Sm,sizeof(double));

  double max = 0.;
  for(int i=0;i<this->Sm;i++){
    if( errors[i] > max ){
      max = errors[i];
    }
  }
  double min = 100.;
  for(int i=0;i<this->Sm;i++){
    if( errors[i] < min ){
      min = errors[i];
    }
  }
  for(int i=0;i<this->Sm;i++){
    ee[i] = (errors[i]-min)/(max-min);
  }

  //Calculate the Voronoi graph
  Voronoi voronoi;
  std::vector<Point> points;
  for(int i=0;i<this->Sm;i++){
    points.push_back( Point(this->x[i],this->y[i]) );
  }
  voronoi.insert(points.begin(),points.end());

  //Locate the centres of the Voronoi cells and then get the vertices of the surrounding face
  //This code is based on an online example from the CGAL voronoi package
  FILE* fh = fopen(filename.c_str(),"w");
  for(int i=0;i<this->Sm;i++){

    if( this->opposite_edges_per_vertex[i].size() != 0 ){
      fprintf(fh,"%12.5f",ee[i]);

      Locate_result f   = voronoi.locate(Point(this->x[i],this->y[i]));
      Face_handle* face = boost::get<Face_handle>(&f);
    
      Ccb_halfedge_circulator ec_start = (*face)->ccb();
      Ccb_halfedge_circulator ec = ec_start;
      do {
	Point p = ec->source()->point();
	fprintf(fh,"%12.5f %12.5f",p.x(),p.y());
      } while( ++ec != ec_start );
      fprintf(fh,"\n");
    }

  }
  fclose(fh);

  free(ee);
}






















/*
//This is an implementation of convex-hull using the boost library
void AdaptiveSource::boundPolygon(){
  typedef boost::tuple<double, double> mypoint;
  typedef boost::geometry::model::ring<mypoint > myring;
  
  mypoint point;
  myring ring;

  for(int i=0;i<this->Sm;i++){
    point = boost::make_tuple(this->x[i],this->y[i]);
    boost::geometry::append(ring,point);
  }
  
  myring hull;
  boost::geometry::convex_hull(ring,hull);

  this->bound_x.resize(hull.size());
  this->bound_y.resize(hull.size());

  for(int i=0;i<hull.size();i++){
    point = hull[i];
    this->bound_x[i] = boost::geometry::get<0>(point);
    this->bound_y[i] = boost::geometry::get<1>(point);
  }
}
*/


