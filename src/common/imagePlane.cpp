#include "imagePlane.hpp"

#include <cmath>
#include <fstream>
#include <cstdlib>
#include <CCfits/CCfits>

#include "tableDefinition.hpp"


//ImagePlane class implementation
//============================================================================================
ImagePlane::ImagePlane(const std::string filepath,int j,int i,double w,double h){
  Ni     = i;
  Nj     = j;
  Nm     = Ni*Nj;
  width  = w;
  height = h;
  S.Ti = Nm;
  S.Tj = Nm;
  C.Ti = Nm;
  C.Tj = Nm;
  B.Ti = Nm;
  B.Tj = Nm;

  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  image.readAllKeys();
  //Check if sizes agree
  std::valarray<float> contents(image.axis(0)*image.axis(1));
  this->readFits(filepath,contents);

  img    = (double*) calloc(Nm,sizeof(double));
  x      = (double*) calloc(Nm,sizeof(double));
  y      = (double*) calloc(Nm,sizeof(double));
  defl_x = (double*) calloc(Nm,sizeof(double));
  defl_y = (double*) calloc(Nm,sizeof(double));
  active = (int*) calloc(Nm,sizeof(int));
  cells  = (SourceCell**) calloc(Nm,sizeof(SourceCell*));
  
  int i0    = floor(Ni/2);
  int j0    = floor(Nj/2);
  double di = width/(Ni);
  double dj = height/(Nj);

  for(int ii=0;ii<Ni;ii++){
    for(int jj=0;jj<Nj;jj++){
      x[ii*Nj+jj]   =  (jj-j0)*di;
      y[ii*Nj+jj]   = -(ii-i0)*dj;//reflect y-axis
      img[ii*Nj+jj] = contents[ii*Nj+jj];
    }
  }
}

ImagePlane::ImagePlane(int j,int i,double w,double h){
  Ni     = i;
  Nj     = j;
  Nm     = Ni*Nj;
  width  = w;
  height = h;
  S.Ti = Nm;
  S.Tj = Nm;
  C.Ti = Nm;
  C.Tj = Nm;
  B.Ti = Nm;
  B.Tj = Nm;

  img    = (double*) calloc(Nm,sizeof(double));
  x      = (double*) calloc(Nm,sizeof(double));
  y      = (double*) calloc(Nm,sizeof(double));
  defl_x = (double*) calloc(Nm,sizeof(double));
  defl_y = (double*) calloc(Nm,sizeof(double));
  active = (int*) calloc(Nm,sizeof(int));
  cells  = (SourceCell**) calloc(Nm,sizeof(SourceCell*));

  int i0    = floor(i/2);
  int j0    = floor(j/2);
  double di = width/(i);
  double dj = height/(j);

  for(int ii=0;ii<Ni;ii++){
    for(int jj=0;jj<Nj;jj++){
      x[ii*Nj+jj] =  (jj-j0)*dj;
      y[ii*Nj+jj] = -(ii-i0)*di;//reflect y-axis
    }
  }
}

ImagePlane::ImagePlane(int i,double w,double h){
  Ni     = 0;
  Nj     = 0;
  Nm     = i;
  width  = w;
  height = h;
  S.Ti = Nm;
  S.Tj = Nm;
  C.Ti = Nm;
  C.Tj = Nm;
  B.Ti = Nm;
  B.Tj = Nm;

  img    = (double*) calloc(Nm,sizeof(double));
  x      = (double*) calloc(Nm,sizeof(double));
  y      = (double*) calloc(Nm,sizeof(double));
  defl_x = (double*) calloc(Nm,sizeof(double));
  defl_y = (double*) calloc(Nm,sizeof(double));
  active = (int*) calloc(Nm,sizeof(int));
  cells  = (SourceCell**) calloc(Nm,sizeof(SourceCell*));
}

ImagePlane::ImagePlane(const ImagePlane& image){
  Ni = image.Ni;
  Nj = image.Nj;
  Nm = image.Nm;
  width = image.width;
  height = image.height;
  img    = (double*) calloc(Nm,sizeof(double));
  x      = (double*) calloc(Nm,sizeof(double));
  y      = (double*) calloc(Nm,sizeof(double));
  defl_x = (double*) calloc(Nm,sizeof(double));
  defl_y = (double*) calloc(Nm,sizeof(double));
  active = (int*) calloc(Nm,sizeof(int));
  cells  = (SourceCell**) calloc(Nm,sizeof(SourceCell*));
  for(int i=0;i<Nm;i++){
    img[i] = image.img[i];
    x[i] = image.x[i];
    y[i] = image.y[i];
  }
}

ImagePlane::~ImagePlane(){
  free(img);
  free(x);
  free(y);
  free(defl_x);
  free(defl_y);
  free(active);
  free(cells);
}


void ImagePlane::writeBin(const std::string filename){
  std::ofstream out(filename,std::ios::out|std::ios::binary);
  for(int i=0;i<this->Nm;i++){
    out.write((const char*) (&this->img[i]),sizeof(double));
  }
  out.close();
}

void ImagePlane::writeImage(const std::string filename){  
  long naxis    = 2;
  long naxes[2] = {(long) this->Ni,(long) this->Nj};
  long Ntot = (long) this->Nm;
  
  //  std::unique_ptr<CCfits::FITS> pFits(nullptr);
  std::auto_ptr<CCfits::FITS> pFits(0);
  pFits.reset( new CCfits::FITS("!"+filename,FLOAT_IMG,naxis,naxes) );
  
  std::vector<long> extAx(2,(long) this->Ni);
  CCfits::ExtHDU* imageExt = pFits->addImage("NEW-EXTENSION",FLOAT_IMG,extAx);
  
  //Need Ni and Nj as index counters to flip image
  std::valarray<float> array(Ntot);
  long count = 0;
  for(int j=0;j<this->Nj;j++){
    for(int i=0;i<this->Ni;i++){
      array[(this->Nj-1-j)*this->Ni+i] = (float) (this->img[count]);
      count++;
    }
  }
  
  long fpixel(1);
  imageExt->write(fpixel,Ntot,array);
  //  pFits->pHDU().addKey("EXPOSURE",13,"Total Exposure Time"); 
  pFits->pHDU().write(fpixel,Ntot,array); 
  //  std::cout << pFits->pHDU() << std::endl;
}

void ImagePlane::readFits(const std::string filepath,std::valarray<float>& contents){
  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  
  std::valarray<float> tmp;
  image.readAllKeys();
  image.read(tmp);
  
  int Nj = image.axis(0);
  int Ni = image.axis(1);

  //convert FITS standard (bottom to top) to the one used in this code (top to bottom)
  for(int j=0;j<Nj;j++){
    for(int i=0;i<Ni;i++){
      contents[j*Ni+i] = tmp[(Nj-j-1)*Ni+i];
    }
  }
}

void ImagePlane::readS(const std::string filepath){
  if( filepath == "0" ){
    
    for(int i=0;i<this->Nm;i++){
      this->S.tri.push_back({i,i,1.});
      this->lookup[i] = i;
    }

  } else {

    std::valarray<float> contents(this->Nm);
    this->readFits(filepath,contents);
    
    int count = 0;
    for(int i=0;i<this->Nm;i++){
      if( contents[i] == 1 ){
	this->S.tri.push_back({i,i,1.});
	this->lookup[i] = count;
	count++;
      } else {
	this->S.tri.push_back({i,i,0});
      }
    }

  }
}

void ImagePlane::maskData(std::map<int,int> lookup,ImagePlane* masked){
  typedef std::map<int,int>::iterator it;
  for(it myiterator=lookup.begin();myiterator!=lookup.end();myiterator++){
    masked->img[ myiterator->second ] = this->img[ myiterator->first ];
    masked->x[ myiterator->second ]   = this->x[ myiterator->first ];
    masked->y[ myiterator->second ]   = this->y[ myiterator->first ];
  }
}

/*
void ImagePlane::setMaskedC(mytable* Cout,mytable* S,mytable* C){
  //this function does the same as multiplying algebraically tables S and C

  int* Sfull = (int*) calloc(this->Nm,sizeof(int));
  for(int k=0;k<this->S->tri.size();k++){
    Sfull[ this->S->tri[k].i ] = 1;
  }

  int i0,j0;
  double val;
  for(int k=0;k<C->tri.size();k++){
    i0  = C->tri[k].i;
    j0  = C->tri[k].j;
    val = C->tri[k].v;

    if( Sfull[i0] == 1 && Sfull[j0] == 1 ){
      Cout->tri.push_back({this->lookup[i0],this->lookup[j0],val});
    }
    
  }

  free(Sfull);
}
*/

void ImagePlane::readC(const std::string flag,const std::string filepath){
  double value;

  if( flag == "uniform" ){

    //There is only one element in the file, repeated along the diagonal
    std::ifstream infile(filepath);
    infile >> value;
    for(int i=0;i<this->Nm;i++){
      this->C.tri.push_back({i,i,1.0/pow(value,2)});
    }
    infile.close();

  } else if( flag == "map" ){

    //There are exactly Ni x Nj (diagonal) elements in the file and in table C
    std::string extension = filepath.substr(filepath.find(".")+1,std::string::npos);

    if( extension == "dat" ){

      std::ifstream infile(filepath);
      int i;
      while( true ){
	infile >> i >> i >> value;
	if( infile.eof() ) break;
	this->C.tri.push_back({i,i,1.0/pow(value,2)});
      }
      infile.close();

    } else if( extension == "fits" ){
      
      std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
      CCfits::PHDU& image = pInfile->pHDU();
      image.readAllKeys();
      std::valarray<float> contents(image.axis(0)*image.axis(1));
      this->readFits(filepath,contents);
      
      for(int i=0;i<this->Ni*this->Nj;i++){
	this->C.tri.push_back({i,i,1.0/pow(contents[i],2)});
	//	  std::cout << 1./pow(contents[i*this->Nj+j],2) << std::endl;
      }

    }


  } else if( flag == "correlated" ){

    //There are at least Ni x Nj (diagonal plus off-diagonal) elements in the file and in table C
    std::ifstream infile(filepath);
    int i,j;
    while( true ){
      infile >> i >> j >> value;
      if( infile.eof() ) break;
      this->C.tri.push_back({i,j,1.0/pow(value,2)});
    }
    infile.close();

  } else {

  }

}

void ImagePlane::readB(const std::string filepath,int i,int j,int ci,int cj){
  if( filepath == "0" ){

    //Create a diagonal matrix of the given dimensions.
    for(int i=0;i<this->B.Ti;i++){
      this->B.tri.push_back({i,i,1.});
    }

  } else {

    int Pi     = i;
    int Pj     = j;
    int Ncropx = ci;
    int Ncropy = cj;

    if( Pi%2 == 0 ){
      if( ci%2 != 0 ){
	std::cout << std::endl << std::endl << "PROBLEM!!! PSF AND CROPPED PSF DIMENSIONS MUST BE BOTH EVEN" << std::endl << std::endl;
      }
    } else {
      if( ci%2 == 0 ){
	std::cout << std::endl << std::endl << "PROBLEM!!! PSF AND CROPPED PSF DIMENSIONS MUST BE BOTH ODD" << std::endl << std::endl;
      }
    }
    if( Pj%2 == 0 ){
      if( cj%2 != 0 ){
	std::cout << std::endl << std::endl << "PROBLEM!!! PSF AND CROPPED PSF DIMENSIONS MUST BE BOTH EVEN" << std::endl << std::endl;
      }
    } else {
      if( cj%2 == 0 ){
	std::cout << std::endl << std::endl << "PROBLEM!!! PSF AND CROPPED PSF DIMENSIONS MUST BE BOTH ODD" << std::endl << std::endl;
      }
    }



    // Read PSF from file
    std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
    CCfits::PHDU& image = pInfile->pHDU();
    image.readAllKeys();
    //Check if sizes agree
    std::valarray<float> contents(image.axis(0)*image.axis(1));
    this->readFits(filepath,contents);

    if( Pj != image.axis(0) ){
      std::cout << std::endl << std::endl << "PROBLEM!!! PSF DIMENSIONS DON'T MATCH" << std::endl << std::endl;
    }
    if( Pi != image.axis(1) ){
      std::cout << std::endl << std::endl << "PROBLEM!!! PSF DIMENSIONS DON'T MATCH" << std::endl << std::endl;
    }



    // Report on the peak location of the PSF
    double vmax = contents[0];
    double imax = 0;
    for(int i=1;i<Pj*Pi;i++){
      if( contents[i] > vmax ){
	vmax = contents[i];
	imax = i;
      }
    }

    if( Pi%2 == 0 && Pj%2 == 0 ){
      int Pi2 = Pi/2;
      int Pj2 = Pj/2;
      int indices[] = {(Pi2-1)*Pj+Pj2-1,(Pi2-1)*Pj+Pj2,Pi2*Pj+Pj2-1,Pi2*Pj+Pj2};
      int loc[]     = {0,0,0,0};
      bool flag     =  true;
      for(int i=0;i<4;i++){
	if( indices[i] == imax ){
	  loc[i] = 1;
	  flag = false;
	}
      }
      if( flag ){
	std::cout << "Particular form of the PSF: it is even-even but the peak lies outside the four central pixels" << std::endl;
      } else {
	std::cout << std::endl << "PSF is even-even with the peak located at: " << std::endl;
	printf("%2d%2d\n%2d%2d\n",loc[0],loc[1],loc[2],loc[3]);
      }
    } else if( Pi%2 != 0 && Pj%2 != 0 ){
      int Pi2 = (Pi-1)/2;
      int Pj2 = (Pj-1)/2;
      int centre = Pi2*Pj + Pj2;
      if( centre != imax ){
	std::cout << "Particular form of the PSF: it is odd-odd but the peak is not the central pixel" << std::endl;	
      }
    } else {
      std::cout << "Particular form of the PSF: it it neither even-even nor odd-odd" << std::endl;
    }




    // Set cropped PSF and normalize
    double* blur = (double*) calloc(Ncropx*Ncropy,sizeof(double));
    int offset_y = (Pi - Ncropy)/2;
    int offset_x = (Pj - Ncropx)/2;
    int offset_tot = (offset_y)*Pj + offset_x;
    double sum = 0.0;
    for(int i=0;i<Ncropy;i++){
      for (int j=0;j<Ncropx;j++){
	blur[i*Ncropx+j] = contents[offset_tot+i*Pj+j];
	sum += blur[i*Ncropx+j];
      }
    }
    double fac = 1.0/sum;
    for(int i=0;i<Ncropx*Ncropy;i++){
      blur[i] *= fac;
    }    
    



    // Set quad-kernel and odd/even functions
    int Nquadx,Nquady;
    void (ImagePlane::*xfunc)(int,int,int,int,int&,int&,int&);
    void (ImagePlane::*yfunc)(int,int,int,int,int&,int&,int&);
    if( Ncropx%2 == 0){
      Nquadx = Ncropx/2;
      xfunc = &ImagePlane::setCroppedLimitsEven;
    } else {
      Nquadx = ceil(Ncropx/2.0);
      xfunc = &ImagePlane::setCroppedLimitsOdd;
    }
    if( Ncropy%2 == 0){
      Nquady = Ncropy/2;
      yfunc = &ImagePlane::setCroppedLimitsEven;
    } else {
      Nquady = ceil(Ncropy/2.0);
      yfunc = &ImagePlane::setCroppedLimitsOdd;
    }

    // Create the blurring matrix
    int Nleft,Nright,Ntop,Nbottom,crop_offsetx,crop_offsety,crop_offset;
    int ic,jc;
    double val;
    for(int i=0;i<this->Ni;i++){
      for(int j=0;j<this->Nj;j++){

	(this->*(xfunc))(j,Ncropx,this->Nj,Nquadx,Nleft,Nright,crop_offsetx);
	(this->*(yfunc))(i,Ncropy,this->Ni,Nquady,Ntop,Nbottom,crop_offsety);
	crop_offset = crop_offsety*Ncropx + crop_offsetx;

	for(int ii=i-Ntop;ii<i+Nbottom;ii++){
	  ic = ii - i + Ntop;
	  for(int jj=j-Nleft;jj<j+Nright;jj++){
	    jc = jj - j + Nleft;
  
	    val = blur[crop_offset + ic*Ncropx + jc];

	    this->B.tri.push_back({i*this->Nj+j,     ii*this->Nj+jj,     val });
	  }
	}

      }
    }

    free(blur);

  }
}
  
void ImagePlane::setCroppedLimitsEven(int k,int Ncrop,int Nimg,int Nquad,int &Npre,int &Npost,int& offset){
  if( k < Nquad ){
    Npre   = k;
    Npost  = Nquad;
    offset = Nquad - k;
  } else if( k > (Nimg - Nquad - 1) ){
    Npre   = Nquad;
    Npost  = Nimg - k;
    offset = 0;
  } else {
    Npre   = Nquad;
    Npost  = Nquad;
    offset = 0;
  }
}

void ImagePlane::setCroppedLimitsOdd(int k,int Ncrop,int Nimg,int Nquad,int &Npre,int &Npost,int& offset){
  if( k < (Nquad - 1) ){
    Npre   = k;
    Npost  = Nquad;
    offset = Nquad - 1 - k;
  } else if( k > (Nimg - Nquad - 1) ){
    Npre   = Nquad - 1;
    Npost  = Nimg - k;
    offset = 0;
  } else {
    Npre   = Nquad-1;
    Npost  = Nquad;
    offset = 0;
  }
}

