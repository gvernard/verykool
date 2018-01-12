#include "imagePlane.hpp"

#include <cmath>
#include <fstream>
#include <cstdlib>
#include <CCfits/CCfits>

#include "tableDefinition.hpp"


//ImagePlane class implementation
//============================================================================================
ImagePlane::ImagePlane(const std::string filepath,int i,int j,double w,double h){
  Nj     = i;
  Ni     = j;
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

  img = (double*) calloc(Nm,sizeof(double));
  x   = (double*) calloc(Nm,sizeof(double));
  y   = (double*) calloc(Nm,sizeof(double));
  active = (int*) calloc(Nm,sizeof(int));
  
  int i0    = floor(Ni/2);
  int j0    = floor(Nj/2);
  double di = width/(Ni-1);
  double dj = height/(Nj-1);
  //  double di = width/Ni;
  //  double dj = height/Nj;

  for(int ii=0;ii<Ni;ii++){
    for(int jj=0;jj<Nj;jj++){
      x[ii*Nj+jj]   =  (jj-j0)*di;
      y[ii*Nj+jj]   = -(ii-i0)*dj;//reflect y-axis
      img[ii*Nj+jj] = contents[ii*Nj+jj];
    }
  }
}

ImagePlane::ImagePlane(int i,int j,double w,double h){
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

  img = (double*) calloc(Nm,sizeof(double));
  x   = (double*) calloc(Nm,sizeof(double));
  y   = (double*) calloc(Nm,sizeof(double));
  active = (int*) calloc(Nm,sizeof(int));

  int i0    = floor(i/2);
  int j0    = floor(j/2);
  double di = width/(i-1);
  double dj = height/(j-1);

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

  img = (double*) calloc(Nm,sizeof(double));
  x   = (double*) calloc(Nm,sizeof(double));
  y   = (double*) calloc(Nm,sizeof(double));
  active = (int*) calloc(Nm,sizeof(int));
}

ImagePlane::~ImagePlane(){
  free(img);
  free(x);
  free(y);
  free(active);
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
      this->C.tri.push_back({i,i,1./pow(value,2)});
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
	this->C.tri.push_back({i,i,1./pow(value,2)});
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
      this->C.tri.push_back({i,j,1./pow(value,2)});
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
    double* blur = (double*) calloc(Ncropx*Ncropy,sizeof(double));


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


    int loffx,roffx,toffy,boffy;

    loffx = floor(Ncropx/2.);
    if( Ncropx % 2 == 0 ){
      roffx = floor(Ncropx/2.);
    } else {
      roffx = floor(Ncropx/2.) + 1;
    }
    toffy = floor(Ncropy/2.);
    if( Ncropy % 2 == 0 ){
      boffy = floor(Ncropy/2.);
    } else {
      boffy = floor(Ncropy/2.) + 1;
    }



    int offset = (floor(Pi/2.)-toffy)*Pi + (floor(Pj/2.)-loffx);
    for(int i=0;i<Ncropy;i++){
      for (int j=0;j<Ncropx;j++){
	blur[i*Ncropx+j] = contents[offset+i*Pi+j];
      }
    }

    /*
    for(int i=0;i<Ncropy;i++){
      for(int j=0;j<Ncropx;j++){
	std::cout << blur[i*Ncropx+j] << " ";
      }
      std::cout << std::endl;
    }
    */






    for(int i=0;i<toffy;i++){

      for(int j=0;j<loffx;j++){

 	for(int ii=-i;ii<boffy;ii++){
	  for(int jj=-j;jj<roffx;jj++){
	    this->B.tri.push_back({i*this->Nj+j,     (i+ii)*this->Nj+(j+jj),     blur[(toffy+ii)*Ncropx+(loffx+jj)]});
	  }
	}

      }

      for(int j=loffx;j<this->Nj-roffx;j++){

 	for(int ii=-i;ii<boffy;ii++){
	  for(int jj=-loffx;jj<roffx;jj++){
	    this->B.tri.push_back({  i*this->Nj+j,   (i+ii)*this->Nj+(j+jj),     blur[(toffy+ii)*Ncropx+(loffx+jj)]});//row-major
	  }
	}

     }

      for(int j=this->Nj-roffx;j<this->Nj;j++){

 	for(int ii=-i;ii<boffy;ii++){
	  for(int jj=-loffx;jj<this->Nj-j;jj++){
	    this->B.tri.push_back({i*this->Nj+j,     (i+ii)*this->Nj+(j+jj),     blur[(toffy+ii)*Ncropx+(loffx+jj)]});
	  }
	}

      }

    }
    


    
    for(int i=toffy;i<this->Ni-boffy;i++){

      for(int j=0;j<loffx;j++){

	for(int ii=-toffy;ii<boffy;ii++){
	  for(int jj=-j;jj<roffx;jj++){
	    this->B.tri.push_back({i*this->Nj+j,     (i+ii)*this->Nj+(j+jj),     blur[(toffy+ii)*Ncropx+(loffx+jj)]});
	  }
	}

      }

      for(int j=loffx;j<this->Nj-roffx;j++){
	
	//Here I assume that the order of the row (i) and column (j) indices is such to create a row-major sparse matrix
	for(int ii=-toffy;ii<boffy;ii++){
	  for(int jj=-loffx;jj<roffx;jj++){
	    this->B.tri.push_back({  i*this->Nj+j,   (i+ii)*this->Nj+(j+jj),     blur[(toffy+ii)*Ncropx+(loffx+jj)]});//row-major
	  }
	}
	
      }

      for(int j=this->Nj-roffx;j<this->Nj;j++){

	for(int ii=-toffy;ii<boffy;ii++){
	  for(int jj=-loffx;jj<this->Nj-j;jj++){
	    this->B.tri.push_back({i*this->Nj+j,     (i+ii)*this->Nj+(j+jj),     blur[(toffy+ii)*Ncropx+(loffx+jj)]});
	  }
	}

      }

    }




    for(int i=this->Ni-boffy;i<this->Ni;i++){

      for(int j=0;j<loffx;j++){
   
	for(int ii=-toffy;ii<this->Ni-i;ii++){
	  for(int jj=-j;jj<loffx;jj++){
	    this->B.tri.push_back({i*this->Nj+j,     (i+ii)*this->Nj+(j+jj),     blur[(toffy+ii)*Ncropx+(loffx+jj)]});
	  }
	}

      }

      for(int j=loffx;j<this->Nj-roffx;j++){
  
	for(int ii=-toffy;ii<this->Ni-i;ii++){
	  for(int jj=-loffx;jj<roffx;jj++){
	    this->B.tri.push_back({  i*this->Nj+j,   (i+ii)*this->Nj+(j+jj),     blur[(toffy+ii)*Ncropx+(loffx+jj)]});//row-major
	  }
	}

      }

      for(int j=this->Nj-roffx;j<this->Nj;j++){
  
	for(int ii=-toffy;ii<this->Ni-i;ii++){
	  for(int jj=-loffx;jj<this->Nj-j;jj++){
	    this->B.tri.push_back({i*this->Nj+j,     (i+ii)*this->Nj+(j+jj),     blur[(toffy+ii)*Ncropx+(loffx+jj)]});
	  }
	}

      }

    }

    free(blur);
    
  }
}
  
