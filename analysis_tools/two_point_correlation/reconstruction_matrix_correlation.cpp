#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "sourcePlane.hpp"
#include "sourceProfile.hpp"


int main(int argc,char* argv[]){
  std::string src_path  = argv[1];


  
  // Read in a Delaunay source
  myDelaunay* delaunay = new myDelaunay(src_path);
  AdaptiveSource* ada_source = new AdaptiveSource(delaunay->N,"curvature"); // the regularization scheme is dummy argument here!!!

  // Set the parameters of the source to the input one
  for(int i=0;i<delaunay->N;i++){
    ada_source->x[i]     = delaunay->x[i];
    ada_source->y[i]     = delaunay->y[i];
    ada_source->src[i]   = delaunay->src[i];
  }




  // Claculating and binning by r the covariance matrix of the reconstructed source

  double Nbins = 100;
  double rmax = 2.0; // arcsec
  double rmin = 0.0; // arcsec
  double dr = (rmax-rmin)/Nbins;
  double* bins = (double*) calloc(Nbins,sizeof(double));
  for(int i=0;i<Nbins;i++){
    bins[i] = rmin + dr*i;
  }

  double* vals_s = (double*) calloc(Nbins,sizeof(double));
  int* counts_s  = (int*) calloc(Nbins,sizeof(int));
  for(int k=0;k<ada_source->Sm;k++){
    for(int q=k+1;q<ada_source->Sm;q++){
      double r = hypot(ada_source->x[q]-ada_source->x[k],ada_source->y[q]-ada_source->y[k]);
      if( r < rmax ){
	int index = (int) floor(r/dr);
	vals_s[index] += ada_source->src[q]*ada_source->src[k];
	counts_s[index]++; 
      }  
    }
  }
  for(int i=0;i<Nbins;i++){
    if( counts_s[i] > 0 ){
      vals_s[i] /= counts_s[i];
    }
  }
  free(counts_s);


  
  FILE* fh = fopen("corr_reconstruction.dat","w");
  for(int i=0;i<Nbins;i++){
    fprintf(fh,"%5.3f %8.3f\n",bins[i],vals_s[i]);
  }
  fclose(fh);
  

  
  free(bins);
  free(vals_s);
  
  delete(delaunay);
  delete(ada_source);

  return 0;
}
