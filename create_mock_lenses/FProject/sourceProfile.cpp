#include "sourceProfile.hpp"

#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <fstream>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

#include <CGAL/Polygon_2_algorithms.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel            K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int,K>    Vb;
typedef CGAL::Triangulation_face_base_with_info_2<unsigned int,K>      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                    Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>                          Delaunay;
typedef K::Point_2                                                     Point;

//Abstract class: BaseProfile
//===============================================================================================================
void BaseProfile::profile(int Sj,int Si,double* sx,double* sy,double* s){
  for(int i=0;i<Si;i++){
    for(int j=0;j<Sj;j++){
      s[i*Sj+j] += this->value(sx[i*Sj+j],sy[i*Sj+j]);
    }
  }
}

void BaseProfile::ellipticalContour(){
  double dt     = (2.*this->pi)/(this->npoly);
  double cosphi = cos(this->pars["pa"]*this->fac);//in rad
  double sinphi = sin(this->pars["pa"]*this->fac);//in rad
  double size   = this->pars["r_eff"]/2.;
  
  for(int i=0;i<this->npoly;i++){
    this->x[i] = ( size*cos(i*dt)*cosphi - size*this->pars["q"]*sin(i*dt)*sinphi ) + this->pars["x0"];
    this->y[i] = ( size*cos(i*dt)*sinphi + size*this->pars["q"]*sin(i*dt)*cosphi ) + this->pars["y0"];
  }
}

//Derived class from BaseProfile: Sersic
//===============================================================================================================
Sersic::Sersic(std::map<std::string,Nlpar*> nlpars){
  this->pars["n"]     = nlpars["n"]->val;
  this->pars["r_eff"] = nlpars["r_eff"]->val;
  this->pars["i_eff"] = nlpars["i_eff"]->val;
  this->pars["q"]     = nlpars["q"]->val;
  this->pars["x0"]    = nlpars["x0"]->val;
  this->pars["y0"]    = nlpars["y0"]->val;
  this->pars["pa"]    = nlpars["pa"]->val;
  this->x = (double*) calloc(this->npoly,sizeof(double));
  this->y = (double*) calloc(this->npoly,sizeof(double));
  ellipticalContour();
}

double Sersic::value(double x,double y){
  double bn = 1.9992*this->pars["n"] - 0.3271;//From Capaccioli 1989
  double u,v,r,fac2;
  double cosphi = cos(this->pars["pa"]*this->fac);
  double sinphi = sin(this->pars["pa"]*this->fac);
  
  u =   (x - this->pars["x0"])*cosphi + (y - this->pars["y0"])*sinphi;
  v = - (x - this->pars["x0"])*sinphi + (y - this->pars["y0"])*cosphi;
  r = sqrt(this->pars["q"]*this->pars["q"]*u*u + v*v);
  fac2 = pow(r/this->pars["r_eff"],1./this->pars["n"]);
  return this->pars["i_eff"]*exp(-bn*fac2 - 1);
}

//Derived class from BaseProfile: ProGauss
//===============================================================================================================
ProGauss::ProGauss(std::map<std::string,Nlpar*> nlpars){
  this->pars["r_eff"] = nlpars["r_eff"]->val;
  this->pars["i_eff"] = nlpars["i_eff"]->val;
  this->pars["q"]     = nlpars["q"]->val;
  this->pars["x0"]    = nlpars["x0"]->val;
  this->pars["y0"]    = nlpars["y0"]->val;
  this->pars["pa"]    = nlpars["pa"]->val;
  this->x = (double*) calloc(this->npoly,sizeof(double));
  this->y = (double*) calloc(this->npoly,sizeof(double));
  ellipticalContour();
}

double ProGauss::value(double x,double y){
  double u,v,r2;
  double cosphi = cos(this->pars["pa"]*this->fac);
  double sinphi = sin(this->pars["pa"]*this->fac);  
  double sdev   = 2*this->pars["r_eff"]*this->pars["r_eff"];
  
  u =   (x - this->pars["x0"])*cosphi + (y - this->pars["y0"])*sinphi;
  v = - (x - this->pars["x0"])*sinphi + (y - this->pars["y0"])*cosphi;
  //    u =   x*cosphi + y*sinphi;
  //    v = - x*sinphi + y*cosphi;
  r2 = (this->pars["q"]*this->pars["q"]*u*u + v*v)/sdev;
  //    return (this->ieff*exp(-r2)/(sqrt(sdev*3.14159)));
  return this->pars["i_eff"]*exp(-r2);
}

//Derived class from BaseProfile: fromFITS
//===============================================================================================================
fromFITS::fromFITS(std::string filename,int Ni,int Nj,double height,double width,double x0,double y0){
  this->mySource = new ImagePlane(filename,Ni,Nj,height,width);
  // ImagePlane sets the coordinate origin in the center of the image, so I can re-position it here.
  for(int i=0;i<this->mySource->Nm;i++){
    this->mySource->x[i] -= x0;
    this->mySource->y[i] -= y0;
  }
  this->mySource->xmin -= x0;
  this->mySource->xmax -= x0;
  this->mySource->ymin -= y0;
  this->mySource->ymax -= y0;
}

double fromFITS::value(double x,double y){
  if( this->mySource->xmin < x && x < this->mySource->xmax && this->mySource->ymin < y && y < this->mySource->ymax ){
    // Source and Image grids MUST be the same, therefore I just need to match the right pixels (no interpolation)
    int i = (int) floor((this->mySource->ymax - y)*this->mySource->Ni/this->mySource->height); // y-axis is reflected
    int j = (int) floor((x - this->mySource->xmin)*this->mySource->Nj/this->mySource->width);
    return this->mySource->img[i*this->mySource->Nj + j];
  } else {
    return 0;
  }
}

//Derived class from BaseProfile: Voronoi
//===============================================================================================================
myVoronoi::myVoronoi(std::string filename){
  this->x = (double*) calloc(this->npoly,sizeof(double));
  this->y = (double*) calloc(this->npoly,sizeof(double));

  std::ifstream file(filename);
  std::string line;
  double val,x,y;

  while( std::getline(file,line) ){
    std::istringstream ss(line);

    ss >> val;
    values.push_back(val);

    std::vector<double> xvec;
    std::vector<double> yvec;
    while( !ss.eof() ){
      ss >> x >> y;
      //      std::cout << x << " " << y << "   ";
      xvec.push_back(x);
      yvec.push_back(y);
    }
    //    std::cout << std::endl;
    sizes.push_back(xvec.size());

    Point* cell = (Point*) malloc(xvec.size()*sizeof(Point));
    for(int i=0;i<xvec.size();i++){
      cell[i] = Point(xvec[i],yvec[i]);
    }
    cells.push_back(cell);
  }
  Ncells = cells.size();


  // check if the polygons are simple.
  for(int i=0;i<Ncells;i++){
    if( !CGAL::is_simple_2(cells[i],cells[i]+sizes[i],K()) ){
      std::cout << i << " polygon is not simple. Something is wrong in the provided list of Voronoi cells." << std::endl;
    }
  }
}

double myVoronoi::value(double x,double y){
  double val = 0.0;

  for(int i=0;i<Ncells;i++){
    if( CGAL::bounded_side_2(cells[i],cells[i]+sizes[i],Point(x,y),K()) == CGAL::ON_UNBOUNDED_SIDE ){
      continue;
    } else {
      val = values[i];
      break;
    }
  }

  return val;
}


//Derived class from BaseProfile: Delaunay
//===============================================================================================================
myDelaunay::myDelaunay(std::string filename){

  // Read v,x,y from file to the class variables
  std::ifstream file(filename);
  std::string line;
  double x,y,v;
  std::vector<double> xvec;
  std::vector<double> yvec;
  std::vector<double> vvec;

  while( std::getline(file,line) ){
    std::istringstream ss(line);
    ss >> v >> x >> y;
    xvec.push_back(x);
    yvec.push_back(y);
    vvec.push_back(v);
  }

  int N = xvec.size();
  this->N = N;
  this->x   = (double*) calloc(N,sizeof(double));
  this->y   = (double*) calloc(N,sizeof(double));
  this->src = (double*) calloc(N,sizeof(double));
  for(int i=0;i<N;i++){
    this->x[i]   = xvec[i];
    this->y[i]   = yvec[i];
    this->src[i] = vvec[i];
  }


  // Create the Dealaunay triangulation
  std::vector< std::pair<Point,int> > points;
  for(int i=0;i<N;i++){
    points.push_back( std::make_pair(Point(this->x[i],this->y[i]),i) );
  }

  Delaunay triangulation;
  triangulation.insert(points.begin(),points.end());


  //Get each Delaunay triangle in my own struct, and number them
  //[constructing this->triangles]
  Delaunay::Finite_faces_iterator fit;
  Delaunay::Face_handle face;
  atriangle triangle;
  int index = 0;
  this->triangles.resize( triangulation.number_of_faces() );
  for(fit=triangulation.finite_faces_begin();fit!=triangulation.finite_faces_end();fit++){
    face = fit;
    face->info() = index;

    triangle.a = (int) face->vertex(0)->info();
    triangle.b = (int) face->vertex(1)->info();
    triangle.c = (int) face->vertex(2)->info();
    
    this->triangles[index] = triangle;
    index++;
  }

  /*
  FILE* fh = fopen("triangles.dat","w");
  for(int q=0;q<this->triangles.size();q++){
    fprintf(fh,"%10.5f%10.5f%10.5f",this->x[this->triangles[q].a],this->x[this->triangles[q].b],this->x[this->triangles[q].c]);
    fprintf(fh,"%10.5f%10.5f%10.5f",this->y[this->triangles[q].a],this->y[this->triangles[q].b],this->y[this->triangles[q].c]);
    fprintf(fh,"\n");
  }
  fclose(fh);
  */

  // Get the convex hull of the triangulation
  Delaunay::Vertex_circulator vc = triangulation.incident_vertices(triangulation.infinite_vertex());
  Delaunay::Vertex_circulator done(vc);
  std::vector<Point> dum;
  do{
    dum.push_back(vc->point());
  }while( ++vc != done );

  this->ch_size = dum.size();
  this->convex_hull = (Point*) malloc(this->ch_size*sizeof(Point));
  for(int i=0;i<this->ch_size;i++){
    this->convex_hull[i] = dum[i];
  }

}


double myDelaunay::value(double x,double y){
  double val;


  // If the point is outside the convex hull of the triangulation set its value to zero, else interpolate within the triangle it in.
  if( CGAL::bounded_side_2(this->convex_hull,this->convex_hull+this->ch_size,Point(x,y),K()) == CGAL::ON_UNBOUNDED_SIDE ){

    val = 0.0;

  } else {
    double wa,wb,wc;
    double ybc,xac,xcb,yac,xxc,yyc,den;
    atriangle triangle;

    for(int j=0;j<this->triangles.size();j++){
      triangle = this->triangles[j];
      
      ybc = this->y[triangle.b] - this->y[triangle.c];//(yb-yc)
      xac = this->x[triangle.a] - this->x[triangle.c];//(xa-xc)
      xcb = this->x[triangle.c] - this->x[triangle.b];//(xc-xb)
      yac = this->y[triangle.a] - this->y[triangle.c];//(ya-yc)
      xxc = x                   - this->x[triangle.c];//(x -xc)
      yyc = y                   - this->y[triangle.c];//(y -yc)
      den = ybc*xac + xcb*yac;
      
      wa = ( ybc*xxc+xcb*yyc)/den;
      wb = (-yac*xxc+xac*yyc)/den;
      wc = 1.0 - wa - wb;
      
      if( 0.0 <= wa && wa <= 1.0 && 0.0 <= wb && wb <= 1.0 && 0.0 <= wc && wc <= 1.0 ){
	val = wa*this->src[triangle.a] + wb*this->src[triangle.b] + wc*this->src[triangle.c];
	break;
      }
    }

  }


  return val;
}
