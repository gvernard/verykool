#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <ctime>

#include "json/json.h"

#include "polygons.hpp"
#include "polyclip.hpp"
#include "initFuncs.hpp"
#include "sourceProfile.hpp"
#include "fitsHeader.hpp"


#include "nonLinearPars.hpp"
#include "massModels.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"


using std::cout;
using std::endl;




double myRandomNumber(double min,double max){
  double r = (double) rand() / (double)RAND_MAX;
  return min + r * (max - min);
}


int main(int argc,char* argv[]){

  //=============== BEGIN:PARSE INPUT =======================
  //  params p;
  //  if ( setParams(argc,argv,p) ){//if successfull, setParams has already set the parameters.
  //    return -1;
  //  }

  //Read in the non-linear parameters
  Json::Value root;
  if( argc > 1 ){//read json object from file
    std::ifstream fin(argv[1]);
    fin >> root;
  } else {//read json object from standard input
    std::cin >> root;
  }

  //Read output location
  std::string output = root["output"].asString();
  int flag_print_all = root["print_all"].asInt();

  //Read information on the triangles that will be used
  const Json::Value jtriangles = root["triangles"];

  //Read the source plane properties as a map of strings
  std::map<std::string,std::string> splane;
  int flag_splane = 0;
  if( root.isMember("splane") ){
    const Json::Value::Members jmembers = root["splane"].getMemberNames();
    for(int i=0;i<jmembers.size();i++){
      splane[jmembers[i]] = root["splane"][jmembers[i]].asString();
    }
    splane["reg"] = "identity"; // dum argument
    flag_splane = 1;
  }

  //Read the image properties
  const Json::Value jiplane = root["iplane"];

  //Read the properties of the source(s)
  const Json::Value jsources = root["sources"];

  //Read the properties of the lenses
  const Json::Value jlenses = root["lenses"];

  //Read the general physical properties (g,phi)
  const Json::Value jphysical = root["physical"];
  //================= END:PARSE INPUT =======================




  //=============== BEGIN:INITIALIZATION =======================
  //Initialize image plane
  ImagePlane mysim(jiplane["inf_x"].asInt(),jiplane["inf_y"].asInt(),jiplane["width"].asDouble(),jiplane["height"].asDouble());

  //Initialize mass model physical parameters
  std::vector<Nlpar*> ext_nlpars;
  for(int k=0;k<jphysical.size();k++){
    ext_nlpars.push_back( new Nlpar(jphysical[k]["nam"].asString(),0,0,jphysical[k]["val"].asDouble(),0,0,0) ); // only nam and val have meaning in this call
  }
  CollectionMassModels* mycollection = new CollectionMassModels(ext_nlpars);



  for(int i=0;i<ext_nlpars.size();i++){
    delete(ext_nlpars[i]);
  }
  ext_nlpars.clear();

  //Initialize mass model parameters for each lens
  mycollection->models.resize(jlenses.size());
  for(int k=0;k<jlenses.size();k++){
    std::string mmodel = jlenses[k]["subtype"].asString();

    if( mmodel == "pert" ){

      const Json::Value::Members jmembers = jlenses[k]["pars"].getMemberNames();
      std::map<std::string,std::string> pars;
      for(int i=0;i<jmembers.size();i++){
	pars[jmembers[i]] = jlenses[k]["pars"][jmembers[i]].asString();
      }
      mycollection->models[k] = FactoryMassModel::getInstance()->createMassModel(mmodel,pars);


      /*
      // This is to print the gradients of the potential
      Pert* pert = dynamic_cast<Pert*>(mycollection->models[k]);
      ImagePlane dumx(pert->Ni,pert->Nj,pert->width,pert->height);
      ImagePlane dumy(pert->Ni,pert->Nj,pert->width,pert->height);
      ImagePlane dumsai(pert->Ni,pert->Nj,pert->width,pert->height);
      for(int i=0;i<pert->Ni*pert->Nj;i++){
	dumx.img[i] = pert->dpsidx[i];
	dumy.img[i] = pert->dpsidy[i];
	dumsai.img[i] = hypot(pert->dpsidx[i],pert->dpsidy[i]);
      }
      dumx.writeImage("dpsi_x.fits");
      dumy.writeImage("dpsi_y.fits");
      dumsai.writeImage("dpsi_sai.fits");
      */



    } else {

      const Json::Value jnlpars = jlenses[k]["nlpars"];
      std::vector<Nlpar*> nlpars;
      for(int i=0;i<jnlpars.size();i++){
	nlpars.push_back( new Nlpar(jnlpars[i]["nam"].asString(),0,0,jnlpars[i]["val"].asDouble(),0,0,0) ); // only nam and val have meaning in this call
      }
      mycollection->models[k] = FactoryMassModel::getInstance()->createMassModel(mmodel,nlpars);
      
      for(int i=0;i<nlpars.size();i++){
	delete(nlpars[i]);
      }
      nlpars.clear();

    }
  }
  //================= END:INITIALIZATION =======================





  //=============== BEGIN:CREATE THE SOURCES =======================
  std::vector<BaseProfile*> sources(jsources.size());
  for(int k=0;k<jsources.size();k++){
    std::string pmodel = jsources[k]["subtype"].asString();
    
    if( pmodel == "voronoi" || pmodel == "delaunay" || pmodel == "fromfits" ){

      const Json::Value::Members jmembers = jsources[k]["pars"].getMemberNames();
      std::map<std::string,std::string> pars;
      for(int i=0;i<jmembers.size();i++){
	pars[jmembers[i]] = jsources[k]["pars"][jmembers[i]].asString();
      }
      sources[k] = FactoryProfile::getInstance()->createProfile(pmodel,pars);

    } else {

      //Create source profile objects
      const Json::Value jnlpars = jsources[k]["nlpars"];
      std::map<std::string,Nlpar*> nlpars;
      for(int i=0;i<jnlpars.size();i++){
	nlpars[jnlpars[i]["nam"].asString()] = new Nlpar(jnlpars[i]["nam"].asString(),0,0,jnlpars[i]["val"].asDouble(),0,0,0); // only nam and val have meaning in this call
      }
      if( nlpars["randomize"]->val ){
	nlpars["x0"]->val = myRandomNumber(-1.5,1.5);
	nlpars["y0"]->val = myRandomNumber(-1.5,1.5);
      }
      sources[k] = FactoryProfile::getInstance()->createProfile(pmodel,nlpars);

    }
  }
  //================= END:CREATE THE SOURCES =======================






  /*
  //=============== BEGIN:PLACE THE SOURCE =======================
  //Set triangles on the image plane
  //  std::vector<triangle> triangles = srcTriangles(mysim.Ni,mysim.Nj,mysim.width,mysim.height,mysim.x,mysim.y,mymodel);
  //  std::cout << triangles.size() << std::endl;

  std::vector<triangle> triangles = createTriangles(jtriangles["Nx"].asInt(),jtriangles["Ny"].asInt(),mysim.width,mysim.height);
  std::cout << triangles.size() << std::endl;

  triangles = selectTriangles(triangles,jtriangles["sel_radius"].asDouble());
  std::cout << triangles.size() << std::endl;

  deflectTriangles(triangles,mycollection);



  //Place the source on the source plane (or a contour of given source brightness)
  srand((unsigned)time(0));
  std::vector<BaseProfile*> sources(jsources.size());
  std::vector<triangle> insiders;
  std::vector<triangle> overlaps;
  int flag = 1;
  std::string myclass = "";
  poly_t sub;
  poly_t cli;
  poly res;


  for(int k=0;k<jsources.size();k++){
    int index = 0;
    while( flag ){
      std::cout << "skata" << std::endl;

      //Create source profile objects
      const Json::Value jnlpars = jsources[k]["nlpars"];
      std::map<std::string,Nlpar*> nlpars;
      for(int i=0;i<jnlpars.size();i++){
	nlpars[jnlpars[i]["nam"].asString()] = new Nlpar(jnlpars[i]["nam"].asString(),0,0,jnlpars[i]["val"].asDouble(),0,0,0); // only nam and val have meaning in this call
      }
      if( nlpars["randomize"]->val ){
	nlpars["x0"]->val = myRandomNumber(-1.5,1.5);
	nlpars["y0"]->val = myRandomNumber(-1.5,1.5);
      }

      std::string pmodel = jsources[k]["subtype"].asString();
      sources[k] = FactoryProfile::getInstance()->createProfile(pmodel,nlpars);
      

      //Find all the triangles within this area
      point  p;
      insiders.clear();
      vec_t s[sources[k]->npoly];
      for(int i=0;i<sources[k]->npoly;i++){
	s[i] = {sources[k]->x[i],sources[k]->y[i]};
      }

      sub = {sources[k]->npoly,0,s};
      for(int i=0;i<triangles.size();i++){
	vec_t c[]  = {{triangles[i].a.x,triangles[i].a.y},{triangles[i].b.x,triangles[i].b.y},{triangles[i].c.x,triangles[i].c.y}};
	cli = {3,0,c};
	
	res = poly_clip(&sub,&cli);
	if( res->len > 2 ){
	  insiders.push_back(triangles[i]);
	}
      }
      cout << "Number of triangles within the source is: " << insiders.size() << endl;
      if( insiders.size() == 0 ){
	continue;
      }
      
      
      //Find the clipped polygons between all the triangles
      std::vector<poly> clipped;
      for(int i=0;i<insiders.size()-1;i++){
	vec_t c[]  = {{insiders[i].a.x,insiders[i].a.y},{insiders[i].b.x,insiders[i].b.y},{insiders[i].c.x,insiders[i].c.y}};
	cli = {3,0,c};

	for(int j=i+1;j<insiders.size();j++){
	  vec_t s[]  = {{insiders[j].a.x,insiders[j].a.y},{insiders[j].b.x,insiders[j].b.y},{insiders[j].c.x,insiders[j].c.y}};
	  sub = {3,0,s};
	  res = poly_clip(&sub,&cli);
	  if( res->len > 2 ){
	    clipped.push_back(res);
	  }
	}
      }
      cout << "Number of clipped polygons is: " << clipped.size() << endl;
      if( clipped.size() == 0 ){
	continue;
      }
      
      
      //Check if any of the triangles overlap, and count the maximum number of overlapping triangles
      point center;
      int* counts = (int*) calloc(clipped.size(),sizeof(int));
      double xarr[3]; 
      double yarr[3];
      std::map<int,int> overlaps_map;
      for(int i=0;i<clipped.size();i++){
	
	double* xpoly = (double*) calloc(clipped[i]->len,sizeof(double));
	double* ypoly = (double*) calloc(clipped[i]->len,sizeof(double));
	for(int j=0;j<clipped[i]->len;j++){
	  xpoly[j] = clipped[i]->v[j].x;
	  ypoly[j] = clipped[i]->v[j].y;
	}
	findBarycenter(clipped[i]->len,xpoly,ypoly,center.x,center.y);
	//    fprintf(fh2,"%10.5f %10.5f\n",center.x,center.y);
	
	for(int j=0;j<insiders.size();j++){
	  xarr[0] = insiders[j].a.x;
	  xarr[1] = insiders[j].b.x;
	  xarr[2] = insiders[j].c.x;
	  yarr[0] = insiders[j].a.y;
	  yarr[1] = insiders[j].b.y;
	  yarr[2] = insiders[j].c.y;
	  
	  if( pnpoly(3,xarr,yarr,center.x,center.y) ){
	    counts[i]++;
	    overlaps_map[j] = 1;
	  }
	}

	free(xpoly);
	free(ypoly);
      }

      overlaps.clear();
      for(std::map<int,int>::iterator iterator=overlaps_map.begin();iterator!=overlaps_map.end();iterator++){
	overlaps.push_back( insiders[iterator->first] );
      }


      std::map<int,int> hist;
      for(int i=0;i<clipped.size();i++){
	if( counts[i] > 0 ){
	  if( hist.count(counts[i]) ){
	    hist[counts[i]]++;
	  } else {
	    hist[counts[i]] = 1;
	  }
	}
      }
      //      for(it iterator=hist.begin();iterator!=hist.end();iterator++){
      //	cout << "images: " << iterator->first << "   locations (triangles): " << iterator->second << endl;
      //      }
      cout << index << " -------------------------------------------" << endl;
      cout << " (x0,y0) = " << sources[k]->pars["x0"] << " , " << sources[k]->pars["y0"] << endl;
      cout << "images: 3   locations (triangles): " << hist[3] << endl;
      cout << "images: 5   locations (triangles): " << hist[5] << endl;



      //produce lensed image for each source
      //      double xdefl,ydefl;
      //      for(int i=0;i<mysim.Ni*mysim.Nj;i++){
      //	mymodel->defl(mysim.x[i],mysim.y[i],xdefl,ydefl);
      //	mysim.img[i] = 0;
      //	for(int j=0;j<sources.size();j++){
      //	  mysim.img[i] += sources[j]->value(xdefl,ydefl);
      //	}
      //      }
      //      char name [3];
      //      sprintf(name,"%03d",index);
      //      //    mysim.writeFits(output+name+".fits");
      //      writeArrayPngPP(output+name+".png",mysim.Nj,mysim.Ni,mysim.img);


      //Stopping criterion if source positions are given
      std::cout << nlpars["randomize"]->val << std::endl;
      if( nlpars["randomize"]->val == 0 ){
	cout << "given positions" << "    " << nlpars["randomize"]->val << endl;
	flag = 0;
      }

      //Stopping criterion for placing sources randomly
      if( hist.count(5) && hist[5] > 0 ){
	if( hist[5] > 2*hist[3] ){
	  myclass = "quad";
	} else {
	  myclass = "arc";
	}
	flag = 0;
      }

      index++;
      free(counts);
      for(int i=0;i<clipped.size();i++){
	poly_free(clipped[i]);
      }

      for(std::map<std::string,Nlpar*>::iterator iterator=nlpars.begin();iterator!=nlpars.end();iterator++){
	delete nlpars[iterator->first];
      }
      nlpars.clear();


      //      cout << "Number of overlapping triangles is: " << overlaps.size() << endl;
      //      if( overlaps.size() > insiders.size()/2. ){
      //	flag = 0;
      //      }


    }

    flag = 1;
    //    printf("%10.5f %10.5f\n",sources[k]->x0,sources[k]->y0);
  }
  //================= END:PLACE THE SOURCE =======================
  */

    


  //=============== BEGIN:CREATE SOURCE PLANE GRID =======================
  //  if( flag_splane ){
  //Initialize source plane grid (it is a map<string,string>)
  if( jsources[0]["subtype"].asString() == "fromfits" ){
    if( jsources[0]["pars"]["Nj"].asInt() > jsources[0]["pars"]["Ni"].asInt() ){
      splane["sx"] = jsources[0]["pars"]["Nj"].asString();
      splane["sy"] = jsources[0]["pars"]["Nj"].asString();
    } else {
      splane["sx"] = jsources[0]["pars"]["Ni"].asString();
      splane["sy"] = jsources[0]["pars"]["Ni"].asString();
    }
    if( jsources[0]["pars"]["width"].asDouble() > jsources[0]["pars"]["height"].asDouble() ){
      splane["size"] = jsources[0]["pars"]["width"].asString();
    } else {
      splane["size"] = jsources[0]["pars"]["height"].asString();
    }
  }
  BaseSourcePlane* mysource = FactorySourcePlane::getInstance()->createSourcePlane(splane);
  
  //  std::map<std::string,std::string> pars;
  //  pars["size"] = splane["size"];
  //  pars["x0"]   = splane["x0"];
  //  pars["y0"]   = splane["y0"];
  //    mysource->setGrid(pars);//Only for floating source plane
  //    mysource->boundPolygon();
  
  //set a fixed source profile
  for(int k=0;k<sources.size();k++){
    sources[k]->profile(mysource->Sj,mysource->Si,mysource->x,mysource->y,mysource->src);
  }
  
  //    mysource->normalize();
  
  // output source profile
  //mysource->outputSource(output+"source.png");
  //  std::cout << output << std::endl;
  mysource->outputSource(output + "vkl_source.fits");
  
  addFitsHeader(output + "vkl_source.fits",splane);
  
  
  delete(mysource);
  //  }
  //================= END:CREATE SOURCE PLANE GRID =======================




  
    /*    
  //=============== BEGIN:PRODUCE IMAGE USING MATRICES =======================
  //Construct lensing matrix
  mytable L;
  L.Ti = mysim.Ni * mysim.Nj;
  L.Tj = mysource->Si * mysource->Sj;
  mysource->constructL(&mysim,mymodel,&L);
  
  //convert sparse L to dense
  double* Lfull = (double*) calloc(mysim.Ni*mysim.Nj*mysource->Si*mysource->Sj,sizeof(double));
  for(int k=0;k<L.tri.size();k++){
    Lfull[L.tri[k].i*(mysource->Si*mysource->Sj) + L.tri[k].j] = L.tri[k].v;
  }

  //Multiply L with source to get the image
  for(int k=0;k<mysim.Ni*mysim.Nj;k++){
        mysim.img[k] = 0;
        for(int i=0;i<mysource->Si*mysource->Sj;i++){
          mysim.img[k] += Lfull[k*(mysource->Si*mysource->Sj) + i] * mysource->src[i];
        }
    //            for(int i=0;i<mysource->Si;i++){
    //              for(int j=0;j<mysource->Sj;j++){
    //		mysim.img[k] += Lfull[k*(mysource->Si*mysource->Sj) + i*mysource->Sj+j] * mysource->src[i*mysource->Sj+j];
    //		//mysim.img[k] += Lfull[k*(mysource->Si*mysource->Sj) + i*mysource->Sj+j] * mysource->src[(mysource->Si-1-i)*mysource->Sj+j];
    //	      }
    //	    }
  }
  //================= END:PRODUCE IMAGE USING MATRICES =======================
  mysim.writeFits(output+"matrix.fits");
    */


  //  double xdefl,ydefl;
  //  mycollection->all_defl(mysim.x[0],mysim.y[0],xdefl,ydefl);
  //  std::cout << sources[0]->value(xdefl,ydefl) << std::endl;



  //=============== BEGIN:PRODUCE IMAGE USING RAY-SHOOTING =======================
  // Set the derivatives of the potential perturbations directly
  /*
  Pert* pert = dynamic_cast<Pert*>(mycollection->models.back());
  std::ifstream input("/net/argo/data/users/gvernard/RESULTS/VKL_MODELS/test_modelling_perturbations/test_2_nopert/base_run/output/linear_perturbations_derivatives.dat");
  for(int i=0;i<pert->dpsi->Sm;i++){
    //input >> pert->dpsi_dx[i];
    pert->dpsi_dx[i] = 0.2;
  }
  for(int i=0;i<pert->dpsi->Sm;i++){
    //input >> pert->dpsi_dy[i];
    pert->dpsi_dy[i] = 0.2;
  }

  std::cout << mysim.xmin << " " << mysim.xmax << std::endl;
  std::cout << mysim.ymin << " " << mysim.ymax << std::endl;
  std::cout << pert->dpsi->xmin << " " << pert->dpsi->xmax << std::endl;
  std::cout << pert->dpsi->ymin << " " << pert->dpsi->ymax << std::endl;
  */


  double xdefl,ydefl;
  for(int i=0;i<mysim.Ni*mysim.Nj;i++){
    mycollection->all_defl(mysim.x[i],mysim.y[i],xdefl,ydefl);
    mysim.img[i] = 0;
    for(int j=0;j<sources.size();j++){
      mysim.img[i] += sources[j]->value(xdefl,ydefl);
    }
    //  std::cout << mysim.img[i] << std::endl;
  }
  //================= END:PRODUCE IMAGE USING RAY-SHOOTING =======================
  //  mysim.writeImage(output+"anal.fits");





    /*
  //=============== BEGIN:CALCULATE MAGNIFICATION ================================
  std::map<std::string,std::string> newpars;
  newpars["type"] = splane["type"];
  newpars["sx"]   = "100";
  newpars["sy"]   = "100";
  newpars["size"] = splane["size"];
  newpars["x0"]   = splane["x0"];
  newpars["y0"]   = splane["y0"];

  BaseSourcePlane* newsource = FactorySourcePlane::getInstance()->createSourcePlane(newpars);
  
  //set a fixed source profile
  for(int k=0;k<sources.size();k++){
    sources[k]->profile(newsource->Sj,newsource->Si,newsource->x,newsource->y,newsource->src);
  }

  double Is = 0.;
  double dss = (newsource->x[newsource->Sj-1]-newsource->x[0])*(newsource->y[0]-newsource->y[(newsource->Si-1)*newsource->Sj])/(newsource->Si*newsource->Sj);
  //  std::cout << "dx " << newsource->x[newsource->Sj-1]-newsource->x[0] << std::endl;
  //  std::cout << "dy " << newsource->y[0]-newsource->y[(newsource->Si-1)*newsource->Sj] << std::endl;
  //  double dss = (1.*1.)/(newsource->Si*newsource->Sj);
  //  std::cout << dss << std::endl;
  for(long i=0;i<newsource->Si*newsource->Sj;i++){
    Is += newsource->src[i] * dss;
  }
  double Ii = 0.;
  double dsi = mysim.width*mysim.height/(mysim.Ni*mysim.Nj);
  //  std::cout << dsi << std::endl;
  for(long i=0;i<mysim.Ni*mysim.Nj;i++){
    Ii += mysim.img[i] * dsi;
  }

  double mag = Ii/Is;
  //================= END:CALCULATE MAGNIFICATION ================================
  */






  //=============== BEGIN:OUTPUT =======================
  //Image of the source (have to set up source plane grid first)

  //Lensed image
  //  writeArrayPngPP(output+"image.png",mysim.Nj,mysim.Ni,mysim.img);
  mysim.writeImage(output+"image_clean.fits");

  //json output
  //  printf("{\"class\": \"%s\", \"mag\": %10.5f, \"x0\": %10.5f, \"y0\": %10.5f}\n",myclass.c_str(),mag,sources[0]->pars["x0"],sources[0]->pars["y0"]);


  if( flag_print_all ){
    FILE* fh;
    std::string filename;
    
    //    filename = output+"triangles.dat";
    //    fh = fopen(filename.c_str(),"w");
    //    for(int i=0;i<triangles.size();i++){
    //      fprintf(fh,"%10.5f %10.5f %10.5f ",triangles[i].a.x,triangles[i].b.x,triangles[i].c.x);
    //      fprintf(fh,"%10.5f %10.5f %10.5f\n",triangles[i].a.y,triangles[i].b.y,triangles[i].c.y);
    //    }
    //    fclose(fh);

    //    filename = output+"insiders.dat";
    //    fh = fopen(filename.c_str(),"w");
    //    for(int i=0;i<insiders.size();i++){
    //      fprintf(fh,"%10.5f %10.5f %10.5f ",insiders[i].a.x,insiders[i].b.x,insiders[i].c.x);
    //      fprintf(fh,"%10.5f %10.5f %10.5f\n",insiders[i].a.y,insiders[i].b.y,insiders[i].c.y);
    //    }
    //    fclose(fh);
    
    //    filename = output+"overlaps.dat";
    //    fh = fopen(filename.c_str(),"w");
    //    for(int i=0;i<overlaps.size();i++){
    //      fprintf(fh,"%10.5f %10.5f %10.5f ",overlaps[i].a.x,overlaps[i].b.x,overlaps[i].c.x);
    //      fprintf(fh,"%10.5f %10.5f %10.5f\n",overlaps[i].a.y,overlaps[i].b.y,overlaps[i].c.y);
    //    }
    //    fclose(fh);
    
    filename = output+"ellipses.dat";
    fh = fopen(filename.c_str(),"w");
    for(int j=0;j<sources.size();j++){
      for(int i=0;i<sources[j]->npoly;i++){
	fprintf(fh," %10.5f",sources[j]->x[i]);
      }
      fprintf(fh,"\n");
      for(int i=0;i<sources[j]->npoly;i++){
	fprintf(fh," %10.5f",sources[j]->y[i]);
      }
      fprintf(fh,"\n");
    }
    fclose(fh);
  }

  //Write lens parameters
  
//  std::ofstream myfile (output+"params.dat",std::ios::out);
//  for(it_type iterator=nlpars.begin();iterator!=nlpars.end();iterator++){
//    myfile << nlpars[iterator->first]->nam << " " << nlpars[iterator->first]->val << std::endl;
//    //    std::cout << nlpars[iterator->first]->nam << " " << nlpars[iterator->first]->val << std::endl;
//  }
//  myfile.close();
  
  //================= END:OUTPUT =======================





  //=============== BEGIN:CLEAN UP =======================
  delete(mycollection);

  for(int j=0;j<sources.size();j++){
    delete(sources[j]);
  }
  //================= END:CLEAN UP =======================


  return 0;
}
