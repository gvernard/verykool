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


using std::cout;
using std::endl;




double myRandomNumber(double min,double max){
  double r = (double) rand() / (double)RAND_MAX;
  return min + r * (max - min);
}


int main(int argc,char* argv[]){

  //=============== BEGIN:PARSE INPUT =======================
  //Read in the json parameters
  Json::Value root;
  std::ifstream fin(argv[1]);
  fin >> root;

  std::string output = root["output"].asString();
  int flag_print_all = root["print_all"].asInt();

  const Json::Value jiplane = root["iplane"];
  const Json::Value jsource = root["source"];
  const Json::Value jlenses = root["lenses"];
  const Json::Value jphysical = root["physical"];
  //================= END:PARSE INPUT =======================




  //=============== BEGIN:CREATE THE LENSES ====================
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
  //================= END:CREATE THE LENSES ====================





  //=============== BEGIN:CREATE THE SOURCES =======================
  BaseProfile* mysource = NULL;
  std::string smodel = jsource["type"].asString();
  if( smodel == "analytic" ){

    std::vector<std::string> names;
    std::vector<std::map<std::string,double> > all_pars;
    for(int i=0;i<jsource["pars"].size();i++){
      std::string name = jsource["pars"][i]["type"].asString();
      names.push_back(name);

      std::map<std::string,double> pars;
      const Json::Value::Members jpars = jsource["pars"][i].getMemberNames();
      for(int j=0;j<jpars.size();j++){
	if( jpars[j] != "type" ){
	  pars[jpars[j]] = jsource["pars"][i][jpars[j]].asDouble();
	}
      }
      all_pars.push_back(pars);
    }
    mysource = new Analytic(names,all_pars);

  } else if( smodel == "delaunay" ){

    std::string filename = jsource["filename"].asString();
    mysource = new myDelaunay(filename);

  } else if( smodel == "fromfits" ){

    std::string filename = jsource["filename"].asString();
    int Ni               = jsource["Ni"].asInt();
    int Nj               = jsource["Nj"].asInt();
    double height        = jsource["height"].asDouble();
    double width         = jsource["width"].asDouble();
    double x0            = jsource["x0"].asDouble();
    double y0            = jsource["y0"].asDouble();
    mysource = new fromFITS(filename,Ni,Nj,height,width,x0,y0);

  } else {

    std::cout << "Unknown source profile type" << std::endl;
    return 1;

  }

  mysource->outputProfile(output + "vkl_source.fits");
  //================= END:CREATE THE SOURCES =======================




  
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

  //Initialize image plane
  ImagePlane mysim(jiplane["inf_x"].asInt(),jiplane["inf_y"].asInt(),jiplane["width"].asDouble(),jiplane["height"].asDouble());

  double xdefl,ydefl;
  for(int i=0;i<mysim.Ni*mysim.Nj;i++){
    mycollection->all_defl(mysim.x[i],mysim.y[i],xdefl,ydefl);
    mysim.img[i] = mysource->value(xdefl,ydefl);
    //  std::cout << mysim.img[i] << std::endl;
  }
  //================= END:PRODUCE IMAGE USING RAY-SHOOTING =======================





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
  //Lensed image
  mysim.writeImage(output+"image_clean.fits");

  //json output
  //  printf("{\"class\": \"%s\", \"mag\": %10.5f, \"x0\": %10.5f, \"y0\": %10.5f}\n",myclass.c_str(),mag,sources[0]->pars["x0"],sources[0]->pars["y0"]);


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
  delete(mysource);
  //================= END:CLEAN UP =======================


  return 0;
}
