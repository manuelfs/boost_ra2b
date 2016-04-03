#include <iostream>
#include <vector>

#include "TString.h"

#include "styles.hpp"
#include "utilities.hpp"

namespace {
  TString luminosity="2.3";
  TString plot_type=".pdf";
  TString plot_style="RA4";
}

using namespace std;

int main(){ 
  time_t begtime, endtime;
  time(&begtime);

  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  
  TString folder(bfolder+"/cms2r0/treemaker/trees/");


  vector<TString> s_t5hnc;
  s_t5hnc.push_back(folder+"*T5HH*1200_*200_*");
  vector<TString> s_t5hc;
  s_t5hc.push_back(folder+"*T5HH*1200_*950_*");


  // Reading ntuples
  vector<sfeats> Samples; 
  Samples.push_back(sfeats(s_t5hnc, "T5HH(1200,200)", 2, 1));
  Samples.push_back(sfeats(s_t5hc,  "T5HH(1200,950)", 2, 2));

  vector<int> ra2b_sam;
  unsigned nsam(Samples.size());
  for(unsigned sam(0); sam < nsam; sam++){
    ra2b_sam.push_back(sam);
  } // Loop over samples



  vector<hfeats> vars;

  // Number of Higgs tags
  vars.push_back(hfeats("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.4&&fjets_csv1>.605)",3,-0.49,2.5, 
			ra2b_sam, "Higgs tags (90<m_{J}<140, #tau_{2}/#tau_{1}<0.4, n_{b}^{L} #geq 1)",
			"Sum$(abs(mc_id)==5&&mc_mom==25)>=3&&Sum$(mc.Pt()>300&&mc_id==25)>=1"));
  vars.back().whichPlots = "3"; vars.back().normalize = true; 

  vars.push_back(hfeats("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.4&&fjets_csv1>.605&&fjets_csv2>.605)",3,-0.49,2.5, 
			ra2b_sam, "Higgs tags (90<m_{J}<140, #tau_{2}/#tau_{1}<0.4, n_{b}^{L} #geq 2)",
			"Sum$(abs(mc_id)==5&&mc_mom==25)>=3&&Sum$(mc.Pt()>300&&mc_id==25)>=1"));
  vars.back().whichPlots = "3"; vars.back().normalize = true; 

  vars.push_back(hfeats("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.4&&fjets_csv1>.605)",3,-0.49,2.5, 
			ra2b_sam, "Higgs tags (90<m_{J}<140, #tau_{2}/#tau_{1}<0.4, n_{b}^{L} #geq 1)","1"));
  vars.back().whichPlots = "3"; vars.back().normalize = true; 

  // Higgs mass
  vars.push_back(hfeats("fjets_pm",100,0,300, 
			ra2b_sam, "m_{J} [GeV]","fjets_tau21<0.4&&fjets_csv1>0.89"));
  vars.back().whichPlots = "3"; vars.back().normalize = true; 

  vars.push_back(hfeats("fjets_pm",100,0,300, 
			ra2b_sam, "m_{J} [GeV]","fjets_tau21<0.4&&fjets_csv1>0.89&&fjets_csv2>0.89"));
  vars.back().whichPlots = "3"; vars.back().normalize = true; 

  // Higgs pT
  vars.push_back(hfeats("mc.Pt()",55,0,1100, ra2b_sam, "Higgs p_{T} [GeV]","mc_id==25"));
  vars.back().whichPlots = "3"; vars.back().normalize = true; 

  plot_distributions(Samples, vars, luminosity, plot_type, plot_style, "nocuts_higgs",false,true);

  time(&endtime); 
  cout<<endl<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
