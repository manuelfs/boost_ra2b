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
  
  TString folder(bfolder+"/cms2r0/treemaker/skim/");


  vector<TString> s_t5hnc;
  s_t5hnc.push_back(folder+"*T5HH*1200_*200*");
  vector<TString> s_t5hc;
  s_t5hc.push_back(folder+"*T5HH*1200_*950*");
  vector<TString> s_tt;
  s_tt.push_back(folder+"*_TTJets*Lept*");
  vector<TString> s_wjets;
  s_wjets.push_back(folder+"*_WJetsToLNu*");
  vector<TString> s_znunu;
  s_znunu.push_back(folder+"*ZJetsToNuNu_HT*");
  vector<TString> s_qcd;
  s_qcd.push_back(folder+"*_QCD_HT*");
  vector<TString> s_other;
  s_other.push_back(folder+"*DYJetsToLL_M-50_HT*");
  s_other.push_back(folder+"*_GJets_HT*");
  s_other.push_back(folder+"*_ST_*");
  s_other.push_back(folder+"*_TTGJets*");
  s_other.push_back(folder+"*_TTZTo*");


  // Reading ntuples
  vector<sfeats> Samples; 
  Samples.push_back(sfeats(s_t5hnc, "T5HH(1200,200)",	 ra2b::c_t5hh, 1)); Samples.back().isSig = true;
  Samples.push_back(sfeats(s_t5hc,  "T5HH(1200,950)",	 ra2b::c_t5hh, 2)); Samples.back().isSig = true;
  Samples.push_back(sfeats(s_tt, "t#bar{t}",		 ra2b::c_tt_1l));
  Samples.push_back(sfeats(s_znunu, "Z#rightarrow#nu#nu",ra2b::c_znunu));
  Samples.push_back(sfeats(s_wjets, "W+jets",		 ra2b::c_wjets));
  Samples.push_back(sfeats(s_qcd, "QCD",		 ra2b::c_qcd));
  Samples.push_back(sfeats(s_other, "Other",		 ra2b::c_other)); 

  vector<int> ra2b_sam;
  unsigned nsam(Samples.size());
  for(unsigned sam(0); sam < nsam; sam++){
    ra2b_sam.push_back(sam);
  } // Loop over samples



  vector<hfeats> vars;

  vars.push_back(hfeats("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.4&&fjets_csv1>.605)",3,-0.49,2.5, 
			ra2b_sam, "Higgs tags (90<m_{J}<140, #tau_{2}/#tau_{1}<0.4, n_{b}^{L} #geq 1)",
			"nleps==0&&ht>500&&mht>400&&njets>=6&&dphi1>0.5&&ntrks==0&&nbm>=1"));
  vars.back().whichPlots = "12"; 

  vars.push_back(hfeats("ht",29,500,3400,ra2b_sam, 
			"H_{T} [GeV]","nleps==0&&ht>500&&mht>200&&njets>=4&&dphi1>0.5&&ntrks==0"));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("mht",18,200,1100,ra2b_sam, 
			"H_{T}^{miss} [GeV]","nleps==0&&ht>500&&mht>200&&njets>=4&&dphi1>0.5&&ntrks==0"));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("njets",10,3.5,13.5,ra2b_sam, 
			"N_{jet}","nleps==0&&ht>500&&mht>200&&njets>=4&&dphi1>0.5&&ntrks==0"));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("nbm",7,-0.5,6.5,ra2b_sam, 
			"N_{b-jet}","nleps==0&&ht>500&&mht>200&&njets>=4&&dphi1>0.5&&ntrks==0"));
  vars.back().whichPlots = "12"; 

  plot_distributions(Samples, vars, luminosity, plot_type, plot_style, "",false,true);

  time(&endtime); 
  cout<<endl<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
