#include <iostream>
#include <vector>

#include "TString.h"

#include "styles.hpp"
#include "utilities.hpp"

namespace {
  TString luminosity="2.3";
  TString plot_type=".pdf";
  TString plot_style="CMSPaper";
}

using namespace std;
using std::cout;
using std::endl;

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

  vars.push_back(hfeats("MHT",55,0,1100, ra2b_sam, "H_{T}^{miss} [GeV]","2",200));
  vars.back().whichPlots = "34"; vars.back().normalize = true; 

  plot_distributions(Samples, vars, luminosity, plot_type, plot_style, "paper",false);

  time(&endtime); 
  cout<<endl<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
