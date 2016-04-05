// plot_t5tttt: Compares kinematic distributions of T1tttt and T5tttt

#include <iostream>
#include <vector>
#include <ctime>

#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TGraph.h"
#include "TLorentzVector.h"

#include "styles.hh"
#include "utilities.hh"
#include "baby_basic.hh"
#include "xcone_utils.hh"

#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace fastjet;
namespace {
  int nEvents=1;
}

int main(){ 
  time_t begtime, endtime;
  time(&begtime);
  styles style("RA4"); style.setDefaultStyle();
  TCanvas can;

  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  TString foldert2(bfolder+"/cms2r0/aovcharova/treemaker/babies/");

  TString masses;
  vector<TString> ntuples;
  vector<TString> mlsp = {"950"};//, "200"};

  cout<<setprecision(5);
  for(size_t ilsp=0; ilsp<mlsp.size(); ilsp++){
      masses = "*_"+mlsp[ilsp]+"_*";
      ntuples = dirlist(foldert2, masses);
      cout<<"Number of ntuples found: "<<ntuples.size()<<endl;
      if(ntuples.size()!=1) continue;

      baby_basic baby(foldert2+ntuples[0]);
      cout<<"Doing "<<foldert2+ntuples[0]<<" with "<<baby.GetEntries()<<" entries"<<endl;
      long nent_max = baby.GetEntries();
      if (nEvents > 0) nent_max = nEvents;
      for(long entry(0); entry<nent_max; entry++){
        cout<<"-----------> Event: "<<baby.event()<<endl;
        baby.GetEntry(entry);

        vector<PseudoJet> pfcands(0);
        for(size_t ipfc = 0; ipfc<baby.pfcands_pt().size(); ++ipfc){

          if(is_nan(baby.pfcands_pt().at(ipfc)) || is_nan(baby.pfcands_eta().at(ipfc))
             || is_nan(baby.pfcands_phi().at(ipfc)) || is_nan(baby.pfcands_m().at(ipfc))) continue;

          if (baby.pfcands_frompv().at(ipfc)==0) continue;
          TLorentzVector this_lv;
          this_lv.SetPtEtaPhiM(baby.pfcands_pt().at(ipfc),baby.pfcands_eta().at(ipfc),baby.pfcands_phi().at(ipfc),baby.pfcands_m().at(ipfc));
          const PseudoJet this_pj(this_lv.Px(), this_lv.Py(), this_lv.Pz(), this_lv.E());

          pfcands.push_back(this_pj);
        }

        cout<<"==================   Raw anti-kT jets with pT > 25 GeV  ===================="<<endl;
        JetDefinition jet_def(antikt_algorithm, 0.4);
        ClusterSequence cs(pfcands, jet_def);
        vector<PseudoJet> diyjets = sorted_by_pt(cs.inclusive_jets(25.));

        for(size_t ijet = 0; ijet < diyjets.size() ; ++ijet){
          const PseudoJet &pj = diyjets.at(ijet);
          cout<<"DIY  jet "<<ijet<<": Pt  "<<setw(10)<<pj.pt()
                            <<" Eta "<<setw(10)<<pj.eta()
                            <<" Phi "<<setw(10)<<pj.phi_std()
                            <<" M "<<setw(10)<<pj.m()<<endl;
        }
        cout<<"==================   XCone example output   ===================="<<endl;

        xcone_utils::nsubj(pfcands);
        xcone_utils::xjets(pfcands);

      } // Loop over baby entries
    
  } // Loop over mlsp
  
  time(&endtime); 
  cout<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;  

}
