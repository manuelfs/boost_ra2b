#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>

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
#include "baby_boost.hh"
#include "xcone_utils.hh"

#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace fastjet;
namespace {
  bool verb = false;
}

int main(){ 
  time_t begtime, endtime;
  time(&begtime);

  styles style("RA4"); style.setDefaultStyle();

  // create baby
  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  TString foldert2(bfolder+"/cms2r0/boost/pfcands/");
  vector<TString> ntuples;
  ntuples = dirlist(foldert2, "*1500_*-100_*_74_V9.root"); 
  if (ntuples.size()==0) cout<<" ERROR: Could not find ntuples"<<endl;
  baby_boost baby(foldert2+ntuples[0]);
  if (verb) cout<<"Doing "<<foldert2+ntuples[0]<<" with "<<baby.GetEntries()<<" entries"<<endl;

  // histogram definitions
  TCanvas can("can","can",500, 500);
  TH1D h_nja("h_nja",";Number of objects;Events",20, 0., 20.);
  TH1D h_njx("h_njx",";Number of objects;Events",20, 0., 20.);
  TH1D h_npart("h_npart",";Number of objects;Events",20, 0., 20.);

  // loopin'
  long nent_max = baby.GetEntries();
  for(long entry(0); entry<nent_max; entry++){
    baby.GetEntry(entry);

    // only look at all-hadronic events with no ISR for now
    if (baby.ntruleps()>0 || baby.isr_tru_pt()>10) continue;

    //   Get the number of partons within acceptance
    //-------------------------------------------------
    unsigned npart(0);
    for (unsigned imc=0; imc<baby.mc_pt().size(); imc++){
      const unsigned pid(abs(baby.mc_id()[imc])), momid(abs(baby.mc_mom()[imc]));
      const bool goodpt(baby.mc_pt()[imc]>30.), goodeta(abs(baby.mc_eta()[imc])<5.), isme(baby.mc_status()[imc]==23);
      if (!goodpt || !goodeta || !isme) continue;
      if (pid==5 && momid==6) npart++;
      else if (pid<6 && momid==24) npart++;
    }

    //   Get constituents for jet building
    //---------------------------------------
    vector<PseudoJet> pfcands(0);
    for(size_t ipfc = 0; ipfc<baby.pfcands_pt().size(); ++ipfc){
      if(is_nan(baby.pfcands_pt().at(ipfc)) || is_nan(baby.pfcands_eta().at(ipfc))
         || is_nan(baby.pfcands_phi().at(ipfc)) || is_nan(baby.pfcands_m().at(ipfc))) continue;
      if (baby.pfcands_frompv().at(ipfc)==0) continue;
      TLorentzVector this_lv;
      this_lv.SetPtEtaPhiM(baby.pfcands_pt().at(ipfc),baby.pfcands_eta().at(ipfc),
                           baby.pfcands_phi().at(ipfc),baby.pfcands_m().at(ipfc));
      const PseudoJet this_pj(this_lv.Px(), this_lv.Py(), this_lv.Pz(), this_lv.E());
      pfcands.push_back(this_pj);
    }

    //   Get the number anti-kT jets with pT > 30
    //----------------------------------------------
    JetDefinition jet_def(antikt_algorithm, 0.4);
    ClusterSequence cs(pfcands, jet_def);
    vector<PseudoJet> aktjets = sorted_by_pt(cs.inclusive_jets(30.));
    int nja(aktjets.size());

    // only use events where the njets we found is the same as the njets from DIY
    if (baby.njets()!=nja) continue;

    //   Get the max number jets with pT > 30 that XCone can find
    //--------------------------------------------------------------
    int njx(nja), njx_temp(njx+1);
    while (njx_temp > njx) {
      XConePlugin xcone_plugin(njx_temp, 0.4, 2.0);
      JetDefinition xcone_jetDef(&xcone_plugin);
      ClusterSequence xcone_seq(pfcands, xcone_jetDef);
      vector<PseudoJet> xcjets = sorted_by_pt(xcone_seq.inclusive_jets(30.));
      if (xcjets.size()==unsigned(njx_temp)) {
        njx++; njx_temp++;
      } else {
        njx_temp = xcjets.size();
      }
    }

    // fill results
    h_njx.Fill(njx);
    h_nja.Fill(nja);
    h_npart.Fill(npart);

  } // Loop over baby entries
    
  // draw results  
  h_npart.Scale(1/h_npart.Integral()); h_npart.SetLineWidth(3); h_npart.SetLineColor(kAzure+1); h_npart.Draw();
  h_njx.Scale(1/h_njx.Integral()); h_njx.SetLineWidth(3);   h_njx.SetLineColor(kGreen+1);   h_njx.Draw("same");
  h_nja.Scale(1/h_nja.Integral()); h_nja.SetLineWidth(3);   h_nja.SetLineColor(kRed+1);     h_nja.Draw("same");
  can.Print("npartons.pdf");

  //print info
  cout<<setprecision(4);
  cout<<"Anti-kT jets - Mean: "<<setw(5)<<h_nja.GetMean()<<" RMS: "<<setw(5)<<h_nja.GetRMS()<<endl;
  cout<<"XCone jets   - Mean: "<<setw(5)<<h_njx.GetMean()<<" RMS: "<<setw(5)<<h_njx.GetRMS()<<endl;
  cout<<"N partons    - Mean: "<<setw(5)<<h_npart.GetMean()<<" RMS: "<<setw(5)<<h_npart.GetRMS()<<endl;

  time(&endtime); 
  cout<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;  
}
