// plot_t5tttt: Compares kinematic distributions of T1tttt and T5tttt

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
  unsigned nEventsMax=10;
  bool verb = true;
}

int main(){ 
  time_t begtime, endtime;
  time(&begtime);
  styles style("RA4"); style.setDefaultStyle();

  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  TString foldert2(bfolder+"/cms2r0/boost/pfcands/");

  TString masses;
  vector<TString> ntuples;
  vector<TString> mlsp = {"950"};//, "200"};

  cout<<setprecision(5);
  for(size_t ilsp=0; ilsp<mlsp.size(); ilsp++){
      TCanvas can("can","can",500, 500);

      masses = "*_"+mlsp[ilsp]+"_*";
      ntuples = dirlist(foldert2, masses); 
      if (verb) cout<<"Number of ntuples found: "<<ntuples.size()<<endl;
      if(ntuples.size()!=1) continue;

      baby_boost baby(foldert2+ntuples[0]);
      if (verb) cout<<"Doing "<<foldert2+ntuples[0]<<" with "<<baby.GetEntries()<<" entries"<<endl;
      unsigned nEvents = 0;
      long nent_max = baby.GetEntries();
      for(long entry(0); entry<nent_max; entry++){
        baby.GetEntry(entry);
        if (verb) cout<<"-----------> Event: "<<entry<<endl;

        //require at least one boosted higgs
        bool boosted05 = false;
        int nb_fromH =0;
        if (baby.isr_tru_pt()>10) continue;
        for (unsigned imc=0; imc<baby.mc_pt().size(); imc++){
          if (abs(baby.mc_id()[imc])==5 && baby.mc_mom()[imc]==25) nb_fromH++;
          if (baby.mc_id()[imc]==25 && baby.mc_pt()[imc]>400.){
            boosted05 = true;
          }
        }
        cout<<nb_fromH<<endl;
        if (!boosted05 || nb_fromH!=4) continue;
        cout<<"Found event"<<endl;
        for (unsigned imc=0; imc<baby.mc_pt().size(); imc++){
          if (baby.mc_id()[imc]==25 || (abs(baby.mc_id()[imc])==5 && baby.mc_mom()[imc]==25))
            cout<<"Pdg id = "<<setw(10)<<baby.mc_id()[imc]
                <<"   pt = "<<setw(10)<<baby.mc_pt()[imc]
                <<"   eta = "<<setw(10)<<baby.mc_eta()[imc]
                <<"   phi = "<<setw(10)<<baby.mc_phi()[imc]<<endl;
        }

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

        // if (verb) cout<<"==================   Raw anti-kT jets with pT > 25 GeV  ===================="<<endl;
        JetDefinition jet_def(antikt_algorithm, 0.4);
        ClusterSequence cs(pfcands, jet_def);
        vector<PseudoJet> diyjets = sorted_by_pt(cs.inclusive_jets(30.));

        for(size_t ijet = 0; ijet < diyjets.size() ; ++ijet){
          const PseudoJet &pj = diyjets.at(ijet);
          if (verb) cout<<"DIY  jet "<<ijet<<": Pt  "<<setw(10)<<pj.pt()
                            <<" Eta "<<setw(10)<<pj.eta()
                            <<" Phi "<<setw(10)<<pj.phi_std()
                            <<" M "<<setw(10)<<pj.m()<<endl;
        }

        // only use events where the njets we found is the same as the njets from DIY
        if (unsigned(baby.njets())!=diyjets.size() && baby.nleps()==0) continue;



        unsigned njets = diyjets.size();
        cout<<"Events has "<<njets<<" akt4 jets"<<endl;

        TH1D htau(TString::Format("tau_%i",int(entry)),TString::Format("N_{jets}(anti-k_{T}) = %i",njets),15, 0.5, 15.5);
        htau.GetXaxis()->SetTitle("Number of requested jets");
        htau.GetYaxis()->SetTitle("N-jettiness");
        TH1D hdtau(TString::Format("dtau_%i",int(entry)),TString::Format("N_{jets}(anti-k_{T}) = %i",njets),15, 0.5, 15.5);
        hdtau.GetXaxis()->SetTitle("Number of requested jets");
        hdtau.GetYaxis()->SetTitle("N-jettiness - (N-1)-jettiness");
        TH1D hrtau(TString::Format("rtau_%i",int(entry)),TString::Format("N_{jets}(anti-k_{T}) = %i",njets),15, 0.5, 15.5);
        hrtau.GetXaxis()->SetTitle("Number of requested jets");
        hrtau.GetYaxis()->SetTitle("N-jettiness / (N-1)-jettiness");
        double last_tau = 0;
        for (unsigned inj=1; inj<16; inj++){
          // if (verb) cout<<"==================  Running XCone N jet = "<<inj<<"   ===================="<<endl;
          XConePlugin xcone_plugin(inj, 0.4, 2.0);
          JetDefinition xcone_jetDef(&xcone_plugin);
          ClusterSequence xcone_seq(pfcands, xcone_jetDef);
          vector<PseudoJet> xcjets = xcone_seq.inclusive_jets();
          if (xcjets.size()==0) continue;
          const NjettinessExtras * extras = njettiness_extras(xcjets[0]);
          // double tau_beam = extras->beamTau();
          double tau_tot = extras->totalTau();
          htau.SetBinContent(inj,tau_tot);
          if (inj>1){
            hdtau.SetBinContent(inj,tau_tot-last_tau);
            hrtau.SetBinContent(inj,tau_tot/last_tau);
          }
          last_tau = tau_tot;

          if (inj==8) xcone_utils::PrintXConeJets(xcjets);
          // double tau_cl = 0;
          // unsigned nxcjets_pt30 = 0;
          // for (unsigned ixcjet=0; ixcjet<xcjets.size(); ixcjet++){
          //   // if (xcjets[ixcjet].perp() < 30) continue;
          //   tau_cl += max(extras->subTau(xcjets[ixcjet]),0.0);
          //   nxcjets_pt30++;
          // }
          // cout<<njets<<" "<<nxcjets_pt30<<endl;
          // if (nxcjets_pt30>njets) {
          //   xcone_utils::PrintXConeJets(xcjets);
          //   if (verb) cout<<"Tau clustered = "<<setw(10)<<tau_cl
          //           <<"    Tau beam = "<<setw(10)<<tau_beam
          //           <<"    Tau tot = "<<setw(10)<<tau_tot<<endl;
          // } else continue;
        }
        htau.Draw();
        TString ttl = htau.GetName();
        ttl+=".pdf";
        can.Print(ttl);
        hdtau.Draw();
        ttl = hdtau.GetName();
        ttl+=".pdf";
        can.Print(ttl);
        hrtau.Draw();
        ttl = hrtau.GetName();
        ttl+=".pdf";
        can.Print(ttl);
        nEvents++;
        if (nEvents>=nEventsMax) break;
      } // Loop over baby entries
    
  } // Loop over mlsp
  
  time(&endtime); 
  cout<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;  

}
