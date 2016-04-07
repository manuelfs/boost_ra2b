// plot_cuts: Macro that plots significance and other S,B quantities

#include <iostream>
#include <vector>

#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TColor.h"
#include "RooStats/NumberCountingUtils.h"

#include "styles.hh"
#include "timer.hh"
#include "alias_ra2b.hh"
#include "utilities.hh"

namespace  {
  TString luminosity="10";
  TString plot_type=".pdf";
  TString plot_style="RA4";
}

using namespace std;
using std::cout;
using std::endl;

int main(){ 
  styles style(plot_style); style.setDefaultStyle();
  vector<hfeats> vars;
  TCanvas can;

  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  
  TString folder(bfolder+"/cms2r0/treemaker/skim/");


  vector<TString> s_t1tnc;
  s_t1tnc.push_back(folder+"*T1tttt*1500_*100*");
  vector<TString> s_t1tc;
  s_t1tc.push_back(folder+"*T1tttt*1200_*800*");
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


  vector<sfeats> Samples; 
  Samples.push_back(sfeats(s_tt, "t#bar{t}",		 ra2b::c_tt_1l));
  Samples.push_back(sfeats(s_znunu, "Z#rightarrow#nu#nu",ra2b::c_znunu));
  Samples.push_back(sfeats(s_wjets, "W+jets",		 ra2b::c_wjets));
  Samples.push_back(sfeats(s_qcd, "QCD",		 ra2b::c_qcd));
  Samples.push_back(sfeats(s_other, "Other",		 ra2b::c_other)); 
  Samples.push_back(sfeats(s_t5hnc, "T5HH(1200,200)",	 ra2b::c_t5hh, 1)); Samples.back().isSig = true;
  Samples.push_back(sfeats(s_t5hc,  "T5HH(1200,950)",	 ra2b::c_t5hh, 2)); Samples.back().isSig = true;
  // Samples.push_back(sfeats(s_t1tnc, "T1tttt(1500,100)",	 ra2b::c_t1t, 1)); Samples.back().isSig = true;
  // Samples.push_back(sfeats(s_t1tc, "T1tttt(1200,800)",	 ra2b::c_t1t, 2)); Samples.back().isSig = true;


  // Reading ntuples
  vector<TChain *> chain;
  for(unsigned sam(0); sam < Samples.size(); sam++){
    TString treeName = getTreeName(Samples[sam].file[0]);
    chain.push_back(new TChain(treeName));
    for(unsigned insam(0); insam < Samples[sam].file.size(); insam++){
      // cout<<"Reading "<<Samples[sam].file[insam]<<endl;
      chain[sam]->Add(Samples[sam].file[insam]);
      // cout<<"Entries "<<chain[sam]->GetEntries()<<endl;
    }
    setAliasRa2b(chain[sam]);
  }

  vector<int> ra2b_sam;
  unsigned nsam(Samples.size());
  for(unsigned sam(0); sam < nsam; sam++){
    ra2b_sam.push_back(sam);
  } // Loop over samples

  //const int scanbins(100);
  // vars.push_back(hfeats("ht",scanbins,500,3000, ra2b_sam, "Cut on H_{T} [GeV]",
  //                       "nleps==0&&mht>200&&njets>=6&&dphi1>0.5&&ntrks==0&&nbm>=1"));
  // vars.push_back(hfeats("mht",scanbins,200,1200, ra2b_sam, "Cut on H_{T}^{miss} [GeV]",
  // 			"nleps==0&&ht>750&&njets>=6&&dphi1>0.5&&ntrks==0&&nbm>=1",500));
  // vars.push_back(hfeats("njets",14,3.5,17.5, ra2b_sam, "Cut on N_{jet}",
  // 			"nleps==0&&ht>750&&mht>200&&nbm>=1&&dphi1>0.5&&ntrks==0",6));

  vars.push_back(hfeats("nbm",7,-0.5,6.5, ra2b_sam, "Cut on N_{b-jet}",
  			"nleps==0&&ht>750&&mht>200&&njets>=6&&dphi1>0.5&&ntrks==0",2));
  vars.push_back(hfeats("nbm",7,-0.5,6.5, ra2b_sam, "Cut on N_{b-jet}",
  			"nleps==0&&ht>750&&mht>200&&njets>=6&&dphi1>0.5&&ntrks==0&&nhiggs==0",2));
  vars.push_back(hfeats("nbm",7,-0.5,6.5, ra2b_sam, "Cut on N_{b-jet}",
  			"nleps==0&&ht>750&&mht>200&&njets>=6&&dphi1>0.5&&ntrks==0&&nhiggs==1",2));
  vars.push_back(hfeats("nbm",7,-0.5,6.5, ra2b_sam, "Cut on N_{b-jet}",
  			"nleps==0&&ht>750&&mht>200&&njets>=6&&dphi1>0.5&&ntrks==0&&nhiggs>=2",2));




  TString plot_tag("_lumi"+luminosity+plot_type);
  double Syserr(pow(0.3,2));
  double legX = 0.66, legY = 0.9, legSingle = 0.061;
  double legW = 0.12, legH = legSingle*3;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(0.052); leg.SetFillColor(0); leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(132);

  TLine line; line.SetLineColor(28); line.SetLineWidth(4); line.SetLineStyle(2);
  const unsigned Nhis(4);
  vector< vector<TH1D*> > histo[Nhis];
  vector<TH1D*> varhisto;
  vector<float> nentries;
  TString hname, pname, variable, leghisto, totCut;
  for(unsigned var(0); var<vars.size(); var++){
    cout<<endl;
    // Generating vector of histograms
    for(unsigned his(0); his < Nhis; his++){
      varhisto.resize(0);
      for(unsigned sam(0); sam < vars[var].samples.size(); sam++){
	int isam = vars[var].samples[sam];
	hname = "histo"; hname += var; hname += his; hname += sam;
	varhisto.push_back(new TH1D(hname, cuts2title(vars[var].cuts), vars[var].nbins, 
				    vars[var].minx, vars[var].maxx));
	varhisto[sam]->SetXTitle(vars[var].title);
	varhisto[sam]->SetLineColor(Samples[isam].color);
	varhisto[sam]->SetLineStyle(Samples[isam].style);
	varhisto[sam]->SetLineWidth(3);
      }
      histo[his].push_back(varhisto);
    }

    //// Plotting Zbi in histo[0], S/B in histo[1], B in histo[2], S in hisot[3] ///
    leg.Clear();
    nentries.resize(0);
    variable = vars[var].varname;
    int bkgind(-1);
    for(unsigned sam(0); sam < vars[var].samples.size(); sam++){
      int isam = vars[var].samples[sam];
      bool isSig = Samples[isam].isSig;
      totCut = Samples[isam].factor+"*"+luminosity+"*weight*("+vars[var].cuts+"&&"+Samples[isam].cut+")"; 
      //cout<<totCut<<endl;
      histo[0][var][sam]->Sumw2();
      chain[isam]->Project(histo[0][var][sam]->GetName(), variable, totCut);
      histo[0][var][sam]->SetBinContent(vars[var].nbins,
					  histo[0][var][sam]->GetBinContent(vars[var].nbins)+
					  histo[0][var][sam]->GetBinContent(vars[var].nbins+1));
      histo[0][var][sam]->SetBinError(vars[var].nbins,
				      sqrt(pow(histo[0][var][sam]->GetBinError(vars[var].nbins),2)+
					   pow(histo[0][var][sam]->GetBinError(vars[var].nbins+1),2)));
      nentries.push_back(histo[0][var][sam]->Integral(1,vars[var].nbins));
      histo[0][var][sam]->SetYTitle("Z_{bi} (#sigma)");
      histo[1][var][sam]->SetYTitle("S/#sqrt{B} (#sigma)");
      histo[2][var][sam]->SetYTitle("Events");
      histo[2][var][sam]->SetLineColor(1);
      if(!isSig){ // Adding previous bkg histos
	for(int bsam(sam-1); bsam >= 0; bsam--){
	  histo[0][var][sam]->Add(histo[0][var][bsam]);
	  break;
	}
	bkgind = sam;
      } 
    } // First loop over samples

    for(int sam(vars[var].samples.size()-1); sam >= 0; sam--){
      int isam = vars[var].samples[sam];
      if(!(Samples[isam].isSig)) continue;

      double maxhisto[Nhis];
      for(unsigned his(0); his<Nhis; his++) maxhisto[his] = -1;
      for(int bin(1); bin<=vars[var].nbins; bin++){
	double Nerr(0);
	const double Nsig(histo[0][var][sam]->Integral(bin,vars[var].nbins+1));
	const double Nbkg(histo[0][var][bkgind]->IntegralAndError(bin,vars[var].nbins+1,Nerr)); 
	histo[2][var][sam]->SetBinContent(bin,Nbkg);
	histo[3][var][sam]->SetBinContent(bin,Nsig);
	if(Nbkg==0){
	  histo[0][var][sam]->SetBinContent(bin,0);
	  histo[1][var][sam]->SetBinContent(bin,0);
	} else {
	  const double Zbi(RooStats::NumberCountingUtils::BinomialExpZ(Nsig, Nbkg, sqrt(pow(Nerr/Nbkg,2)+Syserr)));
	  histo[0][var][sam]->SetBinContent(bin,Zbi>0?Zbi:0);
	  histo[1][var][sam]->SetBinContent(bin,Nsig/sqrt(Nbkg));
	}
	for(unsigned his(0); his<Nhis; his++) 
	  if(maxhisto[his] < histo[his][var][sam]->GetBinContent(bin)) 
	    maxhisto[his] = histo[his][var][sam]->GetBinContent(bin);
      }
    } // Loop over samples

    unsigned sam_nc(vars[var].samples.size()-2);
    unsigned isam_nc = vars[var].samples[sam_nc];
    unsigned sam_c(vars[var].samples.size()-1);
    unsigned isam_c = vars[var].samples[sam_c];

    float maxleg = 1.35;
    leg.Clear();
    legH = legSingle*3;
    leg.SetY1NDC(legY-legH);
    leg.SetHeader("#font[22]{   L = "+luminosity+" fb^{-1}}");
    leg.AddEntry(histo[0][var][sam_nc], Samples[isam_nc].label);
    leg.AddEntry(histo[0][var][sam_c], Samples[isam_c].label);
    float hismax = max(histo[0][var][sam_nc]->GetMaximum(), histo[0][var][sam_c]->GetMaximum())*maxleg;
    hismax = 8.5;
    histo[0][var][sam_nc]->SetMaximum(hismax);
    //histo[0][var][sam_nc]->SetMaximum(10);
    histo[0][var][sam_nc]->SetMinimum(0);
    histo[0][var][sam_nc]->Draw("l hist");
    histo[0][var][sam_c]->Draw("l hist same");
    style.moveYAxisLabel(histo[0][var][sam_nc], hismax, false);
    if(vars[var].cut>0) line.DrawLine(vars[var].cut, 0, vars[var].cut, hismax);
    leg.Draw();
    pname = "plots/zbi_"+vars[var].tag+plot_tag;
    can.SaveAs(pname);
    histo[1][var][sam_nc]->SetMaximum(max(histo[1][var][sam_nc]->GetMaximum(),
      histo[1][var][sam_c]->GetMaximum())*maxleg);
    histo[1][var][sam_nc]->SetMinimum(0);
    histo[1][var][sam_nc]->Draw("l hist");
    histo[1][var][sam_c]->Draw("l hist same");
    if(vars[var].cut>0) line.DrawLine(vars[var].cut, 0, vars[var].cut, histo[1][var][sam_nc]->GetMaximum()*1.05);
    pname = "plots/s_sqrtb_"+vars[var].tag+plot_tag;
    leg.Draw();
    can.SaveAs(pname);
    hismax = max(histo[3][var][sam_nc]->GetMaximum(), histo[3][var][sam_c]->GetMaximum())*maxleg;
    histo[2][var][sam_nc]->SetMaximum(hismax);
    histo[2][var][sam_nc]->Draw("l hist");
    histo[3][var][sam_c]->Draw("l hist same");
    histo[3][var][sam_nc]->Draw("l hist same");
    if(vars[var].cut>0) line.DrawLine(vars[var].cut, 0, vars[var].cut, hismax);
    style.moveYAxisLabel(histo[0][var][sam_nc], hismax, false);
    leg.AddEntry(histo[2][var][sam_nc], "Total bkg.");
    legH = legSingle*4;
    leg.SetY1NDC(legY-legH);
    leg.Draw();
    pname = "plots/sb_"+vars[var].tag+plot_tag;
    can.SaveAs(pname);
    
  }// Loop over variables

  for(unsigned his(0); his < Nhis; his++){
    for(unsigned var(0); var<vars.size(); var++){
      for(unsigned sam(0); sam < vars[var].samples.size(); sam++)
	if(histo[his][var][sam]) histo[his][var][sam]->Delete();
    }
  }
}

