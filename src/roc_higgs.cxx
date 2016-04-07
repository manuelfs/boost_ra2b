// roc_higgs: Macro that plots ROC curves

#include <stdexcept>
#include <iostream>

#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TDirectory.h"
#include "TMarker.h"
#include "TStyle.h"

#include "styles.hh"
#include "utilities.hh"
#include "roc_higgs.hh"
#include "alias_ra2b.hh"

using namespace std;

int main(){
  time_t begtime, endtime;
  time(&begtime);

  styles style("RA4"); style.setDefaultStyle();
  gStyle->SetPadTickX(1); // Tickmarks on the top
  gStyle->SetPadTickY(1); // Tickmarks on the right

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

  vector<TString> s_bkg;
  s_bkg.push_back(folder+"*_TTJets*Lept*");
  s_bkg.push_back(folder+"*_WJetsToLNu*");
  s_bkg.push_back(folder+"*ZJetsToNuNu_HT*");
  s_bkg.push_back(folder+"*_QCD_HT*");
  s_bkg.push_back(folder+"*DYJetsToLL_M-50_HT*");
  s_bkg.push_back(folder+"*_GJets_HT*");
  s_bkg.push_back(folder+"*_ST_*");
  s_bkg.push_back(folder+"*_TTGJets*");
  s_bkg.push_back(folder+"*_TTZTo*");


  ////////////////////////// SAMPLES for the axes /////////////////////////
  vector<sample_class> tt_t5hnc; 
  tt_t5hnc.push_back(sample_class("T5HH(1200,200)", s_t5hnc));
  tt_t5hnc.push_back(sample_class("t#bar{t}", s_tt));

  vector<sample_class> tt_t5hc; 
  tt_t5hc.push_back(sample_class("T5HH(1200,950)", s_t5hc));
  tt_t5hc.push_back(sample_class("t#bar{t}", s_tt));

  vector<sample_class> bkg_t5hc; 
  bkg_t5hc.push_back(sample_class("T5HH(1200,950)", s_t5hc));
  bkg_t5hc.push_back(sample_class("All bkg.", s_bkg));



  ///////////////////// VARIABLES for each ROC /////////////////////
  int ht_col(kGreen+1), hig_col(2);
  int mj_style(8); float mj_size(2);
  vector<marker_class> nb_points, ht_points, mht_points, nj_points, hig_points;
  ht_points.push_back(marker_class(1000,  mj_size, ht_col, mj_style));
  ht_points.push_back(marker_class(1500,  4, ht_col, 29));
  ht_points.push_back(marker_class(2000,  mj_size, ht_col, mj_style));
  ht_points.push_back(marker_class(2500,  mj_size, ht_col, mj_style));

  mht_points.push_back(marker_class(400,  mj_size, 1, mj_style));
  mht_points.push_back(marker_class(600,  4, 1, 29));
  mht_points.push_back(marker_class(800,  mj_size, 1, mj_style));
  mht_points.push_back(marker_class(1000,  mj_size, 1, mj_style));

  nj_points.push_back(marker_class(6,  mj_size, 28, mj_style));
  nj_points.push_back(marker_class(8,  4, 28, 29));
  nj_points.push_back(marker_class(10,  mj_size, 28, mj_style));
  nj_points.push_back(marker_class(12,  mj_size, 28, mj_style));

  nb_points.push_back(marker_class(1,  mj_size, 4, mj_style));
  nb_points.push_back(marker_class(2,  4, 4, 29));
  nb_points.push_back(marker_class(3,  mj_size, 4, mj_style));
  nb_points.push_back(marker_class(4,  mj_size, 4, mj_style));

  hig_points.push_back(marker_class(1,  4, hig_col, 29));
  hig_points.push_back(marker_class(2,  mj_size, hig_col, mj_style));
  hig_points.push_back(marker_class(3,  mj_size, hig_col, mj_style));

  vector<marker_class> hig_green, hig_black, hig_brown, hig_blue;
  int c_green(kGreen+1), c_black(1), c_brown(28), c_blue(4);
  hig_green.push_back(marker_class(1,  4, c_green, 29));
  hig_green.push_back(marker_class(2,  mj_size, c_green, mj_style));
  hig_green.push_back(marker_class(3,  mj_size, c_green, mj_style));

  hig_black.push_back(marker_class(1,  4, c_black, 29));
  hig_black.push_back(marker_class(2,  mj_size, c_black, mj_style));
  hig_black.push_back(marker_class(3,  mj_size, c_black, mj_style));

  hig_brown.push_back(marker_class(1,  4, c_brown, 29));
  hig_brown.push_back(marker_class(2,  mj_size, c_brown, mj_style));
  hig_brown.push_back(marker_class(3,  mj_size, c_brown, mj_style));

  hig_blue.push_back(marker_class(1,  4, c_blue, 29));
  hig_blue.push_back(marker_class(2,  mj_size, c_blue, mj_style));
  hig_blue.push_back(marker_class(3,  mj_size, c_blue, mj_style));

  vector<var_class> general;
  general.push_back(var_class("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.6&&fjets_csv1>.605&&fjets_csv2>.605)",3,0,"N_{H-jet} (#tau_{2}/#tau_{1}<0.6, n_{b}^{L} #geq 2)",hig_col,1,hig_points));
  general.push_back(var_class("ht",4000,0,"H_{T}",ht_col,1,ht_points));
  general.push_back(var_class("mht",1500,0,"H_{T}^{miss}",1,1,mht_points));
  general.push_back(var_class("njets",15,0,"N_{jet}",28,1,nj_points));
  general.push_back(var_class("nbm",7,0,"N_{b-jet}",4,1,nb_points));

  TString nbm("N#lower[-.01]{_{b}}#kern[-0.6]{#lower[.1]{^{M}}}"), nbl("N#lower[-.01]{_{b}}#kern[-0.6]{#lower[.1]{^{L}}}");
  vector<var_class> htags;
  htags.push_back(var_class("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.4&&fjets_csv1>.605)",
			    3,0,"N_{H-jet} (#tau_{2}/#tau_{1}<0.4, "+nbl+" #geq 1)",4,1,hig_blue));
  htags.push_back(var_class("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.4&&fjets_csv1>.89)",3,0,
			    "N_{H-jet} (#tau_{2}/#tau_{1}<0.4, "+nbm+" #geq 1)",ht_col,1,hig_green));
  htags.push_back(var_class("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.6&&fjets_csv1>.89)",
			    3,0,"N_{H-jet} (#tau_{2}/#tau_{1}<0.6, "+nbm+" #geq 1)",1,1,hig_black));
  htags.push_back(var_class("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.3&&fjets_csv1>.89)",
			    3,0,"N_{H-jet} (#tau_{2}/#tau_{1}<0.3, "+nbm+" #geq 1)",28,1,hig_brown));
  htags.push_back(var_class("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.6&&fjets_csv1>.89&&fjets_csv2>.89)",
			    3,0,"N_{H-jet} (#tau_{2}/#tau_{1}<0.6, "+nbm+" #geq 2)",hig_col,1,hig_points));

  vector<var_class> htags2;
  htags2.push_back(var_class("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.4&&fjets_csv1>.89)",
			     3,0,"N_{H-jet} (#tau_{2}/#tau_{1}<0.4,"+nbm+"  #geq 1)",ht_col,1,hig_green));
  htags2.push_back(var_class("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.4&&fjets_csv1>.605&&fjets_csv2>.605)",
			     3,0,"N_{H-jet} (#tau_{2}/#tau_{1}<0.4, "+nbl+" #geq 2)",1,1,hig_black));
  htags2.push_back(var_class("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.3&&fjets_csv1>.605&&fjets_csv2>.605)",
			     3,0,"N_{H-jet} (#tau_{2}/#tau_{1}<0.3, "+nbl+" #geq 2)",28,1,hig_brown));
  htags2.push_back(var_class("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.6&&fjets_csv1>.605&&fjets_csv2>.605)",
			     3,0,"N_{H-jet} (#tau_{2}/#tau_{1}<0.6, "+nbl+" #geq 2)",4,1,hig_blue));
  htags2.push_back(var_class("Sum$(fjets_pm>90&&fjets_pm<140&&fjets_tau21<0.6&&fjets_csv1>.89&&fjets_csv2>.89)",
			     3,0,"N_{H-jet} (#tau_{2}/#tau_{1}<0.6, "+nbm+" #geq 2)",hig_col,1,hig_points));



  ///////////////////// ROCs to plot /////////////////////

  //DrawROC(bkg_t5hc, general, "nleps==0&&ht>750&&mht>400&&njets>=6&&dphi1>0.5&&ntrks==0", "general");
  DrawROC(bkg_t5hc, htags2, "nleps==0&&ht>750&&mht>400&&njets>=6&&dphi1>0.5&&ntrks==0", "htags2");


  // DrawROC(tt_t5hnc, general, "nleps==0&&ht>500&&mht>200&&njets>=4&&dphi1>0.5&&ntrks==0", "NC");
  // DrawROC(tt_t5hc, general, "nleps==0&&ht>500&&mht>200&&njets>=4&&dphi1>0.5&&ntrks==0", "C");


  time(&endtime); 
  cout<<endl<<"ROCs took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}



void DrawROC(vector<sample_class> samples, vector<var_class> vars, TString cuts, TString tag){
  TCanvas can;
  const int nbins(1000);
  vector<vector<TH1D> > histos;
  TString hname, totcut;
  TChain *chain[2];

  for(unsigned sam(0); sam < samples.size(); sam++){
    // Loading chains
    TString treeName = getTreeName(samples[sam].files[0]);
    chain[sam] = new TChain(treeName);
    for(unsigned isam(0); isam < samples[sam].files.size(); isam++){
      int nfiles = chain[sam]->Add(samples[sam].files[isam]);
      if(nfiles==0) cout<<samples[sam].files[isam]<<" not found"<<endl;
    }
    setAliasRa2b(chain[sam]);
    histos.push_back(vector<TH1D>());

    // Projecting variables
    for(unsigned var(0); var<vars.size(); var++){
      float minh(vars[var].minx), maxh(vars[var].maxx);
      if(minh > maxh){
	minh = maxh;
	maxh = vars[var].minx;
      }
      hname = "histo"; hname += sam; hname += var;
      totcut = "weight*("+cuts+"&&"+samples[sam].cut+")";
      histos[sam].push_back(TH1D(hname,"",nbins,minh,maxh));
      chain[sam]->Project(hname, vars[var].varname, totcut);
    } // Loop over variables
  } // Loop over samples

  TH1D base_histo("base",cuts2title(cuts),1,0.005,1.0);
  base_histo.SetXTitle(samples[0].label+" efficiency");
  base_histo.SetYTitle(samples[1].label+" efficiency");
  base_histo.SetMinimum(0.0);
  base_histo.SetMaximum(1.0);
  base_histo.SetDirectory(0);
  base_histo.Draw();

  // Legend
  double legX = 0.2, legY = 0.9, legSingle = 0.067;
  double legW = 0.2, legH = legSingle*vars.size();
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(0.058); leg.SetFillColor(0); leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(132);

  // Making individual graphs
  TGraph graphs[100]; // Had to make it an array because a vector<TGraph> kept crashing
  for(unsigned var(0); var<vars.size(); var++){
    graphs[var] = MakeROC(histos[0][var], histos[1][var], vars[var].minx < vars[var].maxx, vars[var].cuts);
    graphs[var].SetLineColor(vars[var].color);
    graphs[var].SetLineStyle(vars[var].style);
    graphs[var].SetLineWidth(4);
    leg.AddEntry(&(graphs[var]), vars[var].title, "l");
    graphs[var].Draw("lsame");
  } // Loop over variables
  leg.Draw();

  TString pname("plots/roc_"+tag+"_"+format_tag(cuts)+".pdf");  
  can.Print(pname);
  can.SetLogx(1);
  can.SetLogy(1);
  pname.ReplaceAll("roc_","log_roc_");
  base_histo.SetMinimum(1e-4);
  can.Print(pname);

  for(unsigned sam(0); sam < samples.size(); sam++)
    chain[sam]->Delete();
}

TGraph MakeROC(TH1D &good, TH1D &bad, const bool less_is_better, vector<marker_class> cuts){
  const int nbins = good.GetNbinsX();
  if(bad.GetNbinsX() != nbins) throw logic_error("Mismatched number of bins");

  TMarker marker;  

  TGraph graph(0);
  const double good_tot = good.Integral(0, nbins+1);
  const double bad_tot = bad.Integral(0, nbins+1);
  int inibin(0), endbin(nbins+1), dbin(1); unsigned icut(0);
  if(less_is_better){
    inibin = nbins+1;
    endbin = 0;
    dbin = -1;
  }
  for(int bin = inibin; bin*dbin<=endbin*dbin; bin+=dbin){
    const double good_pass = good.Integral(min(endbin,bin), max(endbin,bin));
    const double bad_pass = bad.Integral(min(endbin,bin), max(endbin,bin));
    const double x = good_pass/good_tot;
    const double y = bad_pass/bad_tot;
    graph.SetPoint(graph.GetN(), x, y);
 
    // Plotting the stars
    if(icut<cuts.size()){
      float edge(good.GetXaxis()->GetBinUpEdge(bin));
      if((edge>=cuts[icut].cut&&!less_is_better) || (edge<=cuts[icut].cut&&less_is_better)){
	marker.SetMarkerStyle(cuts[icut].style);marker.SetMarkerColor(cuts[icut].color);
	marker.SetMarkerSize(cuts[icut].size); 
	if(x>0.01 && y>0.0001) marker.DrawMarker(x,y);
	icut++;
      }
    }
  }
  TString name(good.GetName());
  name += "graph";
  graph.SetName(name);
  graph.SetTitle(name);

  return graph;
}

var_class::var_class(TString ivarname, float iminx, float imaxx, TString ititle, int icolor, 
	    int istyle, vector<marker_class> icuts){
  varname = ivarname; minx = iminx; maxx = imaxx; title = ititle;
  cuts = icuts; 
  color = icolor; style = istyle;
}

sample_class::sample_class(TString ilabel, vector<TString> ifiles, TString icut){
  files = ifiles; label = ilabel; cut = icut;
}

marker_class::marker_class(float icut, float isize, int icolor, int istyle){
  cut=icut; size=isize; color=icolor; style=istyle;
}

