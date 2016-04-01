//----------------------------------------------------------------------------
// utilities - Various functions used accross the code
//----------------------------------------------------------------------------

#ifndef H_UTILITIES
#define H_UTILITIES

#include <cstddef>
#include <cstdio>
#include <cmath>

#include <string>
#include <vector>

#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TChain.h"

typedef std::pair<int,double> int_double;
typedef std::pair<double,double> double_double;
const long double PI = acos(-1.L);

namespace ra2b{
  // Had to define the TColor objects in the cpp
  enum {
    c_tt_1l    = 2012, // ucsb_blue
    c_tt_2l    = 2011, // tar_heel_blue
    c_wjets    = 2018, // ucsb_gold
    c_singlet  = 2015,
    c_qcd      = 2020,
    c_ttv      = 2002, // penn_red
    c_other    = 2019
  };

}


struct pfeats{
  pfeats(const std::vector<int> &isamples, const TString &icut = "1", const TString &itagname="");

  std::vector<int> samples;
  TString cut, tagname;
};

class hfeats {
public:
  hfeats(TString ivarname, int inbins, float iminx, float imaxx, std::vector<int> isamples,
         TString ititle="", TString icuts="1", float icut=-1, TString itagname="",bool iskiplog=false, 
         std::vector<double> inevents= std::vector<double>(1,-1.));
  hfeats(TString ivarname, int inbins, float* ibinning, std::vector<int> isamples,
         TString ititle="", TString icuts="1", float icut=-1, TString itagname="",bool iskiplog=false, 
         std::vector<double> inevents= std::vector<double>(1,-1.));
  hfeats(TString ivarnamex, TString ivarnamey, int inbinsx, float iminx, float imaxx, int inbinsy, float iminy, float imaxy,  std::vector<int> isamples,
         TString ititlex, TString ititley, TString icuts, float icutx, float icuty, TString itagname);
  TString title, titlex, titley, varname, varnamex, varnamey, tag, cuts, unit;
  int nbins, nbinsx, nbinsy;
  float *binning;
  float minx, maxx, miny, maxy, cut, cutx, cuty,  maxYaxis, maxRatio;
  std::vector<int> samples;
  TString tagname;
  void format_tag();
  std::vector<double> nevents; //Added for track veto study. Useful to display number of events when hist is filled N times per event
  bool skiplog;
  TString whichPlots; // String that determines which of the [log_]lumi and [log_]shapes plots to make
  bool normalize; //normalizes isData to sum of histograms
  bool PU_reweight; 
  float moveRLegend; // Moves the right legend by this amount
  
};

class sfeats {
public:
  sfeats(std::vector<TString> ifile, TString ilabel, int icolor=1, int istyle=1, TString icut="1",
	 TString samVariable="noPlot");
  std::vector<TString> file;
  TString label, cut, factor,tag;
  int color, style;
  bool isSig, doStack, isData, mcerr, doBand;
  TString samVariable; // Used to plot different variables in the same histogram
};

class sysfeats {
public:
  sysfeats(TString iname, TString ititle);
  TString name;
  TString title;
  std::vector<TString> bincuts;
  std::vector<double> weights;
  void push_back(TString bincut, double weight);
  TString bincut(unsigned i);
  double weight(unsigned i);
  unsigned size();
};


void calc_chi2_diff(TH1D *histo1, TH1D *histo2, float &chi2, int &ndof, float &pvalue, float *average);
void calc_chi2(TH1D *histo, float &chi2, int &ndof, float &pvalue, float &average);
long getYieldErr(TChain& tree, TString cut, double& yield, double& uncertainty);

void plot_distributions(std::vector<sfeats> Samples, std::vector<hfeats> vars, TString luminosity="10", 
			TString filetype=".eps", TString namestyle="LargeLabels", TString dir = "1d", 
			bool doRatio=false, bool showcuts=false);
void plot_2D_distributions(std::vector<sfeats> Samples, std::vector<hfeats> vars, TString luminosity,
                           TString filetype, TString namestyle, TString dir);
TString cuts2title(TString title);
TString invertcut(TString cut);
TString format_tag(TString tag);
double gsl_ran_gamma (const double a, const double b, TRandom3 &rand);
double intGaus(double mean, double sigma, double minX, double maxX);
// yields[Nobs][Nsam] has the entries for each sample for each observable going into kappa
// weights[Nobs][Nsam] has the average weight of each observable for each sample
// powers[Nobs] defines kappa = Product_obs{ Sum_sam{yields[sam][obs]*weights[sam][obs]}^powers[obs] }
double calcKappa(std::vector<std::vector<float> > &entries, std::vector<std::vector<float> > &weights,
		 std::vector<float> &powers, float &mSigma, float &pSigma, bool do_data=false, 
		 bool verbose=false, bool do_plot=false, int nrep=100000);
float Efficiency(float den, float num, float &errup, float &errdown);









float cross_section(const TString &file);
std::vector<TString> dirlist(const TString &folder,
                             const TString &inname="dir",
                             const TString &tag="");
bool eigen2x2(float matrix[2][2], float &eig1, float &eig2);
bool id_big2small(const int_double& left, const int_double& right);
bool dd_big2small(const double_double& left, const double_double& right);
bool dd_small2big(const double_double& left, const double_double& right);
long double DeltaPhi(long double phi1, long double phi2);
long double SignedDeltaPhi(long double phi1, long double phi2);
float dR(float eta1, float eta2, float phi1, float phi2);
TString RoundNumber(double num, int decimals, double denom=1.);
long double AddInQuadrature(long double x, long double y);
long double GetMass(long double e, long double px, long double py, long double pz);
long double GetMT(long double m1, long double pt1, long double phi1,
                  long double m2, long double pt2, long double phi2);
long double GetMT(long double pt1, long double phi1,
                  long double pt2, long double phi2);
bool Contains(const std::string& text, const std::string& pattern);

std::vector<std::string> Tokenize(const std::string& input,
                                  const std::string& tokens=" ");
void AddPoint(TGraph& graph, const double x, const double y);

template<class T>
bool is_nan(const T &x){return x!=x;}

template<class T>
short Sign(T val){
  return (T(0) < val) - (val < T(0));
}

std::string execute(const std::string &cmd);
std::string RemoveTrailingNewlines(std::string str);

std::vector<double> LinearSpacing(size_t npts, double low, double high);

#endif
