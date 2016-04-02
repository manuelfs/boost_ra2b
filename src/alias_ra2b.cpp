//----------------------------------------------------------------------------
// alias_ra2b - Sets aliases for the #$@%-named branches in TreeMaker
//----------------------------------------------------------------------------

#ifndef ROOT_VERSION
#include "alias_ra2b.hpp"
#endif

void setAliasRa2b(TTree *tree){
  tree->SetAlias("ht"	 ,"HT");
  tree->SetAlias("mht"	 ,"MHT");
  tree->SetAlias("njets" ,"NJets");
  tree->SetAlias("nels"	 ,"@Electrons.size()");
  tree->SetAlias("nmus"	 ,"@Muons.size()");
  tree->SetAlias("nleps" ,"@Electrons.size()+@Muons.size()");
  tree->SetAlias("dphi1" ,"DeltaPhi1");
  tree->SetAlias("dphi2" ,"DeltaPhi2");
  tree->SetAlias("dphi3" ,"DeltaPhi3");
  tree->SetAlias("dphi4" ,"DeltaPhi4");
  tree->SetAlias("ntrks" ,"isoMuonTracks+isoElectronTracks+isoPionTracks");
  tree->SetAlias("mc"    ,"GenParticles");
  tree->SetAlias("mc_id" ,"GenParticles_PdgId");
  tree->SetAlias("mc_mom" ,"GenParticles_ParentId");
  tree->SetAlias("fjets_csv1" ,"JetsAK8_bDiscriminatorSubjet1CSV");
  tree->SetAlias("fjets_csv2" ,"JetsAK8_bDiscriminatorSubjet2CSV");
  tree->SetAlias("fjets_tau1" ,"JetsAK8_NsubjettinessTau1");
  tree->SetAlias("fjets_tau2" ,"JetsAK8_NsubjettinessTau2");
  tree->SetAlias("fjets_tau3" ,"JetsAK8_NsubjettinessTau3");
  tree->SetAlias("fjets_tau21" ,"JetsAK8_NsubjettinessTau2/JetsAK8_NsubjettinessTau1");
  tree->SetAlias("fjets_tau32" ,"JetsAK8_NsubjettinessTau3/JetsAK8_NsubjettinessTau2");
  tree->SetAlias("fjets_pm" ,"JetsAK8_prunedMass");
  tree->SetAlias("fjets" ,"JetsAK8");
}
