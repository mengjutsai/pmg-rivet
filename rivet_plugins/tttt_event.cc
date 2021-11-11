#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/AnalysisLoader.hh"

#include <numeric>
#include <functional>

/*
 * Author : Rohin Narayan (narayan@cern.ch)
 *
 *
 * This rivet can be compared to a  simple phasespace equivalent to the the region definitions in the ttH-ML analysis.
 * The histograms need to be normalized to appropriate cross-section and total event weights. In an ATLAS environment
 * these histograms gets converted to ROOT format and the normalizations and handled by a subsequent script outside rivet.
 *
 *
 */

namespace Rivet {
    bool debug = false;
    int countProngs(Particle mother) {
        int n_prongs = 0;
        for(Particle p : mother.children())
            if (p.charge3()!=0) ++n_prongs;
        return n_prongs;
    }



  class ttw_ttH: public Analysis {
  public:

    /// Minimal constructor
    ttw_ttH() : Analysis("ttw_ttH")
    {
    }

    /// Set up projections and book histograms
    void init() {
      FinalState lepfs;
      //Projection to find prompt electrons
      IdentifiedFinalState el_id(lepfs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      electrons.acceptTauDecays(true);
      declare(electrons,"electrons");

      declare(UnstableParticles(),"UFS");

      //Projection to find prompt muons
      IdentifiedFinalState mu_id(lepfs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(true);
      declare(muons,"muons");


      TauFinder tauhadronic(TauFinder::DecayMode::HADRONIC);
      declare(tauhadronic,"TauHadronic");

      declare(HeavyHadrons(Cuts::abseta < 5 && Cuts::pT > 5*GeV), "BCHadrons");

      const FinalState fs(Cuts::abseta < 2.5);
      declare(fs, "FS");
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "Jets");
     // declare(FastJets(FinalState(Cuts::abseta < 2.5 && Cuts::pT>25*GeV),FastJets::ANTIKT,0.4),"Jets");
      //declare(FastJets(FinalState(), FastJets::ANTIKT, 0.4), "Jets");
      declare(MissingMomentum(FinalState(Cuts::abseta < 5 && Cuts::pT >0*GeV)),"MissingET");


      //Histogramming

      // Inclusive region
      book(_h["Inclusive_nJets"],"Inclusive_nJets",11,-0.5,10.5);
      book(_h["Inclusive_HT"],"Inclusive_HT",20,0,1200);
      book(_h["Inclusive_HT_jets"],"Inclusive_HT_jets",20,0,1200);
      book(_h["Inclusive_nBjets"],"Inclusive_nBjets",5,-0.5,4.5);
      book(_h["Inclusive_jet0Pt"],"Inclusive_jet0Pt",20,0,600);
      book(_h["Inclusive_jet0Eta"],"Inclusive_jet0Eta",10,-2.5,2.5);

      // Check lepton number
      book(_h["Inclusive_nLep_noHadTau"],"Inclusive_nLep_noHadTau",4,-0.5,3.5);
      book(_h["Inclusive_nLep_withHadTau"],"Inclusive_nLep_withHadTau",4,-0.5,3.5);
      // book(_h["LepSelection_nLep_noHadTau"],"LepSelection_nLep_noHadTau",4,-0.5,3.5);
      // book(_h["LepSelection_nLep_withHadTau"],"LepSelection_nLep_withHadTau",4,-0.5,3.5);


      // book(_h["Njge4_nJets"],"Njge4_nJets",9,3.5,12.5);
      // book(_h["Njge4_jet0Pt"],"Njge4_jet0Pt",20,0,600);

      book(_h["sumOfWeights"],"sumOfWeights",2,0,2);
      book(_h["Inclusive_sumW"],"Inclusive_sumW",1,-10000000,10000000);
      // book(_h["2lSS0tau_region1_sumW"],"2lSS0tau_region1_sumW",1,-10000000,10000000);
      // book(_h["2lSS0tau_region2_sumW"],"2lSS0tau_region2_sumW",1,-10000000,10000000);
      // book(_h["3l_region_sumW"],"3l_region_sumW",1,-10000000,10000000);


     // book(_h["2tau_region_DR_tau01"],"2tau_region_DR_tau01",10,0,2*M_PI);
     // book(_h["2tau_region_DEta_tau01"],"2tau_region_DEta_tau01",20,-2*M_PI,2*M_PI);
     // book(_h["2tau_region_DPhi_tau01"],"2tau_region_DPhi_tau01",20,-2*M_PI,2*M_PI);
     //
     // book(_h["2lSS0tau_region1_nJets"],"2lSS0tau_region1_nJets",11,-0.5,10.5);
     // book(_h["2lSS0tau_region1_HT"],"2lSS0tau_region1_HT",20,0,1200);
     // book(_h["2lSS0tau_region1_HT_jets"],"2lSS0tau_region1_HT_jets",20,0,1200);
     // book(_h["2lSS0tau_region1_nBjets"],"2lSS0tau_region1_nBjets",5,-0.5,4.5);
     // book(_h["2lSS0tau_region1_bjet0Pt"],"2lSS0tau_region1_bjet0Pt",20,0,400);
     // book(_h["2lSS0tau_region1_lep0Pt"],"2lSS0tau_region1_lep0Pt",20,0,300);
     // book(_h["2lSS0tau_region1_lep0Eta"],"2lSS0tau_region1_lep0Eta",20,-2.5,2.5);
     // book(_h["2lSS0tau_region1_lepJetMinDR"],"2lSS0tau_region1_lepJetMinDR",5,0,0.5);
     // book(_h["2lSS0tau_region1_DR_lep01"],"2lSS0tau_region1_DR_lep01",10,0,2*M_PI);
     // book(_h["2lSS0tau_region1_DEta_lep01"],"2lSS0tau_region1_DEta_lep01",10,0,5);
     // book(_h["2lSS0tau_region1_DPhi_lep01"],"2lSS0tau_region1_DPhi_lep01",10,0,4);
     // book(_h["2lSS0tau_region1_jet0Pt"],"2lSS0tau_region1_jet0Pt",20,0,600);
     // book(_h["2lSS0tau_region1_jet0Eta"],"2lSS0tau_region1_jet0Eta",10,-2.5,2.5);
     //
     // book(_h["2lSS0tau_region2_nJets"],"2lSS0tau_region2_nJets",11,-0.5,10.5);
     // book(_h["2lSS0tau_region2_HT"],"2lSS0tau_region2_HT",20,0,1200);
     // book(_h["2lSS0tau_region2_HT_jets"],"2lSS0tau_region2_HT_jets",20,0,1200);
     // book(_h["2lSS0tau_region2_nBjets"],"2lSS0tau_region2_nBjets",5,-0.5,4.5);
     // book(_h["2lSS0tau_region2_bjet0Pt"],"2lSS0tau_region2_bjet0Pt",20,0,400);
     // book(_h["2lSS0tau_region2_lep0Pt"],"2lSS0tau_region2_lep0Pt",20,0,300);
     // book(_h["2lSS0tau_region2_lep0Eta"],"2lSS0tau_region2_lep0Eta",20,-2.5,2.5);
     // book(_h["2lSS0tau_region2_lepJetMinDR"],"2lSS0tau_region2_lepJetMinDR",5,0,0.5);
     // book(_h["2lSS0tau_region2_DR_lep01"],"2lSS0tau_region2_DR_lep01",10,0,2*M_PI);
     // book(_h["2lSS0tau_region2_DEta_lep01"],"2lSS0tau_region2_DEta_lep01",10,0,5);
     // book(_h["2lSS0tau_region2_DPhi_lep01"],"2lSS0tau_region2_DPhi_lep01",10,0,4);
     // book(_h["2lSS0tau_region2_jet0Pt"],"2lSS0tau_region2_jet0Pt",20,0,600);
     // book(_h["2lSS0tau_region2_jet0Eta"],"2lSS0tau_region2_jet0Eta",10,-2.5,2.5);
     //
     //
     // book(_h["3l_region_nJets"],"3l_region_nJets",11,-0.5,10.5);
     // book(_h["3l_region_HT"],"3l_region_HT",20,0,1200);
     // book(_h["3l_region_HT_jets"],"3l_region_HT_jets",20,0,1200);
     // book(_h["3l_region_nBjets"],"3l_region_nBjets",5,-0.5,4.5);
     // book(_h["3l_region_bjet0Pt"],"3l_region_bjet0Pt",20,0,400);
     // book(_h["3l_region_lep0Pt"],"3l_region_lep0Pt",20,0,300);
     // book(_h["3l_region_lep0Eta"],"3l_region_lep0Eta",20,-2.5,2.5);
     // book(_h["3l_region_lepJetMinDR"],"3l_region_lepJetMinDR",5,0,0.5);
     // book(_h["3l_region_DR_lep01"],"3l_region_DR_lep01",10,0,2*M_PI);
     // book(_h["3l_region_DEta_lep01"],"3l_region_DEta_lep01",10,0,5);
     // book(_h["3l_region_DPhi_lep01"],"3l_region_DPhi_lep01",10,0,4);
     // book(_h["3l_region_jet0Pt"],"3l_region_jet0Pt",20,0,600);
     // book(_h["3l_region_jet0Eta"],"3l_region_jet0Eta",10,-2.5,2.5);



//      book(_h["2lSS0tau_region3_nJets"],"2lSS0tau_region3_nJets",10,-0.5,9.5);
//      book(_h["2lSS0tau_region3_HT"],"2lSS0tau_region3_HT",100,0,2000);
//      book(_h["2lSS0tau_region3_nBjets"],"2lSS0tau_region3_nBjets",4,-0.5,3.5);
//      book(_h["2lSS0tau_region3_bjet0pT"],"2lSS0tau_region3_bjet0Pt",50,0,500);
//      book(_h["2lSS0tau_region3_lep0pT"],"2lSS0tau_region3_lep0Pt",40,0,800);
//      book(_h["2lSS0tau_region3_lepJetMinDR"],"2lSS0tau_region3_lepJetminDR",5,0,0.5);
//      book(_h["2lSS0tau_region3_DR_lep01"],"2lSS0tau_region3_DR_lep01",10,0,2*M_PI);
//      book(_h["2lSS0tau_region3_DEta_lep01"],"2lSS0tau_region3_DEta_lep01",10,0,2*M_PI);
//      book(_h["2lSS0tau_region3_DPhi_lep01"],"2lSS0tau_region3_DPhi_lep01",10,0,2*M_PI);
//
//
//      book(_h["2lSS0tau_region4_nJets"],"2lSS0tau_region4_nJets",10,-0.5,9.5);
//      book(_h["2lSS0tau_region4_HT"],"2lSS0tau_region4_HT",100,0,2000);
//      book(_h["2lSS0tau_region4_nBjets"],"2lSS0tau_region4_nBjets",4,-0.5,3.5);
//      book(_h["2lSS0tau_region4_bjet0pT"],"2lSS0tau_region4_bjet0Pt",50,0,500);
//      book(_h["2lSS0tau_region4_lep0pT"],"2lSS0tau_region4_lep0Pt",40,0,800);
//      book(_h["2lSS0tau_region4_lepJetMinDR"],"2lSS0tau_region4_lepJetminDR",5,0,0.5);
//      book(_h["2lSS0tau_region4_DR_lep01"],"2lSS0tau_region4_DR_lep01",10,0,2*M_PI);
//      book(_h["2lSS0tau_region4_DEta_lep01"],"2lSS0tau_region4_DEta_lep01",10,0,2*M_PI);
//      book(_h["2lSS0tau_region4_DPhi_lep01"],"2lSS0tau_region4_DPhi_lep01",10,0,2*M_PI);
//
//      book(_h["2lSS0tau_region5_nJets"],"2lSS0tau_region5_nJets",10,-0.5,9.5);
//      book(_h["2lSS0tau_region5_HT"],"2lSS0tau_region5_HT",100,0,2000);
//      book(_h["2lSS0tau_region5_nBjets"],"2lSS0tau_region5_nBjets",4,-0.5,3.5);
//      book(_h["2lSS0tau_region5_bjet0pT"],"2lSS0tau_region5_bjet0Pt",50,0,500);
//      book(_h["2lSS0tau_region5_lep0pT"],"2lSS0tau_region5_lep0Pt",40,0,800);
//      book(_h["2lSS0tau_region5_lepJetMinDR"],"2lSS0tau_region5_lepJetminDR",5,0,0.5);
//      book(_h["2lSS0tau_region5_DR_lep01"],"2lSS0tau_region5_DR_lep01",10,0,2*M_PI);
//      book(_h["2lSS0tau_region5_DEta_lep01"],"2lSS0tau_region5_DEta_lep01",10,0,2*M_PI);
//      book(_h["2lSS0tau_region5_DPhi_lep01"],"2lSS0tau_region5_DPhi_lep01",10,0,2*M_PI);
//
//      book(_h["2lSS1tau_region6_nJets"],"2lSS1tau_region6_nJets",10,-0.5,9.5);
//      book(_h["2lSS1tau_region6_HT"],"2lSS1tau_region6_HT",100,0,2000);
//      book(_h["2lSS1tau_region6_nBjets"],"2lSS1tau_region6_nBjets",4,-0.5,3.5);
//      book(_h["2lSS1tau_region6_bjet0pT"],"2lSS1tau_region6_bjet0Pt",50,0,500);
//      book(_h["2lSS1tau_region6_lep0pT"],"2lSS1tau_region6_lep0Pt",40,0,800);
//      book(_h["2lSS1tau_region6_lepJetMinDR"],"2lSS1tau_region6_lepJetminDR",5,0,0.5);
//      book(_h["2lSS1tau_region6_DR_lep01"],"2lSS1tau_region6_DR_lep01",10,0,2*M_PI);
//      book(_h["2lSS1tau_region6_DEta_lep01"],"2lSS1tau_region6_DEta_lep01",10,0,2*M_PI);
//      book(_h["2lSS1tau_region6_DPhi_lep01"],"2lSS1tau_region6_DPhi_lep01",10,0,2*M_PI);
//      book(_h["2lSS1tau_region6_tauPt"],"2lSS1tau_region6_tauPt",50,0,500);
//
//
//      book(_h["2lSS1tau_region7_nJets"],"2lSS1tau_region7_nJets",10,-0.5,9.5);
//      book(_h["2lSS1tau_region7_HT"],"2lSS1tau_region7_HT",100,0,2000);
//      book(_h["2lSS1tau_region7_nBjets"],"2lSS1tau_region7_nBjets",4,-0.5,3.5);
//      book(_h["2lSS1tau_region7_bjet0pT"],"2lSS1tau_region7_bjet0Pt",50,0,500);
//      book(_h["2lSS1tau_region7_lep0pT"],"2lSS1tau_region7_lep0Pt",40,0,800);
//      book(_h["2lSS1tau_region7_lepJetMinDR"],"2lSS1tau_region7_lepJetminDR",5,0,0.5);
//      book(_h["2lSS1tau_region7_DR_lep01"],"2lSS1tau_region7_DR_lep01",10,0,2*M_PI);
//      book(_h["2lSS1tau_region7_DEta_lep01"],"2lSS1tau_region7_DEta_lep01",10,0,2*M_PI);
//      book(_h["2lSS1tau_region7_DPhi_lep01"],"2lSS1tau_region7_DPhi_lep01",10,0,2*M_PI);
//      book(_h["2lSS1tau_region7_tauPt"],"2lSS1tau_region7_tauPt",50,0,500);



    }


    void analyze(const Event& event) {
      // Use the "LFS" projection to require at least one hard charged
      // lepton. This is an experimental signature for the leptonically decaying
      // W. This helps to reduce pure QCD backgrounds.

      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MissingET");
      const double event_met	 = met.vectorEt().mod();

      /*if(zeeFinder.bosons().size()==0 && zmumuFinder.bosons().size()==0)
      {
          MSG_INFO("ZeeFinder size: "<<zeeFinder.size());
          MSG_INFO("ZmumuFinder size: "<<zmumuFinder.size());
          MSG_INFO("Veto Event");
          vetoEvent;
      }*/

      //
      Particles eMinusFromTaus, ePlusFromTaus, muonsFromTaus, antiMuonsFromTaus;
      for(const Particle& tau : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::TAU))
      {
          for(const Particle & p : tau.children())
          {
              if (p.pid()  == PID::EMINUS)
              {
                 eMinusFromTaus.push_back(p);
              }
              else if (p.pid() == PID::EPLUS)
              {
                 ePlusFromTaus.push_back(p);
              }
              else if (p.pid() == PID::MUON)
              {
                 muonsFromTaus.push_back(p);
              }
              else if (p.pid() == PID::ANTIMUON)
              {
                 antiMuonsFromTaus.push_back(p);
              }
          }
      }

 //     for(const Particle &p: eMinusFromTaus)
 //     {
 //         _h["eMinusFromTaus_pt"]->fill(p.pT()/GeV);
 //         _h["eMinusFromTaus_Eta"]->fill(p.eta());
 //         _h["eMinusFromTaus_phi"]->fill(p.phi());
 //     }

 //     for(const Particle &p: ePlusFromTaus)
 //     {
 //         _h["ePlusFromTaus_pt"]->fill(p.pT()/GeV);
 //         _h["ePlusFromTaus_Eta"]->fill(p.eta());
 //         _h["ePlusFromTaus_phi"]->fill(p.phi());
 //     }

 //     for(const Particle &p: muonsFromTaus)
 //     {
 //         _h["MuonsFromTaus_pt"]->fill(p.pT()/GeV);
 //         _h["MuonsFromTaus_Eta"]->fill(p.eta());
 //         _h["MuonsFromTaus_phi"]->fill(p.phi());
 //     }

 //     for(const Particle &p: antiMuonsFromTaus)
 //     {
 //         _h["AntiMuonsFromTaus_pt"]->fill(p.pT()/GeV);
 //         _h["AntiMuonsFromTaus_Eta"]->fill(p.eta());
 //         _h["AntiMuonsFromTaus_phi"]->fill(p.phi());
 //     }

      Particles elVec,muVec,tauVec,lepVec;
      Particles Inclusive_elVec,Inclusive_muVec, Inclusive_tauVec, Inclusive_allVec;
      //Count the total number of leptons
      //
      //

      for (const Particle & el: applyProjection<PromptFinalState>(event,"electrons").particlesByPt())
      {
          Inclusive_elVec.push_back(el);
          Inclusive_allVec.push_back(el);
          if(el.pT()/GeV > 10 && fabs(el.eta()) <2.5)
          {
              elVec.push_back(el);
      	      lepVec.push_back(el);
          }
      }
      for(const Particle &mu: applyProjection<PromptFinalState>(event,"muons").particlesByPt())
      {
          Inclusive_muVec.push_back(mu);
          Inclusive_allVec.push_back(mu);
          if(mu.pT()/GeV >10 && fabs(mu.eta()) <2.5)
          {
              muVec.push_back(mu);
       	      lepVec.push_back(mu);
          }
      }
      const TauFinder &tauhad = applyProjection<TauFinder>(event,"TauHadronic");
      for(const Particle &tau: tauhad.taus())
      {
          Inclusive_tauVec.push_back(tau);
          Inclusive_allVec.push_back(tau);
          if(tau.pT()/GeV >25 )
          {
            int nProng = countProngs(tau);
            if(nProng ==2 || nProng ==3)
            {
              tauVec.push_back(tau);
            }
          }
      }

      elVec = sortByPt(elVec);
      muVec = sortByPt(muVec);
      tauVec= sortByPt(tauVec);
      lepVec= sortByPt(lepVec);

      Inclusive_elVec  = sortByPt(Inclusive_elVec);
      Inclusive_muVec  = sortByPt(Inclusive_muVec);
      Inclusive_tauVec = sortByPt(Inclusive_tauVec);
      Inclusive_allVec = sortByPt(Inclusive_allVec);

      int nLep = lepVec.size();
      int elqsum=0;
      int muqsum=0;

      for(const Particle& el: elVec)
      {
          elqsum += el.charge();
      }
      for(const Particle &mu: muVec)
      {
          muqsum += mu.charge();
      }

      Jets alljets;
      //for(const Jet &jet : applyProjection<FastJets>(event, "Jets").jetsByPt(25*GeV))

      for(const Jet &jet : applyProjection<FastJets>(event, "Jets").jetsByPt(Cuts::pT> 25*GeV))
      {

         // alljets.push_back(jet);
          if(fabs(jet.eta()) < 2.5 )
          {
              alljets.push_back(jet);
          }
      }


      double ht_jets = 0.0;
      double ht = 0.0;
      for(const Jet& j: alljets) {
        ht_jets += j.pT();
        ht += j.pT();
      }
      for(const Particle & lep: lepVec){
        ht += lep.pT();
      }
      for(const Particle & taulep: tauVec){
        ht += taulep.pT();
      }


      // Identify b-jets
      const Particles bhadrons = sortByPt(applyProjection<HeavyHadrons>(event, "BCHadrons").bHadrons());

      Jets bjets, ljets;
      for(const Jet& jet: alljets)
      {
          if(jet.bTagged())
          {
              bjets.push_back(jet);
          }
          else
          {
              ljets.push_back(jet);
          }
     }
     alljets    =   sortByPt(alljets);
     bjets      =   sortByPt(bjets);
     ljets      =   sortByPt(ljets);

     float min_lj_deltaR=100;
     for(const Jet& jet: alljets){
         for(const Particle & part: lepVec){
    	     if(min_lj_deltaR > fabs(deltaR(jet,part))) {min_lj_deltaR = fabs(deltaR(jet,part)); }
      	 }
     }

     if (debug){ // check the lepton from tau are right
       Inclusive_elVec = sortByPt(Inclusive_elVec);
       Inclusive_muVec = sortByPt(Inclusive_muVec);
       eMinusFromTaus = sortByPt(eMinusFromTaus);
       ePlusFromTaus = sortByPt(ePlusFromTaus);

       MSG_INFO("All electron = " << Inclusive_elVec.size() << ", All muon = "<< Inclusive_muVec.size());
       MSG_INFO("electron from tau = " << eMinusFromTaus.size()+ePlusFromTaus.size() << ", muon from tau  = "<< muonsFromTaus.size()+antiMuonsFromTaus.size());
       MSG_INFO("lepton numer = " << lepVec.size() << " , electron number = " << elVec.size() << " , muon number " << muVec.size() << " , tau number = " << tauVec.size());

       if ( (Inclusive_elVec.size() == eMinusFromTaus.size()+ePlusFromTaus.size()) && (Inclusive_elVec.size()>0) ){
         float ele0 = Inclusive_elVec.at(0).pT()/GeV;

         float tmp1 = 0;
         float tmp2 = 0;
         float ele_tau0 = 0;

         if (eMinusFromTaus.size() > 0) {
            tmp1 = eMinusFromTaus.at(0).pT()/GeV;
         }

         if (ePlusFromTaus.size() > 0) {
            tmp2 = ePlusFromTaus.at(0).pT()/GeV;
         }

         if (tmp1 > tmp2){
          ele_tau0 = tmp1;
         }
         else if (tmp2 > tmp1){
          ele_tau0 = tmp2;
         }

         MSG_INFO("electron 0 pt = " << ele0 << ", electron 0 from tau pt =  "<< ele_tau0);
       }
     }

     _h["Inclusive_nLep_noHadTau"]->fill(Inclusive_elVec.size()+Inclusive_muVec.size());
     _h["Inclusive_nLep_withHadTau"]->fill(Inclusive_allVec.size());
     // _h["LepSelection_nLep_noHadTau"]->fill(nLep);
     // _h["LepSelection_nLep_withHadTau"]->fill(nLep+tauVec.size());


      _h["Inclusive_sumW"]->fill(1);
      _h["Inclusive_nJets"]->fill(alljets.size());
      _h["Inclusive_HT"]->fill(ht);
      _h["Inclusive_HT_jets"]->fill(ht_jets);
      _h["Inclusive_nBjets"]->fill(bjets.size());

      if (alljets.size() >= 1){
        _h["Inclusive_jet0Pt"]->fill(alljets.at(0).pT()/GeV);
        _h["Inclusive_jet0Eta"]->fill(alljets.at(0).eta());
      }
      else{
        _h["Inclusive_jet0Pt"]->fill(-99);
        _h["Inclusive_jet0Eta"]->fill(-99);
      }

    // if(alljets.size() > 4)  {
    //   _h["Njge4_nJets"]->fill(alljets.size());
    //   _h["Njge4_jet0Pt"]->fill(alljets.at(0).pT()/GeV);
    // }


    // if(tauVec.size()>=2){
    //     _h["2tau_region_DR_tau01"]->fill(fabs(deltaR(tauVec.at(0),tauVec.at(1))));
    //     _h["2tau_region_DEta_tau01"]->fill(fabs(deltaEta(tauVec.at(0),tauVec.at(1))));
    //     _h["2tau_region_DPhi_tau01"]->fill(fabs(deltaPhi(tauVec.at(0),tauVec.at(1))));
    // }

    //two light-leptons
    // if(nLep==2)
    // {
    //     //same sign + lepton pT
    //     // Region-1
    //     if(lepVec.at(0).charge()*lepVec.at(1).charge() >0 && lepVec.at(0).pT()/GeV >15 && lepVec.at(1).pT()/GeV > 15) // subleading 10 -> 15 GeV
    //     {
    //          _h["2lSS0tau_region1_nJets"]->fill(alljets.size());
    //          _h["2lSS0tau_region1_HT"]->fill(ht);
    //          _h["2lSS0tau_region1_HT_jets"]->fill(ht_jets);
    //          _h["2lSS0tau_region1_nBjets"]->fill(bjets.size());
    //          _h["2lSS0tau_region1_lep0Pt"]->fill(lepVec.at(0).pT()/GeV);
    //          _h["2lSS0tau_region1_lep0Eta"]->fill(lepVec.at(0).eta());
    //          _h["2lSS0tau_region1_lepJetMinDR"]->fill(min_lj_deltaR);
    //          _h["2lSS0tau_region1_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
    //          _h["2lSS0tau_region1_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
    //          _h["2lSS0tau_region1_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
    //
    //          _h["2lSS0tau_region1_sumW"]->fill(1);
    //
    //          if (alljets.size() >= 1){
    //            _h["2lSS0tau_region1_jet0Pt"]->fill(alljets.at(0).pT()/GeV);
    //            _h["2lSS0tau_region1_jet0Eta"]->fill(alljets.at(0).eta());
    //          }
    //          if (bjets.size()>=1){
    //            _h["2lSS0tau_region1_bjet0Pt"]->fill(bjets.at(0).pT()/GeV);
    //          }

	     // 0-tau
           //  if(tauVec.size()==0)
           //  {
       		 // // _h["2lSS0tau_MET"]->fill(event_met/GeV);
           //      // "Region-2"
           //      if(bjets.size()>=1 && alljets.size() >= 4)
           //      {
           // 		     _h["2lSS0tau_region2_nJets"]->fill(alljets.size());
           // 		     _h["2lSS0tau_region2_HT"]->fill(ht);
           // 		     _h["2lSS0tau_region2_HT_jets"]->fill(ht_jets);
           // 		     _h["2lSS0tau_region2_nBjets"]->fill(bjets.size());
           // 		     _h["2lSS0tau_region2_bjet0Pt"]->fill(bjets.at(0).pT()/GeV);
           // 		     _h["2lSS0tau_region2_lep0Pt"]->fill(lepVec.at(0).pT()/GeV);
           //         _h["2lSS0tau_region2_lep0Eta"]->fill(lepVec.at(0).eta());
           // 		     _h["2lSS0tau_region2_lepJetMinDR"]->fill(min_lj_deltaR);
           // 		     _h["2lSS0tau_region2_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
           // 		     _h["2lSS0tau_region2_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
           // 		     _h["2lSS0tau_region2_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
           //         _h["2lSS0tau_region2_jet0Pt"]->fill(alljets.at(0).pT()/GeV);
           //         _h["2lSS0tau_region2_jet0Eta"]->fill(alljets.at(0).eta());
           //
           //         _h["2lSS0tau_region2_sumW"]->fill(1);
                }
         //        // "Region-3"
         //        if(bjets.size()>=2 && alljets.size() >= 4)
         //        {
     		 //     _h["2lSS0tau_region3_nJets"]->fill(alljets.size());
     		 //     _h["2lSS0tau_region3_HT"]->fill(ht);
		     // _h["2lSS0tau_region3_nBjets"]->fill(bjets.size());
		     // _h["2lSS0tau_region3_bjet0pT"]->fill(bjets.at(0).pT()/GeV);
		     // _h["2lSS0tau_region3_lep0pT"]->fill(lepVec.at(0).pT()/GeV);
		     // _h["2lSS0tau_region3_lepJetMinDR"]->fill(min_lj_deltaR);
		     // _h["2lSS0tau_region3_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
		     // _h["2lSS0tau_region3_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
		     // _h["2lSS0tau_region3_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
         //
         //        }
         //        // "Region-4"
         //        if(bjets.size()==1 && alljets.size() >= 3)
         //        {
     		 //     _h["2lSS0tau_region4_nJets"]->fill(alljets.size());
     		 //     _h["2lSS0tau_region4_HT"]->fill(ht);
		     // _h["2lSS0tau_region4_nBjets"]->fill(bjets.size());
		     // _h["2lSS0tau_region4_bjet0pT"]->fill(bjets.at(0).pT()/GeV);
		     // _h["2lSS0tau_region4_lep0pT"]->fill(lepVec.at(0).pT()/GeV);
		     // _h["2lSS0tau_region4_lepJetMinDR"]->fill(min_lj_deltaR);
		     // _h["2lSS0tau_region4_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
		     // _h["2lSS0tau_region4_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
		     // _h["2lSS0tau_region4_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
         //
         //        }
         //        // "Region-5"
         //        if(bjets.size()>=2 && alljets.size() >= 3)
         //        {
     		 //     _h["2lSS0tau_region5_nJets"]->fill(alljets.size());
     		 //     _h["2lSS0tau_region5_HT"]->fill(ht);
		     // _h["2lSS0tau_region5_nBjets"]->fill(bjets.size());
		     // _h["2lSS0tau_region5_bjet0pT"]->fill(bjets.at(0).pT()/GeV);
		     // _h["2lSS0tau_region5_lep0pT"]->fill(lepVec.at(0).pT()/GeV);
		     // _h["2lSS0tau_region5_lepJetMinDR"]->fill(min_lj_deltaR);
		     // _h["2lSS0tau_region5_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
		     // _h["2lSS0tau_region5_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
		     // _h["2lSS0tau_region5_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
         //        }
            // }
	     // 1-tau
     //        else if (tauVec.size()>=1)
     //        {
	   //       _h["2lSS1tau_MET"]->fill(event_met/GeV);
		 // // "Region-6"
		 // _h["2lSS1tau_region6_nJets"]->fill(alljets.size());
     // 	         _h["2lSS1tau_region6_HT"]->fill(ht);
	   //       _h["2lSS1tau_region6_nBjets"]->fill(bjets.size());
	   //       _h["2lSS1tau_region6_lep0pT"]->fill(lepVec.at(0).pT()/GeV);
	   //       _h["2lSS1tau_region6_lepJetMinDR"]->fill(min_lj_deltaR);
	   //       _h["2lSS1tau_region6_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
	   //       _h["2lSS1tau_region6_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
	   //       _h["2lSS1tau_region6_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
		 // _h["2lSS1tau_region6_tauPt"]->fill(tauVec.at(0).pT()/GeV);
     //
	   //       // "Region-7"
     //            if(bjets.size() >= 1 && alljets.size() >= 3)
     //            {
	   //           _h["2lSS1tau_region7_nJets"]->fill(alljets.size());
     // 	             _h["2lSS1tau_region7_HT"]->fill(ht);
	   //           _h["2lSS1tau_region7_nBjets"]->fill(bjets.size());
	   //           _h["2lSS1tau_region7_bjet0pT"]->fill(bjets.at(0).pT()/GeV);
	   //           _h["2lSS1tau_region7_lep0pT"]->fill(lepVec.at(0).pT()/GeV);
	   //           _h["2lSS1tau_region7_lepJetMinDR"]->fill(min_lj_deltaR);
	   //           _h["2lSS1tau_region7_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
	   //           _h["2lSS1tau_region7_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
	   //           _h["2lSS1tau_region7_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
		 //     _h["2lSS1tau_region7_tauPt"]->fill(tauVec.at(0).pT()/GeV);
	   //        }
     //       }
     //     }
     // }
    //  else if (nLep==3){
    //    if( abs(elqsum + muqsum) ==1){
    //      if(lepVec.at(0).pT()/GeV >15 && lepVec.at(1).pT()/GeV > 15 && lepVec.at(2).pT()/GeV > 15){
    //        if(tauVec.size()==0){
    //          if(bjets.size()>=1 && alljets.size() >= 2){
    //        		     _h["3l_region_nJets"]->fill(alljets.size());
    //        		     _h["3l_region_HT"]->fill(ht);
    //                _h["3l_region_HT_jets"]->fill(ht_jets);
    //        		     _h["3l_region_nBjets"]->fill(bjets.size());
    //        		     _h["3l_region_bjet0Pt"]->fill(bjets.at(0).pT()/GeV);
    //        		     _h["3l_region_lep0Pt"]->fill(lepVec.at(0).pT()/GeV);
    //                _h["3l_region_lep0Eta"]->fill(lepVec.at(0).eta());
    //        		     _h["3l_region_lepJetMinDR"]->fill(min_lj_deltaR);
    //        		     _h["3l_region_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
    //        		     _h["3l_region_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
    //        		     _h["3l_region_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
    //                _h["3l_region_jet0Pt"]->fill(alljets.at(0).pT()/GeV);
    //                _h["3l_region_jet0Eta"]->fill(alljets.at(0).eta());
    //                _h["3l_region_sumW"]->fill(1);
    //
    //            }
    //
    //          }
    //      }
    //    }
    // }

  }


    void finalize() {
        //
        // For Powheg
        //
        // MSG_INFO("CROSS SSECTION:"<<crossSection());
        // auto xsec = isnan(crossSection()) ? 1 : crossSection();
        // MSG_INFO("xsec = "<< xsec);
        // MSG_INFO("Sum of weights:"<<sumOfWeights());
        // const double sf = xsec / sumOfWeights();
        // _h["sumOfWeights"]->fill(xsec); // histograms are scaled to xs

        // For Sherpa/MG
        MSG_INFO("CROSS SSECTION:"<<crossSection());
        MSG_INFO("Sum of weights:"<<sumOfWeights());
        const double sf = crossSection() / sumOfWeights();
        for (auto hist : _h) { scale(hist.second, sf); }

        _h["sumOfWeights"]->fill(1);



    }

    //@}


  private:

    // @name Histogram data members
    //@{
    //
    map<std::string,Histo1DPtr> _h;
    // Histo1DPtr _h_nosel_nJets_;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ttw_ttH);
}
