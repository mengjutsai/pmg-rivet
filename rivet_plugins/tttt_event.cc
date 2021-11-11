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



  class tttt_event: public Analysis {
  public:

    /// Minimal constructor
    tttt_event() : Analysis("tttt_event")
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
      book(_h["sumOfWeights"],"sumOfWeights",2,0,2);


      book(_h["Inclusive_nJets"],"Inclusive_nJets",15,-0.5,14.5);
      book(_h["Inclusive_HT"],"Inclusive_HT",15,0,3000);
      book(_h["Inclusive_HT_jets"],"Inclusive_HT_jets",15,0,3000);
      book(_h["Inclusive_nBjets"],"Inclusive_nBjets",7,-0.5,6.5);
      book(_h["Inclusive_nHadTau"],"Inclusive_nHadTau",4,-0.5,3.5);
      book(_h["Inclusive_jet0Pt"],"Inclusive_jet0Pt",20,0,600);
      book(_h["Inclusive_jet0Eta"],"Inclusive_jet0Eta",10,-2.5,2.5);

      book(_h["Inclusive_lep0Pt"],"Inclusive_lep0Pt",20,0,600);
      book(_h["Inclusive_lep0Eta"],"Inclusive_lep0Eta",10,-2.5,2.5);
      book(_h["Inclusive_lep1Pt"],"Inclusive_lep1Pt",20,0,600);
      book(_h["Inclusive_lep1Eta"],"Inclusive_lep1Eta",10,-2.5,2.5);
      book(_h["Inclusive_entries"],"Inclusive_entries",3,-0.5,3.5);


     book(_h["2lSS3l_nJets"],"2lSS3l_nJets",15,-0.5,14.5);
     book(_h["2lSS3l_HT"],"2lSS3l_HT",15,0,3000);
     book(_h["2lSS3l_HT_jets"],"2lSS3l_HT_jets",15,0,3000);
     book(_h["2lSS3l_nBjets"],"2lSS3l_nBjets",7,-0.5,6.5);
     book(_h["2lSS3l_lep0Pt"],"2lSS3l_lep0Pt",20,0,300);
     book(_h["2lSS3l_lep0Eta"],"2lSS3l_lep0Eta",20,-2.5,2.5);
     book(_h["2lSS3l_lepJetMinDR"],"2lSS3l_lepJetMinDR",5,0,0.5);
     book(_h["2lSS3l_DR_lep01"],"2lSS3l_DR_lep01",10,0,2*M_PI);
     book(_h["2lSS3l_DEta_lep01"],"2lSS3l_DEta_lep01",10,0,5);
     book(_h["2lSS3l_DPhi_lep01"],"2lSS3l_DPhi_lep01",10,0,4);
     book(_h["2lSS3l_jet0Pt"],"2lSS3l_jet0Pt",20,0,600);
     book(_h["2lSS3l_jet0Eta"],"2lSS3l_jet0Eta",10,-2.5,2.5);
     book(_h["2lSS3l_bjet0Pt"],"2lSS3l_bjet0Pt",20,0,400);
     book(_h["2lSS3l_entries"],"2lSS3l_entries",3,-0.5,3.5);




     book(_h["2lSS_nJets"],"2lSS_nJets",15,-0.5,14.5);
     book(_h["2lSS_HT"],"2lSS_HT",15,0,3000);
     book(_h["2lSS_HT_jets"],"2lSS_HT_jets",15,0,3000);
     book(_h["2lSS_nBjets"],"2lSS_nBjets",7,-0.5,6.5);
     book(_h["2lSS_lep0Pt"],"2lSS_lep0Pt",20,0,300);
     book(_h["2lSS_lep0Eta"],"2lSS_lep0Eta",20,-2.5,2.5);
     book(_h["2lSS_lepJetMinDR"],"2lSS_lepJetMinDR",5,0,0.5);
     book(_h["2lSS_DR_lep01"],"2lSS_DR_lep01",10,0,2*M_PI);
     book(_h["2lSS_DEta_lep01"],"2lSS_DEta_lep01",10,0,5);
     book(_h["2lSS_DPhi_lep01"],"2lSS_DPhi_lep01",10,0,4);
     book(_h["2lSS_jet0Pt"],"2lSS_jet0Pt",20,0,600);
     book(_h["2lSS_jet0Eta"],"2lSS_jet0Eta",10,-2.5,2.5);
     book(_h["2lSS_bjet0Pt"],"2lSS_bjet0Pt",20,0,400);
     book(_h["2lSS_entries"],"2lSS_entries",3,-0.5,3.5);

     book(_h["2lSS_1b4j_nJets"],"2lSS_1b4j_nJets",15,-0.5,14.5);
     book(_h["2lSS_1b4j_HT"],"2lSS_1b4j_HT",15,0,3000);
     book(_h["2lSS_1b4j_HT_jets"],"2lSS_1b4j_HT_jets",15,0,3000);
     book(_h["2lSS_1b4j_nBjets"],"2lSS_1b4j_nBjets",7,-0.5,6.5);
     book(_h["2lSS_1b4j_bjet0Pt"],"2lSS_1b4j_bjet0Pt",20,0,400);
     book(_h["2lSS_1b4j_lep0Pt"],"2lSS_1b4j_lep0Pt",20,0,300);
     book(_h["2lSS_1b4j_lep0Eta"],"2lSS_1b4j_lep0Eta",20,-2.5,2.5);
     book(_h["2lSS_1b4j_lepJetMinDR"],"2lSS_1b4j_lepJetMinDR",5,0,0.5);
     book(_h["2lSS_1b4j_DR_lep01"],"2lSS_1b4j_DR_lep01",10,0,2*M_PI);
     book(_h["2lSS_1b4j_DEta_lep01"],"2lSS_1b4j_DEta_lep01",10,0,5);
     book(_h["2lSS_1b4j_DPhi_lep01"],"2lSS_1b4j_DPhi_lep01",10,0,4);
     book(_h["2lSS_1b4j_jet0Pt"],"2lSS_1b4j_jet0Pt",20,0,600);
     book(_h["2lSS_1b4j_jet0Eta"],"2lSS_1b4j_jet0Eta",10,-2.5,2.5);
     book(_h["2lSS_1b4j_entries"],"2lSS_1b4j_entries",3,-0.5,3.5);


     book(_h["3l_nJets"],"3l_nJets",15,-0.5,14.5);
     book(_h["3l_HT"],"3l_HT",15,0,3000);
     book(_h["3l_HT_jets"],"3l_HT_jets",15,0,3000);
     book(_h["3l_nBjets"],"3l_nBjets",7,-0.5,6.5);
     book(_h["3l_lep0Pt"],"3l_lep0Pt",20,0,300);
     book(_h["3l_lep0Eta"],"3l_lep0Eta",20,-2.5,2.5);
     book(_h["3l_lepJetMinDR"],"3l_lepJetMinDR",5,0,0.5);
     book(_h["3l_DR_lep01"],"3l_DR_lep01",10,0,2*M_PI);
     book(_h["3l_DEta_lep01"],"3l_DEta_lep01",10,0,5);
     book(_h["3l_DPhi_lep01"],"3l_DPhi_lep01",10,0,4);
     book(_h["3l_jet0Pt"],"3l_jet0Pt",20,0,600);
     book(_h["3l_jet0Eta"],"3l_jet0Eta",10,-2.5,2.5);
     book(_h["3l_bjet0Pt"],"3l_bjet0Pt",20,0,400);
     book(_h["3l_entries"],"3l_entries",3,-0.5,3.5);


     book(_h["3l_1b2j_nJets"],"3l_1b2j_nJets",15,-0.5,14.5);
     book(_h["3l_1b2j_HT"],"3l_1b2j_HT",15,0,3000);
     book(_h["3l_1b2j_HT_jets"],"3l_1b2j_HT_jets",15,0,3000);
     book(_h["3l_1b2j_nBjets"],"3l_1b2j_nBjets",7,-0.5,6.5);
     book(_h["3l_1b2j_bjet0Pt"],"3l_1b2j_bjet0Pt",20,0,400);
     book(_h["3l_1b2j_lep0Pt"],"3l_1b2j_lep0Pt",20,0,300);
     book(_h["3l_1b2j_lep0Eta"],"3l_1b2j_lep0Eta",20,-2.5,2.5);
     book(_h["3l_1b2j_lepJetMinDR"],"3l_1b2j_lepJetMinDR",5,0,0.5);
     book(_h["3l_1b2j_DR_lep01"],"3l_1b2j_DR_lep01",10,0,2*M_PI);
     book(_h["3l_1b2j_DEta_lep01"],"3l_1b2j_DEta_lep01",10,0,5);
     book(_h["3l_1b2j_DPhi_lep01"],"3l_1b2j_DPhi_lep01",10,0,4);
     book(_h["3l_1b2j_jet0Pt"],"3l_1b2j_jet0Pt",20,0,600);
     book(_h["3l_1b2j_jet0Eta"],"3l_1b2j_jet0Eta",10,-2.5,2.5);
     book(_h["3l_1b2j_entries"],"3l_1b2j_entries",3,-0.5,3.5);


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

      elVec = sortByPt(elVec);
      muVec = sortByPt(muVec);
      lepVec= sortByPt(lepVec);


      Inclusive_elVec  = sortByPt(Inclusive_elVec);
      Inclusive_muVec  = sortByPt(Inclusive_muVec);
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

       // Particles test_vec;
       // std::vector<FourMomentum> test_vec;
       // for(const Jet& jet: alljets)
       // {
       //   test_vec.push_back(jet);
       //      // MSG_INFO("test_vec = " << test_vec.size() << ", alljets = "<< alljets.size() );
       //
       // }

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
      // MSG_INFO("alljets = " << alljets.size() << ", alljets_validation = "<< alljets_validation.size() );


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

     // MSG_INFO("bjets = " << bjets.size() << ", bjets_validation = "<< bjets_validation.size() );
     // MSG_INFO("ljets = " << ljets.size() << ", ljets_validation = "<< ljets_validation.size() );


     alljets    =   sortByPt(alljets);
     bjets      =   sortByPt(bjets);
     ljets      =   sortByPt(ljets);


     // Include Hadronic tau in the jet collection
     std::vector<FourMomentum> alljets_withHadTau, ljets_withHadTau;
     alljets_withHadTau = alljets;
     ljets_withHadTau = ljets;

      const TauFinder &tauhad = applyProjection<TauFinder>(event,"TauHadronic");
      for(const Particle &tau: tauhad.taus())
      {
          Inclusive_tauVec.push_back(tau);
          Inclusive_allVec.push_back(tau);
          if(tau.pT()/GeV >25 && fabs(tau.eta()) < 2.5 )
          {
            int nProng = countProngs(tau);
            if(nProng ==2 || nProng ==3)
            {
              tauVec.push_back(tau);
              alljets_withHadTau.push_back(tau);
              ljets_withHadTau.push_back(tau);
            }
          }
      }

      tauVec= sortByPt(tauVec);
      Inclusive_tauVec = sortByPt(Inclusive_tauVec);

     alljets_withHadTau  =   sortByPt(alljets_withHadTau);
     ljets_withHadTau    =   sortByPt(ljets_withHadTau);

     // MSG_INFO("tauVec = " << tauVec.size()  );
     // MSG_INFO("alljets = " << alljets.size() << ", alljets_withHadTau = "<< alljets_withHadTau.size() );
     // MSG_INFO("ljets = " << ljets.size() << ", ljets_withHadTau = "<< ljets_withHadTau.size() );


      double ht_jets = 0.0;
      double ht = 0.0;
      for(const FourMomentum& j: alljets_withHadTau) {
        ht_jets += j.pT();
        ht += j.pT();
      }
      for(const Particle & lep: lepVec){
        ht += lep.pT();
      }


      float min_lj_deltaR=100;
      for(const FourMomentum& jet: alljets_withHadTau){ // Use alljets+Hadronic Tau
          for(const Particle & part: lepVec){
     	     if(min_lj_deltaR > fabs(deltaR(jet,part))) {min_lj_deltaR = fabs(deltaR(jet,part)); }
       	 }
      }


     if (debug){ // to check the lepton from tau is included in the inclusive electron/muon collections -
       // conclusion: yes, ele/mu collections include the e/m from leptonic tau
       Inclusive_elVec = sortByPt(Inclusive_elVec);
       Inclusive_muVec = sortByPt(Inclusive_muVec);
       eMinusFromTaus = sortByPt(eMinusFromTaus);
       ePlusFromTaus = sortByPt(ePlusFromTaus);

       MSG_INFO("All electron = " << Inclusive_elVec.size() << ", All muon = "<< Inclusive_muVec.size());
       MSG_INFO("electron from tau = " << eMinusFromTaus.size()+ePlusFromTaus.size() << ", muon from tau  = "<< muonsFromTaus.size()+antiMuonsFromTaus.size());
       MSG_INFO("After 10GeV, abs(eta) < 2.5");
       MSG_INFO(">> lepton numer = " << lepVec.size() << " , electron number = " << elVec.size() << " , muon number " << muVec.size() << " , tau number = " << tauVec.size());

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

      _h["Inclusive_nJets"]->fill(alljets_withHadTau.size());
      _h["Inclusive_HT"]->fill(ht);
      _h["Inclusive_HT_jets"]->fill(ht_jets);
      _h["Inclusive_nBjets"]->fill(bjets.size());
      _h["Inclusive_nHadTau"]->fill(tauVec.size());
      _h["Inclusive_entries"]->fill(1,1);

      if (alljets_withHadTau.size() >= 1){
        _h["Inclusive_jet0Pt"]->fill(alljets_withHadTau.at(0).pT()/GeV);
        _h["Inclusive_jet0Eta"]->fill(alljets_withHadTau.at(0).eta());
      }
      else{
        _h["Inclusive_jet0Pt"]->fill(-99);
        _h["Inclusive_jet0Eta"]->fill(-99);
      }

      if (lepVec.size() >= 1){
        _h["Inclusive_lep0Pt"]->fill(lepVec.at(0).pT()/GeV);
        _h["Inclusive_lep0Eta"]->fill(lepVec.at(0).eta());
      }

      if (lepVec.size() >= 2){
        _h["Inclusive_lep1Pt"]->fill(lepVec.at(1).pT()/GeV);
        _h["Inclusive_lep1Eta"]->fill(lepVec.at(1).eta());
      }


    // 2LSS+3L region
    if(  ( nLep==2 && (lepVec.at(0).charge()*lepVec.at(1).charge() >0 && lepVec.at(0).pT()/GeV >15 && lepVec.at(1).pT()/GeV > 15) ) || ( nLep==3 &&  abs(elqsum + muqsum) ==1 && (lepVec.at(0).pT()/GeV >15 && lepVec.at(1).pT()/GeV > 15 && lepVec.at(2).pT()/GeV > 15) ) ){

       _h["2lSS3l_nJets"]->fill(alljets_withHadTau.size());
       _h["2lSS3l_HT"]->fill(ht);
       _h["2lSS3l_HT_jets"]->fill(ht_jets);
       _h["2lSS3l_nBjets"]->fill(bjets.size());
       _h["2lSS3l_lep0Pt"]->fill(lepVec.at(0).pT()/GeV);
       _h["2lSS3l_lep0Eta"]->fill(lepVec.at(0).eta());
       _h["2lSS3l_lepJetMinDR"]->fill(min_lj_deltaR);
       _h["2lSS3l_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
       _h["2lSS3l_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
       _h["2lSS3l_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
       _h["2lSS3l_entries"]->fill(1,1);

       if (alljets_withHadTau.size() >= 1){
         _h["2lSS3l_jet0Pt"]->fill(alljets_withHadTau.at(0).pT()/GeV);
         _h["2lSS3l_jet0Eta"]->fill(alljets_withHadTau.at(0).eta());
       }
       if (bjets.size()>=1){
         _h["2lSS3l_bjet0Pt"]->fill(bjets.at(0).pT()/GeV);
       }





    }

    //two light-leptons
    if(nLep==2)
    {
        //same sign + lepton pT
        // Region-1
        if(lepVec.at(0).charge()*lepVec.at(1).charge() >0 && lepVec.at(0).pT()/GeV >15 && lepVec.at(1).pT()/GeV > 15) // subleading 10 -> 15 GeV
        {
             _h["2lSS_nJets"]->fill(alljets_withHadTau.size());
             _h["2lSS_HT"]->fill(ht);
             _h["2lSS_HT_jets"]->fill(ht_jets);
             _h["2lSS_nBjets"]->fill(bjets.size());
             _h["2lSS_lep0Pt"]->fill(lepVec.at(0).pT()/GeV);
             _h["2lSS_lep0Eta"]->fill(lepVec.at(0).eta());
             _h["2lSS_lepJetMinDR"]->fill(min_lj_deltaR);
             _h["2lSS_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
             _h["2lSS_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
             _h["2lSS_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
             _h["2lSS_entries"]->fill(1,1);

             // _h["2lSS_sumW"]->fill(1);

             if (alljets_withHadTau.size() >= 1){
               _h["2lSS_jet0Pt"]->fill(alljets_withHadTau.at(0).pT()/GeV);
               _h["2lSS_jet0Eta"]->fill(alljets_withHadTau.at(0).eta());
             }
             if (bjets.size()>=1){
               _h["2lSS_bjet0Pt"]->fill(bjets.at(0).pT()/GeV);
             }

                // 1b+4j
                if(bjets.size()>=1 && alljets_withHadTau.size() >= 4)
                {
           		     _h["2lSS_1b4j_nJets"]->fill(alljets_withHadTau.size());
           		     _h["2lSS_1b4j_HT"]->fill(ht);
           		     _h["2lSS_1b4j_HT_jets"]->fill(ht_jets);
           		     _h["2lSS_1b4j_nBjets"]->fill(bjets.size());
           		     _h["2lSS_1b4j_bjet0Pt"]->fill(bjets.at(0).pT()/GeV);
           		     _h["2lSS_1b4j_lep0Pt"]->fill(lepVec.at(0).pT()/GeV);
                   _h["2lSS_1b4j_lep0Eta"]->fill(lepVec.at(0).eta());
           		     _h["2lSS_1b4j_lepJetMinDR"]->fill(min_lj_deltaR);
           		     _h["2lSS_1b4j_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
           		     _h["2lSS_1b4j_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
           		     _h["2lSS_1b4j_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
                   _h["2lSS_1b4j_jet0Pt"]->fill(alljets_withHadTau.at(0).pT()/GeV);
                   _h["2lSS_1b4j_jet0Eta"]->fill(alljets_withHadTau.at(0).eta());
                   // _h["2lSS_1b4j_sumW"]->fill(1);
                   _h["2lSS_1b4j_entries"]->fill(1,1);

                }
         }
     }
     else if (nLep==3){
       if( abs(elqsum + muqsum) ==1){
         if(lepVec.at(0).pT()/GeV >15 && lepVec.at(1).pT()/GeV > 15 && lepVec.at(2).pT()/GeV > 15){
         		     _h["3l_nJets"]->fill(alljets_withHadTau.size());
         		     _h["3l_HT"]->fill(ht);
                 _h["3l_HT_jets"]->fill(ht_jets);
         		     _h["3l_nBjets"]->fill(bjets.size());
         		     _h["3l_lep0Pt"]->fill(lepVec.at(0).pT()/GeV);
                 _h["3l_lep0Eta"]->fill(lepVec.at(0).eta());
         		     _h["3l_lepJetMinDR"]->fill(min_lj_deltaR);
         		     _h["3l_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
         		     _h["3l_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
         		     _h["3l_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
                 // _h["3l_sumW"]->fill(1);

                 if (alljets_withHadTau.size() >= 1){
                   _h["3l_jet0Pt"]->fill(alljets_withHadTau.at(0).pT()/GeV);
                   _h["3l_jet0Eta"]->fill(alljets_withHadTau.at(0).eta());
                 }
                 if (bjets.size()>=1){
          		     _h["3l_bjet0Pt"]->fill(bjets.at(0).pT()/GeV);

                 }
                  _h["3l_entries"]->fill(1,1);


             if(bjets.size()>=1 && alljets_withHadTau.size() >= 2){
           		     _h["3l_1b2j_nJets"]->fill(alljets_withHadTau.size());
           		     _h["3l_1b2j_HT"]->fill(ht);
                   _h["3l_1b2j_HT_jets"]->fill(ht_jets);
           		     _h["3l_1b2j_nBjets"]->fill(bjets.size());
           		     _h["3l_1b2j_bjet0Pt"]->fill(bjets.at(0).pT()/GeV);
           		     _h["3l_1b2j_lep0Pt"]->fill(lepVec.at(0).pT()/GeV);
                   _h["3l_1b2j_lep0Eta"]->fill(lepVec.at(0).eta());
           		     _h["3l_1b2j_lepJetMinDR"]->fill(min_lj_deltaR);
           		     _h["3l_1b2j_DR_lep01"]->fill(fabs(deltaR(lepVec.at(0),lepVec.at(1))));
           		     _h["3l_1b2j_DEta_lep01"]->fill(fabs(deltaEta(lepVec.at(0),lepVec.at(1))));
           		     _h["3l_1b2j_DPhi_lep01"]->fill(fabs(deltaPhi(lepVec.at(0),lepVec.at(1))));
                   _h["3l_1b2j_jet0Pt"]->fill(alljets_withHadTau.at(0).pT()/GeV);
                   _h["3l_1b2j_jet0Eta"]->fill(alljets_withHadTau.at(0).eta());
                   // _h["3l_1b2j_sumW"]->fill(1);
                   _h["3l_1b2j_entries"]->fill(1,1);


               }
         }
       }
    }

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
  DECLARE_RIVET_PLUGIN(tttt_event);
}
