// -*- C++ -*-
// #include "Root/TH1F.h"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/HadronicFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Tools/RivetMT2.hh"


//MARIA (ToDo)
// RIVETanalysis: add histos for ttW_pt/eta/rap, tt_eta/rap (pt already exists)
// For the analysis, will spit into two ttW channels: 3l (with dilepton tt: TopDil_Lep1, TopDil_Lep2) and 2lSS (with lep+jets tt: TopLJets_Lep, TopLJets_Had)


/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...




namespace Rivet {

  using namespace Cuts;
  bool debug = false;

  bool hasChild(const HepMC::GenParticle *ptcl, int pdgID) { //If the parent particle decay to some child particles
      for (auto child:Particle(*ptcl).children())
          if (child.pid()==pdgID) return true;
      return false;
  }
 bool hasParent(const HepMC::GenParticle *ptcl, int pdgID) {
      for (auto parent:Rivet::HepMCUtils::particles(ptcl->production_vertex(),HepMC::parents))
        if (parent->pdg_id()==pdgID) return true;
      return false;
  }

  Particle getLastInstance(Particle ptcl) {
    if ( ptcl.genParticle()->end_vertex() ) { // check if has end vertex
      if ( !hasChild(ptcl.genParticle(),ptcl.pid()) ) return ptcl;
      else return getLastInstance(ptcl.children()[0]);
    }
    return ptcl;
  }

  class tttt_parton: public Analysis {
  public:

    /// Constructor
    tttt_parton()
      : Analysis("tttt_parton")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

        declare(FastJets(FinalState(Cuts::abseta < 2.5 && Cuts::pT>25*GeV),FastJets::ANTIKT,0.4),"Jets");
        declare(MissingMomentum(FinalState(Cuts::abseta < 4 && Cuts::pT >0*GeV)),"MissingET");

        book(_h["sumOfWeights"],"sumOfWeights",2,0,2);

        book(_h["Inclusive_nTop"],"Inclusive_nTop",5,0.5,5.5);
        book(_h["Inclusive_nH"],"Inclusive_nH",4,-0.5,3.5);
        book(_h["Inclusive_nTop_FromNonH"],"Inclusive_nTop_FromNonH",5,-0.5,4.5);
        book(_h["Inclusive_nTop_FromH"],"Inclusive_nTop_FromH",5,-0.5,4.5);

        book(_h["Inclusive_pT_ttbar_FromH"],"Inclusive_pT_ttbar_FromH",20,0,600);
        book(_h["Inclusive_mTop1_FromH"],"Inclusive_mTop1_FromH",25,150,200);
        book(_h["Inclusive_mTop2_FromH"],"Inclusive_mTop2_FromH",25,150,200);
        book(_h["Inclusive_mTop1_zoomin_FromH"],"Inclusive_mTop1_zoomin_FromH",20,165,175);
        book(_h["Inclusive_mTop2_zoomin_FromH"],"Inclusive_mTop2_zoomin_FromH",20,165,175);
        book(_h["Inclusive_dR_ttbar_FromH"],"Inclusive_dR_ttbar_FromH",20,0,5);

        book(_h["Inclusive_InvM_FromH"],"Inclusive_InvM_FromH",180,300,1200);
        book(_h["Inclusive_InvM_FromH400"],"Inclusive_InvM_FromH400",20,300,500);
        book(_h["Inclusive_InvM_FromH500"],"Inclusive_InvM_FromH500",20,400,600);
        book(_h["Inclusive_InvM_FromH600"],"Inclusive_InvM_FromH600",20,500,700);
        book(_h["Inclusive_InvM_FromH700"],"Inclusive_InvM_FromH700",20,600,800);
        book(_h["Inclusive_InvM_FromH800"],"Inclusive_InvM_FromH800",20,700,900);
        book(_h["Inclusive_InvM_FromH900"],"Inclusive_InvM_FromH900",20,800,1000);
        book(_h["Inclusive_InvM_FromH1000"],"Inclusive_InvM_FromH1000",20,900,1100);

        book(_h["Inclusive_pT_ttbar_FromNonH"],"Inclusive_pT_ttbar_FromNonH",20,0,600);
        book(_h["Inclusive_mTop1_FromNonH"],"Inclusive_mTop1_FromNonH",25,150,200);
        book(_h["Inclusive_mTop2_FromNonH"],"Inclusive_mTop2_FromNonH",25,150,200);
        book(_h["Inclusive_mTop1_zoomin_FromNonH"],"Inclusive_mTop1_zoomin_FromNonH",20,165,175);
        book(_h["Inclusive_mTop2_zoomin_FromNonH"],"Inclusive_mTop2_zoomin_FromNonH",20,165,175);
        book(_h["Inclusive_dR_ttbar_FromNonH"],"Inclusive_dR_ttbar_FromNonH",20,0,5);


        book(_h["Inclusive_mH"],"Inclusive_mH",180,300,1200);
        book(_h["Inclusive_mH_zoomin"],"Inclusive_mH_zoomin",20,350,450);

        book(_h["Inclusive_mH400"],"Inclusive_mH400",20,300,500);
        book(_h["Inclusive_mH500"],"Inclusive_mH500",20,400,600);
        book(_h["Inclusive_mH600"],"Inclusive_mH600",20,500,700);
        book(_h["Inclusive_mH700"],"Inclusive_mH700",20,600,800);
        book(_h["Inclusive_mH800"],"Inclusive_mH800",20,700,900);
        book(_h["Inclusive_mH900"],"Inclusive_mH900",20,800,1000);
        book(_h["Inclusive_mH1000"],"Inclusive_mH1000",20,900,1100);

        book(_h["Inclusive_pT_H"],"Inclusive_pT_H",25,0,800);


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // const double weight = event.weight();

      // MSG_INFO("#----------------Event--------------#");

      //DO PDG analysis
      vector<const HepMC::GenParticle*> genParticles = HepMCUtils::particles(event.genEvent());

      Particles HCands;
      Particles topCands,topCandsFromH,topCandsFromNonH;

      for (ConstGenParticlePtr part : HepMCUtils::particles(event.genEvent()))
      {
              //Select Tops
              if(fabs(part->pdg_id()) == 6)
              {

                  if ( part->production_vertex() && !hasParent(part, part->pdg_id()) && Particle(part).children().size()>0)
                  {
                      if(debug){
                          std::cout <<"production Top: "<< part->pdg_id()
                                     << ", status: " <<  part->status()
                                     << ", barcode: " << part->barcode()
                                    <<std::endl;
                          std::cout << "pt: " <<  part->momentum().perp()
                                    << ", eta: " <<  part->momentum().eta()
                                    << ", phi: " <<  part->momentum().phi()
                                   <<std::endl;
                          std::cout <<"Top parent number: "<< Particle(part).parents().size()
                                    << ", children number: " << Particle(part).children().size()
                                   <<std::endl;
                          std::cout << "------------------"<<std::endl;
                      }

                      // append W candidates
                      topCands.push_back(Particle(part));

                      if (!( part->production_vertex() && ( hasParent(part,35) ) ) ) {
                          topCandsFromNonH.push_back(Particle(part));
                      }


                    }

              }

              //Select heavy higgs bosons
              if( (fabs(part->pdg_id())==35 || fabs(part->pdg_id())==36) &&Particle(part).children().size())
              {

                    if ( part->production_vertex() && !hasParent(part, part->pdg_id()) && Particle(part).children().size()>0 )
                    {
                     if(debug){
                      std::cout <<"production Heavy Higgs: "<< part->pdg_id()
                                   << ", status: " <<  part->status()
                                   << ", barcode: " << part->barcode()
                                  <<std::endl;
                      std::cout << "pt: " <<  part->momentum().perp()
                                << ", eta: " <<  part->momentum().eta()
                                << ", phi: " <<  part->momentum().phi()
                               <<std::endl;
                      std::cout <<"Heavy Higgs parent number: "<< Particle(part).parents().size()
                                << ", children number: " << Particle(part).children().size()
                               <<std::endl;
                      std::cout << "------------------"<<std::endl;
                     }
                      HCands.push_back(Particle(part));
                    }

              }
          // }
      }

      // MSG_INFO("Number of W boson : "   << wCands.size());

      _h["Inclusive_nTop"]->fill(topCands.size());
      _h["Inclusive_nH"]->fill(HCands.size());
      _h["Inclusive_nTop_FromNonH"]->fill(topCandsFromNonH.size());


      for(const Particle &Higgs: HCands){
        if (debug){
            MSG_INFO("HCands id = " << Higgs.pid() << ", HCands barcode = " << Higgs.genParticle()->barcode() << ", HCands status = " << Higgs.genParticle()->status());
        }
        // if(debug){
        //     for(const Particle &top_p: top.parents()){
        //         MSG_INFO ( "topCands parent id = " << top_p.pid() << ", topCands parent barcode = " << top_p.genParticle()->barcode() << ", topCands parent status = " << top_p.genParticle()->status());
        //     }
        // }

        //check top decay chain
        Particle Higgs_to_decay = getLastInstance(Higgs);
        for(const Particle &Higgs_c: Higgs_to_decay.children()){
         if (debug) {
            MSG_INFO ( "HCands child id = " << Higgs_c.pid() << ", HCands child barcode = " << Higgs_c.genParticle()->barcode() << ", HCands child status = " << Higgs_c.genParticle()->status());
          }
          if(fabs(Higgs_c.pid()) == 6)
           {
               topCandsFromH.push_back(Higgs_c);
           }
        }
        if (debug) MSG_INFO ( "\n");

      }


      if ( topCandsFromH.size() == 2 ){
        double pT_ttbar_FromH = (topCandsFromH.at(0).momentum() + topCandsFromH.at(1).momentum()).pT()/GeV;
        double dR_ttbar_FromH = fabs(deltaR(topCandsFromH.at(0),topCandsFromH.at(1)));

        double InvM_ttbar_FromH = ( topCandsFromH.at(0).momentum() + topCandsFromH.at(1).momentum() ).mass()/GeV;

        // std::cout << "pT_ttbar = " << pT_ttbar << std::endl;
        _h["Inclusive_pT_ttbar_FromH"]->fill(pT_ttbar_FromH);
        _h["Inclusive_mTop1_FromH"]->fill(topCandsFromH.at(0).mass());
        _h["Inclusive_mTop2_FromH"]->fill(topCandsFromH.at(1).mass());
        _h["Inclusive_mTop1_zoomin_FromH"]->fill(topCandsFromH.at(0).mass());
        _h["Inclusive_mTop2_zoomin_FromH"]->fill(topCandsFromH.at(1).mass());
        _h["Inclusive_dR_ttbar_FromH"]->fill(dR_ttbar_FromH);
        _h["Inclusive_InvM_FromH"]->fill(InvM_ttbar_FromH);
        _h["Inclusive_InvM_FromH400"]->fill(InvM_ttbar_FromH);
        _h["Inclusive_InvM_FromH500"]->fill(InvM_ttbar_FromH);
        _h["Inclusive_InvM_FromH600"]->fill(InvM_ttbar_FromH);
        _h["Inclusive_InvM_FromH700"]->fill(InvM_ttbar_FromH);
        _h["Inclusive_InvM_FromH800"]->fill(InvM_ttbar_FromH);
        _h["Inclusive_InvM_FromH900"]->fill(InvM_ttbar_FromH);
        _h["Inclusive_InvM_FromH1000"]->fill(InvM_ttbar_FromH);

      }

      if ( topCandsFromNonH.size() == 2 ){
        double pT_ttbar_FromNonH = (topCandsFromNonH.at(0).momentum() + topCandsFromNonH.at(1).momentum()).pT()/GeV;
        double dR_ttbar_FromNonH = fabs(deltaR(topCandsFromNonH.at(0),topCandsFromNonH.at(1)));

        // std::cout << "pT_ttbar = " << pT_ttbar << std::endl;
        _h["Inclusive_pT_ttbar_FromNonH"]->fill(pT_ttbar_FromNonH);
        _h["Inclusive_mTop1_FromNonH"]->fill(topCandsFromNonH.at(0).mass());
        _h["Inclusive_mTop2_FromNonH"]->fill(topCandsFromNonH.at(1).mass());
        _h["Inclusive_mTop1_zoomin_FromNonH"]->fill(topCandsFromNonH.at(0).mass());
        _h["Inclusive_mTop2_zoomin_FromNonH"]->fill(topCandsFromNonH.at(1).mass());
        _h["Inclusive_dR_ttbar_FromNonH"]->fill(dR_ttbar_FromNonH);

      }


      if (HCands.size()==1){
        _h["Inclusive_mH"]->fill(HCands.at(0).mass());
        _h["Inclusive_mH400"]->fill(HCands.at(0).mass());
        _h["Inclusive_mH500"]->fill(HCands.at(0).mass());
        _h["Inclusive_mH600"]->fill(HCands.at(0).mass());
        _h["Inclusive_mH700"]->fill(HCands.at(0).mass());
        _h["Inclusive_mH800"]->fill(HCands.at(0).mass());
        _h["Inclusive_mH900"]->fill(HCands.at(0).mass());
        _h["Inclusive_mH1000"]->fill(HCands.at(0).mass());

        _h["Inclusive_mH_zoomin"]->fill(HCands.at(0).mass());
        _h["Inclusive_pT_H"]->fill(HCands.at(0).pT());
      }

      _h["Inclusive_nTop_FromH"]->fill(topCandsFromH.size());


      // if ( wCandsFromTop.size() == 2 ){
      //   double dR_WW = fabs(deltaR(wCandsFromTop.at(0),wCandsFromTop.at(1)));
      //   _h["Inclusive_dR_WW"]->fill(dR_WW);
      // }

      // select ditau or dilepton

      // for(const Particle &wboson: wCandsFromTop){
      //
      //   Particle wboson_to_decay = getLastInstance(wboson);
      //
      //   for(const Particle &wboson_c: wboson_to_decay.children()){
      //       if(fabs(wboson_c.pid()) == 15)
      //        {
      //            tauCands.push_back(wboson_c);
      //        }
      //       else if((fabs(wboson_c.pid()) == 11) || fabs(wboson_c.pid()) == 13)
      //        {
      //            lepCands.push_back(wboson_c);
      //        }
      //   }
      // }

      // if (debug){
      //   MSG_INFO("tauCands size = " << tauCands.size() << "lepCands size = " << lepCands.size());
      //
      //   if ( tauCands.size() >= 1 ){
      //     MSG_INFO("tauCands.at(0).pT() = " << tauCands.at(0).pT()); // GeV unit
      //   }
      // }


      // if ( tauCands.size() == 2 ){
      //   _h["Inclusive_pT_tau1"]->fill(tauCands.at(0).pT());
      //   _h["Inclusive_pT_tau2"]->fill(tauCands.at(1).pT());
      //   _h["Inclusive_Eta_tau1"]->fill(tauCands.at(0).eta());
      //   _h["Inclusive_Eta_tau2"]->fill(tauCands.at(1).eta());
      //
      //   double dR_diTau = fabs(deltaR(tauCands.at(0),tauCands.at(1)));
      //   double dEta_diTau = fabs(deltaEta(tauCands.at(0),tauCands.at(1)));
      //   double dPhi_diTau = fabs(deltaPhi(tauCands.at(0),tauCands.at(1)));
      //
      //   _h["Inclusive_dR_diTau"]->fill(dR_diTau);
      //   _h["Inclusive_dEta_diTau"]->fill(dEta_diTau);
      //   _h["Inclusive_dPhi_diTau"]->fill(dPhi_diTau);
      //
      //   if ( tauCands.at(0).pT() > 15 && tauCands.at(1).pT() > 15 ){
      //     if ( fabs(tauCands.at(0).eta()) < 2.5 && fabs(tauCands.at(1).eta()) < 2.5 ){
      //
      //       double dR_diTau_Pt15 = fabs(deltaR(tauCands.at(0),tauCands.at(1)));
      //       double dEta_diTau_Pt15 = fabs(deltaEta(tauCands.at(0),tauCands.at(1)));
      //       double dPhi_diTau_Pt15 = fabs(deltaPhi(tauCands.at(0),tauCands.at(1)));
      //       _h["Inclusive_dR_diTau_Pt15"]->fill(dR_diTau_Pt15);
      //       _h["Inclusive_dEta_diTau_Pt15"]->fill(dEta_diTau_Pt15);
      //       _h["Inclusive_dPhi_diTau_Pt15"]->fill(dPhi_diTau_Pt15);
      //
      //     }
      //   }
      // }


      // if ( lepCands.size() == 2 ){
      //   _h["Inclusive_pT_lep1"]->fill(lepCands.at(0).pT());
      //   _h["Inclusive_pT_lep2"]->fill(lepCands.at(1).pT());
      //   _h["Inclusive_Eta_lep1"]->fill(lepCands.at(0).eta());
      //   _h["Inclusive_Eta_lep2"]->fill(lepCands.at(1).eta());
      //
      //   double dR_diLep = fabs(deltaR(lepCands.at(0),lepCands.at(1)));
      //   _h["Inclusive_dR_diLep"]->fill(dR_diLep);
      //
      //   if ( lepCands.at(0).pT() > 15 && lepCands.at(1).pT() > 15 ){
      //     if ( fabs(lepCands.at(0).eta()) < 2.5 && fabs(lepCands.at(1).eta()) < 2.5 ){
      //       double dR_diLep_Pt15 = fabs(deltaR(lepCands.at(0),lepCands.at(1)));
      //       _h["Inclusive_dR_diLep_Pt15"]->fill(dR_diLep_Pt15);
      //     }
      //   }
      // }


      if (debug){
            // for(const Particle &Wboson: wCandsFromNonTop){
            //   MSG_INFO("W id = " << Wboson.pid() << ", W barcode = " << Wboson.genParticle()->barcode() << ", W status = " << Wboson.genParticle()->status());
            //
            //   for(const Particle &part: Wboson.parents()){
            //     MSG_INFO ( "wCandsFromNonTop Wboson parent id = " << part.pid() << ", Wboson parent barcode = " << part.genParticle()->barcode() << ", Wboson parent status = " << part.genParticle()->status());
            //   }
            //   for(const Particle &part: Wboson.children()){
            //     MSG_INFO ( "wCandsFromNonTop Wboson child id = " << part.pid() << ", Wboson child barcode = " << part.genParticle()->barcode() << ", Wboson child status = " << part.genParticle()->status());
            //   }
            //
            //   MSG_INFO ( "\n");
            //
            // }
            //
            // for(const Particle &Wboson: wCandsFromTop){
            //     MSG_INFO("W id = " << Wboson.pid() << ", W barcode = " << Wboson.genParticle()->barcode() << ", W status = " << Wboson.genParticle()->status());
            //
            // for(const Particle &part: Wboson.parents()){
            //     MSG_INFO ( "wCandsFromTop Wboson parent id = " << part.pid() << ", Wboson parent barcode = " << part.genParticle()->barcode() << ", Wboson parent status = " << part.genParticle()->status());
            // }
            // for(const Particle &part: Wboson.children()){
            //     MSG_INFO ( "wCandsFromTop Wboson child id = " << part.pid() << ", Wboson child barcode = " << part.genParticle()->barcode() << ", Wboson child status = " << part.genParticle()->status());
            // }
            //
            //     MSG_INFO ( "\n");
            // }
      MSG_INFO( "Number of top quark : " << topCands.size());
      MSG_INFO( "Number of H boson:" << HCands.size());
      MSG_INFO( "Number of Top quark is from H:" << topCandsFromH.size());
      MSG_INFO( "Number of Top quark is from Non-H:" << topCandsFromNonH.size());
      MSG_INFO("#----------------DONE--------------#");
    }


      // foreach (const Particle &W, wCands)
      // {
      //     bool topToBWFound = false;
      //     foreach (const Particle &b, bCands)
      //     {
      //         GenVertex *bVert      = b.genParticle()->production_vertex();
      //         GenVertex *wVert      = W.genParticle()->production_vertex();
      //         for(GenVertex::particles_in_const_iterator bvIter = bVert->particles_in_const_begin(); bvIter != bVert->particles_in_const_end(); ++bvIter)
      //         {
      //             for(GenVertex::particles_in_const_iterator wvIter = wVert->particles_in_const_begin();wvIter != wVert->particles_in_const_end(); ++wvIter)
      //             {
      //                 GenParticle *bParent  = (*bvIter);
      //                 GenParticle *wParent  = (*wvIter);
      //
      //                 foreach (const Particle &topQ, topCands)
      //                 {
      //                     if(topQ.genParticle()->barcode() == bParent->barcode() && topQ.genParticle()->barcode() == wParent->barcode())
      //                     {
      //                         topToBWFound = true;
      //                         //You found top->b,w vertex
      //                         topQuarks.push_back(topQ);
      //                         bQuarks.push_back(b);
      //                     }
      //                 }
      //             }
      //         }
      //     }
      //     if(topToBWFound)
      //     {
      //         wBosons.push_back(W);
      //     }
      // }
  //
  //     Particles W_enu, W_munu,W_taunu;
  //     foreach (const Particle &W, wCands)
  //     {
  //         bool wlepDecay = false;
	//       bool whadDecay = false; //MARIA
  //         bool noTopParent = false;
  //         //GenVertex *wVert      = W.genParticle()->production_vertex();
  //         GenVertex *wEndVert   = W.genParticle()->end_vertex();
  //
  //         foreach (const Particle &topW, wBosons)
  //         {
  //             if(topW.genParticle() != W.genParticle())
  //             {
  //                 noTopParent = true;
  //             }
  //         }
  //
  //         /*for(GenVertex::particles_in_const_iterator wvIter = wVert->particles_in_const_begin(); wvIter != wVert->particles_in_const_end(); ++wvIter)
  //         {
  //             GenParticle *wParent = (*wvIter);
  //             if(abs(wParent->pdg_id())!=6 || abs(wParent->pdg_id())!= 24);
  //             {
  //                 noTopParent = true;
  //                 break;
  //             }
  //         }*/
  //
  //         for(GenVertex::particles_out_const_iterator witer = wEndVert->particles_out_const_begin(); witer != wEndVert->particles_out_const_end(); ++witer)
  //         {
  //             int wchPdg =  abs((*witer)->pdg_id());
  //             if (wchPdg ==11 || wchPdg ==12 || wchPdg ==13 || wchPdg ==14 || wchPdg ==15 || wchPdg ==16)
  //             {
  //                 wlepDecay = true;
  //
  //                 if (!noTopParent)
  //                 {
  //                     if(wchPdg ==11 || wchPdg ==12)
  //                     {
  //                         W_enu.push_back(*witer);
  //                     }
  //                     if(wchPdg ==13 || wchPdg ==14)
  //                     {
  //                         W_munu.push_back(*witer);
  //                     }
  //                     if(wchPdg ==15 || wchPdg ==16)
  //                     {
  //                         W_taunu.push_back(*witer);
  //                     }
  //                 }
  //             }
	//           else if(wchPdg ==1 || wchPdg ==2 || wchPdg ==3 || wchPdg ==4) { //MARIA (two entries per event)
  //                 whadDecay = true;
	// 	          LightQuarkfromWTop.push_back(*witer);//will match with particle-level jets
  //             }
  //         }
  //
	//       //associated W from ME
  //         if(wlepDecay && noTopParent)
  //         {
  //             wBosonsME.push_back(W);
  //         }
  //
  //     }//end Wcand loop
  //
  //     //This has a potential problem that we are not caring about dileptonic tt-bar system.
  //     MSG_INFO ("Found Electron or electron neutrino: " << W_enu.size() );
  //     if (W_enu.size() ==2) LepfromWTop = W_enu;
  //     if (W_munu.size() ==2 ) LepfromWTop = W_munu;
  //     if (W_taunu.size() ==2 ) LepfromWTop = W_taunu;
  //
  //
  //
  //     sortByPt(topQuarks);
  //     sortByPt(bQuarks);
  //     sortByPt(wBosons);
  //     sortByPt(wBosonsME);
  //     sortByPt(LepfromWTop);//MARIA
  //     sortByPt(LightQuarkfromWTop);//MARIA
  //
  //
  //
  //     const MissingMomentum &met= applyProjection<MissingMomentum>(event,"MissingET");
  //     const double event_met    = met.vectorEt().mod()*GeV;
  //
  //
  //     const FastJets & jetProj  = applyProjection<FastJets>(event,"Jets");
  //     const Jets alljets        = jetProj.jetsByPt(25*GeV);
  //     double event_ht           = 0;
  //     foreach(const Jet &j, alljets){ event_ht += j.pT()*GeV;}
  //
  //     _h_nJets->fill(alljets.size(),weight);
  //
  //     // dR (quarks & jets) -- bQuarks,
  //     Jets _lightjetsMatchedToLightQ;
  //     foreach (const Jet& _jet, alljets) {
	// if (fabs(_jet.eta()) < 2.5) {
	//   // jet B-hadron matching
	//   foreach (const Particle &quark, LightQuarkfromWTop) {
	//     FourMomentum _quarkfm = quark.momentum();
	//     double _quark_jet_dR = deltaR(_jet.momentum(), _quarkfm);
	//     if(_quark_jet_dR < 0.4) {
	//       _lightjetsMatchedToLightQ.push_back(_jet);
	//       _h_matchedjets_pT->fill(_jet.momentum().pT()/GeV, weight);
	//     }
	//     continue;
	//   }
	// }
  //     }
  //     _h_matchedjets_N->fill(_lightjetsMatchedToLightQ.size(),weight);
  //
  //     if(topQuarks.size()==2 && wBosons.size()==2 && bQuarks.size()==2 && wBosonsME.size()==1 )//&&  LepfromWTop.size() >=1)// && LightQuarkfromWTop.size() >=2)
  //     {
  //         MSG_INFO("WCANDS size: "<< wCands.size());
  //         MSG_INFO("WBOSONSME Size: "<<wBosonsME.size());
  //         MSG_INFO("WBOSONS Size: "<<wBosons.size());
  //         MSG_INFO("bQuarks Size: "<<wBosons.size());
  //         MSG_INFO("topQuarks Size: "<<topQuarks.size());
	//       MSG_INFO("LepfromWTop Size: "<<LepfromWTop.size());
	//       MSG_INFO("LightQuarkfromWTop Size: "<<LightQuarkfromWTop.size());
  //     }
  //
  //     if(topQuarks.size()==2 && wBosons.size()==2 && bQuarks.size()==2 && wBosonsME.size()==1 ) // && (LepfromWTop.size() >=2 && LightQuarkfromWTop.size() >=2) )
  //     {
  //         _h_top_dEta->fill(deltaEta(topQuarks[0],topQuarks[1]),weight);
  //         _h_top_dPhi->fill(deltaPhi(topQuarks[0],topQuarks[1]),weight);
  //         _h_top_dR->fill(deltaR(topQuarks[0],topQuarks[1]),weight);
  //
  //         _h_ttbar_pt->fill( (topQuarks[0].momentum() + topQuarks[1].momentum()).perp()*GeV,weight);
  //         _h_t1_mass->fill(topQuarks[0].momentum().mass()*GeV,weight);
  //         _h_t1_pt->fill(topQuarks[0].momentum().perp()*GeV,weight);
  //         _h_t1_eta->fill(topQuarks[0].momentum().eta(),weight);
  //         _h_t1_phi->fill(topQuarks[0].momentum().phi(),weight);
  //         _h_t2_mass->fill(topQuarks[1].momentum().mass()*GeV,weight);
  //         _h_t2_pt->fill(topQuarks[1].momentum().perp()*GeV,weight);
  //         _h_t2_eta->fill(topQuarks[1].momentum().eta(),weight);
  //         _h_t2_phi->fill(topQuarks[1].momentum().phi(),weight);
  //
  //
  //         _h_b_dEta->fill(deltaEta(bQuarks[0],bQuarks[1]),weight);
  //         _h_b_dPhi->fill(deltaPhi(bQuarks[0],bQuarks[1]),weight);
  //         _h_b_dR->fill(deltaR(bQuarks[0],bQuarks[1]),weight);
  //         _h_b1_pt->fill(bQuarks[0].momentum().perp()*GeV,weight);
  //         _h_b1_eta->fill(bQuarks[0].momentum().eta(),weight);
  //         _h_b1_phi->fill(bQuarks[0].momentum().phi(),weight);
  //         _h_b2_pt->fill(bQuarks[1].momentum().perp()*GeV,weight);
  //         _h_b2_eta->fill(bQuarks[1].momentum().eta(),weight);
  //         _h_b2_phi->fill(bQuarks[1].momentum().phi(),weight);
  //
  //         _h_w_dEta->fill(deltaEta(wBosons[0],wBosons[1]),weight);
  //         _h_w_dPhi->fill(deltaPhi(wBosons[0],wBosons[1]),weight);
  //         _h_w_dR->fill(deltaR(wBosons[0],wBosons[1]),weight);
  //         _h_w1_mass->fill(wBosons[0].momentum().mass()*GeV,weight);
  //         _h_w1_pt->fill(wBosons[0].momentum().perp()*GeV,weight);
  //         _h_w1_eta->fill(wBosons[0].momentum().eta(),weight);
  //         _h_w1_phi->fill(wBosons[0].momentum().phi(),weight);
  //         _h_w2_mass->fill(wBosons[1].momentum().mass()*GeV,weight);
  //         _h_w2_pt->fill(wBosons[1].momentum().perp()*GeV,weight);
  //         _h_w2_eta->fill(wBosons[1].momentum().eta(),weight);
  //         _h_w2_phi->fill(wBosons[1].momentum().phi(),weight);
  //
  //
  //
  //         if(wBosonsME[0].genParticle()->pdg_id()==24)
  //         {
  //             _h_wPlusME_pt->fill(wBosonsME[0].momentum().perp()*GeV,weight);
  //             _h_wPlusME_eta->fill(wBosonsME[0].momentum().eta(),weight);
  //             _h_wPlusME_phi->fill(wBosonsME[0].momentum().phi(),weight);
  //
  //             _h_wPlus_nJets->fill(alljets.size(),weight);
  //             _h_wPlus_MET->fill(event_met,weight);
  //             _h_wPlus_HT->fill(event_ht,weight);
  //         }
  //         else if (wBosonsME[0].genParticle()->pdg_id()==-24)
  //         {
  //             _h_wMinusME_pt->fill(wBosonsME[0].momentum().perp()*GeV,weight);
  //             _h_wMinusME_eta->fill(wBosonsME[0].momentum().eta(),weight);
  //             _h_wMinusME_phi->fill(wBosonsME[0].momentum().phi(),weight);
  //
  //             _h_wMinus_nJets->fill(alljets.size(),weight);
  //             _h_wMinus_MET->fill(event_met,weight);
  //             _h_wMinus_HT->fill(event_ht,weight);
  //         }
  //     }
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
        //
        // For Powheg
        //
        // MSG_INFO("CROSS SSECTION:"<<crossSection());
        // auto xsec = isnan(crossSection()) ? 1: crossSection();
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

    map<std::string,Histo1DPtr> _h;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(tttt_parton);

}
