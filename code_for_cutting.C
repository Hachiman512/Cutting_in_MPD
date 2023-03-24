//C++, ROOT and MiniDst Headers
#include <vector>
#include <cmath>
#include <tgmath.h>
#include "Rtypes.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TPad.h"
#include <TMath.h>

R__ADD_INCLUDE_PATH(/home/anil/PhD/JINR/MPD_Calculations/MiniDst/)  //include the path for the mpdloadlibs.C file
#include "MpdMiniDstReader.h"
#include "MpdMiniDst.h"
#include "MpdMiniEvent.h"
#include "MpdMiniTrack.h"
#include "MpdMiniMcTrack.h"
#include "MpdMiniHelix.h"
#include "MpdMiniPhysicalHelix.h" 
//MPD useful packages
#include "mpdloadlibs.C"

 TH1F *hCHtrack 	      = new TH1F("hCHtrack","Charged tracks Distribution (Reco)",150,0,500); // Charged particles distribution histogram 
 TH1F *hRecTrack_hits   = new TH1F("hRecTrack_hits","Plot for Reconstructed Track_hits",60,0,70); //Reconstructed hits
 TH2D *hRecTrack_VTXxy  = new TH2D("hRecTrack_VTXxy"," 2D Plot for reconstructed trackhits and Primary vertex X and Y (Reco)",100,-0.1,0.1, 100, -0.1, 0.1); 
 TH1F *hRecTrack_VTXz   = new TH1F("hRecTrack_VTXz"," Plot for reconstructed trackhits and Primary vertex Z (Reco)",1000,-160,160); 
 TH1F *hRecTrack_gDCA   = new TH1F("hRecTrack_gDCA"," Plot for reconstructed trackhits and DCA (Reco)",100,-0.1,3); 
 TH2D *hdEdxP           = new TH2D("hdEdxP","Plot for dE/dx and total momentum(P)",100,-2.5,2.5,100,0,10);
 TH1F *hdEdx            = new TH1F("hdEdx", "Plot of energy loss for all the particles",100, 0, 10);
 TH1F *hRecTrack_pT     = new TH1F("hRecTrack_pT","Charged tracks Transverse momentum(pT) Distribution (Reco)",3200,0.01,2.5);
 TH1F *hRecTrack_P 	    = new TH1F("hRecTrack_P","Charged tracks total momentum(P) Distribution (Reco)",3200,0.001,2.5); 
 TH2D *hHitsP           = new TH2D("hHitsP","Plot for Reconstructed hits and total Momentum(P)",500,-0.01,2.5, 50, 0, 60 );   
 TH1F *hRecTrack_eta    = new TH1F("hRecTrack_eta","Charged tracks Rapidity(y) distribution (Reco)",3200,-3,3); 
 TH1F *hpseuRapid 	    = new TH1F("hpseuRapid"," Charged tracks Pseudorapidity(eta) distribution (Reco)",3200,-10,10); 
 TH1F *hRecTrack_phi 	  = new TH1F("hRecTrack_phi"," Charged tracks polar angle(phi) distribution (Reco)",3200,-4,4);  
 TH2D *hetapt           = new TH2D("hetapt","plot of rapidity(y) and transverse momentum(pT)", 500,-2.5,2.5,60,0,2.5);
 TH1F *hetapt1F         = new TH1F("hetapt1F","1D-plot of rapidity(y) and transverse momentum(pT)",500,-2.5,2.5);   
 TH2D *hphiVTXz         = new TH2D("hphiVTXz","Plot for vertexZ and polar angle(phi)",500,-160,160,500,-4,4);   

 //Plotting with cuts
 TH1F *hRecTrack_hitsC   = new TH1F("hRecTrack_hitsC","Plot for Reconstructed Track_hits with cuts",60,0,70); //Reconstructed hits
 TH2D *hRecTrack_VTXxyC  = new TH2D("hRecTrack_VTXxyC"," 2D Plot for reconstructed trackhits and Primary vertex X and Y (Reco) with cuts",100,-0.1,0.1, 100, -0.1, 0.1);
 TH1F *hRecTrack_VTXzC   = new TH1F("hRecTrack_VTXzC"," Plot for reconstructed trackhits and Primary vertex Z (Reco) with cuts",1000,-160,160); 
 TH1F *hRecTrack_PC 	    = new TH1F("hRecTrack_PC","Charged tracks total momentum(P) Distribution (Reco) with cuts",3200,0,2.5); 
 TH2D *hdEdxPC           = new TH2D("hdEdxPC","Plot for dE/dx and total momentum(P) with cuts",100,-2.5,2.5,100,0,10);
 TH1F *hRecTrack_pTC     = new TH1F("hRecTrack_pTC","Charged tracks Transverse momentum(pT) Distribution (Reco) with cuts",3200,0.01,2.5);
 TH1F *hRecTrack_gpTC    = new TH1F("hRecTrack_gpTC","Charged global tracks total momentum(GP) Distribution (Reco) with cuts",3200,0,2.5); 
 TH2D *hHitsPC           = new TH2D("hHitsPC","Plot for Reconstructed hits and total Momentum(P) with cuts",500,-0.01,2.5, 50, 0, 60 );   
 TH1F *hRecTrack_etaC    = new TH1F("hRecTrack_etaC","Charged tracks Rapidity(y) distribution (Reco) with cuts",3200,-0.15,0.15); 
 TH1F *hpseuRapidC 	    = new TH1F("hpseuRapidC"," Charged tracks Pseudorapidity(eta) distribution (Reco) with cuts",3200,-10,10); 
 TH2D *hetaptC           = new TH2D("hetaptC","plot of rapidity(y) and transverse momentum(pT) with cuts", 500,-0.15,0.15,60,0,2);
 TH1F *hRecTrack_gDCAC   = new TH1F("hRecTrack_gDCAC"," Plot for reconstructed trackhits and DCA (Reco) with cuts",100,-0.1,3); 
 
 //_______________________________________________________________________________________
 void trial(const Char_t* inFileName, int myevents, TString outFileName ) {

  //Total time for running this code
  TStopwatch timer;
  timer.Start();  
  
  // Begin new class to read inFileName
  MpdMiniDstReader* miniDstReader = new MpdMiniDstReader(inFileName);
  TFile *fo = new TFile(outFileName.Data(),"RECREATE");  

  //Initialize MiniDst Files
  miniDstReader->Init();

  miniDstReader->SetStatus("*", 0);  
  miniDstReader->SetStatus("Event*", 1);          
  miniDstReader->SetStatus("Track*", 1);
  miniDstReader->SetStatus("McEvent*", 1);
  miniDstReader->SetStatus("McTrack*", 1);
  
  // For mpd reconstruced tracks
  
  Float_t nRecoTracks, mpdTrack_px,mpdTrack_gpx, mpdTrack_py,mpdTrack_gpy,mpdTrack_gpz, mpdTrack_pz, mpdTrack_P,mpdTrack_gP, mpdTrack_Eta, mpdTrack_AbsEta, mpdTrack_pT,mpdTrack_gpT,mpdTrack_phi, mpdTrack_Pneg;
  Float_t mpdTrack_originX,mpdTrack_originY, mpdTrack_originZ,mpdTrack_torigin, mpdTrack_gDCAX,mpdTrack_gDCAY, mpdTrack_gDCAZ, mpdTrack_tDCA, pseuRapid, ratio;
  Int_t  mpdTrack_charge, mpdtrack_pdg;

  //==================================================//
  //            Begin loop of events                    //
  //==================================================//

  //Get number of events of file
  Long64_t events2read = miniDstReader->chain()->GetEntries();

  for (Long64_t i = 0; i < myevents; i++) {
    Bool_t  isOk = miniDstReader->readMiniEvent( i );

    // Retrieve current miniDst (from the given .MiniDst.root file)
    MpdMiniDst *dst = miniDstReader->miniDst();

    // Get MiniEvent information
    MpdMiniEvent *event = dst->event(); //For reconstructed events
    MpdMiniMcEvent *mcEvent = dst->mcEvent(); // For the McEvents
    
    Long_t n_tracks_mpd     = 0;
    //print number of events
    cout << "No. of event: " << i << endl;

    
    //==================================================//
    //		Begin loop of RECO TRACKS		//
    //==================================================//

    // Get number of reconstructed tracks
    Int_t nRecoTracks = dst->numberOfTracks();
    Float_t primaryTrackCounter = 0;

    //print number of reconstructed tracks in this event
    cout << "Total number of reco tracks: " << nRecoTracks << endl;
    
    for (Int_t j = 0; j < nRecoTracks; j++) {
      
      // Get jth reconstructed track
      MpdMiniTrack *miniTrack = dst->track(j);
			MpdMiniMcTrack *mcTrack = dst->mcTrack(j);
      
      Bool_t hasmc 		= miniTrack->hasMcTrack();
      Float_t mcIndex  	= miniTrack->mcTrackIndex();
      Bool_t isprimarytrack 	= miniTrack->isPrimary();
      
      if (mcIndex < 1) continue;
      if (hasmc == false) continue;

      Int_t nHits       = miniTrack->nHits();
      mpdTrack_charge	= miniTrack->charge();
      Float_t dEdx = miniTrack->dEdx();
    
      TVector3 mpdTrack_ptot 	= miniTrack->pMom(); // momentum for primary tracks 
      mpdTrack_px	= mpdTrack_ptot.X();
      mpdTrack_py	= mpdTrack_ptot.Y();
      mpdTrack_pz	= mpdTrack_ptot.Z();
      mpdTrack_P 	= miniTrack->pMom().Mag();  //Get P magnitude
      mpdTrack_Pneg = -1.0*mpdTrack_P;
      mpdTrack_pT	= miniTrack->pPt(); //          (p for primaries, g for global tracks) 
    
      TVector3 mpdTrack_gptot 	= miniTrack->gMom(); /// Return momentum (GeV/c) of the global tracks
      mpdTrack_gpx	= mpdTrack_gptot.X();
      mpdTrack_gpy	= mpdTrack_gptot.Y();
      mpdTrack_gpz	= mpdTrack_gptot.Z();    
      mpdTrack_gP 	= miniTrack->gMom().Mag();  //Get momentum magnitude of global tracks
      mpdTrack_gpT	= miniTrack->gPt(); //          (p for primaries, g for global tracks) 
    
      mpdTrack_Eta 	= 0.5*TMath::Log((mpdTrack_P+mpdTrack_pz)/(mpdTrack_P-mpdTrack_pz));
      mpdTrack_AbsEta 	= TMath::Abs(mpdTrack_Eta);
      mpdTrack_phi 	= TMath::ATan2(mpdTrack_py,mpdTrack_px);     
      pseuRapid = -1.0*TMath::Log(tan(mpdTrack_phi/2));

      // Get primary vertex z-position

      TVector3 primaryVertex  = event->primaryVertex();
      Float_t vertexX = event->primaryVertex().X();
      Float_t vertexY = event->primaryVertex().Y();
      Float_t vertexZ = event->primaryVertex().Z();
    
      TVector3 mpdTrack_origin  = miniTrack->origin();
      mpdTrack_originX	= mpdTrack_origin.X();
      mpdTrack_originY	= mpdTrack_origin.Y();
      mpdTrack_originZ	= mpdTrack_origin.Z();    
      mpdTrack_torigin 	= miniTrack->origin().Mag();  //Get momentum magnitude of global tracks
      
      mpdTrack_gDCAX	= mpdTrack_originX - vertexX;
      mpdTrack_gDCAY	= mpdTrack_originY - vertexY;
      mpdTrack_gDCAZ	= mpdTrack_originZ - vertexZ;
      mpdTrack_tDCA 	= TMath::Sqrt(mpdTrack_gDCAX*mpdTrack_gDCAX + mpdTrack_gDCAY*mpdTrack_gDCAY + mpdTrack_gDCAZ*mpdTrack_gDCAZ);  //Get momentum magnitude of global tracks
      
      Int_t refMultPos = event->refMultPos();
      Int_t refMultNeg = event->refMultNeg();
      Int_t refMult = refMultPos + refMultNeg;

      Int_t refMult2PosEast = event->refMult2PosEast();
      Int_t refMult2NegEast = event->refMult2NegEast();
      Int_t refMult2East = refMult2PosEast + refMult2NegEast;
      Int_t refMult2PosWest = event->refMult2PosWest();
      Int_t refMult2NegWest = event->refMult2NegWest();
      Int_t refMult2West = refMult2PosWest + refMult2NegWest;
      Int_t refMult2 = refMult2West + refMult2East;

      Int_t refMultHalfPosEast = event->refMultHalfPosEast();
      Int_t refMultHalfNegEast = event->refMultHalfNegEast();
      Int_t refMultHalfEast = refMultHalfPosEast + refMultHalfNegEast;
      Int_t refMultHalfPosWest = event->refMultHalfPosWest();
      Int_t refMultHalfNegWest = event->refMultHalfNegWest();
      Int_t refMultHalfWest = refMultHalfPosWest + refMultHalfNegWest;
      Int_t refMult3 = refMultHalfWest + refMultHalfEast;
      Int_t RefMultiTot = refMult;

      hphiVTXz->Fill(vertexZ, mpdTrack_phi);
      hRecTrack_VTXz->Fill(vertexZ); 
      hRecTrack_VTXxy->Fill(vertexX, vertexY);
      hRecTrack_hits->Fill(nHits);


      hRecTrack_eta->Fill(mpdTrack_Eta);
      hetapt->Fill(mpdTrack_Eta,mpdTrack_pT);
      hetapt1F->Fill(mpdTrack_Eta,mpdTrack_pT);
      hRecTrack_pT->Fill(mpdTrack_pT);
      hRecTrack_P->Fill(mpdTrack_P);
      hHitsP->Fill(mpdTrack_P,nHits); 
      hpseuRapid->Fill(pseuRapid);
      hRecTrack_phi->Fill(mpdTrack_phi);
      hdEdx->Fill(dEdx);
      hdEdxP->Fill(mpdTrack_P,dEdx);
      hdEdxP->Fill(mpdTrack_Pneg,dEdx);
      hRecTrack_gDCA->Fill(mpdTrack_tDCA);      
              
      // Filling histograms for charged particles and with constraints like nhits;
      if (nHits<10) continue;
      if (isprimarytrack != 1) continue;
      if (mpdTrack_charge == 0) continue; 
      if( TMath::Abs(mpdTrack_Eta)>0.1) continue; 
      if( TMath::Abs(vertexZ)>50) continue; 
      if( mpdTrack_pT<0.2) continue; 
      if(TMath::Abs(mpdTrack_tDCA)>3) continue;
      if(mpdTrack_P<0.1) continue;
    
      hRecTrack_hitsC->Fill(nHits);
      hRecTrack_etaC->Fill(mpdTrack_Eta);
      hetaptC->Fill(mpdTrack_Eta,mpdTrack_pT);
      hRecTrack_pTC->Fill(mpdTrack_pT);
      hRecTrack_PC->Fill(mpdTrack_P);
      hHitsPC->Fill(mpdTrack_P,nHits); 
      hpseuRapidC->Fill(pseuRapid);
      hdEdxPC->Fill(mpdTrack_P,dEdx);
      hdEdxPC->Fill(mpdTrack_Pneg,dEdx);
      hRecTrack_gDCAC->Fill(mpdTrack_tDCA); 
      hRecTrack_VTXzC->Fill(vertexZ);    
      hRecTrack_VTXxyC->Fill(vertexX, vertexY);  
     
      n_tracks_mpd++;
     
    } //end of reco tracks loop
  
  hCHtrack->Fill(n_tracks_mpd);  
  
  } // end jentry loops


  fo->cd();
 
  hCHtrack->GetXaxis()->SetTitle("Channel Track hits");
  hCHtrack->GetXaxis()->CenterTitle(true);
  hCHtrack->Write();

  hRecTrack_hits->GetXaxis()->SetTitle("Number of Reconstructed hits");
  hRecTrack_hits->GetXaxis()->CenterTitle(true);
  hRecTrack_hits->GetYaxis()->SetTitle("Events");
  hRecTrack_hits->GetYaxis()->CenterTitle(true);
  hRecTrack_hits->Write();

  hRecTrack_VTXxy->GetXaxis()->SetTitle("Primary Vertex X(cm)");
  hRecTrack_VTXxy->GetXaxis()->CenterTitle(true); 
  hRecTrack_VTXxy->GetYaxis()->SetTitle("Primary Vertex Y(cm)");
  hRecTrack_VTXxy->GetYaxis()->CenterTitle(true); 
  hRecTrack_VTXxy->Write();

  hRecTrack_VTXz->GetXaxis()->SetTitle("Primary Vertex Z(cm)");
  hRecTrack_VTXz->GetXaxis()->CenterTitle(true);
  hRecTrack_VTXz->GetYaxis()->SetTitle("Events");
  hRecTrack_VTXz->GetYaxis()->CenterTitle(true);
  hRecTrack_VTXz->Write();

  hRecTrack_gDCA->GetXaxis()->SetTitle("Distance of closest approach (DCA)in cm");
  hRecTrack_gDCA->GetXaxis()->CenterTitle(true);
  hRecTrack_gDCA->GetYaxis()->SetTitle("Events");
  hRecTrack_gDCA->GetYaxis()->CenterTitle(true);
  hRecTrack_gDCA->Write();
  
  hdEdx->GetXaxis()->SetTitle("Energy loss (dE/dx)");
  hdEdx->GetXaxis()->CenterTitle(true);
  hdEdx->GetYaxis()->SetTitle("Events");
  hdEdx->GetYaxis()->CenterTitle(true);
  hdEdx->Write();

  hdEdxP->GetXaxis()->SetTitle("Total momentum (P)");
  hdEdxP->GetXaxis()->CenterTitle(true);
  hdEdxP->GetYaxis()->SetTitle("Energy loss (dE/dx)");
  hdEdxP->GetYaxis()->CenterTitle(true);
  hdEdxP->Write();
  
  
  hRecTrack_pT->GetXaxis()->SetTitle("Transverse momentum pT with primary hits");
  hRecTrack_pT->GetXaxis()->CenterTitle(true);
  hRecTrack_pT->GetYaxis()->SetTitle("Events");
  hRecTrack_pT->GetYaxis()->CenterTitle(true);
  hRecTrack_pT->Write();

  hRecTrack_P->GetXaxis()->SetTitle("Total momentum P with primary hits");
  hRecTrack_P->GetXaxis()->CenterTitle(true);
  hRecTrack_P->GetYaxis()->SetTitle("Events");
  hRecTrack_P->GetYaxis()->CenterTitle(true);
  hRecTrack_P->Write();


  
  hHitsP->GetXaxis()->SetTitle("Total momentum P");
  hHitsP->GetXaxis()->CenterTitle(true);
  hHitsP->GetYaxis()->SetTitle("Number of Hits");
  hHitsP->GetYaxis()->CenterTitle(true);
  hHitsP->Write();

  hRecTrack_eta->GetXaxis()->SetTitle("Rapidity(y)");
  hRecTrack_eta->GetXaxis()->CenterTitle(true);
  hRecTrack_eta->GetYaxis()->SetTitle("Events");
  hRecTrack_eta->GetYaxis()->CenterTitle(true);
  hRecTrack_eta->Write(); 

  hRecTrack_phi->GetXaxis()->SetTitle("Value of polar angle phi");
  hRecTrack_phi->GetXaxis()->CenterTitle(true);
  hRecTrack_phi->GetYaxis()->SetTitle("Events");
  hRecTrack_phi->GetYaxis()->CenterTitle(true);
  hRecTrack_phi->Write();

  hpseuRapid->GetXaxis()->SetTitle("Pseudorapidity(eta)");
  hpseuRapid->GetXaxis()->CenterTitle(true);  
  hpseuRapid->GetYaxis()->SetTitle("Events");
  hpseuRapid->GetYaxis()->CenterTitle(true);  
  hpseuRapid->Write(); 

  hetapt->GetXaxis()->SetTitle("Rapidity(y)");
  hetapt->GetXaxis()->CenterTitle(true);   
  hetapt->GetYaxis()->SetTitle("Transverse momentum(P)");
  hetapt->GetYaxis()->CenterTitle(true);   
  hetapt->Write();
 
  
  hphiVTXz->GetXaxis()->SetTitle("Primary Vertex Z(cm)");
  hphiVTXz->GetXaxis()->CenterTitle(true);  
  hphiVTXz->GetYaxis()->SetTitle("Value of polar angle phi ");
  hphiVTXz->GetYaxis()->CenterTitle(true);  
  hphiVTXz->Write();
 
//Write for cuts

  hRecTrack_hitsC->GetXaxis()->SetTitle("Number of Reconstructed hits");
  hRecTrack_hitsC->GetXaxis()->CenterTitle(true);
  hRecTrack_hitsC->GetYaxis()->SetTitle("Events");
  hRecTrack_hitsC->GetYaxis()->CenterTitle(true);
  hRecTrack_hitsC->Write();


  hRecTrack_VTXzC->GetXaxis()->SetTitle("Primary Vertex Z(cm)");
  hRecTrack_VTXzC->GetXaxis()->CenterTitle(true);
  hRecTrack_VTXzC->GetYaxis()->SetTitle("Events");
  hRecTrack_VTXzC->GetYaxis()->CenterTitle(true);
  hRecTrack_VTXzC->Write();

  hRecTrack_gDCAC->GetXaxis()->SetTitle("Distance of closest approach (DCA)in cm");
  hRecTrack_gDCAC->GetXaxis()->CenterTitle(true);
  hRecTrack_gDCAC->GetYaxis()->SetTitle("Events");
  hRecTrack_gDCAC->GetYaxis()->CenterTitle(true);
  hRecTrack_gDCAC->Write();

  hdEdxPC->GetXaxis()->SetTitle("Total momentum (P)");
  hdEdxPC->GetXaxis()->CenterTitle(true);
  hdEdxPC->GetYaxis()->SetTitle("Energy loss (dE/dx)");
  hdEdxPC->GetYaxis()->CenterTitle(true);
  hdEdxPC->Write();
  
  
  hRecTrack_pTC->GetXaxis()->SetTitle("Transverse momentum P with primary hits");
  hRecTrack_pTC->GetXaxis()->CenterTitle(true);
  hRecTrack_pTC->GetYaxis()->SetTitle("Events");
  hRecTrack_pTC->GetYaxis()->CenterTitle(true);
  hRecTrack_pTC->Write();

  hRecTrack_PC->GetXaxis()->SetTitle("Total momentum P with primary hits");
  hRecTrack_PC->GetXaxis()->CenterTitle(true);
  hRecTrack_PC->GetYaxis()->SetTitle("Events");
  hRecTrack_PC->GetYaxis()->CenterTitle(true);
  hRecTrack_PC->Write();
  
  hHitsPC->GetXaxis()->SetTitle("Total momentum P");
  hHitsPC->GetXaxis()->CenterTitle(true);
  hHitsPC->GetYaxis()->SetTitle("Number of Hits");
  hHitsPC->GetYaxis()->CenterTitle(true);
  hHitsPC->Write();

  hRecTrack_etaC->GetXaxis()->SetTitle("Rapidity(y)");
  hRecTrack_etaC->GetXaxis()->CenterTitle(true);
  hRecTrack_etaC->GetYaxis()->SetTitle("Events");
  hRecTrack_etaC->GetYaxis()->CenterTitle(true);
  hRecTrack_etaC->Write(); 

  hpseuRapidC->GetXaxis()->SetTitle("Pseudorapidity(eta)");
  hpseuRapidC->GetXaxis()->CenterTitle(true);  
  hpseuRapidC->GetYaxis()->SetTitle("Events");
  hpseuRapidC->GetYaxis()->CenterTitle(true);  
  hpseuRapidC->Write(); 

  hetaptC->GetXaxis()->SetTitle("Rapidity(y)");
  hetaptC->GetXaxis()->CenterTitle(true);   
  hetaptC->GetYaxis()->SetTitle("Transverse momentum(P)");
  hetaptC->GetYaxis()->CenterTitle(true);   
  hetaptC->Write();

  hRecTrack_VTXxyC->GetXaxis()->SetTitle("Primary Vertex X(cm)");
  hRecTrack_VTXxyC->GetXaxis()->CenterTitle(true); 
  hRecTrack_VTXxyC->GetYaxis()->SetTitle("Primary Vertex Y(cm)");
  hRecTrack_VTXxyC->GetYaxis()->CenterTitle(true); 
  hRecTrack_VTXxyC->Write();

 
 


  fo->Close();          
 //+++++++++++++++++++++++++++++++++++++++++++++
 timer.Print();

  cout << " Test passed" << endl;
  cout << " All ok " << endl;
  exit(0);

  // Finalize miniDst reader
}

