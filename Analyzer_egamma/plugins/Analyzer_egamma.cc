// -*- C++ -*-
//
// Package:    Demo_e_gamma/Analyzer_egamma
// Class:      Analyzer_egamma
//
/**\class Analyzer_egamma Analyzer_egamma.cc Demo_e_gamma/Analyzer_egamma/plugins/Analyzer_egamma.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Soumya Sarkar
//         Created:  Mon, 13 Mar 2023 18:15:26 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <TLorentzVector.h>

//Utilities Headers
#include <memory>
#include <string>
#include <map>
#include <iostream>
#include <cmath>
#include <fstream>
#include <TClonesArray.h>
#include <TObject.h>
#include <TChain.h>
#include <TROOT.h>
#include <vector>
#include "TLorentzVector.h"
//#include "TH1.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH2.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//Gsf electrons Include Files 
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

//Photons include files
#include "DataFormats/EgammaCandidates/interface/Photon.h"

//Conversions Include file
#include "DataFormats/EgammaCandidates/interface/Conversion.h"

//DeDx hit info include file

#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"

//Track include files

#include "DataFormats/TrackReco/interface/Track.h"

//Calocluster and supercluster include files
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

//GsfTrack include files
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

//PatCandidates/Electron.h
#include "DataFormats/PatCandidates/interface/Electron.h"

//Histogram Headers
#include "TH1.h"
#include "TH2.h"
#include "TGraph2D.h"

//Vector Headers
#include<vector>
#include "TVector3.h"

//EcalRecHit class

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

// class declaration
//

// This will improve performance in multithreaded jobs.

//Structure for EGamma

struct identity
{
  int conv1;
  int conv2;
};



struct Lepton
{ 
  //Generic Variables
  TLorentzVector v;
  TVector3 v1;
  int id; int ind; int momid;
  float wt;

  //GSF Electron.h variables 
  float scPixCharge;
  bool isGsfCtfScPixChargeConsistent;
  bool isGsfScPixChargeConsistent;
  float shFracInnerHits;
  float eSuperClusterOverP;
  float eSeedClusterOverP;
  float eSeedClusterOverPout;
  float eEleClusterOverPout;
  float deltaEtaSuperClusterTrackAtVtx;
  float deltaEtaSeedClusterTrackAtCalo;
  float deltaPhiSuperClusterTrackAtVtx;
  float deltaPhiSeedClusterTrackAtCalo;
  float sigmaEtaEta;
  float sigmaIetaIeta;
  float sigmaIphiIphi;
  float e1x5;
  float r9;
  float hcalOverEcal; // int dept has to be 0,1,2 
  bool  hcalOverEcalValid;
  int convFlags;
  float convDist;
  float convDcot;
  float convRadius;
  float convVtxFitProb;
  float ecalPFClusterIso;
  float hcalPFClusterIso;
  float trackFbrem;
  float superClusterFbrem;
  int numberOfBrems;
  float fbrem_eta_8; //fabs(eta)<0.8
  float fbrem_eta_144; //0.8<fabs(eta)<1.44
  float fbrem_eta_2; //1.57< fabs(eta)<2
  float fbrem_eta_2_; //eta>2
  float ecalEnergy;
  //Conversion.h variables

  // reco::Vertex& conversionVertex = NULL;
  // reco::CaloClusterPtrVector caloCluster;
  //math::XYZVectorF pairMomentum;
  //math::XYZTLorentzVectorF refittedPair4Momentum;
  //std::vector<reco::CaloClusterPtr> &bcMatchingWithTracks= NULL;
  double EoverP;
  double EoverPrefittedTracks;
  double pairInvariantMass;
  double pairCotThetaSeparation;
  double distOfMinimumApproach;
  double dPhiTracksAtVtx;
  double dPhiTracksAtEcal;
  double dEtaTracksAtEcal;
  // std::vector<math::XYZPointF>& ecalImpactPosition;
  math::XYZVectorF pairMomentum;
  uint8_t nSharedHits;
  // std::vector<uint8_t>& nHitsBeforeVtx;
  //std::vector<math::XYZVectorF>& tracksPout;
  //std::vector<math::XYZVectorF>& tracksPin;
  int status;
  bool isConverted;
  int nTracks;

  //DeDxHitInfo.h variables
  
  float charge;
  float pathlength;
  
  //TrackBase.h variables
  bool isTimeOk;
  double chi2;
  double ndof;
  double normalizedChi2;
  // int charge;
  double qoverp;
  double pt;
  double p;
  double px;
  double py;
  double pz;
  double phi;
  double eta;
  double t0;
  double beta;
  double etaError;
  double phiError;    
  
  //Gsftrack variables     
  int chargeMode;
  double qoverpMode;
  double thetaMode;
  double lambdaMode;
  double pMode;
  double ptMode;
  double pxMode;
  double pyMode;
  double pzMode;
  double phiMode;
  double etaMode;
  double thetaModeError;
  double lambdaModeError;
  double etaModeError;
  double phiModeError;
  //double etaError;
  //double phiError;
 
  //SuperCluster variables
  
  size_t size;
  double energy;
  double correctedEnergy;



  //Random variables
  int n_electrons,n_jet_electrons;
};

using namespace std;
using namespace edm;

std::vector<Lepton> elec;
std::vector<Lepton> gsf_elec;
std::vector<Lepton> conversion;
std::vector<Lepton> DeDx;
std::vector<Lepton> goodConversion;
std::vector<Lepton> goodGsf_elec;
std::vector<Lepton> goodIsoGsf_elec;
std::vector<Lepton> badGsf_elec;
std::vector<Lepton> track_general;
std::vector<Lepton> track_displaced;
std::vector<Lepton> GsfTrack;
std::vector<Lepton> T1_conv;
std::vector<Lepton> T2_conv;
std::vector<Lepton> conv_track;
std::vector<Lepton> sc;
std::vector<Lepton> CaloCluster;
std::vector<Lepton> CaloCluster_EB;
std::vector<Lepton> CaloCluster_EE;
std::vector<identity> convtype;
std::vector<int> nevent_goodIso;
using reco::TrackCollection;

class Analyzer_egamma : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Analyzer_egamma(const edm::ParameterSet&);
  ~Analyzer_egamma();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void sort(int opt);
  float delta_phi(float phi1,float phi2);
  float delta_phi_sign(float phi1,float phi2);
  float deltaR(float eta1,float phi1,float eta2,float phi2);
   
private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  
  int nevent;
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  TH1D *demohisto;
  //------------e-gamma members-----------------------
  edm::EDGetTokenT<vector<pat::Electron> > etoken_;  
  TH1D *eg_histo[50];
  edm::EDGetTokenT<vector<reco::GsfElectron> > gsf_e_token1_;
  edm::EDGetTokenT<vector<reco::GsfElectron> > gsf_e_token2_;
  TH1D *gsf_e_histo[100];
  edm::EDGetTokenT<vector<reco::Conversion> > conversion_token_;
  TH1D *conversion_histo[50];
  TH2D *merged_histo[50];
  edm::EDGetTokenT<vector<reco::DeDxHitInfo> > DeDx_token_;
  TH1D *DeDx_histo[50];
  edm::EDGetTokenT<vector<reco::Track> > track_token1_;
  edm::EDGetTokenT<vector<reco::Track> > track_token2_;
  TH1D *track_histo[50];
  edm::EDGetTokenT<vector<reco::GsfTrack> > GsfTrack_token_;
  TH1D *GsfTrack_histo[50];
  edm::EDGetTokenT<vector<reco::SuperCluster> > sc_token_;
  edm::EDGetTokenT<vector<reco::CaloCluster> > calo_token_;
  edm::EDGetTokenT<vector<reco::CaloCluster> > calo_eb_token_;
  edm::EDGetTokenT<vector<reco::CaloCluster> > calo_ee_token_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> > ecalrechit_;
  edm::EDGetTokenT<vector<reco::Photon> > photon_token_;
  TH1D *photon_histo[50];
  TH1D *calo_histo[50];
  TGraph2D *gr_good,*gr_bad,*gr_gsf,*gr_good_ecal,*gr_bad_ecal;
  TH3D *mustache_histo[10];
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Analyzer_egamma::Analyzer_egamma(const edm::ParameterSet& iConfig) :
  // tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
  // etoken_(consumes<vector<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("e_gamma"))),
  gsf_e_token1_(consumes<vector<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("gsf_electron_ged"))),
  gsf_e_token2_(consumes<vector<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("gsf_electron_unonly"))),
  conversion_token_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversions"))),
  DeDx_token_(consumes<vector<reco::DeDxHitInfo> >(iConfig.getParameter<edm::InputTag>("DeDx_token"))),
  track_token1_(consumes<vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("track_general"))),
  track_token2_(consumes<vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("track_displaced"))),
  GsfTrack_token_(consumes<vector<reco::GsfTrack> >(iConfig.getParameter<edm::InputTag>("GsfTrack_token"))),
  sc_token_(consumes<vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("sc_token"))),
  calo_token_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("calo_token"))),
  calo_eb_token_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("calo_eb_token"))),
  calo_ee_token_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("calo_ee_token"))),
  ecalrechit_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("ecalrechit"))),
  photon_token_(consumes<vector<reco::Photon> >(iConfig.getParameter<edm::InputTag>("photon_token")))
{
  edm::Service<TFileService> fs;
  eg_histo[0] = fs->make<TH1D>("Electron_energy","Electron_energy",100, 0, 1000);
  eg_histo[1] = fs->make<TH1D>("Electron_pz","Electron_pz",100, -1000, 1000);
  eg_histo[2] = fs->make<TH1D>("Eta_electron","Eta_electron",100, -5, 5);
  eg_histo[3] = fs->make<TH1D>("Phi_electron","Phi_electron",100, -6, 6);
  eg_histo[4] = fs->make<TH1D>("No_electrons","No_electrons",100,0,99);
  eg_histo[5] = fs->make<TH1D>("pT_electron","pT_electron",100,0,1000);
    
  gsf_e_histo[0] = fs->make<TH1D>("Sc_pixel_charge ged gsf","Sc_pixel_charge",5,0,5);
  gsf_e_histo[1] = fs->make<TH1D>("No_isgsf_electrons ged gsf","No_isgsf_electrons",10,0,10);
  gsf_e_histo[2] = fs->make<TH1D>("No_gsf_electrons ged gsf","No_gsf_electrons",10,0,10);
  gsf_e_histo[3] = fs->make<TH1D>("PDG_ID_electrons ged gsf","PDG_ID_electrons",60,-30,30);
  gsf_e_histo[4] = fs->make<TH1D>("Sc_pixel_charge uncleaned only","Sc_pixel_charge",5,0,5);
  gsf_e_histo[5] = fs->make<TH1D>("No_isgsf_electrons uncleaned only","No_isgsf_electrons",10,0,10);
  gsf_e_histo[6] = fs->make<TH1D>("No_gsf_electrons uncleaned only","No_gsf_electrons",10,0,10);
  gsf_e_histo[7] = fs->make<TH1D>("PDG_ID_electrons uncleaned only","PDG_ID_electrons",60,-30,30);
  gsf_e_histo[8] = fs->make<TH1D>("Pt of ged gsf electrons","Pt_GSf_ged",300,0,300);
  gsf_e_histo[9] = fs->make<TH1D>("isGsfCtfScPixChargeConsistent ged gsf","isGsfCtfScPixChargeConsistent",5,0,5);
  gsf_e_histo[10] = fs->make<TH1D>("isGsfScPixChargeConsistent ged gsf","isGsfScPixChargeConsistent",5,0,5);
  gsf_e_histo[11] = fs->make<TH1D>("shFracInnerHits ged gsf","shFracInnerHits",1000,0,2);  
  gsf_e_histo[12] = fs->make<TH1D>("eSuperClusterOverP ged gsf","eSuperClusterOverP",1000,0,100);
  gsf_e_histo[13] = fs->make<TH1D>("eSeedClusterOverP ged gsf","eSeedClusterOverP",1000,0,100);
  gsf_e_histo[14] = fs->make<TH1D>("eSeedClusterOverPout ged gsf","eSeedClusterOverPout",1000,0,100);
  gsf_e_histo[15] = fs->make<TH1D>("eEleClusterOverP ged gsf","eEleClusterOverP",1000,0,100);
  gsf_e_histo[16] = fs->make<TH1D>("deltaEtaSuperClusterTrackAtVtx ged gsf","deltaEtaSuperClusterAtVtx",1000,-5,5);
  gsf_e_histo[17] = fs->make<TH1D>("deltaEtaSeedClusterTrackAtCalo ged gsf","deltaEtaSeedClusterAtCalo",1000,-5,5);
  gsf_e_histo[18] = fs->make<TH1D>("deltaPhiSuperClusterTrackAtVtx ged gsf","deltaPhiSuperClusterAtVtx",1000,-5,5);
  gsf_e_histo[19] = fs->make<TH1D>("deltaPhiSeedClusterTrackAtCalo ged gsf","deltaPhiSeedClusterAtCalo",1000,-5,5);
  gsf_e_histo[20] = fs->make<TH1D>("sigmaEtaEta ged gsf","sigmaEtaEta",100,0,0.1);
  gsf_e_histo[21] = fs->make<TH1D>("sigmaIetaIeta ged gsf","sigmaIetaIeta",1000,0,1);
  gsf_e_histo[22] = fs->make<TH1D>("sigmaIphiIphi ged gsf","sigmaIphiIphi",1000,0,1);
  gsf_e_histo[23] = fs->make<TH1D>("e1x5 ged gsf","e1x5",1000,0,1);
  gsf_e_histo[24] = fs->make<TH1D>("hcalOverEcal ged gsf","hcalOverEcal",2000,0,200);
  gsf_e_histo[25] = fs->make<TH1D>("hcalOverEcalValid ged gsf","hcalOverEcalValid",5,0,5);
  gsf_e_histo[26] = fs->make<TH1D>("convDist ged gsf","convDist",1000,0,5);
  gsf_e_histo[27] = fs->make<TH1D>("convDcot ged gsf","convDcot",1000,-10,10);
  gsf_e_histo[28] = fs->make<TH1D>("convRadius ged gsf","convRadius",2000,-10,10);
  gsf_e_histo[29] = fs->make<TH1D>("convVtxFitProb ged gsf","convVtxFitProb",2000,-100,100);
  gsf_e_histo[30] = fs->make<TH1D>("ecalPFClusterIso ged gsf","ecalPFClusterIso",2000,0,200);
  gsf_e_histo[31] = fs->make<TH1D>("hcalPFClusterIso ged gsf","hcalPFClusterIso",2000,0,200);
  gsf_e_histo[32] = fs->make<TH1D>("trackFbrem ged gsf","trackFbrem",1000,-5,5);
  gsf_e_histo[33] = fs->make<TH1D>("superClusterFbrem ged gsf","superClusterFbrem",1000,-5,5);
  gsf_e_histo[34] = fs->make<TH1D>("numberOfBrems ged gsf","numberOfBrems",20,0,20);
  gsf_e_histo[35] = fs->make<TH1D>("fbrem_eta_8 ged gsf","fbrem",600,-3,3);
  gsf_e_histo[36] = fs->make<TH1D>("fbrem_eta_144 ged gsf","fbrem",600,-3,3);
  gsf_e_histo[37] = fs->make<TH1D>("fbrem_eta_2 ged gsf","fbrem",600,-3,3);
  gsf_e_histo[38] = fs->make<TH1D>("fbrem_eta2_ ged gsf","fbrem",600,-3,3);
  gsf_e_histo[39] = fs->make<TH1D>("Leading Pt0 of gedGsf electron","pt0_gsf",1000,0,1000);
  gsf_e_histo[40] = fs->make<TH1D>("sub-Leading Pt1 of gedGsf electron","pt1_gsf",1000,0,1000); 
  gsf_e_histo[41]= fs->make<TH1D>("flag for conversion-rejection ged gsf","flag",200,-100,100);
  gsf_e_histo[42] = fs->make<TH1D>("ecalPFClusterIso/pT ged gsf","ecalPFClusterIso/pT",500,0,5);
  gsf_e_histo[43] = fs->make<TH1D>("r9 ged gsf","r9",1000,0,100);
  gsf_e_histo[44] = fs->make<TH1D>("Pt of Good ged_gsf","Pt good_ ged Gsf",10,0,10);
  gsf_e_histo[45] = fs->make<TH1D>("Pt0 of Good  ged_gsf","Pt0 good_ged Gsf",10,0,10);
  gsf_e_histo[46] = fs->make<TH1D>("Pt of Bad ged_gsf","Pt bad_ ged Gsf",10,0,10);
  gsf_e_histo[47] = fs->make<TH1D>("Pt0 of Bad  ged_gsf","Pt0 bad_ged Gsf",10,0,10);
  gsf_e_histo[48] = fs->make<TH1D>("size of good ged gsf","size_good_gsf",10,0,10);
  gsf_e_histo[49] = fs->make<TH1D>("size of bad ged  gsf","size_bad_gsf",10,0,10);
  gsf_e_histo[50] = fs->make<TH1D>("Pt_ged_gsf_electrons_out","Pt_ged_gsf",10,0,10);
  gsf_e_histo[51] = fs->make<TH1D>("r9_good_ged_gsf","r9_good_gsf",150,0,1.5);
  gsf_e_histo[52] = fs->make<TH1D>("r9_bad_ged_gsf","r9_bad_gsf",150,0,1.5);
  gsf_e_histo[53] = fs->make<TH1D>("SigmaEtaEta_good_ged_gsf","SigmaEtaEta_good_gsf",100,0,0.1);
  gsf_e_histo[54] = fs->make<TH1D>("SigmaEtaEta_bad_ged_gsf","SigmaEtaEta_bad_gsf",100,0,0.1);
  gsf_e_histo[55] = fs->make<TH1D>("ecalClusterPFIso_good_ged_gsf","ecalIso_good_gsf",100,0,10);
  gsf_e_histo[56] = fs->make<TH1D>("ecalClusterPFIso_bad_ged_gsf","ecalIso_bad_gsf",100,0,10);
  gsf_e_histo[57] = fs->make<TH1D>("good_dR_rms","good_dR_rms",1000,0,10);
  gsf_e_histo[58] = fs->make<TH1D>("good_dR2","good_dR2",10000,0,100);
  gsf_e_histo[59] = fs->make<TH1D>("good_Energy_dR2","good_Energy_dR2",10000,0,100);
  gsf_e_histo[60] = fs->make<TH1D>("good_Energy_dR_rms","good_Energy_dR_rms",1000,0,10);
  gsf_e_histo[61] = fs->make<TH1D>("bad_dR_rms","bad_dR_rms",1000,0,10);
  gsf_e_histo[62] = fs->make<TH1D>("bad_dR2","bad_dR2",10000,0,100);
  gsf_e_histo[63] = fs->make<TH1D>("bad_Energy_dR2","bad_Energy_dR2",10000,0,100);
  gsf_e_histo[64] = fs->make<TH1D>("bad_Energy_dR_rms","bad_Energy_dR_rms",1000,0,10);
  gsf_e_histo[65] = fs->make<TH1D>("CaloCluster_Eta","CaloCluster_Eta",200,-10,10);
  gsf_e_histo[66] = fs->make<TH1D>("Size of goodIsoGsf","size_goodIsoGsf",10,0,10);
  gsf_e_histo[67] = fs->make<TH1D>("nevent goodIso","nevent_goodIso",10000,0,10000);
  gsf_e_histo[68] = fs->make<TH1D>("goodIso_eta","goodIso_eta",200,-10,10);
  //Conversions Histogram
  conversion_histo[0] = fs->make<TH1D>("Conersion nTracks","nTracks",10,0,10);
  conversion_histo[1] = fs->make<TH1D>("isConverted","isConverted",5,0,5);
  conversion_histo[2] = fs->make<TH1D>("PairInvm_conversions","PairInvm_conversions",550,-5,50);
  conversion_histo[3] = fs->make<TH1D>("paircotThetaSeparation conversions","pairCotThetaSeparation",1000,-50,50);
  conversion_histo[4] = fs->make<TH1D>("EoverP conversions","EoverP",1000,0,100);
  conversion_histo[5] = fs->make<TH1D>("EoverPrefittedTracks conversions","EoverPrefittedTracks",1000,0,100);
  conversion_histo[6] = fs->make<TH1D>("distOfMinimumApproach conversion","distOfMinimumApproach",1000,0,100);
  conversion_histo[7] = fs->make<TH1D>("dPhiTracksAtVtx conversions","dPhiTracksAtVtx",100,0,50);
  conversion_histo[8] = fs->make<TH1D>("dPhiTracksAtEcal conversions","dPhiTracksAtEcal",100,0,50);
  conversion_histo[9] = fs->make<TH1D>("dEtaTracksAtEcal conversions","dEtaTracksAtEcal",100,0,50);
  conversion_histo[10] = fs->make<TH1D>("nSharedHits conversions","nSharedHits",30,0,30);
  conversion_histo[11] = fs->make<TH1D>("Pt conversions","Pt",100,0,100);
  conversion_histo[12] = fs->make<TH1D>("min dR leading ged_gsf_elec and conversion","min dR",100,0,1);
  conversion_histo[13] = fs->make<TH1D>("Pt of good conversions","pT_good_convert",10,0,10);
  conversion_histo[14] = fs->make<TH1D>("Pt0 of good conversions","pT0_good_convert",10,0,10);
  conversion_histo[15] = fs->make<TH1D>("Mass from refitted momentum","Mass",110,-5,50);
  conversion_histo[16] = fs->make<TH1D>("size of goodConversion","size_goodConversion",10,0,10);
  conversion_histo[17] = fs->make<TH1D>("Pt_conversion_out","Pt_conversion",1000,0,10);
  conversion_histo[18] = fs->make<TH1D>("mindR_ged_gsf_conv","min dR",100,0,10);
  conversion_histo[19] = fs->make<TH1D>("Testing dR","dRmin_test",100,0,1);
  conversion_histo[20] = fs->make<TH1D>("dRmin0 of conversion and gsf","dRmin0_conv_gsf",100,0,1);
  conversion_histo[21] = fs->make<TH1D>("dRmin1 of conversion and gsf","dRmin1_conv_gsf",100,0,1);
  conversion_histo[22] = fs->make<TH1D>("dRmin0 of conversion and track","dRmin0_conv_track",100,0,1);
  conversion_histo[23] = fs->make<TH1D>("dRmin1 of conversion and track","dRmin1_conv_track",100,0,1);
  conversion_histo[24] = fs->make<TH1D>("dRmin1-dRmin0 conv_gsf","dRmin1-dRmin0_conv_gsf",2000,-1,1);
  conversion_histo[25] = fs->make<TH1D>("dRmin1-dRmin0 conv_track","dRmin1-dRmin0_conv_track",2000,-1,1);
  conversion_histo[26] = fs->make<TH1D>("dRmin_T1_gsf","dRmin_T1_gsf",100,0,1);
  conversion_histo[27] = fs->make<TH1D>("dRmin_T2_gsf","dRmin_T2_gsf",100,0,1);
  conversion_histo[28] = fs->make<TH1D>("Pt of T1","pT_T1",100,0,100);
  conversion_histo[29] = fs->make<TH1D>("Pt of T2","pT_T2",100,0,100);
  conversion_histo[30] = fs->make<TH1D>("Pt of T1+T2","pT_T1_T2",100,0,100);
  conversion_histo[31] = fs->make<TH1D>("Pt of converted_track","pT_conv_track",100,0,100);
  conversion_histo[32] = fs->make<TH1D>("Pt of T1+T2 - conv_track","pT_T1+T2-conv_track",20,-10,10);
  conversion_histo[33] = fs->make<TH1D>("Pt resolution w.r.t conversions","Pt_reso_conv",200,-5,5);
  conversion_histo[34] = fs->make<TH1D>("Pt resolution w.r.t T1+T2","Pt_reso_T1+T2",200,-5,5);
  conversion_histo[35] = fs->make<TH1D>("Invariant mass of T1+T2","inv_M_T1+T2",1050,-5,100);
  conversion_histo[36] = fs->make<TH1D>("Pair_invm_matching_conv","pair_invm_matching_conv",1050,-5,100);
  conversion_histo[37] = fs->make<TH1D>("dRmin_sc_T1","dRmin_sc_T1",100,0,10);
  conversion_histo[38] = fs->make<TH1D>("dRmin_sc_T2","dRmin_sc_T2",100,0,10);
  conversion_histo[39] = fs->make<TH1D>("dRmin_sc_conv_track","dRmin_sc_conv_track",100,0,10);
  conversion_histo[40] = fs->make<TH1D>("sc_phi","sc_phi",100,-5,5);
  conversion_histo[41] = fs->make<TH1D>("sc_eta","sc_eta",100,-5,5);


  //Merged electrons info in 2-D histogram
    
  merged_histo[0] = fs->make<TH2D>("dRmin0_vs_dRmin1 conv and gsf","dRmin0_dRmin1_conv_gsf",20,0,1,20,0,1);
  merged_histo[1] = fs->make<TH2D>("dRmin0_vs_dRmin1 conv and track","dRmin0_dRmin1_conv_track",20,0,1,20,0,1);
  merged_histo[2] = fs->make<TH2D>("dRmin0_vs_dRmin1-dRmin0_gsf","dRmin0_vs_delta_dRmin_gsf",20,0,1,20,0,1);
  merged_histo[3] = fs->make<TH2D>("dRmin0_vs_dRmin1-dRmin0_track","dRmin0_vs_delta_dRmin_track",20,0,1,20,0,1);
  merged_histo[4] = fs->make<TH2D>("good_Eta_phi_pT_0_10","good_Eta_phi_pT_0_10",200,-1,1,200,-1,1);
  merged_histo[5] = fs->make<TH2D>("good_Eta_phi_pT_10_20","good_Eta_phi_pT_10_20",200,-1,1,200,-1,1);
  merged_histo[6] = fs->make<TH2D>("good_Eta_phi_pT_20_50","good_Eta_phi_pT_20_50",200,-1,1,200,-1,1);
  merged_histo[7] = fs->make<TH2D>("good_Eta_phi_pT_100","good_Eta_phi_pT_100",200,-1,1,200,-1,1);
  merged_histo[8] = fs->make<TH2D>("bad_Eta_phi_pT_0_10","bad_Eta_phi_pT_0_10",200,-1,1,200,-1,1);
  merged_histo[9] = fs->make<TH2D>("bad_Eta_phi_pT_10_20","bad_Eta_phi_pT_10_20",200,-1,1,200,-1,1);
  merged_histo[10] = fs->make<TH2D>("bad_Eta_phi_pT_20_50","bad_Eta_phi_pT_20_50",200,-1,1,200,-1,1);
  merged_histo[11] = fs->make<TH2D>("bad_Eta_phi_pT_100","bad_Eta_phi_pT_100",200,-1,1,200,-1,1);
  merged_histo[12] = fs->make<TH2D>("dEta_dPhi_gsfISo_anevent","dEta_dPhi_gsfIso",200,-1,1,200,-1,1);


  //DEDxHitINfo Histograms
  DeDx_histo[0] = fs->make<TH1D>("DeDx charge","charge",2000,-1000,1000);
  DeDx_histo[1] = fs->make<TH1D>("DeDx pathlength","pathlength",100,0,10);

  //Track Histograms
  track_histo[0] = fs->make<TH1D>("normalisedChi2 general tracks","normalisedChi2",100,-50,50);
  track_histo[1] = fs->make<TH1D>("Charge General Tracks","charge",100,-50,50);
  track_histo[2] = fs->make<TH1D>("qoverp general tracks","qoverp",200,-100,100);
  track_histo[3] = fs->make<TH1D>("pT general tracks","pT",1000,0,10);
  track_histo[4] = fs->make<TH1D>("Phi genral tracks","Phi",100,-5,5);
  track_histo[5] = fs->make<TH1D>("Eta general tracks","Eta",100,-5,5);
  track_histo[6] = fs->make<TH1D>("t0 general tracks","t0",200,-100,100);
  track_histo[7] = fs->make<TH1D>("beta general tracks","beta",100,0,1);
  track_histo[8] = fs->make<TH1D>("etaError general tracks","etaError",100,-50,50);
  track_histo[9] = fs->make<TH1D>("PhiError general tracks","Phierror",100,-50,50);
  track_histo[10] = fs->make<TH1D>("normalisedChi2 displaced tracks","normalisedChi2",100,-50,50);
  track_histo[11] = fs->make<TH1D>("Charge displaced Tracks","charge",100,-50,50);
  track_histo[12] = fs->make<TH1D>("qoverp displaced tracks","qoverp",200,-100,100);
  track_histo[13] = fs->make<TH1D>("pT displaced tracks","pT",10,0,10);
  track_histo[14] = fs->make<TH1D>("Phi displaced tracks","Phi",100,-5,5);
  track_histo[15] = fs->make<TH1D>("Eta displaced tracks","Eta",100,-5,5);
  track_histo[16] = fs->make<TH1D>("t0 displaced tracks","t0",200,-100,100);
  track_histo[17] = fs->make<TH1D>("beta displaced tracks","beta",100,0,1);
  track_histo[18] = fs->make<TH1D>("etaError displaced tracks","etaError",100,-50,50);
  track_histo[19] = fs->make<TH1D>("PhiError displaced tracks","Phierror",100,-50,50);
  track_histo[20] = fs->make<TH1D>("dRmin of good general tracks","dRmin_good_general",100,0,10);
  track_histo[21] = fs->make<TH1D>("dRmin of bad general tracks","dRmin_bad_general",100,0,10);
  track_histo[22] = fs->make<TH1D>("dRmin of good displaced tracks","dRmin_good_displaced",100,0,10);
  track_histo[23] = fs->make<TH1D>("dRmin of bad displaced tracks","dRmin_bad_displaced",100,0,10);

  //GsfTrack Histograms
  
  GsfTrack_histo[0] = fs->make<TH1D>("chargeMode GsfTrack","chargeMode",100,-50,50);
  GsfTrack_histo[1] = fs->make<TH1D>("qoverpMode GsfTrack","qoverpMode",200,-100,100);
  GsfTrack_histo[2] = fs->make<TH1D>("thetaMode GsfTrack","thetaMode",200,-10,10);
  GsfTrack_histo[3] = fs->make<TH1D>("lambdaMode GsfTrack","lambdaMode",200,-10,10);
  GsfTrack_histo[4] = fs->make<TH1D>("pTMode GsfTrack","pTMode",10,0,10);
  GsfTrack_histo[5] = fs->make<TH1D>("phiMode GsfTrack","phiMode",100,-5,5);
  GsfTrack_histo[6] = fs->make<TH1D>("etaMode GsfTrack","etaMode",100,-5,5);
  GsfTrack_histo[7] = fs->make<TH1D>("thetaModeError GsfTrack","thetaModeError",200,-10,10);
  GsfTrack_histo[8] = fs->make<TH1D>("lambdaModeError GsfTrack","lambdaModeError",200,-10,10);
  GsfTrack_histo[9] = fs->make<TH1D>("etaModeError GsfTrack","etaModeError",200,-10,10);
  GsfTrack_histo[10] = fs->make<TH1D>("phiModeError GsfTrack","phiModeError",200,-10,10);
  GsfTrack_histo[11] = fs->make<TH1D>("dRmin GsfTrack_goodGsf","dRmin_good_GsfTrack",100,0,10);
  GsfTrack_histo[12] = fs->make<TH1D>("dRmin GsfTrack_badGsf","dRmin_bad_GsfTrack",100,0,10);

  //CaloCluster Histogram
  calo_histo[0] = fs->make<TH1D>("Calocluster_EB_Eta","CaloCluster_EB_Eta",200,-10,10);
  calo_histo[1] = fs->make<TH1D>("Calocluster_EE_Eta","CaloCluster_EE_Eta",200,-10,10);
  calo_histo[2] = fs->make<TH1D>("SC_eta","SC_eta",200,-10,10);
  calo_histo[3] = fs->make<TH1D>("no_calo_04","no_calo_04",10,0,10);
  calo_histo[4] = fs->make<TH1D>("dRmin0+dRmin1_goodIso_calo","dRmin0+dRmin1_goodIso_calo",100,0,10);
  calo_histo[5] = fs->make<TH1D>("Energy_weighted_dR_by_E_elec","Energy_dR_by_E_elec",150,0,150);

  //Photon Histogram

  photon_histo[0] = fs->make<TH1D>("Photon seed Eta","Photon_seed_Eta",200,-10,10);
  photon_histo[1] = fs->make<TH1D>("Photon seed Phi","Photon_seed_Phi",100,-5,5);
  photon_histo[2] = fs->make<TH1D>("Converted r9","converted_r9",100,0,1);
  photon_histo[3] = fs->make<TH1D>("unconverted r9","unconverted_r9",100,0,1);

  //TGraph2D 
    
  gr_good = fs->make<TGraph2D>();
  gr_good->SetTitle("GoodGsf;dEta;dPhi;Energy");
  gr_bad = fs->make<TGraph2D>();
  gr_bad->SetTitle("badGsf;dEta;dPhi;Energy");
  gr_gsf = fs->make<TGraph2D>();
  gr_gsf->SetTitle("Gsf;dEta;dPhi;Energy");
  gr_good_ecal = fs->make<TGraph2D>();
  gr_good_ecal->SetTitle("GoodGsf_ECAL;Eta;Phi:Energy");
  gr_bad_ecal = fs->make<TGraph2D>();
  gr_bad_ecal->SetTitle("BadGsf_ECAL;Eta;Phi:Energy");




  //TH3D 

  mustache_histo[0] = fs->make<TH3D>("goodGsf Energy vs eta-phi","energy_eta_phi",500,-5,5,500,-5,5,100,0,100);    
  mustache_histo[1] = fs->make<TH3D>("badGsf Energy vs eta-phi","energy_eta_phi",500,-5,5,500,-5,5,100,0,100);


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

Analyzer_egamma::~Analyzer_egamma() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void Analyzer_egamma::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // using namespace edm;

  //Clean the Container
  elec.clear();
  gsf_elec.clear();
  goodIsoGsf_elec.clear();
  conversion.clear();
  DeDx.clear();
  goodConversion.clear();
  goodGsf_elec.clear();
  badGsf_elec.clear();
  track_general.clear();
  track_displaced.clear();
  GsfTrack.clear();
  convtype.clear();
  T1_conv.clear();
  T2_conv.clear();
  conv_track.clear();
  sc.clear();
  CaloCluster.clear();
  CaloCluster_EB.clear();
  CaloCluster_EE.clear();
  nevent_goodIso.clear();
  //Define the handler,a token and get the information by token
  //Handle<vector<pat::Electron> > myelec;
  //iEvent.getByToken(etoken_,myelec);
  Handle<vector<reco::GsfElectron> > mygsf_elec_ged;
  iEvent.getByToken(gsf_e_token1_,mygsf_elec_ged);
  Handle<vector<reco::GsfElectron> > mygsf_elec_unonly;
  iEvent.getByToken(gsf_e_token2_,mygsf_elec_unonly);
  Handle<vector<reco::Conversion> > myconversion;
  iEvent.getByToken(conversion_token_,myconversion);
  Handle<vector<reco::DeDxHitInfo> > myDeDx;
  iEvent.getByToken(DeDx_token_,myDeDx);
  Handle<vector<reco::Track> > myTrack_general;
  iEvent.getByToken(track_token1_,myTrack_general);
  Handle<vector<reco::Track> > myTrack_displaced;
  iEvent.getByToken(track_token2_,myTrack_displaced);
  Handle<vector<reco::GsfTrack> > myGsfTrack;
  iEvent.getByToken(GsfTrack_token_,myGsfTrack);
  Handle<vector<reco::SuperCluster> > mySc;
  iEvent.getByToken(sc_token_,mySc);
  Handle<vector<reco::CaloCluster> > myCalo;
  iEvent.getByToken(calo_token_,myCalo);
  Handle<vector<reco::CaloCluster> > myCalo_eb;
  iEvent.getByToken(calo_eb_token_,myCalo_eb);
  Handle<vector<reco::CaloCluster> > myCalo_ee;
  iEvent.getByToken(calo_ee_token_,myCalo_ee);
  Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> > myrechit;
  iEvent.getByToken(ecalrechit_,myrechit);
  Handle<vector<reco::Photon> > myphoton;
  iEvent.getByToken(photon_token_,myphoton);


  
  fstream file("ECAL_rechit.txt",ios::app);
  fstream file1("Photon_evt50_rechit.txt",ios::app);
  fstream file2("unconverted_photon_evt53_rechit.txt",ios::app);


  nevent++;
  //Define Temporary Variables

  //EGamma temp;
  //int n_e;

  // Lepton temp;
  int n_isgsf,n_gsf;
 
  n_gsf=0;
  n_isgsf=0;

  //if collection is valid loop over electron events
  /* if(myelec.isValid())
     {
     for(const pat::Electron &ie : *myelec)
     {
     temp.v_e.SetPxPyPzE(ie.px(),ie.py(),ie.pz(),ie.energy());
     temp.eta = temp.v_e.Eta();      
     temp.phi = temp.v_e.Phi();
     temp.pT = temp.v_e.Pt();
     if(ie.isElectron())
     {
     n_e++;
     }

     elec.push_back(temp);

     eg_histo[0]->Fill(temp.v_e.E());     
     eg_histo[1]->Fill(temp.v_e.Pz());
     eg_histo[2]->Fill(temp.eta);
     eg_histo[3]->Fill(temp.phi);
     eg_histo[5]->Fill(temp.pT);
     }
     temp.n_electrons = n_e;
     elec.push_back(temp);
     eg_histo[4]->Fill(n_e);

     //for (const auto& track : iEvent.get(tracksToken_)) {
     // do something with track parameters, e.g, plot the charge.
     // int charge = track.charge();
     }*/
  for(const reco::GsfElectron &ie : *mygsf_elec_ged)
    {
      Lepton temp;
      temp.scPixCharge = ie.scPixCharge();
      temp.id = ie.pdgId();
      temp.v.SetPtEtaPhiM(ie.pt(),ie.eta(),ie.phi(),0.000511);
      temp.isGsfCtfScPixChargeConsistent = ie.isGsfCtfScPixChargeConsistent();
      temp.isGsfScPixChargeConsistent = ie.isGsfScPixChargeConsistent();
      temp.shFracInnerHits = ie.shFracInnerHits();
      temp.eSuperClusterOverP = ie.eSuperClusterOverP();
      temp.eSeedClusterOverP = ie.eSeedClusterOverP();
      temp.eSeedClusterOverPout = ie.eSeedClusterOverPout();
      temp.eEleClusterOverPout = ie.eEleClusterOverPout();
      temp.deltaEtaSuperClusterTrackAtVtx = ie.deltaEtaSuperClusterTrackAtVtx();
      temp.deltaEtaSeedClusterTrackAtCalo = ie.deltaEtaSeedClusterTrackAtCalo();
      temp.deltaPhiSuperClusterTrackAtVtx = ie.deltaPhiSuperClusterTrackAtVtx();
      temp.deltaPhiSeedClusterTrackAtCalo = ie.deltaPhiSeedClusterTrackAtCalo();
      temp.sigmaEtaEta = ie.sigmaEtaEta();
      temp.sigmaIetaIeta = ie.sigmaIetaIeta();
      temp.sigmaIphiIphi = ie.sigmaIphiIphi();
      temp.e1x5 = ie.e1x5();
      temp.hcalOverEcal = ie.hcalOverEcal(0);
      temp.hcalOverEcalValid = ie.hcalOverEcalValid();
      temp.convDist = ie.convDist();
      temp.convDcot = ie.convDcot();
      temp.convRadius = ie.convRadius();
      temp.convRadius = ie.convRadius();
      temp.ecalPFClusterIso = ie.ecalPFClusterIso();
      temp.hcalPFClusterIso = ie.hcalPFClusterIso();
      temp.trackFbrem = ie.trackFbrem();
      temp.superClusterFbrem = ie.superClusterFbrem();
      temp.numberOfBrems = ie.numberOfBrems();
      temp.convFlags = ie.convFlags();
      temp.r9 = ie.r9();
      temp.ecalEnergy = ie.ecalEnergy();
      if(abs(ie.eta()<0.8)) temp.fbrem_eta_8 = ie.fbrem();
      if(abs(ie.eta())>0.8 && abs(ie.eta())<1.44) temp.fbrem_eta_144 = ie.fbrem();
      if(abs(ie.eta())>1.57 && abs(ie.eta())<2) temp.fbrem_eta_2 = ie.fbrem();
      if(abs(ie.eta())>2) temp.fbrem_eta_2_ = ie.fbrem();
	  
      gsf_elec.push_back(temp);

      if(ie.isElectron())
	{
	  n_isgsf++;
	}
      n_gsf++;

      gsf_e_histo[0]->Fill(temp.scPixCharge);
      gsf_e_histo[3]->Fill(temp.id);
      gsf_e_histo[8]->Fill(temp.v.Pt());
      gsf_e_histo[9]->Fill(temp.isGsfCtfScPixChargeConsistent);
      gsf_e_histo[10]->Fill(temp.isGsfScPixChargeConsistent);
      gsf_e_histo[11]->Fill(temp.shFracInnerHits);
      gsf_e_histo[12]->Fill(temp.eSuperClusterOverP);
      gsf_e_histo[13]->Fill(temp.eSeedClusterOverP);
      gsf_e_histo[14]->Fill(temp.eSeedClusterOverPout);
      gsf_e_histo[15]->Fill(temp.eEleClusterOverPout);
      gsf_e_histo[16]->Fill(temp.deltaEtaSuperClusterTrackAtVtx);
      gsf_e_histo[17]->Fill(temp.deltaEtaSeedClusterTrackAtCalo);
      gsf_e_histo[18]->Fill(temp.deltaPhiSuperClusterTrackAtVtx);
      gsf_e_histo[19]->Fill(temp.deltaPhiSeedClusterTrackAtCalo);
      gsf_e_histo[20]->Fill(temp.sigmaEtaEta);
      gsf_e_histo[21]->Fill(temp.sigmaEtaEta);
      gsf_e_histo[22]->Fill(temp.sigmaIetaIeta);
      gsf_e_histo[23]->Fill(temp.sigmaIphiIphi);
      gsf_e_histo[24]->Fill(temp.e1x5);
      gsf_e_histo[25]->Fill(temp.hcalOverEcal);
      gsf_e_histo[26]->Fill(temp.hcalOverEcalValid);
      gsf_e_histo[27]->Fill(temp.convDist);
      gsf_e_histo[28]->Fill(temp.convDcot);
      gsf_e_histo[29]->Fill(temp.convRadius);
      gsf_e_histo[30]->Fill(temp.ecalPFClusterIso);
      gsf_e_histo[31]->Fill(temp.hcalPFClusterIso);
      gsf_e_histo[32]->Fill(temp.trackFbrem);
      gsf_e_histo[33]->Fill(temp.superClusterFbrem);
      gsf_e_histo[34]->Fill(temp.numberOfBrems);
      gsf_e_histo[35]->Fill(temp.fbrem_eta_8);
      gsf_e_histo[36]->Fill(temp.fbrem_eta_144);
      gsf_e_histo[37]->Fill(temp.fbrem_eta_2);
      gsf_e_histo[38]->Fill(temp.fbrem_eta_2_);
      gsf_e_histo[41]->Fill(temp.convFlags);
      gsf_e_histo[42]->Fill(temp.ecalPFClusterIso/temp.v.Pt());
      gsf_e_histo[43]->Fill(temp.r9);
    }
  //temp.n_gsf_e = n_is_gsf_e;
  // gsf_elec.push_back(gsf_temp);
  gsf_e_histo[1]->Fill(n_isgsf);
  gsf_e_histo[2]->Fill(n_gsf);

  n_gsf=0;
  n_isgsf=0;

  //Sorting     
  sort(0);

  for(unsigned int i=0;i<gsf_elec.size();i++) gsf_e_histo[50]->Fill(gsf_elec.at(i).v.Pt());

  if((int)gsf_elec.size()>1){
    gsf_e_histo[39]->Fill(gsf_elec.at(0).v.Pt());
    gsf_e_histo[40]->Fill(gsf_elec.at(1).v.Pt());
  }


  for(const reco::GsfElectron &ie : *mygsf_elec_unonly)
    {
      Lepton temp;
      temp.scPixCharge = ie.scPixCharge();
      temp.id = ie.pdgId();
      temp.v.SetPtEtaPhiM(ie.pt(),ie.eta(),ie.phi(),0.000511);
      //gsf_elec.push_back(temp);


      if(ie.isElectron())
	{
	  n_isgsf++;
	}
      n_gsf++;

      gsf_e_histo[4]->Fill(temp.scPixCharge);
      gsf_e_histo[7]->Fill(ie.pdgId());
    }
  //temp.n_gsf_e = n_is_gsf_e;                                                                                                                                                                          
  // gsf_elec.push_back(gsf_temp);                                                                                                                                                                      
  gsf_e_histo[5]->Fill(n_isgsf);
  gsf_e_histo[6]->Fill(n_gsf);


  //for (const auto& track : iEvent.get(tracksToken_)) {                                                                                                                                                
  // do something with track parameters, e.g, plot the charge.                                                                                                                                          
  // int charge = track.charge();                                                                                                                                                                       
  for(const reco::Conversion &ie : *myconversion)
    {
      Lepton temp;
      Lepton temp_gsf;
      temp.isConverted = ie.isConverted();
      temp.nTracks = ie.nTracks();
      temp.pairInvariantMass = ie.pairInvariantMass();
      temp.pairCotThetaSeparation = ie.pairCotThetaSeparation();
      temp.EoverP = ie.EoverP();
      temp.EoverPrefittedTracks = ie.EoverPrefittedTracks();
      temp.distOfMinimumApproach = ie.distOfMinimumApproach();
      // temp.dPhiTracksAtVtx = ie.dPhiTracksAtVtx();
      //temp.dPhiTracksAtEcal = ie.dPhiTracksAtEcal();
      //temp.dEtaTracksAtEcal = ie.dEtaTracksAtEcal();
      temp.nSharedHits = ie.nSharedHits();
      //std::cout<<"Phi: "<<ie.phi()<<"eta: "<<ie.eta()<<std::endl; 
      temp.v.SetPtEtaPhiM(ie.refittedPair4Momentum().Pt(),ie.refittedPair4Momentum().Eta(),ie.refittedPair4Momentum().Phi(),ie.refittedPair4Momentum().M());
      conversion.push_back(temp);
      //bool good_convert = false;
      //bool good_gsf = false;
      /*for(unsigned int i=0; i<gsf_elec.size();i++){
	if(gsf_elec.at(i).v.DeltaR(temp.v) < 0.1){
	good_convert = true;
	break;
	}
	}

	for(unsigned int i=0; i<gsf_elec.size();i++){
	if(gsf_elec.at(i).v.DeltaR(temp.v) >= 0.4){
	good_gsf = true;
	temp_gsf = gsf_elec.at(i);
	break;
	}
	}*/


      //if(good_convert) goodConversion.push_back(temp);
      //if(good_gsf) goodGsf_elec.push_back(temp_gsf);

      conversion_histo[0]->Fill(temp.nTracks);
      conversion_histo[1]->Fill(temp.isConverted);
      conversion_histo[2]->Fill(temp.pairInvariantMass);
      conversion_histo[3]->Fill(temp.pairCotThetaSeparation);
      conversion_histo[4]->Fill(temp.EoverP);
      conversion_histo[5]->Fill(temp.EoverPrefittedTracks);
      conversion_histo[6]->Fill(temp.distOfMinimumApproach);
      //conversion_histo[7]->Fill(ie.dPhiTracksAtVtx());
      //conversion_histo[8]->Fill(ie.dPhiTracksAtEcal());
      //conversion_histo[9]->Fill(ie.dEtaTracksAtEcal());
      conversion_histo[10]->Fill(temp.nSharedHits);
      conversion_histo[11]->Fill(temp.v.Pt());
      conversion_histo[15]->Fill(temp.v.M());
        
    }

  sort(0);

  for(unsigned int i=0;i<conversion.size();i++) conversion_histo[17]->Fill(conversion.at(i).v.Pt());
  
  float dRmin0,test;

  if(gsf_elec.size()>0){

    dRmin0 = 9999;
    test = 9999;
    for(unsigned int i=0;i<conversion.size();i++){
      float temp0 = gsf_elec.at(0).v.DeltaR(conversion.at(i).v);
      float test0 = deltaR(gsf_elec.at(0).v.Eta(),gsf_elec.at(0).v.Phi(),conversion.at(i).v.Eta(),conversion.at(i).v.Phi());
      if(temp0 < dRmin0) dRmin0 = temp0;
      if(test0 < test ) test = test0;
    }
    conversion_histo[12]->Fill(dRmin0);
    conversion_histo[19]->Fill(test);
  }

  for(unsigned int i=0;i<gsf_elec.size();i++){
    float dRmin = 9999;
    for(unsigned int j=0;j<conversion.size();j++){
      float temp = gsf_elec.at(i).v.DeltaR(conversion.at(j).v);
      if(temp<dRmin) dRmin = temp;
    }
    conversion_histo[18]->Fill(dRmin);
  }
          
   
  // conversion_histo[12]->Fill(dRmin0);


  // bool good_convert = false;
  //bool good_gsf = false;
  //Lepton temp;
  int ind_gsf;
  int ind_conv;
  float dRmin;
    
  for(unsigned int i=0;i<gsf_elec.size();i++){
    dRmin = 9999;
    for(unsigned int j=0;j<conversion.size();j++){
      float temp = gsf_elec.at(i).v.DeltaR(conversion.at(j).v);
      if(temp < dRmin){
	dRmin = temp;
	ind_gsf = i;
	ind_conv = j;
      }
    }
    
    if(dRmin < 0.1 && gsf_elec.size()>0 && conversion.size()>0) {

      badGsf_elec.push_back(gsf_elec.at(ind_gsf));
      goodConversion.push_back(conversion.at(ind_conv));
    }

    if(dRmin > 0.4 && gsf_elec.size()>0 && conversion.size()>0) goodGsf_elec.push_back(gsf_elec.at(ind_gsf));
  }

  /* if(dRmin < 0.1 && gsf_elec.size()>0 && conversion.size()>0) {

     badGsf_elec.push_back(gsf_elec.at(ind_gsf));
     goodConversion.push_back(conversion.at(ind_conv));
     }
  
     if(dRmin > 0.4 && gsf_elec.size()>0 && conversion.size()>0) goodGsf_elec.push_back(gsf_elec.at(ind_gsf));
  */   
  
  sort(0);
        
  //Conversion

  for(unsigned int i=0;i<goodConversion.size();i++){

    conversion_histo[13]->Fill(goodConversion.at(i).v.Pt());
  }
     
  if(goodConversion.size()>0)  conversion_histo[14]->Fill(goodConversion.at(0).v.Pt());

  //goodGsf 
  for(unsigned int i=0;i<goodGsf_elec.size();i++){

    gsf_e_histo[44]->Fill(goodGsf_elec.at(i).v.Pt());
  }
    
  if(goodGsf_elec.size()>0) gsf_e_histo[45]->Fill(goodGsf_elec.at(0).v.Pt());

  //badGsf
  for(unsigned int i=0;i<badGsf_elec.size();i++){

    gsf_e_histo[46]->Fill(badGsf_elec.at(i).v.Pt());
  }

  if(badGsf_elec.size()>0) gsf_e_histo[47]->Fill(badGsf_elec.at(0).v.Pt());

  //Storing the size of objects
  gsf_e_histo[48]->Fill(goodGsf_elec.size());
  gsf_e_histo[49]->Fill(badGsf_elec.size());
  conversion_histo[16]->Fill(goodConversion.size());
  
  //Showr shape variables for good and bad gsf

  for(unsigned int i=0;i<goodGsf_elec.size();i++){ 
    gsf_e_histo[51]->Fill(goodGsf_elec.at(i).r9);
    gsf_e_histo[53]->Fill(goodGsf_elec.at(i).sigmaEtaEta);
    gsf_e_histo[55]->Fill(goodGsf_elec.at(i).ecalPFClusterIso);
  }

  for(unsigned int i=0;i<badGsf_elec.size();i++){
    gsf_e_histo[52]->Fill(badGsf_elec.at(i).r9);
    gsf_e_histo[54]->Fill(badGsf_elec.at(i).sigmaEtaEta);
    gsf_e_histo[56]->Fill(badGsf_elec.at(i).ecalPFClusterIso);
  }

  for(unsigned int i=0;i<goodGsf_elec.size();i++){
    if(goodGsf_elec.at(i).ecalPFClusterIso/goodGsf_elec.at(i).v.Pt() <0.1) goodIsoGsf_elec.push_back(goodGsf_elec.at(i));
  }

  sort(0);

  for(const reco::DeDxHitInfo &ie : *myDeDx)
    {
      Lepton temp;
      temp.charge = ie.charge(ie.size());
      temp.pathlength = ie.pathlength(ie.size());
        
      DeDx_histo[0]->Fill(temp.charge);
      DeDx_histo[1]->Fill(temp.pathlength);                                                                                                             
	                                                                                                                                 
	                                                                                                                                       
    }


  for(const reco::Track &ie : *myTrack_general)
    {
      Lepton temp;
      temp.normalizedChi2 = ie.normalizedChi2();
      temp.charge = ie.charge();
      temp.qoverp = ie.qoverp();
      temp.p = ie.p();
      temp.pt = ie.pt();
      temp.px = ie.px();
      temp.py= ie.py();
      temp.pz= ie.pz();
      temp.phi = ie.phi();
      temp.eta = ie.eta();
      temp.v.SetPtEtaPhiM(ie.pt(),ie.eta(),ie.phi(),0.00511);
      temp.t0 = ie.t0();   
      temp.beta = ie.beta();
      temp.etaError = ie.etaError();
      temp.phiError = ie.phiError();
        
      track_general.push_back(temp);

      track_histo[0]->Fill(temp.normalizedChi2);
      track_histo[1]->Fill(temp.charge);
      track_histo[2]->Fill(temp.qoverp);
      track_histo[3]->Fill(temp.pt);
      track_histo[4]->Fill(temp.phi);
      track_histo[5]->Fill(temp.eta);
      track_histo[6]->Fill(temp.t0);
      track_histo[7]->Fill(temp.beta);
      track_histo[8]->Fill(temp.etaError);
      track_histo[9]->Fill(temp.phiError);

    }

  for(const reco::Track &ie : *myTrack_displaced)
    {
      Lepton temp;
      temp.normalizedChi2 = ie.normalizedChi2();
      temp.charge = ie.charge();
      temp.qoverp = ie.qoverp();
      temp.p = ie.p();
      temp.pt = ie.pt();
      temp.phi = ie.phi();
      temp.eta = ie.eta();
      temp.t0 = ie.t0();
      temp.beta = ie.beta();
      temp.etaError = ie.etaError();
      temp.phiError = ie.phiError();

      track_displaced.push_back(temp);

      track_histo[10]->Fill(temp.normalizedChi2);
      track_histo[11]->Fill(temp.charge);
      track_histo[12]->Fill(temp.qoverp);
      track_histo[13]->Fill(temp.pt);
      track_histo[14]->Fill(temp.phi);
      track_histo[15]->Fill(temp.eta);
      track_histo[16]->Fill(temp.t0);
      track_histo[17]->Fill(temp.beta);
      track_histo[18]->Fill(temp.etaError);
      track_histo[19]->Fill(temp.phiError);

    }

  sort(0);

  for(unsigned int i=0;i<goodGsf_elec.size();i++){
    float dRmin = 9999;
    //int ind;
    for(unsigned int j=0;j<track_general.size();j++){
      float temp = deltaR(track_general.at(j).eta,track_general.at(j).phi,goodGsf_elec.at(i).v.Eta(),goodGsf_elec.at(i).v.Phi());
      if(temp < dRmin){
	dRmin = temp;
	// ind = i;
      }
    }
    track_histo[20]->Fill(dRmin);
  }

  for(unsigned int i=0;i<badGsf_elec.size();i++){
    float dRmin = 9999;
    //int ind;
    for(unsigned int j=0;j<track_general.size();j++){
      float temp = deltaR(track_general.at(j).eta,track_general.at(j).phi,badGsf_elec.at(i).v.Eta(),badGsf_elec.at(i).v.Phi());
      if(temp < dRmin){
	dRmin = temp;
	// ind =i;
      }
    }
    track_histo[21]->Fill(dRmin);
  }

  for(unsigned int i=0;i<goodGsf_elec.size();i++){
    float dRmin = 9999;
    //int ind;
    for(unsigned int j=0;j<track_displaced.size();j++){
      float temp = deltaR(track_displaced.at(j).eta,track_displaced.at(j).phi,goodGsf_elec.at(i).v.Eta(),goodGsf_elec.at(i).v.Phi());
      if(temp < dRmin){
	dRmin = temp;
	// ind =i;
      }
    }
    track_histo[22]->Fill(dRmin);
  }

  for(unsigned int i=0;i<badGsf_elec.size();i++){
    float dRmin = 9999;
    //int ind;
    for(unsigned int j=0;j<track_displaced.size();j++){
      float temp = deltaR(track_displaced.at(j).eta,track_displaced.at(j).phi,badGsf_elec.at(i).v.Eta(),badGsf_elec.at(i).v.Phi());
      if(temp < dRmin){
	dRmin = temp;
	// ind =i;
      }
    }
    track_histo[23]->Fill(dRmin);
  }

  for(const reco::GsfTrack &ie : *myGsfTrack)
    {
      Lepton temp;
      temp.chargeMode = ie.chargeMode();
      temp.qoverpMode = ie.qoverpMode();
      temp.thetaMode = ie.thetaMode();
      temp.lambdaMode = ie.lambdaMode();
      temp.pMode = ie.pMode();
      temp.ptMode = ie.ptMode();
      temp.pxMode = ie.pxMode();
      temp.pyMode = ie.pyMode();
      temp.pzMode = ie.pzMode();
      temp.phiMode = ie.phiMode();
      temp.etaMode = ie.etaMode();
      temp.thetaModeError = ie.thetaModeError();
      temp.lambdaModeError = ie.lambdaModeError();
      temp.etaModeError = ie.etaModeError();
      temp.phiModeError = ie.phiModeError();

      GsfTrack.push_back(temp);
 
      GsfTrack_histo[0]->Fill(temp.chargeMode);
      GsfTrack_histo[1]->Fill(temp.qoverpMode);
      GsfTrack_histo[2]->Fill(temp.thetaMode);
      GsfTrack_histo[3]->Fill(temp.lambdaMode);
      GsfTrack_histo[4]->Fill(temp.ptMode);
      GsfTrack_histo[5]->Fill(temp.phiMode);
      GsfTrack_histo[6]->Fill(temp.etaMode);
      GsfTrack_histo[7]->Fill(temp.thetaModeError);
      GsfTrack_histo[8]->Fill(temp.lambdaModeError);
      GsfTrack_histo[9]->Fill(temp.etaModeError);
      GsfTrack_histo[10]->Fill(temp.phiModeError);
   
    }

  sort(0);

  for(unsigned int i=0;i<goodGsf_elec.size();i++){
    float dRmin = 9999;
    //int ind;                                                                                                                                                                                   

    for(unsigned int j=0;j<GsfTrack.size();j++){
      float temp = deltaR(GsfTrack.at(j).eta,GsfTrack.at(j).phi,goodGsf_elec.at(i).v.Eta(),goodGsf_elec.at(i).v.Phi());
      if(temp < dRmin){
	dRmin = temp;
	// ind =i;      
      }
    }
    GsfTrack_histo[11]->Fill(dRmin);
  }


  for(unsigned int i=0;i<badGsf_elec.size();i++){
    float dRmin = 9999;
    //int ind;                                                                                                                                                                        
                                                                                          
    for(unsigned int j=0;j<GsfTrack.size();j++){
      float temp = deltaR(GsfTrack.at(j).eta,GsfTrack.at(j).phi,badGsf_elec.at(i).v.Eta(),badGsf_elec.at(i).v.Phi());
      if(temp < dRmin){
	dRmin = temp;
	// ind =i;                                                                                                                                                                                  
                                              
      }
    }
    GsfTrack_histo[12]->Fill(dRmin);
  }

    
  for(unsigned int i=0;i<conversion.size();i++){
    float dRmin0_gsf = 9999;
    float dRmin1_gsf = 9999;
    float dRmin0_track = 9999;
    float dRmin1_track = 9999;
    //int ind0_gsf,ind1_gsf;
    int ind0_track,ind1_track;
    for(unsigned int j=0;j<gsf_elec.size();j++){
      float dR = conversion.at(i).v.DeltaR(gsf_elec.at(j).v);
      if(dR < dRmin0_gsf){
	dRmin1_gsf = dRmin0_gsf;
	dRmin0_gsf = dR;
	// ind0_gsf = i;
      }
      if(dR < dRmin1_gsf && dR>dRmin0_gsf){
	dRmin1_gsf = dR;
	// ind1_gsf = i;
      }
    }
    for(unsigned int j=0;j<track_general.size();j++){
      float dR = deltaR(conversion.at(i).v.Eta(),conversion.at(i).v.Phi(),track_general.at(j).eta,track_general.at(j).phi);
      if(dR < dRmin0_track){
	dRmin1_track = dRmin0_track;
	ind1_track = ind0_track;
	dRmin0_track =dR;
	ind0_track = j;
      }
      if(dR < dRmin1_track && dR>dRmin0_track){
	dRmin1_track = dR;
	ind1_track = j;
      }
    }
      
    conversion_histo[20]->Fill(dRmin0_gsf);
    conversion_histo[21]->Fill(dRmin1_gsf);
    conversion_histo[22]->Fill(dRmin0_track);
    conversion_histo[23]->Fill(dRmin1_track);
    conversion_histo[24]->Fill(dRmin1_gsf-dRmin0_gsf);
    conversion_histo[25]->Fill(dRmin1_track-dRmin0_track);

    merged_histo[0]->Fill(dRmin0_gsf,dRmin1_gsf);
    merged_histo[1]->Fill(dRmin0_track,dRmin1_track);
    merged_histo[2]->Fill(dRmin0_gsf,dRmin1_gsf-dRmin0_gsf);
    merged_histo[3]->Fill(dRmin0_track,dRmin1_track-dRmin0_track);

    if(dRmin0_track<0.1 && dRmin1_track<0.1){
      T1_conv.push_back(track_general.at(ind0_track));
      T2_conv.push_back(track_general.at(ind1_track));
      conv_track.push_back(conversion.at(i));
    }

  }

  sort(0);

  for(unsigned int i=0;i<T1_conv.size();i++){
    float dRmin = 9999;
    for(unsigned int j=0;j<gsf_elec.size();j++){
      float temp = deltaR(gsf_elec.at(j).v.Eta(),gsf_elec.at(j).v.Phi(),T1_conv.at(i).eta,T1_conv.at(i).phi);
      if(temp<dRmin) dRmin = temp;
    }
    conversion_histo[26]->Fill(dRmin);
  }

  for(unsigned int i=0;i<T2_conv.size();i++){
    float dRmin = 9999;
    for(unsigned int j=0;j<gsf_elec.size();j++){
      float temp = deltaR(gsf_elec.at(j).v.Eta(),gsf_elec.at(j).v.Phi(),T2_conv.at(i).eta,T2_conv.at(i).phi);
      if(temp<dRmin) dRmin = temp;
    }
    conversion_histo[27]->Fill(dRmin);
  }

  for(unsigned int i=0;i<conv_track.size();i++){
       
    float pT2 = (T1_conv.at(i).px+T2_conv.at(i).px)*(T1_conv.at(i).px+T2_conv.at(i).px)
      + (T1_conv.at(i).py+T2_conv.at(i).py)*(T1_conv.at(i).py+T2_conv.at(i).py);
    float M = (T1_conv.at(i).v+T2_conv.at(i).v).M();
    conversion_histo[28]->Fill(T1_conv.at(i).pt);
    conversion_histo[29]->Fill(T2_conv.at(i).pt);
    conversion_histo[30]->Fill(TMath::Sqrt(pT2));
    conversion_histo[31]->Fill(conv_track.at(i).v.Pt());
    conversion_histo[32]->Fill(TMath::Sqrt(pT2)-conv_track.at(i).v.Pt());
    conversion_histo[33]->Fill((TMath::Sqrt(pT2)-conv_track.at(i).v.Pt())/conv_track.at(i).v.Pt());
    conversion_histo[34]->Fill((TMath::Sqrt(pT2)-conv_track.at(i).v.Pt())/TMath::Sqrt(pT2));
    conversion_histo[35]->Fill(M);
    conversion_histo[36]->Fill(conv_track.at(i).pairInvariantMass);
 
  }


  for(const reco::SuperCluster &ie : *mySc)
    {
      Lepton temp;
      temp.energy = ie.energy();
      temp.eta = ie.eta();
      temp.phi = ie.phi();
      temp.size = ie.size();
          
      sc.push_back(temp);

      conversion_histo[40]->Fill(temp.phi);
      conversion_histo[41]->Fill(temp.eta);

    }

  for(const reco::CaloCluster &ie : *myCalo)
    {
      Lepton temp;
      temp.energy = ie.energy();
      temp.eta = ie.eta();
      temp.phi = ie.phi();
      //temp.size = ie.size();

      CaloCluster.push_back(temp);
      gsf_e_histo[65]->Fill(temp.eta,temp.energy/10);
     
      //conversion_histo[40]->Fill(temp.phi);
      //conversion_histo[41]->Fill(temp.eta);

    }

  for(const reco::CaloCluster &ie : *myCalo_eb)
    {
      Lepton temp;
      temp.energy = ie.energy();
      temp.eta = ie.eta();
      temp.phi = ie.phi();
      //temp.size = ie.size();                                                                                                          

      CaloCluster_EB.push_back(temp);
      calo_histo[0]->Fill(temp.eta);

      //conversion_histo[40]->Fill(temp.phi);                                                                                           
      //conversion_histo[41]->Fill(temp.eta);                                                                                           

    }

  for(const reco::CaloCluster &ie : *myCalo_ee)
    {
      Lepton temp;
      temp.energy = ie.energy();
      temp.eta = ie.eta();
      temp.phi = ie.phi();
      //temp.size = ie.size();              
      CaloCluster_EB.push_back(temp);
      calo_histo[1]->Fill(temp.eta);

      //conversion_histo[40]->Fill(temp.phi);       
      //conversion_histo[41]->Fill(temp.eta);                                                                      
                                                    

    }

  int nrec = 0;
  if(nevent == 50){
    for(const EcalRecHit &ie : *myrechit)
      {
	
	EBDetId a(ie.id());
	file<<ie.energy()<<" "<<a.ieta()<<" "<<a.iphi()<<endl;
	
	if(fabs(a.ieta()-85)<=20 && fabs(a.iphi()-299)<=20){
	  file1<<ie.energy()<<" "<<(a.ieta()-85)<<" "<<(a.iphi()-299)<<endl;
	}
      }
  }
  
  file.close();
  file1.close();
  
  
  if(nevent == 53){
    for(const EcalRecHit &ie : *myrechit)
      {

	EBDetId a(ie.id());

	
        if(fabs(a.ieta()-60)<=20 && fabs(a.iphi()-186)<=20){
          file2<<ie.energy()<<" "<<(a.ieta()-60)<<" "<<(a.iphi()-186)<<endl;
        }
      }
  }
  
  file2.close();
  
  
  for(const reco::Photon &ie : *myphoton)
    {

      
      photon_histo[0]->Fill(ie.caloPosition().Eta());
      photon_histo[1]->Fill(ie.caloPosition().Phi());
      
      if(ie.isEB() && ie.hasConversionTracks()){
	photon_histo[2]->Fill(ie.r9());
	DetId id = (ie.superCluster())->seed()->seed();
	EBDetId a(id);
	cout<<"Has COnversion Tracks: "<<a.ieta()<<" "<<a.iphi()<<" "<<nevent<<endl;
      }
      
      else if(ie.isEB() && !(ie.hasConversionTracks())){
	photon_histo[3]->Fill(ie.r9());
	DetId id = (ie.superCluster())->seed()->seed();
        EBDetId a(id);
	cout<<"Has No COnversion Tracks: "<<a.ieta()<<" "<<a.iphi()<<" "<<nevent<<endl;
      }
    }
  
  
  
  
  
  for(unsigned int i=0;i<gsf_elec.size();i++){
    float dRmin = 9999;
    int ind;
    for(unsigned int j=0;j<CaloCluster.size();j++){
      float temp = deltaR(gsf_elec.at(i).v.Eta(),gsf_elec.at(i).v.Phi(),CaloCluster.at(j).eta,CaloCluster.at(j).phi);
      if(temp < dRmin){
	dRmin = temp;
	ind = j;
      }
    }
    if(dRmin<0.4){
      double eta_calo = CaloCluster.at(ind).eta;
      double eta_gsf = gsf_elec.at(i).v.Eta();
      double phi_calo = CaloCluster.at(ind).phi;
      double phi_gsf = gsf_elec.at(i).v.Phi();
	    
	    
      gr_gsf->SetPoint(gr_gsf->GetN(),eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf),CaloCluster.at(ind).energy);
	 
      // mustache_histo[0]->Fill(CaloCluster.at(ind).eta-goodGsf_elec.at(i).v.Eta(),delta_phi(CaloCluster.at(ind).phi,goodGsf_elec.at(
      // ).v.Phi(),CaloCluster.at(ind).energy);                                     
	 
    }
  }

	
	
  for(unsigned int i=0;i<goodIsoGsf_elec.size();i++){
    float dRmin = 9999;
    int ind;
    for(unsigned int j=0;j<CaloCluster.size();j++){
      float temp = deltaR(goodIsoGsf_elec.at(i).v.Eta(),goodIsoGsf_elec.at(i).v.Phi(),CaloCluster.at(j).eta,CaloCluster.at(j).phi);
      if(temp < dRmin){
	dRmin = temp;
	ind = j;
      }
    }
    if(dRmin<0.4){
      double eta_calo = CaloCluster.at(ind).eta;
      double eta_gsf = goodIsoGsf_elec.at(i).v.Eta();
      double phi_calo = CaloCluster.at(ind).phi;
      double phi_gsf = goodIsoGsf_elec.at(i).v.Phi();
    
         
      gr_good->SetPoint(gr_good->GetN(),eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf),CaloCluster.at(ind).energy);
      if(goodIsoGsf_elec.at(i).v.Pt()<10)  merged_histo[4]->Fill(eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf));
      if(goodIsoGsf_elec.at(i).v.Pt()>10 && goodIsoGsf_elec.at(i).v.Pt()<20)  merged_histo[5]->Fill(eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf));
      if(goodIsoGsf_elec.at(i).v.Pt()>20 && goodIsoGsf_elec.at(i).v.Pt()<50)  merged_histo[6]->Fill(eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf));
      if(goodIsoGsf_elec.at(i).v.Pt()>100)  merged_histo[7]->Fill(eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf));


      // mustache_histo[0]->Fill(CaloCluster.at(ind).eta-goodGsf_elec.at(i).v.Eta(),delta_phi(CaloCluster.at(ind).phi,goodGsf_elec.at(i).v.Phi(),CaloCluster.at(ind).energy);

    }
  }
  
  for(unsigned int i=0;i<badGsf_elec.size();i++){
    float dRmin = 9999;
    int ind;
    for(unsigned int j=0;j<CaloCluster.size();j++){
      float temp = deltaR(badGsf_elec.at(i).v.Eta(),badGsf_elec.at(i).v.Phi(),CaloCluster.at(j).eta,CaloCluster.at(j).phi);
      if(temp < dRmin){
	dRmin = temp;
	ind =j;
      }
    }
    if(dRmin<0.1){
         
      double eta_calo = CaloCluster.at(ind).eta;
      double eta_gsf = badGsf_elec.at(i).v.Eta();
      double phi_calo = CaloCluster.at(ind).phi;
      double phi_gsf = badGsf_elec.at(i).v.Phi();


      gr_bad->SetPoint(gr_bad->GetN(),eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf),CaloCluster.at(ind).energy);
      if(badGsf_elec.at(i).v.Pt()<10)  merged_histo[8]->Fill(eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf));
      if(badGsf_elec.at(i).v.Pt()>10 && badGsf_elec.at(i).v.Pt()<20)  merged_histo[9]->Fill(eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf));
      if(badGsf_elec.at(i).v.Pt()>20 && badGsf_elec.at(i).v.Pt()<50)  merged_histo[10]->Fill(eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf));
      if(badGsf_elec.at(i).v.Pt()>100)  merged_histo[11]->Fill(eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf));


      // mustache_histo[1]->Fill(CaloCluster.at(ind).eta-badGsf_elec.at(i).v.Eta(),CaloCluster.at(ind).phi-badGsf_elec.at(i).v.Phi \
      //			 (),CaloCluster.at(ind).energy);
    }
  }
    
  for(unsigned int i=0;i<goodIsoGsf_elec.size();i++) nevent_goodIso.push_back(nevent);

  for(unsigned int i=0;i<nevent_goodIso.size();i++) gsf_e_histo[67]->Fill(nevent_goodIso.at(i));

  gsf_e_histo[66]->Fill(goodIsoGsf_elec.size());
     
  if(nevent == 237){
    for(unsigned int i=0;i<goodIsoGsf_elec.size();i++){
      // float dRmin = 9999;
      // int ind;
      for(unsigned int j=0;j<CaloCluster.size();j++){
	float dR = deltaR(goodIsoGsf_elec.at(i).v.Eta(),goodIsoGsf_elec.at(i).v.Phi(),CaloCluster.at(j).eta,CaloCluster.at(j).phi);
	if(dR<0.4){
	  double eta_calo = CaloCluster.at(j).eta;
	  double eta_gsf = goodIsoGsf_elec.at(i).v.Eta();
	  double phi_calo = CaloCluster.at(j).phi;
	  double phi_gsf = goodIsoGsf_elec.at(i).v.Phi();
	  
	  merged_histo[12]->Fill(eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf));
	}
      }
    }
  }
  
   
  //Calculating number of CaloCLuster within dR<0.4
  if(goodIsoGsf_elec.size()>0){
    int ncalo=0;
    for(unsigned int i=0;i<CaloCluster.size();i++){
      float dR = deltaR(goodIsoGsf_elec.at(0).v.Eta(),goodIsoGsf_elec.at(0).v.Phi(),CaloCluster.at(i).eta,CaloCluster.at(i).phi);
      if(dR<0.4){
	ncalo++;
	// merged_histo[12]->Fill(eta_calo-eta_gsf,delta_phi_sign(phi_calo,phi_gsf));
      }
    }
    calo_histo[3]->Fill(ncalo);
  }

  if(goodIsoGsf_elec.size()>0){
    float dRmin0 = 9999;
    float dRmin1 = 9999;
    for(unsigned int i=0;i<CaloCluster.size();i++){
      float dR = deltaR(goodIsoGsf_elec.at(0).v.Eta(),goodIsoGsf_elec.at(0).v.Phi(),CaloCluster.at(i).eta,CaloCluster.at(i).phi);
      if(dR < dRmin0){
	dRmin1 = dRmin0;
	dRmin0 = dR;
      }
      if(dR < dRmin1 && dR > dRmin0){
	dRmin1 = dR;
      }
    }
    calo_histo[4]->Fill(dRmin0+dRmin1);
  }

  if(goodIsoGsf_elec.size()>0){
    float sum = 0;
    float E_elec = goodIsoGsf_elec.at(0).v.E();
    for(unsigned int i=0;i<CaloCluster.size();i++){
      float dR = deltaR(goodIsoGsf_elec.at(0).v.Eta(),goodIsoGsf_elec.at(0).v.Phi(),CaloCluster.at(i).eta,CaloCluster.at(i).phi);
      float E = CaloCluster.at(i).energy;
      sum = sum + (E/dR);
    }
       
    calo_histo[5]->Fill(sum/E_elec);
  }
	 
 

  for(unsigned int i=0;i<goodIsoGsf_elec.size();i++){
    gsf_e_histo[68]->Fill(goodIsoGsf_elec.at(i).v.Eta());
    float sum_dR = 0;
    float sum_EdR = 0;
    unsigned int ncalo = 0;
    for(unsigned int j=0;j<CaloCluster.size();j++){
      float eta1 = goodIsoGsf_elec.at(i).v.Eta();
      float phi1 = goodIsoGsf_elec.at(i).v.Phi();
      float eta2 = CaloCluster.at(j).eta;
      float phi2 = CaloCluster.at(j).phi;
      float dR = deltaR(eta1,phi1,eta2,phi2);
      float E = CaloCluster.at(j).energy;
      if(dR<0.2){
	sum_dR = sum_dR + dR*dR;
	sum_EdR = sum_EdR + E*dR*dR;
	ncalo++;
      }
    }
    gsf_e_histo[57]->Fill(sum_dR/ncalo);
    gsf_e_histo[58]->Fill(sum_dR);
    gsf_e_histo[59]->Fill(sum_EdR);
    gsf_e_histo[60]->Fill(sum_EdR/ncalo);
  }


  for(unsigned int i=0;i<badGsf_elec.size();i++){
    float sum_dR = 0;
    float sum_EdR = 0;
    unsigned int ncalo = 0;
    for(unsigned int j=0;j<CaloCluster.size();j++){
      float eta1 = badGsf_elec.at(i).v.Eta();
      float phi1 = badGsf_elec.at(i).v.Phi();
      float eta2 = CaloCluster.at(j).eta;
      float phi2 = CaloCluster.at(j).phi;
      float dR = deltaR(eta1,phi1,eta2,phi2);
      float E = CaloCluster.at(j).energy;
      if(dR<0.2){
	sum_dR = sum_dR + dR*dR;
	sum_EdR = sum_EdR + E*dR*dR;
	ncalo++;
      }
    }
    gsf_e_histo[61]->Fill(sum_dR/ncalo);
    gsf_e_histo[62]->Fill(sum_dR);
    gsf_e_histo[63]->Fill(sum_EdR);
    gsf_e_histo[64]->Fill(sum_EdR/ncalo);
  }

  for(unsigned int i=0;i<goodIsoGsf_elec.size();i++){
    gr_good_ecal->SetPoint(gr_good_ecal->GetN(),goodIsoGsf_elec.at(i).v.Eta(),goodIsoGsf_elec.at(i).v.Phi(),goodIsoGsf_elec.at(i).ecalEnergy);
  }

  for(unsigned int i=0;i<badGsf_elec.size();i++){
    gr_bad_ecal->SetPoint(gr_bad_ecal->GetN(),badGsf_elec.at(i).v.Eta(),badGsf_elec.at(i).v.Phi(),badGsf_elec.at(i).ecalEnergy);
  }
  for(unsigned int i=0;i<T1_conv.size();i++){
    float dRmin = 9999;
    for(unsigned int j=0;j<sc.size();j++){
      float temp = deltaR(T1_conv.at(i).eta,T1_conv.at(i).phi,sc.at(j).eta,sc.at(j).phi);
      if(temp < dRmin) dRmin = temp;
    }
    conversion_histo[37]->Fill(dRmin);
  }

  for(unsigned int i=0;i<T2_conv.size();i++){
    float dRmin = 9999;
    for(unsigned int j=0;j<sc.size();j++){
      float temp = deltaR(T2_conv.at(i).eta,T2_conv.at(i).phi,sc.at(j).eta,sc.at(j).phi);
      if(temp < dRmin) dRmin = temp;
    }
    conversion_histo[38]->Fill(dRmin);
  }

  for(unsigned int i=0;i<conv_track.size();i++){
    float dRmin = 9999;
    for(unsigned int j=0;j<sc.size();j++){
      float temp = deltaR(conv_track.at(i).v.Eta(),conv_track.at(i).v.Phi(),sc.at(j).eta,sc.at(j).phi);
      if(temp < dRmin) dRmin = temp;
    }
    conversion_histo[39]->Fill(dRmin);
  }

     

    



#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void Analyzer_egamma::beginJob() {
  // please remove this method if not needed
  nevent = 0;

}

// ------------ method called once each job just after ending the event loop  ------------
void Analyzer_egamma::endJob() {
  // please remove this method if not needed
  std::cout<<"Total Events:"<<nevent<<std::endl;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Analyzer_egamma::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

void Analyzer_egamma::sort(int opt)
{
  if(opt==0){
    for(int i=0;i<(int)gsf_elec.size()-1;i++){
      for(int j=i+1;j<(int)gsf_elec.size();j++){
	if(gsf_elec[i].v.Pt()<gsf_elec[j].v.Pt()) swap(gsf_elec.at(i),gsf_elec.at(j));
      }
    }
    for(int i=0;i<(int)conversion.size()-1;i++){
      for(int j=i+1;j<(int)conversion.size();j++){
	if(conversion[i].v.Pt()<conversion[j].v.Pt()) swap(conversion.at(i),conversion.at(j));
      }
    }

    for(int i=0;i<(int)goodConversion.size()-1;i++){
      for(int j=i+1;j<(int)goodConversion.size();j++){
	if(goodConversion[i].v.Pt()<goodConversion[j].v.Pt()) swap(goodConversion.at(i),goodConversion.at(j));
      }
    }
    for(int i=0;i<(int)goodGsf_elec.size()-1;i++){
      for(int j=i+1;j<(int)goodGsf_elec.size();j++){
	if(goodGsf_elec[i].v.Pt()<goodGsf_elec[j].v.Pt()) swap(goodGsf_elec.at(i),goodGsf_elec.at(j));
      }
    }
    for(int i=0;i<(int)goodIsoGsf_elec.size()-1;i++){
      for(int j=i+1;j<(int)goodIsoGsf_elec.size();j++){
	if(goodIsoGsf_elec[i].v.Pt()<goodIsoGsf_elec[j].v.Pt()) swap(goodIsoGsf_elec.at(i),goodIsoGsf_elec.at(j));
      }
    }

    for(int i=0;i<(int)badGsf_elec.size()-1;i++){
      for(int j=i+1;j<(int)badGsf_elec.size();j++){
	if(badGsf_elec[i].v.Pt()<badGsf_elec[j].v.Pt()) swap(badGsf_elec.at(i),badGsf_elec.at(j));
      }
    }
    for(int i=0;i<(int)track_general.size()-1;i++){
      for(int j=i+1;j<(int)track_general.size();j++){
	if(track_general[i].pt<track_general[j].pt) swap(track_general.at(i),track_general.at(j));
      }
    }
    for(int i=0;i<(int)track_displaced.size()-1;i++){
      for(int j=i+1;j<(int)track_displaced.size();j++){
	if(track_displaced[i].pt<track_displaced[j].pt) swap(track_displaced.at(i),track_displaced.at(j));
      }
    }

    for(int i=0;i<(int)GsfTrack.size()-1;i++){
      for(int j=i+1;j<(int)GsfTrack.size();j++){
	if(GsfTrack[i].ptMode<GsfTrack[j].ptMode) swap(GsfTrack.at(i),GsfTrack.at(j));
      }
    }
    for(int i=0;i<(int)T1_conv.size()-1;i++){
      for(int j=i+1;j<(int)T1_conv.size();j++){
	if(T1_conv[i].pt<T1_conv[j].pt) swap(T1_conv.at(i),T1_conv.at(j));
      }
    }

    for(int i=0;i<(int)T2_conv.size()-1;i++){
      for(int j=i+1;j<(int)T2_conv.size();j++){
	if(T2_conv[i].pt<T2_conv[j].pt) swap(T2_conv.at(i),T2_conv.at(j));
      }
    }
    for(int i=0;i<(int)conv_track.size()-1;i++){
      for(int j=i+1;j<(int)conv_track.size();j++){
	if(conv_track[i].v.Pt()<conv_track[j].v.Pt()) swap(conv_track.at(i),conv_track.at(j));
      }
    }


  }
}

float Analyzer_egamma::delta_phi(float phi1, float phi2)
{
  //The correct deltaPhi falls in the interval [0 , pi]
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}


float Analyzer_egamma::delta_phi_sign(float phi1, float phi2)
{
  float dphi = phi1-phi2;
  if(dphi > TMath::Pi()) dphi = 2*TMath::Pi()-fabs(dphi);
  else if(dphi < -TMath::Pi()) dphi = fabs(dphi)-2*TMath::Pi();
  return dphi;
}

float Analyzer_egamma::deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float dphi = Analyzer_egamma::delta_phi(phi1,phi2);
  float deta = fabs(eta1-eta2);
  
  float dr = TMath::Sqrt(dphi*dphi+deta*deta);
  return dr;
}
  
  
//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer_egamma);
