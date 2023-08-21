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
#include "TGraph.h"

//Vector Headers
#include<vector>
#include "TVector3.h"

//EcalRecHit class

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

// Genparticle header File

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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
  float E;
  float dR;
  TVector3 v1;
  int id; int ind; int momid;
  int dauid;
  float wt;

  //GSF Electron.h variables 
  float scPixCharge;
  bool isGsfCtfScPixChargeConsistent;
  bool isGsfScPixChargeConsistent;
  float shFracInnerHits;
  float correctedEcalEnergy;
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
  float rawEnergy;
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
  double dxy;
  double dz;
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
  
  //ecalrechit collection

  int ieta;
  int iphi;

  //Photon variables

  int seedieta;
  int seediphi;
  bool hasConversionTracks;
  bool isEB,isEE,isES;


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
std::vector<Lepton> ecalrechit_EB;
std::vector<Lepton> Photon;
std::vector<Lepton> Two_Gsfelec;
std::vector<Lepton> One_Gsfelec;
std::vector<Lepton> genElec;
std::vector<Lepton> genElec_momJPsi;
std::vector<Lepton> genElec_mom0_isJPsi;
std::vector<Lepton> genElec_momn_isJPsi;
std::vector<Lepton> genJPsi;
std::vector<Lepton> genJPsi_dau_iselec;
std::vector<Lepton> recoJPsi;
std::vector<Lepton> genJPsi_to_ee;
std::vector<identity> convtype;
std::vector<int> nevent_goodIso;
std::vector<int> nevent_convertphotons;
std::vector<int> nevent_photons;
std::vector<int> nevent_unconvertphotons;

//Text File Declaratiom
std::fstream file_gsfrechit_EB("Gsf_elec_EBrechit.txt",ios::app);
std::fstream file_photonrechit_unconvert_EB("Photon_unconvert_EBrechit.txt",ios::app);
std::fstream file_photonrechit_convert_EB("Photon_convert_EBrechit.txt",ios::app);

/*std::fstream file_convert_photons("Converted_Photons.txt",ios::app);
std::fstream file_unconvert_photons("Unconverted_Photons.txt",ios::app);
std::fstream file_twogsf_EB("Two_Gsf_rechit_EB.txt",ios::app);
std::fstream file_onegsf_EB("One_Gsf_rechit_EB.txt",ios::app);
std::fstream file_delta_seedietaiphi_twogsf_EB_eta_05("TwoGsf_delta_seedietaiphi_EB_eta_05.txt",ios::app);
std::fstream file_delta_seedietaiphi_twogsf_EB_eta_05_1("TwoGsf_delta_seedietaiphi_EB_05_1.txt",ios::app);
std::fstream file_delta_seedietaiphi_twogsf_EB_eta_1_2("TwoGsf_delta_seedietaiphi_EB_1_2.txt",ios::app);
std::fstream file_delta_seedietaiphi_twogsf_EB_eta_2_("TwoGsf_delta_seedietaiphi_EB_2_.txt",ios::app);
std::fstream file_eff("Efficiency_vs_dR",ios::app);
*/
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
  int den[10]={0};
  int num[10]={0};
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
  edm::EDGetTokenT<vector<reco::GenParticle> > gen_token_;
  TH1D *gen_histo[100];
  TH1D *photon_histo[50];
  TH1D *calo_histo[50];
  TH2D *JPsi_histo[50];
  TGraph2D *gr_good,*gr_bad,*gr_gsf,*gr_good_ecal,*gr_bad_ecal;
  TGraph *gr_JPsi;
  TH3D *mustache_histo[10];
  TH1D *var_histo[10];
  
  
  
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
  photon_token_(consumes<vector<reco::Photon> >(iConfig.getParameter<edm::InputTag>("photon_token"))),
  gen_token_(consumes<vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("gen_token")))
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
  gsf_e_histo[41] = fs->make<TH1D>("flag for conversion-rejection ged gsf","flag",200,-100,100);
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
  gsf_e_histo[69] = fs->make<TH1D>("ecal_rechit_gsf","ecal_rechit_gsf",2560,0,256);
  gsf_e_histo[70] = fs->make<TH1D>("di-electron invariant mass","di_electron_invm",1000,0,100);
  gsf_e_histo[71] = fs->make<TH1D>("dR of ee from recoJPsi","dR_JPsi",1000,0,1);
  gsf_e_histo[72] = fs->make<TH1D>("pT of recoJPsi if 2.6<Mee<3.4","pT of recoJPsi if 2.6<Mee<3.4",1000,0,1000);
  gsf_e_histo[73] = fs->make<TH1D>("JPsi_M_dR005","JPsi_M_dR005",200,2,4);
  gsf_e_histo[74] = fs->make<TH1D>("JPsi_M_dR005to01","JPsi_M_dR005to01",200,2,4);
  gsf_e_histo[75] = fs->make<TH1D>("JPsi_M_dR01to015","JPsi_M_dR01to015",200,2,4);
  gsf_e_histo[76] = fs->make<TH1D>("JPsi_M_dR015to02","JPsi_M_dR015to02",200,2,4);
  gsf_e_histo[77] = fs->make<TH1D>("JPsi_M_dR02toinf","JPsi_M_dR02toinf",200,2,4);
  gsf_e_histo[78] = fs->make<TH1D>("JPsi_M_dR0p03","JPsi_M_dR0p03",200,2,4);
  gsf_e_histo[79] = fs->make<TH1D>("dxy of gsf electrons","dxy of gsf electrons",200,-1,1);
  gsf_e_histo[80] = fs->make<TH1D>("dz of gsf electrons","dz of gsf electrons",2000,-10,10);
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
  photon_histo[4] = fs->make<TH1D>("nevent converted photons","nevent_converted_photons",10000,0,10000);
  photon_histo[5] = fs->make<TH1D>("nevent photons","nevent_photons",10000,0,10000);
  photon_histo[6] = fs->make<TH1D>("nevent unconverted photons","nevent_unconverted_photons",10000,0,10000);

  // 2D histograms for JPsi analysis

  JPsi_histo[0] = fs->make<TH2D>("dR_vs_pT0ofGsf","dR_vs_pT0ofGsf",200,0,200,1000,0,1); 
  JPsi_histo[1] = fs->make<TH2D>("dR_vs_pTofJPsi","dR_vs_pT0ofJPsi",200,0,200,1000,0,1);
  JPsi_histo[2] = fs->make<TH2D>("Gen_dR_vs_pT0","Gen_dR_vs_pT0",200,0,200,1000,0,1);
  JPsi_histo[3] = fs->make<TH2D>("efficiency vs dR for JPsi reconstruction","efficiency vs dR for JPsi reconstruction",100,0,1,100,0,1);
  JPsi_histo[4] = fs->make<TH2D>("dxy0 vs recoJPsi pT","dxy0 vs recoJPsi pT",600,0,600,200,-10,10);
  JPsi_histo[5] = fs->make<TH2D>("dxy1 vs recoJPsi pT","dxy1 vs recoJPsi pT",600,0,600,200,-10,10);
  JPsi_histo[6] = fs->make<TH2D>("dz0 vs recoJPsi pT","dz0 vs recoJPsi pT",600,0,600,200,-10,10);
  JPsi_histo[7] = fs->make<TH2D>("dz1 vs recoJPsi pT","dz1 vs recoJPsi pT",600,0,600,200,-10,10);
  // Histograms for Genparticles

  gen_histo[0] = fs->make<TH1D>("Mee_if_mom_JPsi","Mee_if_mom_isJPsi",1000,0,100);
  gen_histo[1] = fs->make<TH1D>("Mee_genElec","Mee_genElec",1000,0,100);
  gen_histo[2] = fs->make<TH1D>("size of genElec","Size_genElec",100,0,100);
  gen_histo[3] = fs->make<TH1D>("size of genElec_ifMom_isJPsi","Size_genElec_momJPsi",100,0,100);
  gen_histo[4] = fs->make<TH1D>("Number of Mothers","no_of_mothers",100,0,100);
  gen_histo[5] = fs->make<TH1D>("mom_pdgID","mom_pdgID",600,0,600);
  gen_histo[6] = fs->make<TH1D>("last_momID","last_momID",1200,-600,600);
  gen_histo[7] = fs->make<TH1D>("first_momID","first_momID",1200,-600,600);
  gen_histo[8] = fs->make<TH1D>("Mee_if_mom0_isJPsi","Mee_if_mom0_isJPsi",1000,0,100);
  gen_histo[9] = fs->make<TH1D>("Mee_if_momn_isJPsi","Mee_if_momn_isJPsi",1000,0,100);
  gen_histo[10] = fs->make<TH1D>("size of genElec_mom0_isJPsi","size of genElec_mom0_isJPsi",100,0,100);
  gen_histo[11] = fs->make<TH1D>("size of genElec_momn_isJPsi","size of genElec_momn_isJPsi",100,0,100);
  gen_histo[12] = fs->make<TH1D>("JPsi mass","JPsi mass",1000,0,100);
  gen_histo[13] = fs->make<TH1D>("size of genJPsi_dau_iselec","size of genJPsi_dau_iselec",100,0,100);
  gen_histo[14] = fs->make<TH1D>("pdgId of genJsi daughter","pdgId of genJPsi daughter",1200,-600,600);
  gen_histo[15] = fs->make<TH1D>("Mee if dau is JPsi","Mee_if_dau_is_JPsi",1000,0,100);
  gen_histo[16] = fs->make<TH1D>("JPsi mom's Pdg ID","JPsi mom's Pdg ID",1000,0,1000);
  gen_histo[17] = fs->make<TH1D>("JPsi size","JPsi size",100,0,100);
  gen_histo[18] = fs->make<TH1D>("dR gen electrons","dR b/w gen electrons",1000,0,10);
  gen_histo[19] = fs->make<TH1D>("Pt of JPsi","Pt of JPsi",100,0,100);
  gen_histo[20] = fs->make<TH1D>("dRmin_gsfElec_genElec","dRmin_gsfElec_genElec",1000,0,1);
  gen_histo[21] = fs->make<TH1D>("genElec_dR-gsf_elec_dR","genElec_dR-gsf_elec_dR",1000,0,1);
  gen_histo[22] = fs->make<TH1D>("gsf_elec_dR_if_ dR<0.1","gsf_elec_dR_if dR<0.1",1000,0,1);
  gen_histo[23] = fs->make<TH1D>("Gen_elec_dR_if dR<0.1","Gen_elec_dR_if dR<0.1",1000,0,1);
  gen_histo[24] = fs->make<TH1D>("efficiency of JPsi reconstruction","efficiency of JPsi reconstruction",10,0,1);          
  gen_histo[25] = fs->make<TH1D>("dR of ee from JPSi","dR of ee from JPsi",1000,0,1);
  gen_histo[26] = fs->make<TH1D>("gen JPsi to ee size","gen JPsi to ee size",10,0,10);
  gen_histo[27] = fs->make<TH1D>("recoJPsi size","recoJPsi size",10,0,10);
  gen_histo[28] = fs->make<TH1D>("genElec dR if 2.6<Mee<3.4","genElec dR if 2.6<Mee<3.4",100,0,1);
  gen_histo[29] = fs->make<TH1D>("genElec dR if reco matched genElec ee","genElec dR if reco matched genElec ee",100,0,1);
  gen_histo[30] = fs->make<TH1D>("dR ee if genElec","dR b/w ee if genElec_momJPsi",1000,0,10);
  gen_histo[31] = fs->make<TH1D>("dR ee genElec","dR ee if genElec",200,0,2);
  gen_histo[32] = fs->make<TH1D>("pT of genJPsi from genElec_momJPsi","pT of genJPsi from genElec_momJPsi",200,0,200);
  gen_histo[33] = fs->make<TH1D>("pT of recoJPsi","pT of recoJPsi",1000,0,1000);
  gen_histo[34] = fs->make<TH1D>("size of genElec","size of genElec",10,0,10);
  gen_histo[35] = fs->make<TH1D>("Mom PDGID of all Gen electrons","Mom PDGID of all Gen electrons",1200,-600,600);
  gen_histo[36] = fs->make<TH1D>("Mom PDGID of all Gen JPsi","Mom PDGID of all Gen JPsi",1200,-600,600);
  gen_histo[37] = fs->make<TH1D>("pT of gen JPsi","pT of Gen JPsi",1000,0,1000);
  gen_histo[38] = fs->make<TH1D>("dR of Gen JPsi","dR of Gen JPsi",100,0,10);
  gen_histo[39] = fs->make<TH1D>("dphi of Gen JPsi","dphi of Gen JPsi",100,-5,5);
  gen_histo[40] = fs->make<TH1D>("phi of genElec","phi of genElec",100,-5,5);
  gen_histo[41] = fs->make<TH1D>("eta of genElec","eta of genElec",100,-5,5);
  gen_histo[42] = fs->make<TH1D>("phi of genJPsi","phi of genJPsi",100,-5,5);
  gen_histo[43] = fs->make<TH1D>("eta of genJPsi","eta of genJPsi",100,-5,5);
  gen_histo[44] = fs->make<TH1D>("pT of gen elec","pT of gen elec",500,0,500);
  gen_histo[45] = fs->make<TH1D>("pT0 of gen Elec 0.2<dR<0.3","pT0 of gen Elec 0.2<dR<0.3",5000,0,500);
  gen_histo[46] = fs->make<TH1D>("pT1 of gen Elec 0.2<dR<0.3","pT1 of gen Elec 0.2<dR<0.3",5000,0,500);
  gen_histo[47] = fs->make<TH1D>("dR of genElec if matched to gsf","dR of genElec if matched to gsf",1000,0,10);
  gen_histo[48] = fs->make<TH1D>("dR of gsf_elec if matched to genElec","dR gsf_elec if matched to genElec",1000,0,10);
  gen_histo[49] = fs->make<TH1D>("deltaE/Egen for matched gsf to gen","deltaE/Egen for matched gsf to gen",100,0,10);
  gen_histo[50] = fs->make<TH1D>("deltaE_raw/Egen for matched gsf to gen","deltaE_raw/Egen for matched gsf to gen",100,0,10);
  //gen_histo[51] = fs->make<TH1D>("dR genee if mom is Z","dR gen_momZ ee matched to gsf",1000,0,10);
  //Histogram with variable Binning
  
  float pTbins0[7] = {0,5,10,15,20,30,2000};

  var_histo[0] = fs->make<TH1D>("pT of ee from genElec_momJPsi var bin","pT of ee from genElec_momJPsi var bin",6,pTbins0);
  var_histo[1] = fs->make<TH1D>("pT of recoJPsi var bin","pT of recJPsi var bin",6,pTbins0);


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
  gr_JPsi = fs->make<TGraph>();
  gr_JPsi->SetTitle("Efficiency_vs_dR;dR;eff");



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
  ecalrechit_EB.clear();
  Photon.clear();
  nevent_convertphotons.clear();
  nevent_photons.clear();
  nevent_unconvertphotons.clear();
  Two_Gsfelec.clear();
  One_Gsfelec.clear();
  genElec_momJPsi.clear();
  genElec.clear();
  genElec_mom0_isJPsi.clear();
  genElec_momn_isJPsi.clear();
  genJPsi.clear();
  genJPsi_dau_iselec.clear();
  recoJPsi.clear();
  genJPsi_to_ee.clear();
  //genZ_to_ee.clear();
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
  Handle<vector<reco::GenParticle> > mygenparticle;
  iEvent.getByToken(gen_token_,mygenparticle);

  

  nevent++;
  //Define Temporary Variables

  //EGamma temp;
  //int n_e;

  // Lepton temp;
  int n_isgsf,n_gsf;
  
  n_gsf=0;
  n_isgsf=0;
  
  for(const reco::GsfElectron &ie : *mygsf_elec_ged)
    {
      Lepton temp;
      temp.scPixCharge = ie.scPixCharge();
      temp.id = ie.pdgId();
      temp.v.SetPtEtaPhiM(ie.pt(),ie.eta(),ie.phi(),0.000511);
      temp.correctedEcalEnergy = ie.correctedEcalEnergy();
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
      temp.isEB = ie.isEB();
      temp.isEE = ie.isEE();
 
      DetId id = (ie.superCluster())->seed()->seed();
      EBDetId a(id);
      temp.seedieta = a.ieta();
      temp.seediphi = a.iphi();

      temp.rawEnergy = (ie.superCluster())->rawEnergy();

      temp.dxy = (ie.gsfTrack())->dxy();
      temp.dz  = (ie.gsfTrack())->dz();      
      if(abs(ie.eta()<0.8)) temp.fbrem_eta_8 = ie.fbrem();
      if(abs(ie.eta())>0.8 && abs(ie.eta())<1.44) temp.fbrem_eta_144 = ie.fbrem();
      if(abs(ie.eta())>1.57 && abs(ie.eta())<2) temp.fbrem_eta_2 = ie.fbrem();
      if(abs(ie.eta())>2) temp.fbrem_eta_2_ = ie.fbrem();
      if(ie.pt()>3) gsf_elec.push_back(temp);
      
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

  if((int)gsf_elec.size()>1){
    float Mee = (gsf_elec.at(0).v+gsf_elec.at(1).v).M();
    gsf_e_histo[70]->Fill(Mee);
  }

  for(const reco::GsfElectron &ie : *mygsf_elec_unonly)
    {
      Lepton temp;
      temp.scPixCharge = ie.scPixCharge();
      temp.id = ie.pdgId();
      temp.v.SetPtEtaPhiM(ie.pt(),ie.eta(),ie.phi(),0.000511);


      if(ie.isElectron())
	{
	  n_isgsf++;
	}
      n_gsf++;

      gsf_e_histo[4]->Fill(temp.scPixCharge);
      gsf_e_histo[7]->Fill(ie.pdgId());
    }                                                                                                                                                                      
  gsf_e_histo[5]->Fill(n_isgsf);
  gsf_e_histo[6]->Fill(n_gsf);
                                                                                                                                                                       
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

  
  sort(0);
        
  /********************************************************************************************Conversion*********************************************************************************/

  for(unsigned int i=0;i<goodConversion.size();i++){

    conversion_histo[13]->Fill(goodConversion.at(i).v.Pt());
  }
     
  if(goodConversion.size()>0)  conversion_histo[14]->Fill(goodConversion.at(0).v.Pt());

  /***********************************************************************************************goodGsf***********************************************************************************/ 
  for(unsigned int i=0;i<goodGsf_elec.size();i++){

    gsf_e_histo[44]->Fill(goodGsf_elec.at(i).v.Pt());
  }
    
  if(goodGsf_elec.size()>0) gsf_e_histo[45]->Fill(goodGsf_elec.at(0).v.Pt());

  /***********************************************************************************************badGsf**************************************************************************************/
  for(unsigned int i=0;i<badGsf_elec.size();i++){

    gsf_e_histo[46]->Fill(badGsf_elec.at(i).v.Pt());
  }

  if(badGsf_elec.size()>0) gsf_e_histo[47]->Fill(badGsf_elec.at(0).v.Pt());

  //Storing the size of objects
  gsf_e_histo[48]->Fill(goodGsf_elec.size());
  gsf_e_histo[49]->Fill(badGsf_elec.size());
  conversion_histo[16]->Fill(goodConversion.size());
  
  /*************************************************************************Shower shape variables for good and bad gsf*****************************************************************************/

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


  /****************************************************************************************************************************************************************************************************

                                                                                GoodGsf and Track matching. Identifying T1 and T2
										

  ******************************************************************************************************************************************************************************************************/



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
 
    }


  /******************************************************************************************************************************************************************************************************

                                                                   ECAL Rechit and matching it With Photon and Gsf Electron Collection
										   
  ******************************************************************************************************************************************************************************************************/



  //Converted and Unconverted loop starts
  
  for(const EcalRecHit &ie : *myrechit)
    {
      Lepton temp;
      EBDetId a(ie.id());
      temp.ieta = a.ieta();
      temp.iphi = a.iphi();
      temp.energy = ie.energy();
      ecalrechit_EB.push_back(temp);
      
    }


  for(const reco::Photon &ie : *myphoton)
    {
      Lepton temp;
      temp.hasConversionTracks = ie.hasConversionTracks();
      temp.isEB = ie.isEB();
      temp.isEE = ie.isEE();
      temp.r9 = ie.r9();
      
      DetId id = (ie.superCluster())->seed()->seed();
      EBDetId a(id);
      temp.seedieta = a.ieta();
      temp.seediphi = a.iphi();
      
      Photon.push_back(temp);
      
      photon_histo[0]->Fill(ie.caloPosition().Eta());
      photon_histo[1]->Fill(ie.caloPosition().Phi());
    }

  
  for(unsigned int i=0;i<Photon.size();i++){
    if(Photon.at(i).isEB && Photon.at(i).hasConversionTracks){
      nevent_convertphotons.push_back(nevent);
      photon_histo[4]->Fill(nevent);
      photon_histo[2]->Fill(Photon.at(i).r9); 
    }
  }
  
  for(unsigned int i=0;i<Photon.size();i++){
    if(Photon.at(i).isEB && !(Photon.at(i).hasConversionTracks)){
      nevent_unconvertphotons.push_back(nevent);
      photon_histo[6]->Fill(nevent);
      photon_histo[3]->Fill(Photon.at(i).r9);
    }
  }
  
  for(unsigned int i=0;i<Photon.size();i++){
    if(Photon.at(i).isEB){
      nevent_photons.push_back(nevent);
      photon_histo[5]->Fill(nevent);
    }
  }

  for(unsigned int i=0;i<Photon.size();i++){
    int seedieta = Photon.at(i).seedieta;
    int seediphi = Photon.at(i).seediphi;
    for(unsigned int j=0;j<ecalrechit_EB.size();j++){
      int ieta = ecalrechit_EB.at(j).ieta;
      int iphi = ecalrechit_EB.at(j).iphi;
      int dieta = ieta-seedieta;
      int diphi = iphi-seediphi;
      float energy = ecalrechit_EB.at(j).energy;
      bool isConversion = fabs(diphi)<=20 && fabs(dieta)<=20 && Photon.at(i).isEB && Photon.at(i).hasConversionTracks;
      bool isUnconversion = fabs(diphi)<=20 && fabs(dieta)<=20 && Photon.at(i).isEB && !(Photon.at(i).hasConversionTracks);
      
      //File Format : Evt number__Photon pos(i)__Energy__dieta__diphi  

      if(isConversion) file_photonrechit_convert_EB<<nevent<<" "<<i<<" "<<energy<<" "<<dieta<<" "<<diphi<<endl;
      if(isUnconversion) file_photonrechit_unconvert_EB<<nevent<<" "<<i<<" "<<energy<<" "<<dieta<<" "<<diphi<<endl;
    }
  }


  for(unsigned int i=0;i<gsf_elec.size();i++){
    int seedieta = gsf_elec.at(i).seedieta;
    int seediphi = gsf_elec.at(i).seediphi;
    for(unsigned int j=0;j<ecalrechit_EB.size();j++){
      int ieta = ecalrechit_EB.at(j).ieta;
      int iphi = ecalrechit_EB.at(j).iphi;
      int dieta = ieta-seedieta;
      int diphi = iphi-seediphi;
      float energy = ecalrechit_EB.at(j).energy;
      bool isGood = fabs(dieta)<=20 && fabs(diphi)<=20;
      
      //File Format : Evt number__Gsf pos(i)__Energy__dieta__diphi  
      if(isGood) file_gsfrechit_EB<<nevent<<" "<<i<<" "<<energy<<" "<<dieta<<" "<<diphi<<endl;

      //if(isConversion) file_photonrechit_convert_EB<<nevent<<" "<<i<<" "<<energy<<" "<<dieta<<" "<<diphi<<endl;
      //if(isUnconversion) file_photonrechit_unconvert_EB<<nevent<<" "<<i<<" "<<energy<<" "<<dieta<<" "<<diphi<<endl;
    }
  }


     
  /****************************************************************************************************************************************************************************************************
 
                                                                    ECAL rechit of 2 GSF electron and 1 GSF electron event and dR(e1,e2)<1

  ****************************************************************************************************************************************************************************************************/  
  
  sort(1); //Sorting based on Energy in ecal

  if(gsf_elec.size()==2){
    float dR_e0e1 = gsf_elec.at(0).v.DeltaR(gsf_elec.at(1).v);
    bool is_two_electron = (dR_e0e1<1);
    bool is_one_electron = (dR_e0e1>1);
    if(is_two_electron){
      Two_Gsfelec.push_back(gsf_elec.at(0));
      Two_Gsfelec.push_back(gsf_elec.at(1));
    }
    if(is_one_electron) One_Gsfelec.push_back(gsf_elec.at(0));
  }
  
  if(gsf_elec.size()==1) One_Gsfelec.push_back(gsf_elec.at(0));

  if(Two_Gsfelec.size()>0){
    for(unsigned int i=0;i<ecalrechit_EB.size();i++){
      int seedieta = Two_Gsfelec.at(0).seedieta;
      int seediphi = Two_Gsfelec.at(0).seediphi;
      int ieta = ecalrechit_EB.at(i).ieta;
      int iphi = ecalrechit_EB.at(i).iphi;
      int dieta = ieta-seedieta;
      int diphi = iphi-seediphi;
      float energy = ecalrechit_EB.at(i).energy;
      
      //if(fabs(dieta)<=20 && fabs(diphi)<=20 && Two_Gsfelec.at(0).isEB) file_twogsf_EB<<nevent<<" "<<energy<<" "<<dieta<<" "<<diphi<<endl;
  }
  }
  
  if(One_Gsfelec.size()>0){
    
    for(unsigned int i=0;i<ecalrechit_EB.size();i++){
      int seedieta = One_Gsfelec.at(0).seedieta;
      int seediphi = One_Gsfelec.at(0).seediphi;
      int ieta = ecalrechit_EB.at(i).ieta;
      int iphi = ecalrechit_EB.at(i).iphi;
      int dieta = ieta-seedieta;
      int diphi = iphi-seediphi;
      float energy = ecalrechit_EB.at(i).energy;
      
      //if(fabs(dieta)<=20 && fabs(diphi)<=20 && One_Gsfelec.at(0).isEB) file_onegsf_EB<<nevent<<" "<<energy<<" "<<dieta<<" "<<diphi<<endl;
    }
  }

  /**************************************************************************** Rechits for all gsf electrons **************************************************************************************/

  for(unsigned int i=0;i<gsf_elec.size();i++){
    for(unsigned int j=0;j<ecalrechit_EB.size();j++){ 
      if(gsf_elec.at(i).isEB) gsf_e_histo[69]->Fill(ecalrechit_EB.at(j).energy);
    }
  }
 
  
  /***************************************************************************** delta(seedieta),delta(seediphi),dR *******************************************************************************/
  
  
  if(Two_Gsfelec.size()>0){
    float dR = Two_Gsfelec.at(0).v.DeltaR(Two_Gsfelec.at(1).v);
    int ieta0 = Two_Gsfelec.at(0).seedieta;
    int ieta1 = Two_Gsfelec.at(1).seedieta;
    int iphi0 = Two_Gsfelec.at(0).seediphi;
    int iphi1 = Two_Gsfelec.at(1).seediphi;
    int dieta = fabs(ieta0-ieta1);
    int diphi = fabs(iphi0-iphi1);
    
    //if(Two_Gsfelec.at(0).v.Eta()<0.5) file_delta_seedietaiphi_twogsf_EB_eta_05<<nevent<<" "<<dieta<<" "<<diphi<<" "<<dR<<endl;
    //if(Two_Gsfelec.at(0).v.Eta()>0.5 && Two_Gsfelec.at(0).v.Eta()<1) file_delta_seedietaiphi_twogsf_EB_eta_05_1<<nevent<<" "<<dieta<<" "<<diphi<<" "<<dR<<endl;
    //if(Two_Gsfelec.at(0).v.Eta()>1 && Two_Gsfelec.at(0).v.Eta()<2) file_delta_seedietaiphi_twogsf_EB_eta_1_2<<nevent<<" "<<dieta<<" "<<diphi<<" "<<dR<<endl;
    //if(Two_Gsfelec.at(0).v.Eta()>2) file_delta_seedietaiphi_twogsf_EB_eta_2_<<nevent<<" "<<dieta<<" "<<diphi<<" "<<dR<<endl;
  }
 
  
  /*************************************************************************************************************************************************************************************************  

                                                                     GSF Electron, Track  and CaloCluster variables Calculation and matching with CaloClusters
								    
  ***************************************************************************************************************************************************************************************************/
  
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
      // (),CaloCluster.at(ind).energy);
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

  /*************************************************************************************************************************************************************************************************
                                
                                                                    Analysis of two Gsf electrons from JPsi Starts to see if we have some really close by electrons
								   
  **************************************************************************************************************************************************************************************************/

  sort(0);

  if(gsf_elec.size()>1){
    Lepton temp;
    float dR = gsf_elec.at(0).v.DeltaR(gsf_elec.at(1).v);
    float M = (gsf_elec.at(0).v+gsf_elec.at(1).v).M();
    float pT0 = gsf_elec.at(0).v.Pt();
    float pT1 = gsf_elec.at(1).v.Pt();
    float pT = (gsf_elec.at(0).v+gsf_elec.at(1).v).Pt();
    float eta = (gsf_elec.at(0).v+gsf_elec.at(1).v).Eta();
    float phi = (gsf_elec.at(0).v+gsf_elec.at(1).v).Phi();
    temp.v.SetPtEtaPhiM(pT,eta,phi,M);
    temp.dR = dR;
    //recoJPsi.push_back(temp);
    bool isOS = (gsf_elec.at(0).id)*(gsf_elec.at(1).id) < 0;
    bool isJPsi = isOS && M>2.6 && M<3.4;
    if(isJPsi){
      gsf_e_histo[71]->Fill(dR);
      gsf_e_histo[72]->Fill(pT);
      JPsi_histo[0]->Fill(pT0,dR);
      JPsi_histo[1]->Fill(pT,dR);
      recoJPsi.push_back(temp);
    }
    if(isJPsi && dR<0.03) gsf_e_histo[78]->Fill(M);
    if(isJPsi && dR<0.05) gsf_e_histo[73]->Fill(M);
    if(isJPsi && dR>0.05 && dR<0.1) gsf_e_histo[74]->Fill(M);
    if(isJPsi && dR>0.1 && dR<0.15) gsf_e_histo[75]->Fill(M);
    if(isJPsi && dR>0.15 && dR<0.2) gsf_e_histo[76]->Fill(M);
    if(isJPsi && dR>0.2) gsf_e_histo[77]->Fill(M);
  }


  /***************************************************************************************************************************************************************************************************

                                                                     Analysis With Generator level info Starts. Checking for electrons with mom JPsi
								     
  ****************************************************************************************************************************************************************************************************/

  
  for(const reco::GenParticle &ie : *mygenparticle)
    {
      Lepton temp;
      temp.v.SetPtEtaPhiM(ie.pt(),ie.eta(),ie.phi(),ie.mass());
      temp.E = ie.energy();
      temp.id = ie.pdgId();
      int id = ie.pdgId();
      if(fabs(ie.pdgId()) == 11 && ie.status()==1 && ie.pt()>3) genElec.push_back(temp);
      if(fabs(id) == 443){
	gen_histo[12]->Fill(ie.mass());
	genJPsi.push_back(temp);
      }
      size_t n = ie.numberOfMothers();
      size_t ndau = ie.numberOfDaughters();
      if(fabs(id) == 443) gen_histo[37]->Fill(ie.pt());

      if(n>0){
	const  reco::Candidate *lastmom = ie.mother(n-1);
	const  reco::Candidate *firstmom = ie.mother(0);
	
	int lastmomid = lastmom->pdgId();
	int firstmomid = firstmom->pdgId();
	if(fabs(firstmomid)==443 && fabs(id)==11 && ie.status()==1) genElec_mom0_isJPsi.push_back(temp);
	if(fabs(lastmomid)==443 && fabs(id)==11 && ie.status()==1) genElec_momn_isJPsi.push_back(temp);
	gen_histo[6]->Fill(lastmomid);
	gen_histo[7]->Fill(firstmomid);
      }
      
      gen_histo[4]->Fill(n);      

      //Daughter Loop
      for(size_t i=0; i<ndau; i++){
	Lepton temp2;
	const reco::Candidate *dau = ie.daughter(i);
	int dauid = dau->pdgId();
	temp2.id = dauid;
	temp2.v.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
	if(dau->status()==1) gen_histo[14]->Fill(dauid);
	if(fabs(id) == 443 && fabs(dauid) == 11 && dau->status() ==1) genJPsi_dau_iselec.push_back(temp2);
      }

      //Mother Loops
      for(size_t i=0; i<n; i++){
	const reco::Candidate *mom = ie.mother(i);
	int momid = fabs(mom->pdgId());
	temp.momid = momid;
	if(fabs(id) == 443) gen_histo[16]->Fill(momid);
	if(fabs(id) == 11) gen_histo[35]->Fill(momid);
	if(fabs(id) == 443) gen_histo[36]->Fill(momid);
	if(momid == 443 && ie.status() == 1 && fabs(ie.pdgId())==11 && ie.pt()>3){
	  genElec_momJPsi.push_back(temp);
	  //break;
	}
      }

      for(size_t i=0;i<n; i++){
	const reco::Candidate *mom = ie.mother(i);
	int momid = fabs(mom->pdgId());
	temp.momid = momid;
	if(fabs(ie.pdgId())==11 && ie.status()==1) gen_histo[5]->Fill(momid);
      }
    }
  
  sort(0);

  gen_histo[34]->Fill(genElec.size());

  /****************************************************************************************************************************************************************************************************
 
                                     Checking Properties of genElec whose mom and daughters are found using different methods. Not much useful!just use for cross-checkings

 *****************************************************************************************************************************************************************************************************/

  if(genElec_momJPsi.size()>1){
    float Mee = (genElec_momJPsi.at(0).v + genElec_momJPsi.at(1).v).M();
    bool isOS = (genElec_momJPsi.at(0).id)*(genElec_momJPsi.at(1).id) < 0;
    if(isOS) gen_histo[0]->Fill(Mee);
  }


  if(genElec.size()>1){
    float Mee = (genElec.at(0).v+genElec.at(1).v).M();
    float dR = genElec.at(0).v.DeltaR(genElec.at(1).v);
    float pT0 = genElec.at(0).v.Pt();
    bool isOS = (genElec.at(0).id)*(genElec.at(1).id) < 0;
    gen_histo[18]->Fill(dR);
    JPsi_histo[2]->Fill(pT0,dR);
    if(isOS) gen_histo[1]->Fill(Mee);
  }

  if(genElec_mom0_isJPsi.size()>1){
    float Mee = (genElec_mom0_isJPsi.at(0).v + genElec_mom0_isJPsi.at(1).v).M();
    bool isOS = (genElec_mom0_isJPsi.at(0).id)*(genElec_mom0_isJPsi.at(1).id) < 0;
    if(isOS) gen_histo[8]->Fill(Mee);
  }
  
  
  if(genElec_momn_isJPsi.size()>1){
    float Mee = (genElec_momn_isJPsi.at(0).v + genElec_momn_isJPsi.at(1).v).M();
    bool isOS = (genElec_momn_isJPsi.at(0).id)*(genElec_momn_isJPsi.at(1).id) < 0;
    if(isOS) gen_histo[9]->Fill(Mee);
  }
 
  if(genJPsi_dau_iselec.size()>1 ){
    float Mee = (genJPsi_dau_iselec.at(0).v + genJPsi_dau_iselec.at(1).v).M();
    bool isOS = (genJPsi_dau_iselec.at(0).id)*(genJPsi_dau_iselec.at(1).id) < 0;
    float dR = genJPsi_dau_iselec.at(0).v.DeltaR(genJPsi_dau_iselec.at(1).v);
    if(isOS) gen_histo[15]->Fill(Mee);
  }
  if(genJPsi.size()>0) gen_histo[19]->Fill(genJPsi.at(0).v.Pt());
  
  if(genJPsi.size()>1){
    float dR = genJPsi.at(0).v.DeltaR(genJPsi.at(1).v);
    float dphi = delta_phi(genJPsi.at(0).v.Phi(),genJPsi.at(1).v.Phi());
    gen_histo[39]->Fill(dphi);
    gen_histo[38]->Fill(dR);
  }

  for(unsigned int i=0;i<genElec.size();i++){
    gen_histo[40]->Fill(genElec.at(i).v.Phi());
    gen_histo[41]->Fill(genElec.at(i).v.Eta());
  }
  for(unsigned int i=0;i<genJPsi.size();i++){
    gen_histo[42]->Fill(genJPsi.at(i).v.Phi());
    gen_histo[43]->Fill(genJPsi.at(i).v.Eta());
  }

  for(unsigned int i=0;i<genElec.size();i++) gen_histo[44]->Fill(genElec.at(i).v.Pt());

  gen_histo[13]->Fill(genJPsi_dau_iselec.size());
  gen_histo[17]->Fill(genJPsi.size());
  gen_histo[2]->Fill(genElec.size());
  gen_histo[3]->Fill(genElec_momJPsi.size());
  gen_histo[10]->Fill(genElec_mom0_isJPsi.size());
  gen_histo[11]->Fill(genElec_momn_isJPsi.size());
  

  /**************************************************************************** RECO::GSF and GenElec matching and plotting **************************************************************************/

  dRmin = 9999;
  for(unsigned int i=0;i<gsf_elec.size();i++){
    for(unsigned int j=0;j<genElec.size();j++){
      float dR = gsf_elec.at(i).v.DeltaR(genElec.at(j).v);
      if(dR<dRmin) dRmin = dR;
    }
    gen_histo[20]->Fill(dRmin);
    if(genElec.size()>1 && gsf_elec.size()>1){
      float dRgen = genElec.at(0).v.DeltaR(genElec.at(1).v);
      float dRgsf = gsf_elec.at(0).v.DeltaR(gsf_elec.at(1).v);
      if(dRmin<0.1){
	gen_histo[21]->Fill(fabs(dRgen-dRgsf));
	gen_histo[22]->Fill(dRgsf);
	gen_histo[23]->Fill(dRgen);
      }	
      else gen_histo[21]->Fill(-1);
    }  
  }


  sort(0);
  /****************************************************************************** Some properties of reco and gen JPsi *******************************************************************************/
  
  if(genElec.size()>1){
    float dR = genElec.at(0).v.DeltaR(genElec.at(1).v);
    gen_histo[31]->Fill(dR);
  }

  //Small exercise to check pT on gen and pT of reco electrons

  if(genElec_momJPsi.size()>1){
    float pT= (genElec_momJPsi.at(0).v+genElec_momJPsi.at(1).v).Pt();
    gen_histo[32]->Fill(pT);
    var_histo[0]->Fill(pT);
  }
  for(unsigned int i=0;i<recoJPsi.size();i++){
    float pT = recoJPsi.at(i).v.Pt();
    gen_histo[33]->Fill(pT);
    var_histo[1]->Fill(pT);
  }


  /****************************************************************** Efficiency of reconstruction of recoJPsi to genJPsi in dR and pT bins*************************************************/

  sort(0);
  
  if(genElec.size()>1){ 
    float dR = genElec.at(0).v.DeltaR(genElec.at(1).v);
    float Mee = (genElec.at(0).v+genElec.at(1).v).M();
    bool isOS = (genElec.at(0).id)*(genElec.at(1).id)<0;
    if(dR<0.05){
      den[0]++;
      if(recoJPsi.size()>0) num[0]++;
    }

    if(dR>0.05 && dR<0.1){
      den[1]++;
      if(recoJPsi.size()>0) num[1]++;
    }

    if(dR>0.1 && dR<0.2){
      den[2]++;
      if(recoJPsi.size()>0) num[2]++;
    }
    if(dR>0.2 && dR<0.3){
      den[3]++;
      gen_histo[45]->Fill(genElec.at(0).v.Pt());
      gen_histo[46]->Fill(genElec.at(1).v.Pt());
      if(recoJPsi.size()>0) num[3]++;
    }
    if(dR>0.3){
      den[4]++;
      if(recoJPsi.size()>0) num[4]++;
    }
  }
  
  /**************************************************************************Checking for Prompt Gsf electrons as a function of JPsi pT ***************************************************************/

  for(unsigned int i=0;i<gsf_elec.size();i++){
    gsf_e_histo[79]->Fill(gsf_elec.at(i).dxy);
    gsf_e_histo[80]->Fill(gsf_elec.at(i).dz);
  }
  
  if(recoJPsi.size()>0){
    float Mee = (gsf_elec.at(0).v+gsf_elec.at(1).v).M();
    bool isOS = (gsf_elec.at(0).id)*(gsf_elec.at(1).id)<0;
    if(Mee>2.6 && Mee<3.4 && isOS){
      JPsi_histo[4]->Fill(gsf_elec.at(0).dxy,recoJPsi.at(0).v.Pt());
      JPsi_histo[5]->Fill(gsf_elec.at(1).dxy,recoJPsi.at(0).v.Pt());
      JPsi_histo[6]->Fill(gsf_elec.at(0).dz,recoJPsi.at(0).v.Pt());
      JPsi_histo[7]->Fill(gsf_elec.at(1).dz,recoJPsi.at(0).v.Pt());
    }
  }
  
  /****************************************************************************************************************************************************************************************************

                                                                        Checking properties of Gsf Electron matched to genElectron
								       
  ****************************************************************************************************************************************************************************************************/
  
  //Checking dR of reco Gsf and genGsf where both the gen is matched to reco Gsf

  if(gsf_elec.size()>1){
    float dRmin0=999;
    float dRmin1=999;
    int ind0=0;
    int ind1=0;
    for(unsigned int i=0;i<genElec.size();i++){
      float dR0 = gsf_elec.at(0).v.DeltaR(genElec.at(i).v);
      float dR1 = gsf_elec.at(1).v.DeltaR(genElec.at(i).v);
      if(dR0<dRmin0){
	dRmin0 = dR0;
	ind0 = i;
      }
      if(dR1<dRmin1){
	dRmin1 = dR1;
	ind1 = i;
      }
    }
    if(genElec.size()>0){
      float dR_gsf = gsf_elec.at(0).v.DeltaR(gsf_elec.at(1).v);
      float dR_gen = genElec.at(ind0).v.DeltaR(genElec.at(ind1).v);
      if(dRmin0<0.1 && dRmin1<0.1 && ind0!=ind1){
	gen_histo[47]->Fill(dR_gen);
	gen_histo[48]->Fill(dR_gsf);
      }
    }
  }
  
  //Plotting some user-defined variable for gsf and gen eelctron that are matched

  for(unsigned int i=0;i<gsf_elec.size();i++){
    float dRmin = 999;
    int ind=0;
    for(unsigned int j=0;j<genElec.size();j++){
      float dR = gsf_elec.at(i).v.DeltaR(genElec.at(j).v);
      if(dR<dRmin){
	dRmin = dR;
	ind = j;
      }
    }
    float Esup = gsf_elec.at(i).correctedEcalEnergy;
    float Eraw = gsf_elec.at(i).rawEnergy;
    if(genElec.size()>0 && dRmin<0.1){
      float Egen = genElec.at(ind).E;
      gen_histo[49]->Fill(fabs(Esup-Egen)/Egen);
      gen_histo[50]->Fill(fabs(Eraw-Egen)/Egen);
    }
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

  std::cout<<"***********************************************Efficiency with another approach (16th aug 2023) ********************************************"<<endl;
  std::cout<<"                                                                                                                                            "<<endl;
  std::cout<<"Numerator: "<<num[0]<<" "<<"Denominator: "<<den[0]<<"Efficiency if dR<0.05:    "<<(float)num[0]/(float)den[0]<<endl;
  std::cout<<"Numerator: "<<num[1]<<" "<<"Denominator: "<<den[1]<<"Efficiency if0.05<dR<0.1: "<<(float)num[1]/(float)den[1]<<endl;
  std::cout<<"Numerator: "<<num[2]<<" "<<"Denominator: "<<den[2]<<"Efficiency if 0.1<dR<0.2: "<<(float)num[2]/(float)den[2]<<endl;
  std::cout<<"Numerator: "<<num[3]<<" "<<"Denominator: "<<den[3]<<"Efficiency if 0.2<dR<0.3  "<<(float)num[3]/(float)den[3]<<endl;
  std::cout<<"Numerator: "<<num[4]<<" "<<"Denominator: "<<den[4]<<"Efficiency if dR>0.3:     "<<(float)num[4]/(float)den[4]<<endl;
  

  /*file_convert_photons.close();
  file_unconvert_photons.close();
  file_twogsf_EB.close();
  file_onegsf_EB.close();
  file_delta_seedietaiphi_twogsf_EB_eta_05.close();
  file_delta_seedietaiphi_twogsf_EB_eta_05_1.close();
  file_delta_seedietaiphi_twogsf_EB_eta_1_2.close();
  file_delta_seedietaiphi_twogsf_EB_eta_2_.close();
  file_eff.close();*/
  file_gsfrechit_EB.close();
  file_photonrechit_convert_EB.close();
  file_photonrechit_unconvert_EB.close();

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
    for(int i=0;i<(int)genElec.size()-1;i++){
      for(int j=i+1;j<(int)genElec.size();j++){
        if(genElec[i].v.Pt()<genElec[j].v.Pt()) swap(genElec.at(i),genElec.at(j));
      }
    }

    for(int i=0;i<(int)genElec_momJPsi.size()-1;i++){
      for(int j=i+1;j<(int)genElec_momJPsi.size();j++){
        if(genElec_momJPsi[i].v.Pt()<genElec_momJPsi[j].v.Pt()) swap(genElec_momJPsi.at(i),genElec_momJPsi.at(j));
      }
    }
    for(int i=0;i<(int)genJPsi_to_ee.size()-1;i++){
      for(int j=i+1;j<(int)genJPsi_to_ee.size();j++){
	if(genJPsi_to_ee[i].v.Pt()<genJPsi_to_ee[j].v.Pt()) swap(genJPsi_to_ee.at(i),genJPsi_to_ee.at(j));
      }
    }
    
 }
  
  if(opt==1){
    for(int i=0;i<(int)gsf_elec.size()-1;i++){
      for(int j=i+1;j<(int)gsf_elec.size();j++){
        if(gsf_elec[i].v.E()<gsf_elec[j].v.E()) swap(gsf_elec.at(i),gsf_elec.at(j));
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
