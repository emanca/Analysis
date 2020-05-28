// -*- C++ -*-
//
// Package:    TrackAnalysis/HitAnalyzer
// Class:      HitAnalyzer
//
/**\class HitAnalyzer HitAnalyzer.cc TrackAnalysis/HitAnalyzer/plugins/HitAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michail Bachtis
//         Created:  Mon, 21 Mar 2016 14:17:37 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// required for Transient Tracks
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
// required for vtx fitting
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#define _MuMass_ 0.1056583745
#define _MuMassErr_ 0.0000000024

//
// class declaration
//

class MuonCalibAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit MuonCalibAnalyzer(const edm::ParameterSet &);
  ~MuonCalibAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  bool muonIDZ(const pat::Muon &, const reco::Vertex &);
  bool muonIDOnia(const pat::Muon &, const reco::Vertex &);
  bool selectResonance(const pat::Muon &, const pat::Muon &, const reco::Vertex &);
  int findGenParticle(const pat::Muon &, const reco::GenParticleCollection &);

  // ----------member data ---------------------------
  edm::EDGetTokenT<pat::MuonCollection> muons_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;
  edm::EDGetTokenT<pat::METCollection> mets_;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;

  TFile *fout;
  TTree *tree;

  float pt1;
  float c1;
  float sc1;
  float dsc1;
  float gc1;
  float eta1;
  float phi1;
  float dpt1;
  float dxy1;
  float dz1;
  float mcpt1;
  float mceta1;
  float mcphi1;

  float pt2;
  float c2;
  float sc2;
  float dsc2;
  float gc2;
  float eta2;
  float phi2;
  float dpt2;
  float dxy2;
  float dz2;
  float mcpt2;
  float mceta2;
  float mcphi2;

  // kinematic fit
  float ptvtx1;
  float etavtx1;
  float phivtx1;
  float ptErrvtx1;
  float ptvtx2;
  float etavtx2;
  float phivtx2;
  float ptErrvtx2;

  float mass;
  float massErr;
  float genMass;
  float recoil;
  float massvtx;
  float xPV;
  float yPV;
  float zPV;
  float xvtx;
  float yvtx;
  float zvtx;
  float mcxvtx;
  float mcyvtx;
  float mczvtx;
  float mcxhs;
  float mcyhs;
  float mczhs;

  float dxy1_bmsp;
  float dxy2_bmsp;
  float dz1_bmsp;
  float dz2_bmsp;
  float sigmaz_bmsp;
  float widthx_bmsp;
  float widthy_bmsp;

  float dxy1_mcvtx;
  float dxy2_mcvtx;
  float dz1_mcvtx;
  float dz2_mcvtx;

  int run;
  float genWeight;
  bool isOnia_;
};

bool MuonCalibAnalyzer::muonIDZ(const pat::Muon &muon, const reco::Vertex &vertex)
{
  bool acc = muon.pt() > 10 && fabs(muon.eta()) < 2.4;
  bool id = muon::isTightMuon(muon, vertex);
  bool iso = muon.isolationR03().sumPt / muon.pt() < 0.2;
  bool vtx = fabs(muon.innerTrack()->dxy(vertex.position())) < 0.02 && fabs(muon.innerTrack()->dz(vertex.position())) < 0.2;

  return id && vtx && iso && acc;
}

bool MuonCalibAnalyzer::muonIDOnia(const pat::Muon &muon, const reco::Vertex &vertex)
{
  bool acc = muon.pt() > 3.0 && fabs(muon.eta()) < 2.4;
  bool id = muon::isSoftMuon(muon, vertex);
  bool vtx = fabs(muon.innerTrack()->dxy(vertex.position())) < 0.02 && fabs(muon.innerTrack()->dz(vertex.position())) < 0.2;
  return id && vtx && acc;
}

bool MuonCalibAnalyzer::selectResonance(const pat::Muon &m1, const pat::Muon &m2, const reco::Vertex &vertex)
{

  if (!(m1.isTrackerMuon() && m2.isTrackerMuon()))
    return false;

  if (!(m1.charge() + m2.charge() == 0))
    return false;

  if (isOnia_)
  {
    if (!muonIDOnia(m1, vertex))
      return false;
    if (!muonIDOnia(m2, vertex))
      return false;

    float m = (m1.p4() + m2.p4()).M();

    if ((m > 2.7 && m < 3.5) || (m > 9 && m < 10))
      return true;

    return false;
  }
  else
  {
    if (!muonIDZ(m1, vertex))
      return false;
    if (!muonIDZ(m2, vertex))
      return false;

    float m = (m1.p4() + m2.p4()).M();

    if (m > 70 && m < 120)
      return true;

    return false;
  }
}

int MuonCalibAnalyzer::findGenParticle(const pat::Muon &mu, const reco::GenParticleCollection &gen)
{
  float drMin = 10000;
  int index = -1;
  for (unsigned int i = 0; i < gen.size(); ++i)
  {
    const auto &p = gen[i];
    float dr = deltaR(p.eta(), p.phi(), mu.eta(), mu.phi());
    if (dr < 0.1)
    {
      if (dr < drMin)
      {
        drMin = dr;
        index = (int)i;
      }
    }
  }
  return index;
}

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonCalibAnalyzer::MuonCalibAnalyzer(const edm::ParameterSet &iConfig)
{
  //now do what ever initialization is needed
  muons_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  genParticles_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  vertices_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  mets_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("met"));
  genInfoToken_ = consumes<GenEventInfoProduct>(edm::InputTag("generator"));

  isOnia_ = iConfig.getParameter<bool>("isOnia");

  fout = new TFile("muonTree.root", "RECREATE");
  tree = new TTree("tree", "tree");

  tree->Branch("pt1", &pt1, "pt1/F");
  tree->Branch("pt2", &pt2, "pt2/F");
  tree->Branch("mcpt1", &mcpt1, "mcpt1/F");
  tree->Branch("mcpt2", &mcpt2, "mcpt2/F");
  tree->Branch("gc1", &gc1, "gc1/F");
  tree->Branch("gc2", &gc2, "gc2/F");
  tree->Branch("c1", &c1, "c1/F");
  tree->Branch("c2", &c2, "c2/F");
  tree->Branch("sc1", &sc1, "sc1/F");
  tree->Branch("sc2", &sc2, "sc2/F");
  tree->Branch("dsc1", &dsc1, "dsc1/F");
  tree->Branch("dsc2", &dsc2, "dsc2/F");
  tree->Branch("eta1", &eta1, "eta1/F");
  tree->Branch("eta2", &eta2, "eta2/F");
  tree->Branch("mceta1", &mceta1, "mceta1/F");
  tree->Branch("mceta2", &mceta2, "mceta2/F");
  tree->Branch("phi1", &phi1, "phi1/F");
  tree->Branch("phi2", &phi2, "phi2/F");
  tree->Branch("mcphi1", &mcphi1, "mcphi1/F");
  tree->Branch("mcphi2", &mcphi2, "mcphi2/F");
  tree->Branch("cErr1", &dpt1, "cErr1/F");
  tree->Branch("cErr2", &dpt2, "cErr2/F");
  tree->Branch("dxy1", &dxy1, "dxy1/F");
  tree->Branch("dxy2", &dxy2, "dxy2/F");
  tree->Branch("dz1", &dz1, "dz1/F");
  tree->Branch("dz2", &dz2, "dz2/F");
  tree->Branch("mass", &mass, "mass/F");
  tree->Branch("massErr", &massErr, "massErr/F");
  tree->Branch("genMass", &genMass, "genMass/F");
  tree->Branch("run", &run, "run/I");
  tree->Branch("recoil", &recoil, "recoil/F");
  tree->Branch("ptvtx1", &ptvtx1, "ptvtx1/F");
  tree->Branch("etavtx1", &etavtx1, "etavtx1/F");
  tree->Branch("phivtx1", &phivtx1, "phivtx1/F");
  tree->Branch("ptErrvtx1", &ptErrvtx1, "ptErrvtx1/F");
  tree->Branch("ptvtx2", &ptvtx2, "ptvtx2/F");
  tree->Branch("etavtx2", &etavtx2, "etavtx2/F");
  tree->Branch("phivtx2", &phivtx2, "phivtx2/F");
  tree->Branch("ptErrvtx2", &ptErrvtx2, "ptErrvtx2/F");
  tree->Branch("massvtx", &massvtx, "massvtx/F");
  tree->Branch("xPV", &xPV, "xPV/F");
  tree->Branch("yPV", &yPV, "yPV/F");
  tree->Branch("zPV", &zPV, "zPV/F");
  tree->Branch("xvtx", &xvtx, "xvtx/F");
  tree->Branch("yvtx", &yvtx, "yvtx/F");
  tree->Branch("zvtx", &zvtx, "zvtx/F");
  tree->Branch("mcxvtx", &mcxvtx, "mcxvtx/F");
  tree->Branch("mcyvtx", &mcyvtx, "mcyvtx/F");
  tree->Branch("mczvtx", &mczvtx, "mczvtx/F");
  tree->Branch("mcxhs", &mcxhs, "mcxhs/F");
  tree->Branch("mcyhs", &mcyhs, "mcyhs/F");
  tree->Branch("mczhs", &mczhs, "mczhs/F");
  tree->Branch("genWeight", &genWeight, "genWeight/F");
  tree->Branch("dxy1_bmsp", &dxy1_bmsp, "dxy1_bmsp/F");
  tree->Branch("dxy2_bmsp", &dxy2_bmsp, "dxy2_bmsp/F");
  tree->Branch("dz1_bmsp", &dz1_bmsp, "dz1_bmsp/F");
  tree->Branch("dz2_bmsp", &dz2_bmsp, "dz2_bmsp/F");
  tree->Branch("sigmaz_bmsp", &sigmaz_bmsp, "sigmaz_bmsp/F");
  tree->Branch("widthx_bmsp", &widthx_bmsp, "widthx_bmsp/F");
  tree->Branch("widthy_bmsp", &widthy_bmsp, "widthy_bmsp/F");
  tree->Branch("dxy1_mcvtx", &dxy1_mcvtx, "dxy1_mcvtx/F");
  tree->Branch("dxy2_mcvtx", &dxy2_mcvtx, "dxy2_mcvtx/F");
  tree->Branch("dz1_mcvtx", &dz1_mcvtx, "dz1_mcvtx/F");
  tree->Branch("dz2_mcvtx", &dz2_mcvtx, "dz2_mcvtx/F");
}

MuonCalibAnalyzer::~MuonCalibAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void MuonCalibAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  using namespace edm;

  reco::GenParticleCollection genMuons;

  if (iEvent.run() < 50000)
  {
    Handle<reco::GenParticleCollection> genParticlesH;
    iEvent.getByToken(genParticles_, genParticlesH);
    for (unsigned int i = 0; i < genParticlesH->size(); ++i)
    {
      const auto &p = (*genParticlesH)[i];
      if (abs(p.pdgId()) == 13 and p.status() == 1)
        genMuons.push_back(p);
    }
    const auto &proton = (*genParticlesH)[0];
    mcxhs = proton.vx();
    mcyhs = proton.vy();
    mczhs = proton.vz();

    Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(genInfoToken_, genInfo);
    GenEventInfoProduct genInfoP = *(genInfo.product());
    genWeight = genInfoP.weight();
  }

  run = iEvent.run();

  Handle<reco::VertexCollection> vertexH;
  iEvent.getByToken(vertices_, vertexH);

  if (vertexH->size() == 0)
    return;

  const reco::Vertex &vertex = vertexH->at(0);
  xPV = vertex.position().x();
  yPV = vertex.position().y();
  zPV = vertex.position().z();

  Handle<pat::MuonCollection> muonH;
  iEvent.getByToken(muons_, muonH);

  unsigned N = muonH->size();

  if (N < 2)
    return;

  Handle<pat::METCollection> metH;
  iEvent.getByToken(mets_, metH);

  reco::Candidate::LorentzVector metV = (*metH)[0].p4();

  reco::BeamSpot beamSpot;
  Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);

  for (unsigned int i = 0; i < N - 1; ++i)
  {
    for (unsigned int j = i + 1; j < N; ++j)
    {

      const pat::Muon &mu1 = (*muonH)[i];
      const pat::Muon &mu2 = (*muonH)[j];

      if (selectResonance((*muonH)[i], (*muonH)[j], vertex))
      {
        const pat::Muon &pos = mu1.charge() > 0 ? mu1 : mu2;
        const pat::Muon &neg = mu1.charge() < 0 ? mu1 : mu2;

        recoil = (-(pos.p4() + neg.p4() + metV)).pt();

        pt1 = pos.innerTrack()->pt();
        if (pos.isGlobalMuon() && pos.standAloneMuon().isNonnull())
        {
          sc1 = 1.0 / pos.standAloneMuon()->pt();
          dsc1 = pos.standAloneMuon()->ptError() / pos.standAloneMuon()->pt();
        }
        else
        {
          sc1 = -1.0;
          dsc1 = -1.0;
        }

        dpt1 = pos.innerTrack()->ptError() / pt1;
        c1 = 1.0 / pt1;
        eta1 = pos.innerTrack()->eta();
        phi1 = pos.innerTrack()->phi();
        int ptr1 = findGenParticle(pos, genMuons);

        //	   printf("Track length=%f\n",pos.innerTrack()->outerPosition().rho()-pos.innerTrack()->innerPosition().rho());

        if (ptr1 >= 0)
        {
          mcpt1 = genMuons[ptr1].pt();
          mceta1 = genMuons[ptr1].eta();
          mcphi1 = genMuons[ptr1].phi();
          gc1 = 1.0 / mcpt1;
          mcxvtx = genMuons[ptr1].vx();
          mcyvtx = genMuons[ptr1].vy();
          mczvtx = genMuons[ptr1].vz();

          //impact parameters wrt gen vtx (MC only)
          math::XYZPoint point(mcxvtx, mcyvtx, mczvtx);
          dxy1_mcvtx = -1. * pos.innerTrack()->dxy(point);
          dxy2_mcvtx = -1. * neg.innerTrack()->dxy(point);
          dz1_mcvtx = pos.innerTrack()->dz(point);
          dz2_mcvtx = neg.innerTrack()->dz(point);
        }
        else
        {
          mcpt1 = -1.0;
          mceta1 = -1.0;
          mcphi1 = -1.0;
          gc1 = -1.0;
          mcxvtx = -9999.;
          mcyvtx = -9999.;
          mczvtx = -9999.;
          dxy1_mcvtx = -9999.;
          dxy2_mcvtx = -9999.;
          dz1_mcvtx = -9999.;
          dz2_mcvtx = -9999.;
        }
        //impact parameters wrt PV
        dxy1 = pos.innerTrack()->dxy(vertex.position());
        dxy2 = neg.innerTrack()->dxy(vertex.position());
        dz1 = pos.innerTrack()->dz(vertex.position());
        dz2 = neg.innerTrack()->dz(vertex.position());

        //impact parameters wrt beamspot
        if (beamSpotHandle.isValid())
        {
          beamSpot = *beamSpotHandle;
          math::XYZPoint point(beamSpot.x0(), beamSpot.y0(), beamSpot.z0());
          dxy1_bmsp = -1. * pos.innerTrack()->dxy(point);
          dxy2_bmsp = -1. * neg.innerTrack()->dxy(point);
          dz1_bmsp = pos.innerTrack()->dz(point);
          dz2_bmsp = neg.innerTrack()->dz(point);
          sigmaz_bmsp = beamSpot.sigmaZ();
          widthx_bmsp = beamSpot.BeamWidthX();
          widthy_bmsp = beamSpot.BeamWidthY();
        }
        else
        {
          dxy1_bmsp = -9999.;
          dxy2_bmsp = -9999.;
          dz1_bmsp =  -9999.;
          dz2_bmsp =  -9999.;
          sigmaz_bmsp = -9999.;
          widthx_bmsp = -9999.;
          widthy_bmsp = -9999.;
        }
        
        pt2 = neg.innerTrack()->pt();

        if (neg.isGlobalMuon() && neg.standAloneMuon().isNonnull())
        {
          sc2 = 1.0 / neg.standAloneMuon()->pt();
          dsc2 = neg.standAloneMuon()->ptError() / neg.standAloneMuon()->pt();
        }

        dpt2 = neg.innerTrack()->ptError() / pt2;
        c2 = 1.0 / pt2;
        eta2 = neg.innerTrack()->eta();
        phi2 = neg.innerTrack()->phi();
        int ptr2 = findGenParticle(neg, genMuons);

        if (ptr2 >= 0)
        {
          mcpt2 = genMuons[ptr2].pt();
          mceta2 = genMuons[ptr2].eta();
          mcphi2 = genMuons[ptr2].phi();
          gc2 = 1.0 / mcpt2;
        }
        else
        {
          mcpt2 = -1.;
          mceta2 = -1.;
          mcphi2 = -1.;
          gc2 = -1.;
        }

        math::PtEtaPhiMLorentzVector posVec(pt1, eta1, phi1, 0.105658);
        math::PtEtaPhiMLorentzVector negVec(pt2, eta2, phi2, 0.105658);
        mass = (posVec + negVec).M();
        massErr = 0.5 * mass * sqrt(dpt1 * dpt1 + dpt2 * dpt2);

        genMass = -1.0;
        if (iEvent.run() < 50000 && mcpt1 > 0 && mcpt2 > 0)
        {
          math::PtEtaPhiMLorentzVector posVec(mcpt1, mceta1, mcphi1, 0.105658);
          math::PtEtaPhiMLorentzVector negVec(mcpt2, mceta2, mcphi2, 0.105658);
          genMass = (posVec + negVec).M();
        }

        // fit dimuon vertex
        edm::ESHandle<TransientTrackBuilder> TTBuilder;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", TTBuilder);
        reco::TransientTrack m1_tk = TTBuilder->build(pos.muonBestTrack());
        reco::TransientTrack m2_tk = TTBuilder->build(neg.muonBestTrack());
        std::vector<RefCountedKinematicParticle> parts;
        KinematicParticleFactoryFromTransientTrack pFactory;
        
        double chi = 0, ndf = 0;
        float mMu = _MuMass_, dmMu = _MuMassErr_;
        parts.push_back(pFactory.particle(m1_tk, mMu, chi, ndf, dmMu));
        parts.push_back(pFactory.particle(m2_tk, mMu, chi, ndf, dmMu));

        KinematicParticleVertexFitter VtxFitter;
        RefCountedKinematicTree KinTree = VtxFitter.fit(parts);
        
        TLorentzVector mu1_tlv;
        TLorentzVector mu2_tlv;
        RefCountedKinematicVertex dimu_vertex;

        if (KinTree->isEmpty() == 1 || KinTree->isConsistent() == 0)
        {
          std::cout << "Kinematic Fit unsuccessful" << std::endl;
          // assign dummy values
          ptvtx1 = -99.;
          etavtx1 = -99.;
          phivtx1 = -99.;
          ptErrvtx1 = -99.;
          ptvtx2 = -99.;
          etavtx2 = -99.;
          phivtx2 = -99.;
          ptErrvtx2 = -99.;
          massvtx = -99.;
        }

        else
        {
          //accessing the tree components
          KinTree->movePointerToTheTop();
          //We are now at the top of the decay tree getting the dimuon reconstructed KinematicParticle
          RefCountedKinematicParticle dimu_kinfit = KinTree->currentParticle();
          //getting the dimuon decay vertex
          //RefCountedKinematicVertex
          dimu_vertex = KinTree->currentDecayVertex();
          xvtx = dimu_vertex->position().x();
          yvtx = dimu_vertex->position().y();
          zvtx = dimu_vertex->position().z();
          //Now navigating down the tree
          bool child = KinTree->movePointerToTheFirstChild();
          if (child)
          {
            RefCountedKinematicParticle mu1_kinfit = KinTree->currentParticle();
            AlgebraicVector7 mu1_kinfit_par = mu1_kinfit->currentState().kinematicParameters().vector();
            AlgebraicSymMatrix77 mu1_kinfit_cov = mu1_kinfit->currentState().kinematicParametersError().matrix();
            ptErrvtx1 = sqrt(mu1_kinfit_cov(3, 3) + mu1_kinfit_cov(4, 4));
            mu1_tlv.SetXYZM(mu1_kinfit_par.At(3), mu1_kinfit_par.At(4), mu1_kinfit_par.At(5), mu1_kinfit_par.At(6));
          }
          //Now navigating down the tree
          bool nextchild = KinTree->movePointerToTheNextChild();
          if (nextchild)
          {
            RefCountedKinematicParticle mu2_kinfit = KinTree->currentParticle();
            AlgebraicVector7 mu2_kinfit_par = mu2_kinfit->currentState().kinematicParameters().vector();
            AlgebraicSymMatrix77 mu2_kinfit_cov = mu2_kinfit->currentState().kinematicParametersError().matrix();
            ptErrvtx2 = sqrt(mu2_kinfit_cov(3, 3) + mu2_kinfit_cov(4, 4));
            mu2_tlv.SetXYZM(mu2_kinfit_par.At(3), mu2_kinfit_par.At(4), mu2_kinfit_par.At(5), mu2_kinfit_par.At(6));
          }
        } // end else - isEmpty()
        ptvtx1 = mu1_tlv.Pt();
        etavtx1 = mu1_tlv.Eta();
        phivtx1 = mu1_tlv.Phi();

        ptvtx2 = mu2_tlv.Pt();
        etavtx2 = mu2_tlv.Eta();
        phivtx2 = mu2_tlv.Phi();

        massvtx = (mu1_tlv + mu2_tlv).M();

        tree->Fill();
        return;
      }
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void MuonCalibAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void MuonCalibAnalyzer::endJob()
{
  fout->Write();
  fout->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MuonCalibAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MuonCalibAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MuonCalibAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MuonCalibAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonCalibAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonCalibAnalyzer);
