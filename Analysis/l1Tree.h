//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 22 14:15:13 2014 by ROOT version 5.34/10
// from TTree L1Tree/L1Tree
// found on file: /mnt/hadoop/cms/store/user/luck/L1Emulator/L1Tree_MinBiasSkim_v1.root
//////////////////////////////////////////////////////////

#ifndef l1Tree_h
#define l1Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxmaxRCTREG = 1;

class l1Tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //L1Analysis::L1AnalysisEventDataFormat *Event;
   Int_t           run;
   Int_t           event;
   Int_t           lumi;
   Int_t           bx;
   ULong64_t       orbit;
   ULong64_t       time;
   vector<TString> hlt;
   Double_t        puWeight;
 //L1Analysis::L1AnalysisGCTDataFormat *GCT;
   Int_t           IsoEmSize;
   vector<float>   IsoEmEta;
   vector<float>   IsoEmPhi;
   vector<float>   IsoEmRnk;
   vector<int>     IsoEmBx;
   Int_t           NonIsoEmSize;
   vector<float>   NonIsoEmEta;
   vector<float>   NonIsoEmPhi;
   vector<float>   NonIsoEmRnk;
   vector<int>     NonIsoEmBx;
   Int_t           CJetSize;
   vector<float>   CJetEta;
   vector<float>   CJetPhi;
   vector<float>   CJetRnk;
   vector<int>     CJetBx;
   Int_t           FJetSize;
   vector<float>   FJetEta;
   vector<float>   FJetPhi;
   vector<float>   FJetRnk;
   vector<int>     FJetBx;
   Int_t           TJetSize;
   vector<float>   TJetEta;
   vector<float>   TJetPhi;
   vector<float>   TJetRnk;
   vector<int>     TJetBx;
   Int_t           EtMissSize;
   vector<float>   EtMiss;
   vector<float>   EtMissPhi;
   vector<float>   EtMissBX;
   Int_t           HtMissSize;
   vector<float>   HtMiss;
   vector<float>   HtMissPhi;
   vector<float>   HtMissBX;
   Int_t           EtHadSize;
   vector<float>   EtHad;
   vector<float>   EtHadBX;
   Int_t           EtTotSize;
   vector<float>   EtTot;
   vector<float>   EtTotBX;
   Int_t           HFRingEtSumSize;
   vector<float>   HFRingEtSumEta;
   Float_t         HFBitCountsSize;
   vector<float>   HFBitCountsEta;
 //L1Analysis::L1AnalysisGTDataFormat *GT;
   vector<unsigned long long> tw1;
   vector<unsigned long long> tw2;
   vector<unsigned long long> tt;
   ULong_t         partrig_tcs;
   ULong_t         gpsTimehi;
   ULong_t         gpsTimelo;
   ULong_t         bstMasterStatus;
   ULong_t         bstturnCountNumber;
   ULong_t         bstlhcFillNumber;
   ULong_t         bstbeamMode;
   ULong_t         bstparticleTypeBeam1;
   ULong_t         bstparticleTypeBeam2;
   ULong_t         bstbeamMomentum;
   ULong_t         bsttotalIntensityBeam1;
   ULong_t         bsttotalIntensityBeam2;
   Int_t           Nele;
   vector<int>     Bxel;
   vector<float>   Rankel;
   vector<float>   Phiel;
   vector<float>   Etael;
   vector<bool>    Isoel;
   Int_t           Njet;
   vector<int>     Bxjet;
   vector<float>   Rankjet;
   vector<float>   Phijet;
   vector<float>   Etajet;
   vector<bool>    Taujet;
   vector<bool>    Fwdjet;
   Int_t           RankETT;
   Bool_t          OvETT;
   Int_t           RankHTT;
   Bool_t          OvHTT;
   Int_t           RankETM;
   Int_t           PhiETM;
   Bool_t          OvETM;
   Int_t           RankHTM;
   Int_t           PhiHTM;
   Bool_t          OvHTM;
 //L1Analysis::L1AnalysisRCTDataFormat *RCT;
   Int_t           maxRCTREG_;
   Int_t           RegSize;
   vector<float>   RegEta;
   vector<float>   RegPhi;
   vector<float>   RegGEta;
   vector<float>   RegGPhi;
   vector<float>   RegRnk;
   vector<int>     RegVeto;
   vector<int>     RegBx;
   vector<int>     RegOverFlow;
   vector<int>     RegMip;
   vector<int>     RegFGrain;
   Int_t           EmSize;
   vector<int>     IsIsoEm;
   vector<float>   EmEta;
   vector<float>   EmPhi;
   vector<float>   EmRnk;
   vector<int>     EmBx;

   // List of branches
   TBranch        *b_Event_run;   //!
   TBranch        *b_Event_event;   //!
   TBranch        *b_Event_lumi;   //!
   TBranch        *b_Event_bx;   //!
   TBranch        *b_Event_orbit;   //!
   TBranch        *b_Event_time;   //!
   TBranch        *b_Event_hlt;   //!
   TBranch        *b_Event_puWeight;   //!
   TBranch        *b_GCT_IsoEmSize;   //!
   TBranch        *b_GCT_IsoEmEta;   //!
   TBranch        *b_GCT_IsoEmPhi;   //!
   TBranch        *b_GCT_IsoEmRnk;   //!
   TBranch        *b_GCT_IsoEmBx;   //!
   TBranch        *b_GCT_NonIsoEmSize;   //!
   TBranch        *b_GCT_NonIsoEmEta;   //!
   TBranch        *b_GCT_NonIsoEmPhi;   //!
   TBranch        *b_GCT_NonIsoEmRnk;   //!
   TBranch        *b_GCT_NonIsoEmBx;   //!
   TBranch        *b_GCT_CJetSize;   //!
   TBranch        *b_GCT_CJetEta;   //!
   TBranch        *b_GCT_CJetPhi;   //!
   TBranch        *b_GCT_CJetRnk;   //!
   TBranch        *b_GCT_CJetBx;   //!
   TBranch        *b_GCT_FJetSize;   //!
   TBranch        *b_GCT_FJetEta;   //!
   TBranch        *b_GCT_FJetPhi;   //!
   TBranch        *b_GCT_FJetRnk;   //!
   TBranch        *b_GCT_FJetBx;   //!
   TBranch        *b_GCT_TJetSize;   //!
   TBranch        *b_GCT_TJetEta;   //!
   TBranch        *b_GCT_TJetPhi;   //!
   TBranch        *b_GCT_TJetRnk;   //!
   TBranch        *b_GCT_TJetBx;   //!
   TBranch        *b_GCT_EtMissSize;   //!
   TBranch        *b_GCT_EtMiss;   //!
   TBranch        *b_GCT_EtMissPhi;   //!
   TBranch        *b_GCT_EtMissBX;   //!
   TBranch        *b_GCT_HtMissSize;   //!
   TBranch        *b_GCT_HtMiss;   //!
   TBranch        *b_GCT_HtMissPhi;   //!
   TBranch        *b_GCT_HtMissBX;   //!
   TBranch        *b_GCT_EtHadSize;   //!
   TBranch        *b_GCT_EtHad;   //!
   TBranch        *b_GCT_EtHadBX;   //!
   TBranch        *b_GCT_EtTotSize;   //!
   TBranch        *b_GCT_EtTot;   //!
   TBranch        *b_GCT_EtTotBX;   //!
   TBranch        *b_GCT_HFRingEtSumSize;   //!
   TBranch        *b_GCT_HFRingEtSumEta;   //!
   TBranch        *b_GCT_HFBitCountsSize;   //!
   TBranch        *b_GCT_HFBitCountsEta;   //!
   TBranch        *b_GT_tw1;   //!
   TBranch        *b_GT_tw2;   //!
   TBranch        *b_GT_tt;   //!
   TBranch        *b_GT_partrig_tcs;   //!
   TBranch        *b_GT_gpsTimehi;   //!
   TBranch        *b_GT_gpsTimelo;   //!
   TBranch        *b_GT_bstMasterStatus;   //!
   TBranch        *b_GT_bstturnCountNumber;   //!
   TBranch        *b_GT_bstlhcFillNumber;   //!
   TBranch        *b_GT_bstbeamMode;   //!
   TBranch        *b_GT_bstparticleTypeBeam1;   //!
   TBranch        *b_GT_bstparticleTypeBeam2;   //!
   TBranch        *b_GT_bstbeamMomentum;   //!
   TBranch        *b_GT_bsttotalIntensityBeam1;   //!
   TBranch        *b_GT_bsttotalIntensityBeam2;   //!
   TBranch        *b_GT_Nele;   //!
   TBranch        *b_GT_Bxel;   //!
   TBranch        *b_GT_Rankel;   //!
   TBranch        *b_GT_Phiel;   //!
   TBranch        *b_GT_Etael;   //!
   TBranch        *b_GT_Isoel;   //!
   TBranch        *b_GT_Njet;   //!
   TBranch        *b_GT_Bxjet;   //!
   TBranch        *b_GT_Rankjet;   //!
   TBranch        *b_GT_Phijet;   //!
   TBranch        *b_GT_Etajet;   //!
   TBranch        *b_GT_Taujet;   //!
   TBranch        *b_GT_Fwdjet;   //!
   TBranch        *b_GT_RankETT;   //!
   TBranch        *b_GT_OvETT;   //!
   TBranch        *b_GT_RankHTT;   //!
   TBranch        *b_GT_OvHTT;   //!
   TBranch        *b_GT_RankETM;   //!
   TBranch        *b_GT_PhiETM;   //!
   TBranch        *b_GT_OvETM;   //!
   TBranch        *b_GT_RankHTM;   //!
   TBranch        *b_GT_PhiHTM;   //!
   TBranch        *b_GT_OvHTM;   //!
   TBranch        *b_RCT_maxRCTREG_;   //!
   TBranch        *b_RCT_RegSize;   //!
   TBranch        *b_RCT_RegEta;   //!
   TBranch        *b_RCT_RegPhi;   //!
   TBranch        *b_RCT_RegGEta;   //!
   TBranch        *b_RCT_RegGPhi;   //!
   TBranch        *b_RCT_RegRnk;   //!
   TBranch        *b_RCT_RegVeto;   //!
   TBranch        *b_RCT_RegBx;   //!
   TBranch        *b_RCT_RegOverFlow;   //!
   TBranch        *b_RCT_RegMip;   //!
   TBranch        *b_RCT_RegFGrain;   //!
   TBranch        *b_RCT_EmSize;   //!
   TBranch        *b_RCT_IsIsoEm;   //!
   TBranch        *b_RCT_EmEta;   //!
   TBranch        *b_RCT_EmPhi;   //!
   TBranch        *b_RCT_EmRnk;   //!
   TBranch        *b_RCT_EmBx;   //!

   l1Tree(TTree *tree=0);
   virtual ~l1Tree();
   //virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

//#ifdef l1Tree_cxx
l1Tree::l1Tree(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/mnt/hadoop/cms/store/user/luck/L1Emulator/L1Tree_MinBiasSkim_v1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/mnt/hadoop/cms/store/user/luck/L1Emulator/L1Tree_MinBiasSkim_v1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/mnt/hadoop/cms/store/user/luck/L1Emulator/L1Tree_MinBiasSkim_v1.root:/l1NtupleProducer");
      dir->GetObject("L1Tree",tree);

   }
   Init(tree);
}

l1Tree::~l1Tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t l1Tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t l1Tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void l1Tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_Event_run);
   fChain->SetBranchAddress("event", &event, &b_Event_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_Event_lumi);
   fChain->SetBranchAddress("bx", &bx, &b_Event_bx);
   fChain->SetBranchAddress("orbit", &orbit, &b_Event_orbit);
   fChain->SetBranchAddress("time", &time, &b_Event_time);
   fChain->SetBranchAddress("hlt", &hlt, &b_Event_hlt);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_Event_puWeight);
   fChain->SetBranchAddress("IsoEmSize", &IsoEmSize, &b_GCT_IsoEmSize);
   fChain->SetBranchAddress("IsoEmEta", &IsoEmEta, &b_GCT_IsoEmEta);
   fChain->SetBranchAddress("IsoEmPhi", &IsoEmPhi, &b_GCT_IsoEmPhi);
   fChain->SetBranchAddress("IsoEmRnk", &IsoEmRnk, &b_GCT_IsoEmRnk);
   fChain->SetBranchAddress("IsoEmBx", &IsoEmBx, &b_GCT_IsoEmBx);
   fChain->SetBranchAddress("NonIsoEmSize", &NonIsoEmSize, &b_GCT_NonIsoEmSize);
   fChain->SetBranchAddress("NonIsoEmEta", &NonIsoEmEta, &b_GCT_NonIsoEmEta);
   fChain->SetBranchAddress("NonIsoEmPhi", &NonIsoEmPhi, &b_GCT_NonIsoEmPhi);
   fChain->SetBranchAddress("NonIsoEmRnk", &NonIsoEmRnk, &b_GCT_NonIsoEmRnk);
   fChain->SetBranchAddress("NonIsoEmBx", &NonIsoEmBx, &b_GCT_NonIsoEmBx);
   fChain->SetBranchAddress("CJetSize", &CJetSize, &b_GCT_CJetSize);
   fChain->SetBranchAddress("CJetEta", &CJetEta, &b_GCT_CJetEta);
   fChain->SetBranchAddress("CJetPhi", &CJetPhi, &b_GCT_CJetPhi);
   fChain->SetBranchAddress("CJetRnk", &CJetRnk, &b_GCT_CJetRnk);
   fChain->SetBranchAddress("CJetBx", &CJetBx, &b_GCT_CJetBx);
   fChain->SetBranchAddress("FJetSize", &FJetSize, &b_GCT_FJetSize);
   fChain->SetBranchAddress("FJetEta", &FJetEta, &b_GCT_FJetEta);
   fChain->SetBranchAddress("FJetPhi", &FJetPhi, &b_GCT_FJetPhi);
   fChain->SetBranchAddress("FJetRnk", &FJetRnk, &b_GCT_FJetRnk);
   fChain->SetBranchAddress("FJetBx", &FJetBx, &b_GCT_FJetBx);
   fChain->SetBranchAddress("TJetSize", &TJetSize, &b_GCT_TJetSize);
   fChain->SetBranchAddress("TJetEta", &TJetEta, &b_GCT_TJetEta);
   fChain->SetBranchAddress("TJetPhi", &TJetPhi, &b_GCT_TJetPhi);
   fChain->SetBranchAddress("TJetRnk", &TJetRnk, &b_GCT_TJetRnk);
   fChain->SetBranchAddress("TJetBx", &TJetBx, &b_GCT_TJetBx);
   fChain->SetBranchAddress("EtMissSize", &EtMissSize, &b_GCT_EtMissSize);
   fChain->SetBranchAddress("EtMiss", &EtMiss, &b_GCT_EtMiss);
   fChain->SetBranchAddress("EtMissPhi", &EtMissPhi, &b_GCT_EtMissPhi);
   fChain->SetBranchAddress("EtMissBX", &EtMissBX, &b_GCT_EtMissBX);
   fChain->SetBranchAddress("HtMissSize", &HtMissSize, &b_GCT_HtMissSize);
   fChain->SetBranchAddress("HtMiss", &HtMiss, &b_GCT_HtMiss);
   fChain->SetBranchAddress("HtMissPhi", &HtMissPhi, &b_GCT_HtMissPhi);
   fChain->SetBranchAddress("HtMissBX", &HtMissBX, &b_GCT_HtMissBX);
   fChain->SetBranchAddress("EtHadSize", &EtHadSize, &b_GCT_EtHadSize);
   fChain->SetBranchAddress("EtHad", &EtHad, &b_GCT_EtHad);
   fChain->SetBranchAddress("EtHadBX", &EtHadBX, &b_GCT_EtHadBX);
   fChain->SetBranchAddress("EtTotSize", &EtTotSize, &b_GCT_EtTotSize);
   fChain->SetBranchAddress("EtTot", &EtTot, &b_GCT_EtTot);
   fChain->SetBranchAddress("EtTotBX", &EtTotBX, &b_GCT_EtTotBX);
   fChain->SetBranchAddress("HFRingEtSumSize", &HFRingEtSumSize, &b_GCT_HFRingEtSumSize);
   fChain->SetBranchAddress("HFRingEtSumEta", &HFRingEtSumEta, &b_GCT_HFRingEtSumEta);
   fChain->SetBranchAddress("HFBitCountsSize", &HFBitCountsSize, &b_GCT_HFBitCountsSize);
   fChain->SetBranchAddress("HFBitCountsEta", &HFBitCountsEta, &b_GCT_HFBitCountsEta);
   fChain->SetBranchAddress("tw1", &tw1, &b_GT_tw1);
   fChain->SetBranchAddress("tw2", &tw2, &b_GT_tw2);
   fChain->SetBranchAddress("tt", &tt, &b_GT_tt);
   fChain->SetBranchAddress("partrig_tcs", &partrig_tcs, &b_GT_partrig_tcs);
   fChain->SetBranchAddress("gpsTimehi", &gpsTimehi, &b_GT_gpsTimehi);
   fChain->SetBranchAddress("gpsTimelo", &gpsTimelo, &b_GT_gpsTimelo);
   fChain->SetBranchAddress("bstMasterStatus", &bstMasterStatus, &b_GT_bstMasterStatus);
   fChain->SetBranchAddress("bstturnCountNumber", &bstturnCountNumber, &b_GT_bstturnCountNumber);
   fChain->SetBranchAddress("bstlhcFillNumber", &bstlhcFillNumber, &b_GT_bstlhcFillNumber);
   fChain->SetBranchAddress("bstbeamMode", &bstbeamMode, &b_GT_bstbeamMode);
   fChain->SetBranchAddress("bstparticleTypeBeam1", &bstparticleTypeBeam1, &b_GT_bstparticleTypeBeam1);
   fChain->SetBranchAddress("bstparticleTypeBeam2", &bstparticleTypeBeam2, &b_GT_bstparticleTypeBeam2);
   fChain->SetBranchAddress("bstbeamMomentum", &bstbeamMomentum, &b_GT_bstbeamMomentum);
   fChain->SetBranchAddress("bsttotalIntensityBeam1", &bsttotalIntensityBeam1, &b_GT_bsttotalIntensityBeam1);
   fChain->SetBranchAddress("bsttotalIntensityBeam2", &bsttotalIntensityBeam2, &b_GT_bsttotalIntensityBeam2);
   fChain->SetBranchAddress("Nele", &Nele, &b_GT_Nele);
   fChain->SetBranchAddress("Bxel", &Bxel, &b_GT_Bxel);
   fChain->SetBranchAddress("Rankel", &Rankel, &b_GT_Rankel);
   fChain->SetBranchAddress("Phiel", &Phiel, &b_GT_Phiel);
   fChain->SetBranchAddress("Etael", &Etael, &b_GT_Etael);
   fChain->SetBranchAddress("Isoel", &Isoel, &b_GT_Isoel);
   fChain->SetBranchAddress("Njet", &Njet, &b_GT_Njet);
   fChain->SetBranchAddress("Bxjet", &Bxjet, &b_GT_Bxjet);
   fChain->SetBranchAddress("Rankjet", &Rankjet, &b_GT_Rankjet);
   fChain->SetBranchAddress("Phijet", &Phijet, &b_GT_Phijet);
   fChain->SetBranchAddress("Etajet", &Etajet, &b_GT_Etajet);
   fChain->SetBranchAddress("Taujet", &Taujet, &b_GT_Taujet);
   fChain->SetBranchAddress("Fwdjet", &Fwdjet, &b_GT_Fwdjet);
   fChain->SetBranchAddress("RankETT", &RankETT, &b_GT_RankETT);
   fChain->SetBranchAddress("OvETT", &OvETT, &b_GT_OvETT);
   fChain->SetBranchAddress("RankHTT", &RankHTT, &b_GT_RankHTT);
   fChain->SetBranchAddress("OvHTT", &OvHTT, &b_GT_OvHTT);
   fChain->SetBranchAddress("RankETM", &RankETM, &b_GT_RankETM);
   fChain->SetBranchAddress("PhiETM", &PhiETM, &b_GT_PhiETM);
   fChain->SetBranchAddress("OvETM", &OvETM, &b_GT_OvETM);
   fChain->SetBranchAddress("RankHTM", &RankHTM, &b_GT_RankHTM);
   fChain->SetBranchAddress("PhiHTM", &PhiHTM, &b_GT_PhiHTM);
   fChain->SetBranchAddress("OvHTM", &OvHTM, &b_GT_OvHTM);
   fChain->SetBranchAddress("maxRCTREG_", &maxRCTREG_, &b_RCT_maxRCTREG_);
   fChain->SetBranchAddress("RegSize", &RegSize, &b_RCT_RegSize);
   fChain->SetBranchAddress("RegEta", &RegEta, &b_RCT_RegEta);
   fChain->SetBranchAddress("RegPhi", &RegPhi, &b_RCT_RegPhi);
   fChain->SetBranchAddress("RegGEta", &RegGEta, &b_RCT_RegGEta);
   fChain->SetBranchAddress("RegGPhi", &RegGPhi, &b_RCT_RegGPhi);
   fChain->SetBranchAddress("RegRnk", &RegRnk, &b_RCT_RegRnk);
   fChain->SetBranchAddress("RegVeto", &RegVeto, &b_RCT_RegVeto);
   fChain->SetBranchAddress("RegBx", &RegBx, &b_RCT_RegBx);
   fChain->SetBranchAddress("RegOverFlow", &RegOverFlow, &b_RCT_RegOverFlow);
   fChain->SetBranchAddress("RegMip", &RegMip, &b_RCT_RegMip);
   fChain->SetBranchAddress("RegFGrain", &RegFGrain, &b_RCT_RegFGrain);
   fChain->SetBranchAddress("EmSize", &EmSize, &b_RCT_EmSize);
   fChain->SetBranchAddress("IsIsoEm", &IsIsoEm, &b_RCT_IsIsoEm);
   fChain->SetBranchAddress("EmEta", &EmEta, &b_RCT_EmEta);
   fChain->SetBranchAddress("EmPhi", &EmPhi, &b_RCT_EmPhi);
   fChain->SetBranchAddress("EmRnk", &EmRnk, &b_RCT_EmRnk);
   fChain->SetBranchAddress("EmBx", &EmBx, &b_RCT_EmBx);
   Notify();
}

Bool_t l1Tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void l1Tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
/* Int_t l1Tree::Cut(Long64_t entry) */
/* { */
/* // This function may be called from Loop. */
/* // returns  1 if entry is accepted. */
/* // returns -1 otherwise. */
/*    return 1; */
/* } */
//#endif // #ifdef l1Tree_cxx
