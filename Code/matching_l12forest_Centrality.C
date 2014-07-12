#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TString.h>
#include <stdio.h> 
#include <TLegendEntry.h>
#include <TGraphAsymmErrors.h>

#include <vector>
#include <map>

using namespace std;
#include "l1ExtraTree.h"
#include "l1Tree.h"

int limit0=10;
int limit1=20;
int limit2=60;
int limit3=100;
int limit4=140; 

long makeKey(long run, long event){
  return (10000000*run + event);
}

void matching_l12forest_Centrality(){

  const TString l1_input = "/data/ginnocen/L1UpgradeAnalyzer_2July2014_v2.root";
  TFile *lFile = TFile::Open(l1_input);
  TTree *l1t = (TTree*)lFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");
  
  Int_t l1_evt, l1_run, l1_etsumPlus,l1_etsumMinus;
  
  l1t->SetBranchAddress("event",&l1_evt);
  l1t->SetBranchAddress("run",&l1_run);
  l1t->SetBranchAddress("centrality_hwPt",&l1_etsumPlus);
  l1t->SetBranchAddress("centrality_hwEta",&l1_etsumMinus);

  const TString forest_input = "/data/richard/0.root";
  TFile *fFile = TFile::Open(forest_input);
  TTree *fTree = (TTree*)fFile->Get("akPu3CaloJetAnalyzer/t");
  TTree *fEvtTree = (TTree*)fFile->Get("hiEvtAnalyzer/HiTree");
  fTree->AddFriend(fEvtTree);
  TTree *fSkimTree = (TTree*)fFile->Get("skimanalysis/HltTree");
  fTree->AddFriend(fSkimTree);

  Int_t f_evt, f_run, f_lumi, hiBin;
  float hiHF;
  fTree->SetBranchAddress("evt",&f_evt);
  fTree->SetBranchAddress("run",&f_run);
  fTree->SetBranchAddress("lumi",&f_lumi);
  fTree->SetBranchAddress("hiBin",&hiBin);
  fTree->SetBranchAddress("hiHF",&hiHF);
  
  Int_t pcollisionEventSelection, pHBHENoiseFilter;
  fTree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  fTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  map<long, int> kmap;
  map<long, int> kmapcal;

  int l_entries = l1t->GetEntries();
  
  for(long j = 0; j < l_entries; ++j){

    l1t->GetEntry(j);
    long key = makeKey(l1_run, l1_evt);

    pair<long,int> p(key,j);
    kmap.insert(p);
    kmapcal.insert(p);
    //printf("id number %d, run number %d, event number %d\n",j,l1_run,l1_evt);
    //printf("etsumPlus %d\n",l1_etsumPlus);
    //printf("etsumMinus %d\n",l1_etsumMinus);
  
  }
  
  
  TH2D *hcorrl1EtsumPlusVscorrl1EtsumMinus = new TH2D("hcorrl1EtsumPlusVscorrl1EtsumMinus","L1 EtsumPlus vs L1 EtsumMinus; L1 EtsumMinus ; L1 EtsumPlus",100,0,6000,100,0,6000); 
  TH2D *hcorrl1EtsumPlusVscorrl1EtsumMinusNoEvSel = new TH2D("hcorrl1EtsumPlusVscorrl1EtsumMinusNoEvSel","L1 EtsumPlus vs L1 EtsumMinus; L1 EtsumMinus ; L1 EtsumPlus",100,0,6000,100,0,6000);
  TH2D *hcorrl1EtsumVsofflineCentrality = new TH2D("hcorrl1EtsumVsofflineCentrality","L1 Etsum vs offline centrality; offline centrality ; L1 Etsum",100,0,200,100,0,6000);
  TH2D *hcorrofflineCentralityVsl1Etsum = new TH2D("hcorrofflineCentralityVsl1Etsum","Offline centrality vs L1 Etsum; L1 Etsum ; offline centrality",100,0,6000,100,0,200);  
  TH2D *hcorrL1CentralityVsfflineCentrality = new TH2D("hcorrL1CentralityVsfflineCentrality","Online centrality vs offline centrality; L1 Etsum ; offline centrality",200,0,200,200,0,200);  
  TH1D *hofflineCentrality = new TH1D("hofflineCentrality","Offline cent; Offline centrality; Entries",40,0.,200.);
  TH1D *hl1Etsum = new TH1D("hl1Etsum","L1 Etsum ; L1 Etsum; Entries ",40,0.,6000.);
  TH1D *hofflineEtsum = new TH1D("hofflineEtsum","Offline Etsum ; Offline Etsum; Entries ",40,0.,6000.);
  TH2D *hcorrOfflineEtsumVsL1Etsum = new TH2D("hcorrOfflineEtsumVsL1Etsum","Offline Etsum vs L1 Etsum; Offline Etsum; L1 Etsum",100,0.,6000.,100,0.,6000.);
  
  TProfile *profilel1EtsumVsofflineCentrality = new TProfile("profilel1EtsumVsofflineCentrality","L1 Etsum vs offline centrality; offline centrality ; L1 Etsum",100,0,200,0,6000);
  TProfile *profileofflineCentralityVsl1Etsum = new TProfile("profileofflineCentralityVsl1Etsum","Offline centrality vs L1 Etsum; L1 Etsum ; offline centrality",400,0,6000,0,200);
  TProfile *profileofflineCentralityVsl1Etsum_Calibration = new TProfile("profileofflineCentralityVsl1Etsum_Calibration","Offline centrality vs L1 Etsum; L1 Etsum ; offline centrality",200,0,6000,0,200);
  TProfile *profilel1EtsumVsofflineCentrality_Calibration = new TProfile("profilel1EtsumVsofflineCentrality_Calibration","Offline centrality vs L1 Etsum; L1 Etsum ; offline centrality",50,0,200,0,6000);
  TProfile *profilel1CentralityVsofflineCentrality = new TProfile("profilel1CentralityVsofflineCentrality","L1 centrality vs Offline Centrality; Offline Centrality ; L1 centrality",100,0,200,0,200);
  
  TH1D *hOffline_Bin0 = new TH1D("hOffline_Bin0","hOffline_Bin0 ; Offline Centrality; Entries ",260,-30.,230.);
  TH1D *hOffline_Bin1 = new TH1D("hOffline_Bin1","hOffline_Bin1 ; Offline Centrality; Entries ",260,-30.,230.);
  TH1D *hOffline_Bin2 = new TH1D("hOffline_Bin2","hOffline_Bin2 ; Offline Centrality; Entries ",260,-30.,230.);
  TH1D *hOffline_Bin3 = new TH1D("hOffline_Bin3","hOffline_Bin3 ; Offline Centrality; Entries ",260,-30.,230.);
  TH1D *hOffline_Bin4 = new TH1D("hOffline_Bin4","hOffline_Bin4 ; Offline Centrality; Entries ",260,-30.,230.);
  TH1D *hOffline_Bin5 = new TH1D("hOffline_Bin5","hOffline_Bin5 ; Offline Centrality; Entries ",260,-30.,230.);
  
  const int nBins = 100; 
  const double maxCen = 200.;
  const Double_t L1_THRESHOLD[12] = {0.,10.,20.,30.,40.,50.,100.,140.,160.,170.,180.,190.};
  
  TH1D *accepted[12][1];
  accepted[0][0] = new TH1D("accepted_cen_0","",nBins,0,maxCen);
  accepted[1][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_10");
  accepted[2][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_20");
  accepted[3][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_30");
  accepted[4][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_40");
  accepted[5][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_50");
  accepted[6][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_100");
  accepted[7][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_140");
  accepted[8][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_160");
  accepted[9][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_170");
  accepted[10][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_180");
  accepted[11][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_190");
  
  TH1D *fCenOffline[1]; 
  fCenOffline[0] = new TH1D("fCenOffline",";Offline Centrality",nBins,0,maxCen);

  int countCalib = 0;

  int entries = fTree->GetEntries();
    
  for(long j = 1; j < entries; ++j){
    fTree->GetEntry(j);
    long keycal = makeKey(f_run, f_evt);
    
    map<long,int>::const_iterator gotcal = kmapcal.find(keycal);
    
    if (pcollisionEventSelection==0 || pHBHENoiseFilter ==0) continue;
    
    if(gotcal == kmapcal.end()){
      continue;      
    } else {
      l1t->GetEntry(gotcal->second);
      kmapcal.erase(keycal);
      
      int etsum=l1_etsumPlus+l1_etsumMinus;
      //printf("et sum plus + minus  %d\n",etsum);
      //printf("hiBin  %d\n",hiBin);
      profileofflineCentralityVsl1Etsum_Calibration->Fill(etsum,hiBin); 
      profilel1EtsumVsofflineCentrality_Calibration->Fill(hiBin,etsum); 
      countCalib++;
    }  
  }  


  printf("Matching entries: %d\n",countCalib);

  TF1 *fprofileofflineCentralityVsl1Etsum_Calibration = new TF1("fprofileofflineCentralityVsl1Etsum_Calibration","pol9",0,2500);
  profileofflineCentralityVsl1Etsum_Calibration->Fit("fprofileofflineCentralityVsl1Etsum_Calibration");
  TF1 *fprofileofflinel1EtsumVsCentrality_Calibration = new TF1("fprofileofflinel1EtsumVsCentrality_Calibration","pol9",0,200);
  profilel1EtsumVsofflineCentrality_Calibration->Fit("fprofileofflinel1EtsumVsCentrality_Calibration");


  double L1centrality=0.;
      
  for(long j = 1; j < entries; ++j){
  
    if(j % 10000 == 0) printf("%ld / %d\n",j,entries);
    fTree->GetEntry(j);
    if (pcollisionEventSelection==0 || pHBHENoiseFilter ==0) continue;

    long key = makeKey(f_run, f_evt);

    map<long,int>::const_iterator got = kmap.find(key);
    if(got == kmap.end()){
      continue;      
    } else {
      l1t->GetEntry(got->second);  
      kmap.erase(key);
      //printf("id number %d\n",j);
      int etsum=l1_etsumPlus+l1_etsumMinus;
      printf("et value %d\n",etsum);
      printf("hiBin %d\n",hiBin);
      printf("hiHF %f\n",hiHF);
	  
	  //if (pcollisionEventSelection==0 || pHBHENoiseFilter ==0) { hcorrl1EtsumPlusVscorrl1EtsumMinusNoEvSel->Fill(l1_etsumMinus,l1_etsumPlus); continue;  }
	
	  hcorrl1EtsumPlusVscorrl1EtsumMinus->Fill(l1_etsumMinus,l1_etsumPlus);
	  hcorrl1EtsumVsofflineCentrality->Fill(hiBin,etsum);
	  hcorrofflineCentralityVsl1Etsum->Fill(etsum,hiBin);
	  hofflineCentrality->Fill(hiBin);
	  hl1Etsum->Fill(etsum);
	  hofflineEtsum->Fill(hiHF);
	  hcorrOfflineEtsumVsL1Etsum->Fill(hiHF,etsum); 
	  profileofflineCentralityVsl1Etsum->Fill(etsum,hiBin);  
	  profilel1EtsumVsofflineCentrality->Fill(hiBin,etsum);  	  
	  //L1centrality=fprofileofflineCentralityVsl1Etsum_Calibration->Eval(etsum);
	  L1centrality=fprofileofflinel1EtsumVsCentrality_Calibration->GetX(etsum);
	  printf("L1centrality %f\n",L1centrality);
	  hcorrL1CentralityVsfflineCentrality->Fill(hiBin,L1centrality);
	  profilel1CentralityVsofflineCentrality->Fill(hiBin,L1centrality);
	  
	  double Etsum_Bin0=fprofileofflinel1EtsumVsCentrality_Calibration->Eval((double)(limit0));
      double Etsum_Bin1=fprofileofflinel1EtsumVsCentrality_Calibration->Eval((double)(limit1));
      double Etsum_Bin2=fprofileofflinel1EtsumVsCentrality_Calibration->Eval((double)(limit2));
      double Etsum_Bin3=fprofileofflinel1EtsumVsCentrality_Calibration->Eval((double)(limit3));
      double Etsum_Bin4=fprofileofflinel1EtsumVsCentrality_Calibration->Eval((double)(limit4));
      
      if(etsum>Etsum_Bin0) {hOffline_Bin0->Fill(hiBin);}
      if(etsum>Etsum_Bin1 && etsum<Etsum_Bin0) {hOffline_Bin1->Fill(hiBin);}
      if(etsum>Etsum_Bin2 && etsum<Etsum_Bin1) {hOffline_Bin2->Fill(hiBin);}
      if(etsum>Etsum_Bin3 && etsum<Etsum_Bin2) {hOffline_Bin3->Fill(hiBin);}
      if(etsum>Etsum_Bin4 && etsum<Etsum_Bin3) {hOffline_Bin4->Fill(hiBin);}
      if(etsum<Etsum_Bin4) {hOffline_Bin5->Fill(hiBin);}
      
	  fCenOffline[0]->Fill(hiBin);
	  
	  for(int k = 0; k < 12; ++k){
	    if(L1centrality>L1_THRESHOLD[k]){
	      accepted[k][0]->Fill(hiBin);
	    }
	  }
    }  
  }  



  //TFile*fouput=new TFile("CalibrationEtsum_saturation.root","recreate"); 
  TFile*fouput=new TFile(Form("Results/CalibrationEtsum_CentLim_%d_%d_%d_%d_%d_8July2014.root",limit0,limit1,limit2,limit3,limit4),"recreate"); 
  fouput->cd();
  
  hcorrl1EtsumPlusVscorrl1EtsumMinus->Write();
  hcorrl1EtsumPlusVscorrl1EtsumMinusNoEvSel->Write();
  hcorrl1EtsumVsofflineCentrality->Write();
  hcorrofflineCentralityVsl1Etsum->Write();
  hofflineCentrality->Write();
  hl1Etsum->Write();
  hofflineEtsum->Write();
  hcorrOfflineEtsumVsL1Etsum->Write();
  profilel1EtsumVsofflineCentrality->Write();
  profileofflineCentralityVsl1Etsum->Write();
  profileofflineCentralityVsl1Etsum_Calibration->Write();
  fprofileofflineCentralityVsl1Etsum_Calibration->SetName("fprofileofflineCentralityVsl1Etsum_Calibration");
  fprofileofflineCentralityVsl1Etsum_Calibration->Write();
  fprofileofflinel1EtsumVsCentrality_Calibration->SetName("fprofileofflinel1EtsumVsCentrality_Calibration");
  fprofileofflinel1EtsumVsCentrality_Calibration->Write();
  hcorrL1CentralityVsfflineCentrality->Write();
  profilel1CentralityVsofflineCentrality->Write();
  hOffline_Bin0->Write();
  hOffline_Bin1->Write();
  hOffline_Bin2->Write();
  hOffline_Bin3->Write();
  hOffline_Bin4->Write();
  hOffline_Bin5->Write();
  

  TGraphAsymmErrors *a[12][1];
  for(int k = 0; k < 12; ++k){
    for(int l = 0; l < 1; ++l)
    {
      a[k][l] = new TGraphAsymmErrors();
      a[k][l]->BayesDivide(accepted[k][l],fCenOffline[l]);
    }
  }
  a[0][0]->SetName("asymm_cen_0_cen");
  a[1][0]->SetName("asymm_cen_10_cen");
  a[2][0]->SetName("asymm_cen_20_cen");
  a[3][0]->SetName("asymm_cen_30_cen");
  a[4][0]->SetName("asymm_cen_40_cen");
  a[5][0]->SetName("asymm_cen_50_cen");
  a[6][0]->SetName("asymm_cen_100_cen");
  a[7][0]->SetName("asymm_cen_140_cen");
  a[8][0]->SetName("asymm_cen_160_cen");
  a[9][0]->SetName("asymm_cen_170_cen");
  a[10][0]->SetName("asymm_cen_180_cen");
  a[11][0]->SetName("asymm_cen_190_cen");
  
  fCenOffline[0]->Write();
   
  for(int k = 0; k < 12; ++k){
    for(int l = 0; l < 1; ++l)
    {
      accepted[k][l]->Write();
      a[k][l]->Write();
    }
  }
  
  fouput->Close();
  lFile->Close();
  fFile->Close();

}



int main(){
  matching_l12forest_Centrality();
  return 0; 
}
