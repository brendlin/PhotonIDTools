#include <vector>
#include "TString.h"

//
// This struct is designed to save histograms for each variable, in each bin of et and eta.
//
struct photonVariables
{
  photonVariables(const char* name,int n_et,int n_eta,bool isConv) {

    net = n_et;
    neta = n_eta;

    variables.push_back("CutHadLeakage");
    variables.push_back("Reta37"       );
    variables.push_back("Rphi33"       );
    variables.push_back("weta2"        );
    variables.push_back("fracm"        );
    variables.push_back("wtot"         );
    variables.push_back("w1"           );
    variables.push_back("deltae"       );
    variables.push_back("DEmaxs1"      );

    for (int et=0;et<n_et;et++) {
      for (int eta=0;eta<n_eta;eta++) {
        TString nm_tmp = Form("%s_%d_%d",name,et,eta);
        const char* nm = nm_tmp.Data();
        double rphi_lo = (isConv) ? 0.45 : 0.75;
        h_CutHadLeakage.push_back(new TH1F(Form("CutHadLeakage_%s",nm),"CutHadLeakage",100,  -0.03,     0.07));
        h_Reta37       .push_back(new TH1F(Form("Reta37_%s"       ,nm),"Reta37"       ,100,   0.80,    1.010));
        h_Rphi33       .push_back(new TH1F(Form("Rphi33_%s"       ,nm),"Rphi33"       ,100,rphi_lo,     1.00));
        h_weta2        .push_back(new TH1F(Form("weta2_%s"        ,nm),"weta2"        ,100, 0.0065, 0.015990));
        h_fracm        .push_back(new TH1F(Form("fracm_%s"        ,nm),"fracm"        ,100,  0.001,   0.7990));
        h_wtot         .push_back(new TH1F(Form("wtot_%s"         ,nm),"wtot"         ,100,      1,   4.9990));
        h_w1           .push_back(new TH1F(Form("w1_%s"           ,nm),"w1"           ,100,   0.45,    0.850));
        h_deltae       .push_back(new TH1F(Form("deltae_%s"       ,nm),"deltae"       ,100, 0.0001,   799.90));
        h_DEmaxs1      .push_back(new TH1F(Form("DEmaxs1_%s"      ,nm),"DEmaxs1"      ,100, 0.6001,   1.0490));
        h_CutHadLeakage.back()->GetXaxis()->SetTitle("R_{Had}"     );
        h_Reta37       .back()->GetXaxis()->SetTitle("R_{#eta}"    );
        h_Rphi33       .back()->GetXaxis()->SetTitle("R_{#phi}"    );
        h_weta2        .back()->GetXaxis()->SetTitle("W_{#eta^{}2}");
        h_fracm        .back()->GetXaxis()->SetTitle("f_{side}"    );
        h_wtot         .back()->GetXaxis()->SetTitle("W_{stot}"    );
        h_w1           .back()->GetXaxis()->SetTitle("W_{s3}"      );
        h_deltae       .back()->GetXaxis()->SetTitle("#Delta^{}E"  );
        h_DEmaxs1      .back()->GetXaxis()->SetTitle("E_{ratio}"   );
      }
    }
    return;
  }

  std::vector<std::string> variables;
  std::vector<TH1F*> h_CutHadLeakage;
  std::vector<TH1F*> h_Reta37       ;
  std::vector<TH1F*> h_Rphi33       ;
  std::vector<TH1F*> h_weta2        ;
  std::vector<TH1F*> h_fracm        ;
  std::vector<TH1F*> h_wtot         ;
  std::vector<TH1F*> h_w1           ;
  std::vector<TH1F*> h_deltae       ;
  std::vector<TH1F*> h_DEmaxs1      ;

  int net;
  int neta;

};


//
// This struct is designed to save all of the relevant cuts / info needed to evaluate a menu.
//
struct photonID
{
  photonID(int n_et,int n_eta) {
    net = n_et;
    neta = n_eta;
    CutHadLeakage = new Double_t[net*neta];
    Reta37        = new Double_t[net*neta];
    Rphi33        = new Double_t[net*neta];
    weta2         = new Double_t[net*neta];
    fracm         = new Double_t[net*neta];
    wtot          = new Double_t[net*neta];
    w1            = new Double_t[net*neta];
    deltae        = new Double_t[net*neta];
    DEmaxs1       = new Double_t[net*neta];
    isolation     = "";
    EtaBinThresholds       = new Double_t[neta];
    EtBinThresholds        = new Double_t[net-1];

    do_CutHadLeakage = false;
    do_Reta37        = false;
    do_Rphi33        = false;
    do_weta2         = false;
    do_fracm         = false;
    do_wtot          = false;
    do_w1            = false;
    do_deltae        = false;
    do_DEmaxs1       = false;

  }

  Double_t* CutHadLeakage;
  Double_t* Reta37;
  Double_t* Rphi33;
  Double_t* weta2;
  Double_t* fracm;
  Double_t* wtot;
  Double_t* w1;
  Double_t* deltae;
  Double_t* DEmaxs1;

  Double_t* EtaBinThresholds;
  Double_t* EtBinThresholds;

  std::string isolation;

  int net;
  int neta;

  bool do_CutHadLeakage;
  bool do_Reta37       ;
  bool do_Rphi33       ;
  bool do_weta2        ;
  bool do_fracm        ;
  bool do_wtot         ;
  bool do_w1           ;
  bool do_deltae       ;
  bool do_DEmaxs1      ;

  void Set_CutHadLeakage(Double_t* c){ do_CutHadLeakage = true; for (int i=0; i < net*neta; i++) CutHadLeakage[i] = c[i]; }
  void Set_Reta37       (Double_t* c){ do_Reta37        = true; for (int i=0; i < net*neta; i++) Reta37[i]        = c[i]; }
  void Set_Rphi33       (Double_t* c){ do_Rphi33        = true; for (int i=0; i < net*neta; i++) Rphi33[i]        = c[i]; }
  void Set_weta2        (Double_t* c){ do_weta2         = true; for (int i=0; i < net*neta; i++) weta2[i]         = c[i]; }
  void Set_fracm        (Double_t* c){ do_fracm         = true; for (int i=0; i < net*neta; i++) fracm[i]         = c[i]; }
  void Set_wtot         (Double_t* c){ do_wtot          = true; for (int i=0; i < net*neta; i++) wtot[i]          = c[i]; }
  void Set_w1           (Double_t* c){ do_w1            = true; for (int i=0; i < net*neta; i++) w1[i]            = c[i]; }
  void Set_deltae       (Double_t* c){ do_deltae        = true; for (int i=0; i < net*neta; i++) deltae[i]        = c[i]; }
  void Set_DEmaxs1      (Double_t* c){ do_DEmaxs1       = true; for (int i=0; i < net*neta; i++) DEmaxs1[i]       = c[i]; }

  // No lower bin edge (0); Yes upper bin edge (2.47)
  void Set_EtaBinThresholds(Double_t* c){for (int i=0; i < neta; i++) EtaBinThresholds[i] = c[i]; }
  // No lower bin edge (0); No upper bin edge (inf).
  void Set_EtBinThresholds (Double_t* c){for (int i=0; i < net-1 ; i++) EtBinThresholds [i] = c[i]; }
};

void EvaluatePhotonID_InclusivePhoton(TTree* tree, photonID* iddef,bool doConv, TH2* denominator, TH2* numerator,photonVariables* hists = NULL){

  float y_pt, y_eta_cl_s2,y_Reta,y_Rphi,y_weta2,y_fracs1,y_weta1,y_f1,y_wtots1,y_Rhad,y_Rhad1;
  float y_Eratio,y_e277,y_deltae;
  int y_convType;
  bool y_isTruthMatchedPhoton;
  double mcTotWeightNoPU_PIDuse;
  float mc_weight_pu,mc_weight_gen;

  bool isRz = tree->GetListOfBranches()->FindObject("ph.pt");

  tree->SetBranchAddress(isRz ? "ph.pt"      : "y_pt"       ,&y_pt  );
  tree->SetBranchAddress(isRz ? "ph.eta2"    : "y_eta_cl_s2",&y_eta_cl_s2);
  tree->SetBranchAddress(isRz ? "ph.reta"    : "y_Reta"     ,&y_Reta  );
  tree->SetBranchAddress(isRz ? "ph.rphi"    : "y_Rphi"     ,&y_Rphi  );
  tree->SetBranchAddress(isRz ? "ph.weta2"   : "y_weta2"    ,&y_weta2 );
  tree->SetBranchAddress(isRz ? "ph.fside"   : "y_fracs1"   ,&y_fracs1);
  tree->SetBranchAddress(isRz ? "ph.wstot"   : "y_wtots1"   ,&y_wtots1);
  tree->SetBranchAddress(isRz ? "ph.w1"      : "y_weta1"    ,&y_weta1 );
  tree->SetBranchAddress(isRz ? "ph.deltae"  : "y_deltae"   ,&y_deltae);
  tree->SetBranchAddress(isRz ? "ph.eratio"  : "y_Eratio"   ,&y_Eratio);
  tree->SetBranchAddress(isRz ? "ph.rhad1"   : "y_Rhad1"    ,&y_Rhad1 );
  tree->SetBranchAddress(isRz ? "ph.rhad"    : "y_Rhad"     ,&y_Rhad  );
  tree->SetBranchAddress(isRz ? "ph.convFlag": "y_convType" ,&y_convType);
  // f1
  // e277
  if (isRz) {
    tree->SetBranchAddress("mc_weight.pu" ,&mc_weight_pu);
    tree->SetBranchAddress("mc_weight.gen",&mc_weight_gen);
  } else {
    tree->SetBranchAddress("y_isTruthMatchedPhoton",&y_isTruthMatchedPhoton);
    tree->SetBranchAddress("mcTotWeightNoPU_PIDuse",&mcTotWeightNoPU_PIDuse);
  }

  TAxis* xaxis = denominator->GetXaxis();
  TAxis* yaxis = denominator->GetYaxis();

  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
    tree->GetEntry(i);
    if (!(i%1000000)) std::cout << Form("Entry %i million",int(i/1000000.)) << std::endl;

    float the_pt = isRz ? y_pt/1000. : y_pt;
    
    // Binning for the cuts
    int etbin_cuts = 0; // start by assuming 0th bin
    int etabin_cuts = 0; // start by assuming 0th bin
    for (int et=0; et < iddef->net-1; et++) {
      if (the_pt < 0.001 * iddef->EtBinThresholds[et]) break;
      etbin_cuts = et+1;
    }
    if (fabs(y_eta_cl_s2) > iddef->EtaBinThresholds[iddef->neta-1] ) continue;
    for (int eta=0; eta < iddef->neta-1; eta++) {
      if (fabs(y_eta_cl_s2) < iddef->EtaBinThresholds[eta]) break;
      etabin_cuts = eta+1;
    }

    int cutbin_cuts = etbin_cuts*(iddef->neta) + etabin_cuts;

    // Binning for the histograms
    int etbin_hist = 9999;
    int etabin_hist = 9999;
    double etval_hist = -9999999;
    double etaval_hist = -9999999;
    if (the_pt <  xaxis->GetBinLowEdge(1)) continue;
    if (the_pt >= xaxis->GetBinLowEdge(denominator->GetNbinsX()+1)) continue;

    for (int x=1;x<denominator->GetNbinsX()+1;x++) {
      if (xaxis->GetBinLowEdge(x) <= the_pt && the_pt < xaxis->GetBinLowEdge(x+1))
      {
        etbin_hist = x-1; // start at 0
        etval_hist = xaxis->GetBinLowEdge(x);
        continue;
      }
    }

    for (int y=1;y<denominator->GetNbinsY()+1;y++) {
      if (yaxis->GetBinLowEdge(y) <= fabs(y_eta_cl_s2) && fabs(y_eta_cl_s2) < yaxis->GetBinLowEdge(y+1))
      {
        etabin_hist = y-1; // start at 0
        etaval_hist = yaxis->GetBinLowEdge(y);
        continue;
      }
    }

    int cutbin_hist = etbin_hist*(denominator->GetNbinsY()) + etabin_hist;

    // Phase space: converted or unconverted:
    if (doConv == (y_convType == 0)) continue;

    // Fill denominator histogram
    bool passDen = true;
    if (!isRz) passDen = passDen && y_isTruthMatchedPhoton;
    if (!passDen) continue;
    
    double weight = isRz ? mc_weight_pu*mc_weight_gen : mcTotWeightNoPU_PIDuse;

    double the_rhad = y_Rhad1;
    if ( 0.8 < fabs(y_eta_cl_s2) && fabs(y_eta_cl_s2) < 1.37 ) the_rhad = y_Rhad;

    denominator->Fill(etval_hist,etaval_hist,weight);

    if (hists) {
      hists->h_CutHadLeakage[cutbin_hist]->Fill(the_rhad,weight);
      hists->h_Reta37       [cutbin_hist]->Fill(y_Reta  ,weight);
      hists->h_Rphi33       [cutbin_hist]->Fill(y_Rphi  ,weight);
      hists->h_weta2        [cutbin_hist]->Fill(y_weta2 ,weight);
      hists->h_fracm        [cutbin_hist]->Fill(y_fracs1,weight);
      hists->h_wtot         [cutbin_hist]->Fill(y_wtots1,weight);
      hists->h_w1           [cutbin_hist]->Fill(y_weta1 ,weight);
      hists->h_deltae       [cutbin_hist]->Fill(y_deltae,weight);
      hists->h_DEmaxs1      [cutbin_hist]->Fill(y_Eratio,weight);
    }

    bool passNum = true;
    if (iddef->do_CutHadLeakage) passNum = passNum && the_rhad < iddef->CutHadLeakage[cutbin_cuts];
    if (iddef->do_Reta37       ) passNum = passNum && y_Reta   > iddef->Reta37[cutbin_cuts];
    if (iddef->do_Rphi33       ) passNum = passNum && y_Rphi   > iddef->Rphi33[cutbin_cuts];
    if (iddef->do_weta2        ) passNum = passNum && y_weta2  < iddef->weta2[cutbin_cuts];
    if (iddef->do_fracm        ) passNum = passNum && y_fracs1 < iddef->fracm[cutbin_cuts];
    if (iddef->do_wtot         ) passNum = passNum && y_wtots1 < iddef->wtot[cutbin_cuts];
    if (iddef->do_w1           ) passNum = passNum && y_weta1  < iddef->w1[cutbin_cuts];
    if (iddef->do_deltae       ) passNum = passNum && y_deltae < iddef->deltae[cutbin_cuts];
    if (iddef->do_DEmaxs1      ) passNum = passNum && y_Eratio > iddef->DEmaxs1[cutbin_cuts];
    if (!passNum) continue;

    // Fill numerator histogram
    numerator->Fill(etval_hist,etaval_hist,weight);
    
  } // End loop over entries
  
  return;
}
