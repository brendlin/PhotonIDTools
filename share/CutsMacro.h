#include <vector>

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

  int net;
  int neta;

  void Set_CutHadLeakage(Double_t* c){for (int i=0; i < net*neta; i++) CutHadLeakage[i] = c[i]; }
  void Set_Reta37       (Double_t* c){for (int i=0; i < net*neta; i++) Reta37[i]        = c[i]; }
  void Set_Rphi33       (Double_t* c){for (int i=0; i < net*neta; i++) Rphi33[i]        = c[i]; }
  void Set_weta2        (Double_t* c){for (int i=0; i < net*neta; i++) weta2[i]         = c[i]; }
  void Set_fracm        (Double_t* c){for (int i=0; i < net*neta; i++) fracm[i]         = c[i]; }
  void Set_wtot         (Double_t* c){for (int i=0; i < net*neta; i++) wtot[i]          = c[i]; }
  void Set_w1           (Double_t* c){for (int i=0; i < net*neta; i++) w1[i]            = c[i]; }
  void Set_deltae       (Double_t* c){for (int i=0; i < net*neta; i++) deltae[i]        = c[i]; }
  void Set_DEmaxs1      (Double_t* c){for (int i=0; i < net*neta; i++) DEmaxs1[i]       = c[i]; }
};

void EvaluatePhotonID_InclusivePhoton(TTree* tree, photonID* iddef,bool doConv, TH2* denominator, TH2* numerator){
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
    
    int etbin,etabin;
    double etval,etaval;
    if (the_pt < xaxis->GetBinLowEdge(1)) continue;
    if (the_pt > xaxis->GetBinLowEdge(denominator->GetNbinsX()+1)) continue;
    
    for (int x=1;x<denominator->GetNbinsX()+1;x++) {
      if (xaxis->GetBinLowEdge(x) < the_pt && the_pt < xaxis->GetBinLowEdge(x+1))
      {
        etbin = x-1; // start at 0
        etval = xaxis->GetBinLowEdge(x);
        continue;
      }
    }

    for (int y=1;y<denominator->GetNbinsY()+1;y++) {
      if (yaxis->GetBinLowEdge(y) < fabs(y_eta_cl_s2) && fabs(y_eta_cl_s2) < yaxis->GetBinLowEdge(y+1))
      {
        etabin = y-1; // start at 0
        etaval = yaxis->GetBinLowEdge(y);
        continue;
      }
    }

    int cutbin = etbin*(iddef->neta) + etabin;

    // Phase space: converted or unconverted:
    if (doConv == (y_convType == 0)) continue;

    // Fill denominator histogram
    bool passDen = true;
    if (!isRz) passDen = passDen && y_isTruthMatchedPhoton;
    if (!passDen) continue;
    
    double weight = isRz ? mc_weight_pu*mc_weight_gen : mcTotWeightNoPU_PIDuse;

    denominator->Fill(etval,etaval,weight);

    bool passNum = true;
    double the_rhad = y_Rhad1;
    if ( 0.8 < fabs(y_eta_cl_s2) && fabs(y_eta_cl_s2) < 1.37 ) the_rhad = y_Rhad;
    passNum = passNum && the_rhad < iddef->CutHadLeakage[cutbin];
    passNum = passNum && y_Reta   > iddef->Reta37[cutbin];
    passNum = passNum && y_Rphi   > iddef->Rphi33[cutbin];
    passNum = passNum && y_weta2  < iddef->weta2[cutbin];
    passNum = passNum && y_fracs1 < iddef->fracm[cutbin];
    passNum = passNum && y_wtots1 < iddef->wtot[cutbin];
    passNum = passNum && y_weta1  < iddef->w1[cutbin];
    passNum = passNum && y_deltae < iddef->deltae[cutbin];
    passNum = passNum && y_Eratio > iddef->DEmaxs1[cutbin];
    if (!passNum) continue;

    // Fill numerator histogram
    numerator->Fill(etval,etaval,weight);
    
  } // End loop over entries
  
  return;
}
