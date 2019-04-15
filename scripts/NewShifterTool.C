
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>
#include "TTree.h"

// RooFit
#include "RooHistPdf.h"
#include "RooArgList.h"
#include "TH1.h"

// TMVA
#include "TMVA/PDF.h"

double xlow = -5;
double xhigh = 7;
int nEvents = 10000;
int nbins = 50;
TRandom3* rnd = new TRandom3();


Float_t mcVar;

TTree* getMCTree() {
  TTree* mcTree = new TTree("mcTree","mcTree");
  TBranch *mcBranch = mcTree->Branch("mcBranch", &mcVar, "mcVar/F");
  for (int i=0;i<nEvents;++i) {
    mcVar = rnd->Gaus(0.9,1); // mu 
    mcTree->Fill();
  }
  return mcTree;
}

TH1F* getDataHist() {
  TH1F* hist = new TH1F("data","data",nbins,-5,7);
  hist->Sumw2();
  for (int i=0;i<nEvents;++i) {
    hist->Fill(rnd->Gaus(1.3,1.2));
  }
  return hist;
}

TTree* mcTree;
TH1F* dataHist;
TH1F* mcHist_th1;
TH1F* mcHist_orig;
TH1* kde;

double MakeaSingleFit(const double *xx) {

  mcHist_th1->Reset();
  for (int i=0; i<mcTree->GetEntries(); ++i) {
    mcTree->GetEntry(i);
    mcHist_th1->Fill((mcVar - xx[2])*xx[1] + xx[2] + xx[0]);
//     for (int j=0;j<29; ++j) {
//       mcHist_th1->Fill((mcVar*rnd->Gaus(1,0.05) - xx[2])*xx[1] + xx[2] + xx[0],1./30.);
//     }
  }

  //
  // Umm, maybe we need to implement KDE smoothing, so that the convergence goes better.
  // Not necessary if you fit for the shift first...?
  //
//   TMVA::PDF* TMVAPDF = new TMVA::PDF("temp", mcHist_th1
//                                      ,TMVA::KDEKernel::kGauss
//                                      ,TMVA::KDEKernel::kNonadaptiveKDE // kNonadaptiveKDE, kAdaptiveKDE
//                                      ,TMVA::KDEKernel::kSampleMirror
//                                      ,1
//                                      ,0); // normalize

//   TMVAPDF->ValidatePDF(0);
//   kde = TMVAPDF->GetPDFHist();
//   kde->Rebin(kde->GetNbinsX()/float(nbins));

  TH1* ref_hist = mcHist_th1;

//   mcHist_th1->Scale(1/100.);

  double chi2 = 0;
  for (int i=0; i<nbins; ++i){
    double e2 = (ref_hist->GetBinError(i+1)*ref_hist->GetBinError(i+1) + 
                 dataHist->GetBinError(i+1)*dataHist->GetBinError(i+1));
    if (!e2) continue;
    double diff = ref_hist->GetBinContent(i+1) - dataHist->GetBinContent(i+1);
    chi2 += diff*diff / e2 ;
  }
  
  std::cout << "Trying shift = " << xx[0]
            << ", width = " <<  xx[1]
            << ", mu = " << xx[2]
            << " has Chi2: " << chi2 << std::endl;
  return chi2/200.;
}

void NewShifterTool() {

  gROOT->SetBatch(false);

  GU::SetupStyle();

  mcTree = getMCTree();
  dataHist = getDataHist();
  mcHist_th1 = new TH1F("mcHist_th1","mcHist_th1",nbins,xlow,xhigh);
  mcHist_orig = new TH1F("mcHist_th1","mcHist_th1",nbins,xlow,xhigh);

  TCanvas* c = GU::RatioCanvas("my_ratiocanvas","blah",600,500);

  mcHist_th1->Sumw2();
  mcHist_orig->Sumw2();

  for (int i=0; i<mcTree->GetEntries(); ++i) {
    mcTree->GetEntry(i);
    mcHist_orig->Fill(mcVar);
  }
  double mean = mcHist_orig->GetMean();
  
  //MakeaSingleFit();

  const char * minName = "Minuit2";
  const char *algoName = "";

  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

  // set tolerance , etc...
//   minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
//   minimum->SetMaxIterations(10000);  // for GSL
//   minimum->SetTolerance(20.);
  minimum->SetStrategy( 0 );
  minimum->SetPrintLevel(1);

  // create function wrapper for minimizer
  // a IMultiGenFunction type
  ROOT::Math::Functor f(&MakeaSingleFit,3);
  double step[3] = {10.,10.,10.};
  // starting point

  minimum->SetFunction(f);

  std::cout << "Using mean " << mean << std::endl;
  double variable[3] = { 0.0 , 1.0 , mean };
  double starting_points[3] = {variable[0],variable[1],variable[2]};

  // Get a good starting point for the shift!
  if (true)
  {
    minimum->SetVariable(0,"shift" ,variable[0], step[0]);
    minimum->SetFixedVariable(1,"width" ,variable[1]);
    minimum->SetFixedVariable(2,"mu"    ,variable[2]);
    minimum->Minimize();
    // Save the new shift value:
    variable[0] = minimum->X()[0];
  }

  // Get a good starting point for the shift!
  if (true)
  {
    minimum->SetFixedVariable(0,"shift" ,variable[0]);
    minimum->SetVariable(1,"width" ,variable[1], step[1]);
    minimum->SetFixedVariable(2,"mu"    ,variable[2]);
    minimum->Minimize();
    // Save the new shift value:
    variable[1] = minimum->X()[1];
  }

  minimum->SetVariable(0,"shift" ,variable[0], step[0]);
  minimum->SetVariable(1,"width" ,variable[1], step[1]);
  minimum->SetVariable(2,"mu"    ,variable[2], step[2]);

  // do the minimization
  minimum->Minimize();
  minimum->PrintResults();

  const double *xs = minimum->X();
  std::cout << "Starting point: f(" << starting_points[0] << "," << starting_points[1] << "," << starting_points[2] << std::endl;
  std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] << "): chi2 = " << minimum->MinValue()  << std::endl;
  
  // Plot the result:  
  MakeaSingleFit(xs);

  c->cd();
//   mcHist_orig->Rebin(nbins/100);
//   mcHist_th1->Rebin(nbins/100);
//   dataHist->Rebin(nbins/100);

  mcHist_orig->SetLineWidth(3); mcHist_orig->SetLineColor(3); mcHist_orig->SetMarkerColor(3);
  mcHist_th1->SetLineWidth(3); mcHist_th1->SetLineColor(2); mcHist_th1->SetMarkerColor(2);
  GU::AddHistogram(*c,*dataHist,"hist");
  GU::AddRatio(*c,*mcHist_orig,*dataHist,"hist");
  GU::AddRatio(*c,*mcHist_th1,*dataHist,"hist");
  GU::FormatCanvasAxes(*c);
  GU::FullFormatCanvasDefault(*c);
  GU::SetYaxisRanges(*GU::GetBotPad(*c),0,2);
  GU::SetColors(*c);
  //mcHist_orig->Draw("hist");
//   mcHist_th1->Draw("histsame");
//   kde->Draw("histsame");
//   dataHist->Draw("pEsame");
  return 0;
}

// RooRealVar x("x","x",xlow,xhigh) ;

//   RooGaussian g1("g","g",x,RooFit::RooConst(1),RooFit::RooConst(1));
//   x.setBins(100) ;
//   RooDataSet* data = g1.generate(x,10000) ;
  

// RooGaussian gmc("g","g",x,RooFit::RooConst(1),RooFit::RooConst(0.8));
// RooDataSet* mcUnbinned = gmc.generate(x,10000) ;


  // Create a binned dataset with 20 bins and 500 events
//   RooDataHist* datahist = data->binnedClone() ;
//   RooDataHist* mcHist = new RooDataHist("mcHist","mcHist",x,RooFit::Import(*mcHist_th1));

  // Represent data in dh as pdf in x
//   RooHistPdf mcHistPdf("mcHistPdf","mcHistPdf",x,*mcHist,0) ;
//   mcHistPdf.fitTo(*datahist);

  // Plot unbinned data and histogram pdf overlaid
//   RooPlot* frame1 = x.frame(RooFit::Title("Low statistics histogram pdf"),RooFit::Bins(50)) ;
//   data->plotOn(frame1) ;
//   mcHistPdf.plotOn(frame1) ; 

//   c->cd();
//   frame1->Draw();

