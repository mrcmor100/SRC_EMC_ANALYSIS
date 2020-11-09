#include <TMath.h>
// Macro to make plots comparing weighted MC events and data
// Histograms were produced via. mc_reweight.C
// Author: Eric Pooser, pooser@jlab.org

bool doSub = false;

// Input ROOT file created via mc_reweight.C
TFile *inFile;
// Input ROOT file directories
TDirectory *mcWgtDir, *dataDir, *subDataDir, *alumDir;
// Weighted 1D MC histos
TH1D *h_xFocalMCWgt, *h_xpFocalMCWgt, *h_yFocalMCWgt, *h_ypFocalMCWgt;
TH1D *h_yTarMCWgt, *h_xpTarMCWgt, *h_ypTarMCWgt, *h_deltaMCWgt;
TH1D *h_thetaMCWgt, *h_q2MCWgt, *h_w2MCWgt, *h_xbjMCWgt;
// Weighted 2D MC histos
TH2D *h2_xVxpFocalMCWgt, *h2_xVyFocalMCWgt, *h2_xVypFocalMCWgt;
TH2D *h2_xpVyFocalMCWgt, *h2_xpVypFocalMCWgt, *h2_yVypFocalMCWgt;
// Data 1D histos
TH1D *h_xFocalData, *h_xpFocalData, *h_yFocalData, *h_ypFocalData;
TH1D *h_yTarData, *h_xpTarData, *h_ypTarData, *h_deltaData;
TH1D *h_thetaData, *h_q2Data, *h_w2Data, *h_xbjData;
TH1D *h_yTarRatio, *h_xpTarRatio, *h_ypTarRatio, *h_deltaRatio;
TH1D *h_thetaRatio, *h_q2Ratio, *h_w2Ratio, *h_xFocalRatio, *h_yFocalRatio;
TH1D *h_xbjRatio;
//1D Al histos
TH1D *h_xFocalAlData, *h_xpFocalAlData, *h_yFocalAlData, *h_ypFocalAlData;
TH1D *h_yTarAlData, *h_xpTarAlData, *h_ypTarAlData, *h_deltaAlData;
TH1D *h_thetaAlData, *h_q2AlData, *h_w2AlData, *h_xbjAlData;

// Data 2D histos
TH2D *h2_xVxpFocalData, *h2_xVyFocalData, *h2_xVypFocalData;
TH2D *h2_xpVyFocalData, *h2_xpVypFocalData, *h2_yVypFocalData;
// Comparison canvas'
TCanvas *c_tarComp, *c_kinComp, *c_focalComp, *c_tarComp2;
// Legends
TLegend *l_yTarComp, *l_xpTarComp, *l_ypTarComp, *l_deltaComp;
TLegend *l_thetaComp, *l_q2Comp, *l_w2Comp, *l_xbjComp;
TLegend *l_xFocalComp, *l_yFocalComp;
// Lines
TLine *ln_yTarComp, *ln_xpTarComp, *ln_ypTarComp, *ln_deltaComp;
TLine *ln_thetaComp, *ln_q2Comp, *ln_w2Comp, *ln_xbjComp;
TLine *ln_xFocalComp, *ln_yFocalComp;

// Define constants
// Canvas size paramters
static const Double_t canWidth   = 1600.0;
static const Double_t canHeight  = 800.0;
// X-axis limits
static const Double_t yTarXMin   = -7.0;
static const Double_t yTarXMax   = 7.0;
static const Double_t xFocalXMin = -30.0;
static const Double_t xFocalXMax = 45.0;
static const Double_t yFocalXMin = -30.0;
static const Double_t yFocalXMax = 45.0;
static const Double_t xpTarXMin  = -60.0;
static const Double_t xpTarXMax  = 60.0;
static const Double_t ypTarXMin  = -40.0;
static const Double_t ypTarXMax  = 40.0;
static const Double_t deltaXMin  = -10.0;
static const Double_t deltaXMax  = 22.0;
static const Double_t thetaXMin  = 11;
static const Double_t thetaXMax  = 13;
static const Double_t q2XMin     = 2.;
static const Double_t q2XMax     = 6.;
static const Double_t w2XMin     = 9.;
static const Double_t w2XMax     = 14;
static const Double_t xbjXMin    = 0.10;
static const Double_t xbjXMax    = 0.5;
// Y-axis limits
static const Double_t yTarYMin   = 0.0;
static const Double_t yTarYMax   = 45000.0;
static const Double_t xpTarYMin  = 0.0;
static const Double_t xpTarYMax  = 50.0;
static const Double_t ypTarYMin  = 0.0;
static const Double_t ypTarYMax  = 150.0;
static const Double_t deltaYMin  = 0.0;
static const Double_t deltaYMax  = 240.0;
static const Double_t thetaYMin  = 0.0;
static const Double_t thetaYMax  = 110.0;
static const Double_t q2YMin     = 0.0;
static const Double_t q2YMax     = 80.0;
static const Double_t w2YMin     = 0.0;
static const Double_t w2YMax     = 85.0;
static const Double_t xbjYMin    = 0.0;
static const Double_t xbjYMax    = 110.0;
// Ratio limits
static const Double_t idealRatio = 1.0;
static const Double_t ratioMin   = 0.0;
static const Double_t ratioMax   = 2.0;

// Scale factor for MC events to match data
//static const Double_t scaleFactor = 371;
//static const Double_t scaleFactor = 157.5;

//const char* mcComparisonFile = "/volatile/hallc/xem2/cmorean/mc-single-arm/comparison/kin%d_m%dc_%s.root";
//const char* mcComparisonFile = "/group/c-xem2/cmorean/E12-10-008/monteCarlos/comparison/kin%d_m%dc_%s.root";
const char* mcComparisonFile = "./dataToMC/kin%d_m%dc_%s.root";
double d2r    = 3.14159 / 180.;
double beamEnergy, hsec, thetac;


double cal_ep (double momentum, double delta) {
  return (momentum * (1 + delta / 100));
}

double cal_w2 (double Mp, double eb, double thetai, double ep) {
  return (Mp*Mp + 2.* Mp*(eb - ep) - 4.*eb*ep*(TMath::Power(TMath::Sin(0.5* d2r*thetai),2)));
}

double cal_xbj(double eb, double thetai, double Mp, double ep) {
  return ((eb*ep*(1.0 - TMath::Cos(d2r*thetai)))/(Mp*(eb - ep)));
}
double cal_xbj2(double eb, double ep, double thetai, double Mp) {
  return (4*eb*ep * pow(sin(d2r*thetai / 2),2) / (2*Mp*(eb-ep)));
}

std::fstream& GetTblLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

const char* kinFileDir = "/volatile/hallc/xem2/cmorean/runlist/kinSettings.dat";
static const Double_t protonMass    = 0.938272;    // GeV

void compPlots(int kin , const char* targ, int model, Double_t scaleFactor = 1) {

  if(std::strcmp(targ,"H1")==0) {doSub = true;}
  else if(std::strcmp(targ,"D2")==0) {doSub=true;}
  else {doSub = false; }

//============================================
//   Input Kinematics Settings for MC
//============================================
  string kinFile         = kinFileDir;
  float hsec, thetac;
  fstream kinTblFile(kinFile);
  //Kin0 has a 2% offset in hsec.  
  GetTblLine(kinTblFile, kin+1);  //Changes pointer.
  kinTblFile >> beamEnergy >> hsec >> thetac;
  cout << beamEnergy << " " << hsec << " " << thetac << endl;


  const char* csTableDir = "/volatile/hallc/xem2/cmorean/cs_tables/kin%d_%s_m%d.out";
  string csTableFile     = Form(csTableDir,kin, targ, model);

  TGraph2D *born = new TGraph2D();
  TGraph2D *xt = new TGraph2D();
  TGraph2D *cst = new TGraph2D();
  
  //variables to hold elements in each line
  double beamEnergy, ep, theta, x_bj, Q2, amu, radCorrFactor, C_cor;
  double Sig_Born, Sig_Born_In, Sig_Born_QE,
    Sig_Rad, Sig_Rad_EL, Sig_Rad_QE, Sig_Rad_DIS;
  int i = 0;
  string line;
  ifstream crossFile(csTableFile);

  //Open file
  while(getline(crossFile, line)) {
    if(i == 0) {
      i+=1;
      continue;
    }
    istringstream ss(line);
    i+=1;
    ss >> beamEnergy >> ep >> theta  >> x_bj >> Q2 >>
      Sig_Born >> Sig_Born_In >> Sig_Born_QE >> Sig_Rad >>
      Sig_Rad_EL >> Sig_Rad_QE >> Sig_Rad_DIS >> C_cor;
    
    if(i%50 == 0) {
      cout << ep << " " << theta << " " << x_bj << " " << Sig_Rad  << endl;
    }
    if(x_bj < 1.95) {
      born->SetPoint(i, ep, theta, Sig_Born);
      xt->SetPoint(i, ep, theta, x_bj);
      cst->SetPoint(i, ep, theta, Sig_Rad);
    }
  }
  
  //Aluminum Scale Factors
  Bool_t cryo = false;
  // Global ROOT settings
  gStyle->SetOptStat(0);
  
  // Obtain the comparison ROOT file
  inFile = new TFile(Form(mcComparisonFile,kin,model, targ), "UPDATE");
  // Obtain the directories which contain the histograms of interest
  mcWgtDir = dynamic_cast <TDirectory*> (inFile->FindObjectAny("mcWgtDir"));
  dataDir  = dynamic_cast <TDirectory*> (inFile->FindObjectAny("dataDir"));
  subDataDir  = dynamic_cast <TDirectory*> (inFile->FindObjectAny("subDataDir"));
  alumDir  = dynamic_cast <TDirectory*> (inFile->FindObjectAny("alumDir"));
  //dataDir  = dynamic_cast <TDirectory*> (inFile->GetMotherDir(););

  // Obtain the 1D histograms of interest
  h_yTarMCWgt  = dynamic_cast <TH1D*> (mcWgtDir->FindObjectAny("h_yTarMCWgt"));
  h_xFocalMCWgt  = dynamic_cast <TH1D*> (mcWgtDir->FindObjectAny("h_xFocalMCWgt"));
  h_yFocalMCWgt  = dynamic_cast <TH1D*> (mcWgtDir->FindObjectAny("h_yFocalMCWgt"));
  h_xpTarMCWgt = dynamic_cast <TH1D*> (mcWgtDir->FindObjectAny("h_xpTarMCWgt"));
  h_ypTarMCWgt = dynamic_cast <TH1D*> (mcWgtDir->FindObjectAny("h_ypTarMCWgt"));
  h_deltaMCWgt = dynamic_cast <TH1D*> (mcWgtDir->FindObjectAny("h_deltaMCWgt"));
  h_thetaMCWgt = dynamic_cast <TH1D*> (mcWgtDir->FindObjectAny("h_thetaMCWgt"));
  h_q2MCWgt    = dynamic_cast <TH1D*> (mcWgtDir->FindObjectAny("h_q2MCWgt"));
  h_w2MCWgt    = dynamic_cast <TH1D*> (mcWgtDir->FindObjectAny("h_w2MCWgt"));
  h_xbjMCWgt    = dynamic_cast <TH1D*> (mcWgtDir->FindObjectAny("h_xbjMCWgt"));

  if(doSub==false) {
    h_yTarData   = dynamic_cast <TH1D*> (dataDir->FindObjectAny("h_yTarData"));
    h_xFocalData = dynamic_cast <TH1D*> (dataDir->FindObjectAny("h_xFocalData"));
    h_yFocalData = dynamic_cast <TH1D*> (dataDir->FindObjectAny("h_yFocalData"));
    h_xpTarData  = dynamic_cast <TH1D*> (dataDir->FindObjectAny("h_xpTarData"));
    h_ypTarData  = dynamic_cast <TH1D*> (dataDir->FindObjectAny("h_ypTarData"));
    h_deltaData  = dynamic_cast <TH1D*> (dataDir->FindObjectAny("h_deltaData"));
    h_thetaData  = dynamic_cast <TH1D*> (dataDir->FindObjectAny("h_thetaData"));
    h_q2Data     = dynamic_cast <TH1D*> (dataDir->FindObjectAny("h_q2Data"));
    h_w2Data     = dynamic_cast <TH1D*> (dataDir->FindObjectAny("h_w2Data"));
    h_xbjData    = dynamic_cast <TH1D*> (dataDir->FindObjectAny("h_xbjData"));
  } else {
    h_yTarData   = dynamic_cast <TH1D*> (subDataDir->FindObjectAny("h_yTarSubData"));
    h_xFocalData = dynamic_cast <TH1D*> (subDataDir->FindObjectAny("h_xFocalSubData"));
    h_yFocalData = dynamic_cast <TH1D*> (subDataDir->FindObjectAny("h_yFocalSubData"));
    h_xpTarData  = dynamic_cast <TH1D*> (subDataDir->FindObjectAny("h_xpTarSubData"));
    h_ypTarData  = dynamic_cast <TH1D*> (subDataDir->FindObjectAny("h_ypTarSubData"));
    h_deltaData  = dynamic_cast <TH1D*> (subDataDir->FindObjectAny("h_deltaSubData"));
    h_thetaData  = dynamic_cast <TH1D*> (subDataDir->FindObjectAny("h_thetaSubData"));
    h_q2Data     = dynamic_cast <TH1D*> (subDataDir->FindObjectAny("h_q2SubData"));
    h_w2Data     = dynamic_cast <TH1D*> (subDataDir->FindObjectAny("h_w2SubData"));
    h_xbjData    = dynamic_cast <TH1D*> (subDataDir->FindObjectAny("h_xbjSubData"));

    h_yTarAlData   = dynamic_cast <TH1D*> (alumDir->FindObjectAny("h_yTarAlData"));
    h_xFocalAlData = dynamic_cast <TH1D*> (alumDir->FindObjectAny("h_xFocalAlData"));
    h_yFocalAlData = dynamic_cast <TH1D*> (alumDir->FindObjectAny("h_yFocalAlData"));
    h_xpTarAlData  = dynamic_cast <TH1D*> (alumDir->FindObjectAny("h_xpTarAlData"));
    h_ypTarAlData  = dynamic_cast <TH1D*> (alumDir->FindObjectAny("h_ypTarAlData"));
    h_deltaAlData  = dynamic_cast <TH1D*> (alumDir->FindObjectAny("h_deltaAlData"));
    h_thetaAlData  = dynamic_cast <TH1D*> (alumDir->FindObjectAny("h_thetaAlData"));
    h_q2AlData     = dynamic_cast <TH1D*> (alumDir->FindObjectAny("h_q2AlData"));
    h_w2AlData     = dynamic_cast <TH1D*> (alumDir->FindObjectAny("h_w2AlData"));
    h_xbjAlData    = dynamic_cast <TH1D*> (alumDir->FindObjectAny("h_xbjAlData"));

  }
  // Obtain the 2D histos of interest
  h2_xVxpFocalMCWgt  = dynamic_cast <TH2D*> (mcWgtDir->FindObjectAny("h2_xVxpFocalMCWgt"));
  h2_xVyFocalMCWgt   = dynamic_cast <TH2D*> (mcWgtDir->FindObjectAny("h2_xVyFocalMCWgt"));
  h2_xVypFocalMCWgt  = dynamic_cast <TH2D*> (mcWgtDir->FindObjectAny("h2_xVypFocalMCWgt"));
  h2_xpVyFocalMCWgt  = dynamic_cast <TH2D*> (mcWgtDir->FindObjectAny("h2_xpVyFocalMCWgt"));
  h2_xpVypFocalMCWgt = dynamic_cast <TH2D*> (mcWgtDir->FindObjectAny("h2_xpVypFocalMCWgt"));
  h2_yVypFocalMCWgt  = dynamic_cast <TH2D*> (mcWgtDir->FindObjectAny("h2_yVypFocalMCWgt"));
  h2_xVxpFocalData   = dynamic_cast <TH2D*> (dataDir->FindObjectAny("h2_xVxpFocalData"));
  h2_xVyFocalData    = dynamic_cast <TH2D*> (dataDir->FindObjectAny("h2_xVyFocalData"));
  h2_xVypFocalData   = dynamic_cast <TH2D*> (dataDir->FindObjectAny("h2_xVypFocalData"));
  h2_xpVyFocalData   = dynamic_cast <TH2D*> (dataDir->FindObjectAny("h2_xpVyFocalData"));
  h2_xpVypFocalData  = dynamic_cast <TH2D*> (dataDir->FindObjectAny("h2_xpVypFocalData"));
  h2_yVypFocalData   = dynamic_cast <TH2D*> (dataDir->FindObjectAny("h2_yVypFocalData"));
  // Clone histos for calculating ratios
  h_yTarRatio  = dynamic_cast <TH1D*> (h_yTarData->Clone("h_yTarData"));
  h_xFocalRatio= dynamic_cast <TH1D*> (h_xFocalData->Clone("h_xFocalData"));
  h_yFocalRatio= dynamic_cast <TH1D*> (h_yFocalData->Clone("h_yFocalData"));
  h_xpTarRatio = dynamic_cast <TH1D*> (h_xpTarData->Clone("h_xpTarData"));
  h_ypTarRatio = dynamic_cast <TH1D*> (h_ypTarData->Clone("h_ypTarData"));
  h_deltaRatio = dynamic_cast <TH1D*> (h_deltaData->Clone("h_deltaData"));
  h_thetaRatio = dynamic_cast <TH1D*> (h_thetaData->Clone("h_thetaData"));
  h_q2Ratio    = dynamic_cast <TH1D*> (h_q2Data->Clone("h_q2Data"));
  h_w2Ratio    = dynamic_cast <TH1D*> (h_w2Data->Clone("h_w2Data"));
  h_xbjRatio   = dynamic_cast <TH1D*> (h_xbjData->Clone("h_xbjData"));

  c_tarComp = new TCanvas("c_tarComp", "Comparison of Target Variables", canWidth, canHeight); 
  c_tarComp->Divide(3, 2);
  c_tarComp->cd(1);
  h_xpTarMCWgt->Scale(scaleFactor);
  h_xpTarMCWgt->SetOption("HIST");
  h_xpTarMCWgt->SetFillColor(8);
  h_xpTarMCWgt->SetLineColor(8);
  h_xpTarMCWgt->SetFillStyle(3001);
  h_xpTarMCWgt->GetXaxis()->SetRangeUser(xpTarXMin, xpTarXMax);
  //h_xpTarMCWgt->GetYaxis()->SetRangeUser(xpTarYMin, xpTarYMax);
  h_xpTarMCWgt->SetTitle("Weighted MC to Data Comparison: X'_{tar}");
  h_xpTarMCWgt->Draw();

  h_xpTarData->SetMarkerColor(4);
  h_xpTarData->SetMarkerStyle(22);
  h_xpTarData->Draw("E, SAME");
  if(doSub==true) {
    h_xpTarAlData->SetMarkerColor(5);
    h_xpTarAlData->SetMarkerStyle(22);
    h_xpTarAlData->Draw("E, SAME");
  }

  l_xpTarComp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_xpTarComp->AddEntry(h_xpTarMCWgt, "Weighted MC", "F");
  l_xpTarComp->AddEntry(h_xpTarData, "Data", "EP");
  if(doSub==true) {
    l_xpTarComp->AddEntry(h_xpTarAlData, "Data", "EP");
  }
  l_xpTarComp->Draw();

  c_tarComp->cd(2);
  h_ypTarMCWgt->Scale(scaleFactor);
  h_ypTarMCWgt->SetOption("HIST");
  h_ypTarMCWgt->SetFillColor(8);
  h_ypTarMCWgt->SetLineColor(8);
  h_ypTarMCWgt->SetFillStyle(3001);
  h_ypTarMCWgt->GetXaxis()->SetRangeUser(ypTarXMin, ypTarXMax);
  //h_ypTarMCWgt->GetYaxis()->SetRangeUser(ypTarYMin, ypTarYMax);
  h_ypTarMCWgt->SetTitle("Weighted MC to Data Comparison: Y'_{tar}");
  h_ypTarMCWgt->Draw();

  h_ypTarData->SetMarkerColor(4);
  h_ypTarData->SetMarkerStyle(22);
  h_ypTarData->Draw("E, SAME");
  if(doSub==true) {
    h_ypTarAlData->SetMarkerColor(5);
    h_ypTarAlData->SetMarkerStyle(22);
    h_ypTarAlData->Draw("E, SAME");
  }

  l_ypTarComp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_ypTarComp->AddEntry(h_ypTarMCWgt, "Weighted MC", "F");
  l_ypTarComp->AddEntry(h_ypTarData, "Data", "EP");
  if(doSub==true) {
    l_ypTarComp->AddEntry(h_ypTarAlData, "Data", "EP");
  }
  l_ypTarComp->Draw();

  c_tarComp->cd(3);
  h_deltaMCWgt->Scale(scaleFactor);
  h_deltaMCWgt->SetOption("HIST");
  h_deltaMCWgt->SetFillColor(8);
  h_deltaMCWgt->SetLineColor(8);
  h_deltaMCWgt->SetFillStyle(3001);
  h_deltaMCWgt->GetXaxis()->SetRangeUser(deltaXMin, deltaXMax);
  //h_deltaMCWgt->GetYaxis()->SetRangeUser(deltaYMin, deltaYMax);
  h_deltaMCWgt->SetTitle("Weighted MC to Data Comparison: #delta");
  h_deltaMCWgt->Draw();

  h_deltaData->SetMarkerColor(4);
  h_deltaData->SetMarkerStyle(22);
  h_deltaData->Draw("E, SAME");
  if(doSub==true) {
    h_deltaAlData->SetMarkerColor(5);
    h_deltaAlData->SetMarkerStyle(22);
    h_deltaAlData->Draw("E, SAME");
  }

  l_deltaComp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_deltaComp->AddEntry(h_deltaMCWgt, "Weighted MC", "F");
  l_deltaComp->AddEntry(h_deltaData, "Data", "EP");
  if(doSub==true) {
    l_deltaComp->AddEntry(h_deltaAlData, "Data", "EP");
  }
  l_deltaComp->Draw();

  c_tarComp->cd(4);
  h_xpTarRatio->Divide(h_xpTarMCWgt);
  h_xpTarRatio->SetMarkerStyle(22);
  h_xpTarRatio->SetMarkerSize(1.25);
  h_xpTarRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_xpTarRatio->GetXaxis()->SetRangeUser(xpTarXMin, xpTarXMax);
  h_xpTarRatio->GetYaxis()->SetTitle("X'_{tar}^{data} / X'_{tar}^{MC}");
  h_xpTarRatio->SetTitle("Ratio of X'_{tar}^{data} to X'_{tar}^{MC}");
  h_xpTarRatio->Draw();

  ln_xpTarComp = new TLine(xpTarXMin, idealRatio, xpTarXMax, idealRatio);
  ln_xpTarComp->SetLineStyle(9);
  ln_xpTarComp->SetLineWidth(2);
  ln_xpTarComp->SetLineColor(38);
  ln_xpTarComp->Draw();

  c_tarComp->cd(5);
  h_ypTarRatio->Divide(h_ypTarMCWgt);
  h_ypTarRatio->SetMarkerStyle(22);
  h_ypTarRatio->SetMarkerSize(1.25);
  h_ypTarRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_ypTarRatio->GetXaxis()->SetRangeUser(ypTarXMin, ypTarXMax);
  h_ypTarRatio->GetYaxis()->SetTitle("Y'_{tar}^{data} / Y'_{tar}^{MC}");
  h_ypTarRatio->SetTitle("Ratio of Y'_{tar}^{data} to Y'_{tar}^{MC}");
  h_ypTarRatio->Draw();

  ln_ypTarComp = new TLine(ypTarXMin, idealRatio, ypTarXMax, idealRatio);
  ln_ypTarComp->SetLineStyle(9);
  ln_ypTarComp->SetLineWidth(2);
  ln_ypTarComp->SetLineColor(38);
  ln_ypTarComp->Draw();

  c_tarComp->cd(6);
  h_deltaRatio->Divide(h_deltaMCWgt);
  h_deltaRatio->SetMarkerStyle(22);
  h_deltaRatio->SetMarkerSize(1.25);
  h_deltaRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_deltaRatio->GetXaxis()->SetRangeUser(deltaXMin, deltaXMax);
  h_deltaRatio->GetYaxis()->SetTitle("#delta_{data} / #delta_{MC}");
  h_deltaRatio->SetTitle("Ratio of #delta_{data} to #delta_{MC}");
  h_deltaRatio->Draw();
  
  ln_deltaComp = new TLine(deltaXMin, idealRatio, deltaXMax, idealRatio);
  ln_deltaComp->SetLineStyle(9);
  ln_deltaComp->SetLineWidth(2);
  ln_deltaComp->SetLineColor(38);
  ln_deltaComp->Draw();

  c_tarComp2 = new TCanvas("c_tarComp2", "Comparison of Target Variables 2", canWidth, canHeight); 
  c_tarComp2->Divide(3, 2);
  c_tarComp2->cd(1);
  h_yTarMCWgt->Scale(scaleFactor);
  h_yTarMCWgt->SetOption("HIST");
  h_yTarMCWgt->SetFillColor(8);
  h_yTarMCWgt->SetLineColor(8);
  h_yTarMCWgt->SetFillStyle(3001);
  h_yTarMCWgt->GetXaxis()->SetRangeUser(yTarXMin, yTarXMax);
  //h_yTarMCWgt->GetYaxis()->SetRangeUser(yTarYMin, yTarYMax);
  h_yTarMCWgt->SetTitle("Weighted MC to Data Comparison: Y_{tar}");
  h_yTarMCWgt->Draw();

  h_yTarData->SetMarkerColor(4);
  h_yTarData->SetMarkerStyle(22);
  h_yTarData->Draw("E, SAME");
  if(doSub==true) {
    h_yTarAlData->SetMarkerColor(5);
    h_yTarAlData->SetMarkerStyle(22);
    h_yTarAlData->Draw("E, SAME");
  }

  l_yTarComp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_yTarComp->AddEntry(h_yTarMCWgt, "Weighted MC", "F");
  l_yTarComp->AddEntry(h_yTarData, "Data", "EP");
  if(doSub==true) {
    l_yTarComp->AddEntry(h_yTarAlData, "Data", "EP");
  }
  l_yTarComp->Draw();

  c_tarComp2->cd(2);
  h_xFocalMCWgt->Scale(scaleFactor);
  h_xFocalMCWgt->SetOption("HIST");
  h_xFocalMCWgt->SetFillColor(8);
  h_xFocalMCWgt->SetLineColor(8);
  h_xFocalMCWgt->SetFillStyle(3001);
  h_xFocalMCWgt->GetXaxis()->SetRangeUser(xFocalXMin, xFocalXMax);
  //h_xFocalMCWgt->GetYaxis()->SetRangeUser(xFocalYMin, xFocalYMax);
  h_xFocalMCWgt->SetTitle("Weighted MC to Data Comparison: X_{Foc}");
  h_xFocalMCWgt->Draw();

  h_xFocalData->SetMarkerColor(4);
  h_xFocalData->SetMarkerStyle(22);
  h_xFocalData->Draw("E, SAME");
  if(doSub==true) {
    h_xFocalAlData->SetMarkerColor(5);
    h_xFocalAlData->SetMarkerStyle(22);
    h_xFocalAlData->Draw("E, SAME");
  }

  l_xFocalComp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_xFocalComp->AddEntry(h_xFocalMCWgt, "Weighted MC", "F");
  l_xFocalComp->AddEntry(h_xFocalData, "Data", "EP");
  if(doSub==true) {
    l_xFocalComp->AddEntry(h_xFocalAlData, "Data", "EP");
  }
  l_xFocalComp->Draw();

  c_tarComp2->cd(3);
  h_yFocalMCWgt->Scale(scaleFactor);
  h_yFocalMCWgt->SetOption("HIST");
  h_yFocalMCWgt->SetFillColor(8);
  h_yFocalMCWgt->SetLineColor(8);
  h_yFocalMCWgt->SetFillStyle(3001);
  h_yFocalMCWgt->GetXaxis()->SetRangeUser(yFocalXMin, yFocalXMax);
  //h_yFocalMCWgt->GetYaxis()->SetRangeUser(yFocalYMin, yFocalYMax);
  h_yFocalMCWgt->SetTitle("Weighted MC to Data Comparison: Y_{Foc}");
  h_yFocalMCWgt->Draw();

  h_yFocalData->SetMarkerColor(4);
  h_yFocalData->SetMarkerStyle(22);
  h_yFocalData->Draw("E, SAME");
  if(doSub==true) {
    h_yFocalAlData->SetMarkerColor(5);
    h_yFocalAlData->SetMarkerStyle(22);
    h_yFocalAlData->Draw("E, SAME");
  }

  l_yFocalComp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_yFocalComp->AddEntry(h_yFocalMCWgt, "Weighted MC", "F");
  l_yFocalComp->AddEntry(h_yFocalData, "Data", "EP");
  if(doSub==true) {
    l_yFocalComp->AddEntry(h_yFocalAlData, "Data", "EP");
  }
  l_yFocalComp->Draw();

  c_tarComp2->cd(4);
  h_yTarRatio->Divide(h_yTarMCWgt);
  h_yTarRatio->SetMarkerStyle(22);
  h_yTarRatio->SetMarkerSize(1.25);
  h_yTarRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_yTarRatio->GetXaxis()->SetRangeUser(yTarXMin, yTarXMax);
  h_yTarRatio->GetYaxis()->SetTitle("Y_{tar}^{data} / Y_{tar}^{MC}");
  h_yTarRatio->SetTitle("Ratio of Y_{tar}^{data} to Y_{tar}^{MC}");
  h_yTarRatio->Draw();

  ln_yTarComp = new TLine(yTarXMin, idealRatio, yTarXMax, idealRatio);
  ln_yTarComp->SetLineStyle(9);
  ln_yTarComp->SetLineWidth(2);
  ln_yTarComp->SetLineColor(38);
  ln_yTarComp->Draw();

  c_tarComp2->cd(5);
  h_xFocalRatio->Divide(h_xFocalMCWgt);
  h_xFocalRatio->SetMarkerStyle(22);
  h_xFocalRatio->SetMarkerSize(1.25);
  h_xFocalRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_xFocalRatio->GetXaxis()->SetRangeUser(xFocalXMin, xFocalXMax);
  h_xFocalRatio->GetYaxis()->SetTitle("Y_{foc}^{data} / Y_{foc}^{MC}");
  h_xFocalRatio->SetTitle("Ratio of Y_{foc}^{data} to Y_{foc}^{MC}");
  h_xFocalRatio->Draw();

  ln_xFocalComp = new TLine(xFocalXMin, idealRatio, xFocalXMax, idealRatio);
  ln_xFocalComp->SetLineStyle(9);
  ln_xFocalComp->SetLineWidth(2);
  ln_xFocalComp->SetLineColor(38);
  ln_xFocalComp->Draw();

  c_tarComp2->cd(6);
  h_yFocalRatio->Divide(h_yFocalMCWgt);
  h_yFocalRatio->SetMarkerStyle(22);
  h_yFocalRatio->SetMarkerSize(1.25);
  h_yFocalRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_yFocalRatio->GetXaxis()->SetRangeUser(yFocalXMin, yFocalXMax);
  h_yFocalRatio->GetYaxis()->SetTitle("x_{foc}^{data} / x_{foc}^{MC}");
  h_yFocalRatio->SetTitle("Ratio of x_{foc}^{data} to x_{foc}^{MC}");
  h_yFocalRatio->Draw();
  
  ln_yFocalComp = new TLine(yFocalXMin, idealRatio, yFocalXMax, idealRatio);
  ln_yFocalComp->SetLineStyle(9);
  ln_yFocalComp->SetLineWidth(2);
  ln_yFocalComp->SetLineColor(38);
  ln_yFocalComp->Draw();
  //
  //
  c_kinComp = new TCanvas("c_kinComp", "Comparison of Kinematic Variables", canWidth, canHeight); 
  c_kinComp->Divide(3, 2);
  c_kinComp->cd(1);
  /*
  h_thetaMCWgt->Scale(scaleFactor);
  h_thetaMCWgt->SetOption("HIST");
  h_thetaMCWgt->SetFillColor(8);
  h_thetaMCWgt->SetLineColor(8);
  h_thetaMCWgt->SetFillStyle(3001);
  h_thetaMCWgt->GetXaxis()->SetRangeUser(thetaXMin, thetaXMax);
  //h_thetaMCWgt->GetYaxis()->SetRangeUser(thetaYMin, thetaYMax);
  h_thetaMCWgt->SetTitle("Weighted MC to Data Comparison: #theta");
  h_thetaMCWgt->Draw();

  h_thetaData->SetMarkerColor(4);
  h_thetaData->SetMarkerStyle(22);
  h_thetaData->Draw("E, SAME");

  l_thetaComp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_thetaComp->AddEntry(h_thetaMCWgt, "Weighted MC", "F");
  l_thetaComp->AddEntry(h_thetaData, "Data", "EP");
  l_thetaComp->Draw();
  */

  h_xbjMCWgt->Scale(scaleFactor);
  h_xbjMCWgt->SetOption("HIST");
  h_xbjMCWgt->SetFillColor(8);
  h_xbjMCWgt->SetLineColor(8);
  h_xbjMCWgt->SetFillStyle(3001);
  h_xbjMCWgt->GetXaxis()->SetRangeUser(xbjXMin, xbjXMax);
  //h_xbjMCWgt->GetYaxis()->SetRangeUser(xbjYMin, xbjYMax);
  h_xbjMCWgt->SetTitle("Weighted MC to Data Comparison: #X_{bj}");
  h_xbjMCWgt->Draw();

  h_xbjData->SetMarkerColor(4);
  h_xbjData->SetMarkerStyle(22);
  h_xbjData->Draw("E, SAME");
  if(doSub==true) {
    h_xbjAlData->SetMarkerColor(5);
    h_xbjAlData->SetMarkerStyle(22);
    h_xbjAlData->Draw("E, SAME");
  }

  l_xbjComp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_xbjComp->AddEntry(h_xbjMCWgt, "Weighted MC", "F");
  l_xbjComp->AddEntry(h_xbjData, "Data", "EP");
  if(doSub==true) {
    l_xbjComp->AddEntry(h_xbjAlData, "Data", "EP");
  }
  l_xbjComp->Draw();


  c_kinComp->cd(2);
  h_q2MCWgt->Scale(scaleFactor);
  h_q2MCWgt->SetOption("HIST");
  h_q2MCWgt->SetFillColor(8);
  h_q2MCWgt->SetLineColor(8);
  h_q2MCWgt->SetFillStyle(3001);
  h_q2MCWgt->GetXaxis()->SetRangeUser(q2XMin, q2XMax);
  //h_q2MCWgt->GetYaxis()->SetRangeUser(q2YMin, q2YMax);
  h_q2MCWgt->SetTitle("Weighted MC to Data Comparison: Q^{2}");
  h_q2MCWgt->Draw();

  h_q2Data->SetMarkerColor(4);
  h_q2Data->SetMarkerStyle(22);
  h_q2Data->Draw("E, SAME");
  if(doSub==true) {
    h_q2AlData->SetMarkerColor(5);
    h_q2AlData->SetMarkerStyle(22);
    h_q2AlData->Draw("E, SAME");
  }

  l_q2Comp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_q2Comp->AddEntry(h_q2MCWgt, "Weighted MC", "F");
  l_q2Comp->AddEntry(h_q2Data, "Data", "EP");
  if(doSub==true) {
    l_q2Comp->AddEntry(h_q2AlData, "Data", "EP");
  }
  l_q2Comp->Draw();

  c_kinComp->cd(3);
  h_w2MCWgt->Scale(scaleFactor);
  h_w2MCWgt->SetOption("HIST");
  h_w2MCWgt->SetFillColor(8);
  h_w2MCWgt->SetLineColor(8);
  h_w2MCWgt->SetFillStyle(3001);
  h_w2MCWgt->GetXaxis()->SetRangeUser(w2XMin, w2XMax);
  //h_w2MCWgt->GetYaxis()->SetRangeUser(w2YMin, w2YMax);
  h_w2MCWgt->SetTitle("Weighted MC to Data Comparison: W^{2}");
  h_w2MCWgt->Draw();

  h_w2Data->SetMarkerColor(4);
  h_w2Data->SetMarkerStyle(22);
  h_w2Data->Draw("E, SAME");
  if(doSub==true) {
    h_w2AlData->SetMarkerColor(5);
    h_w2AlData->SetMarkerStyle(22);
    h_w2AlData->Draw("E, SAME");
  }

  l_w2Comp = new TLegend(0.1, 0.7, 0.25, 0.9);
  l_w2Comp->AddEntry(h_w2MCWgt, "Weighted MC", "F");
  l_w2Comp->AddEntry(h_w2Data, "Data", "EP");
  if(doSub==true) {
    l_w2Comp->AddEntry(h_w2AlData, "Data", "EP");
  }
  l_w2Comp->Draw();

  c_kinComp->cd(4);
  h_xbjRatio->Divide(h_xbjMCWgt);
  h_xbjRatio->SetMarkerStyle(22);
  h_xbjRatio->SetMarkerSize(1.25);
  h_xbjRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_xbjRatio->GetXaxis()->SetRangeUser(xbjXMin, xbjXMax);
  h_xbjRatio->GetYaxis()->SetTitle("#X_{data} / #X_{MC}");
  h_xbjRatio->SetTitle("Ratio of #X_{data} to #X_{MC}");
  h_xbjRatio->Draw();

  ln_xbjComp = new TLine(xbjXMin, idealRatio, xbjXMax, idealRatio);
  ln_xbjComp->SetLineStyle(9);
  ln_xbjComp->SetLineWidth(2);
  ln_xbjComp->SetLineColor(38);
  ln_xbjComp->Draw();

  c_kinComp->cd(5);
  h_q2Ratio->Divide(h_q2MCWgt);
  h_q2Ratio->SetMarkerStyle(22);
  h_q2Ratio->SetMarkerSize(1.25);
  h_q2Ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_q2Ratio->GetXaxis()->SetRangeUser(q2XMin, q2XMax);
  h_q2Ratio->GetYaxis()->SetTitle("Q^{2}_{data} / Q^{2}_{MC}");
  h_q2Ratio->SetTitle("Ratio of Q^{2}_{data} to Q^{2}_{MC}");
  h_q2Ratio->Draw();

  ln_q2Comp = new TLine(q2XMin, idealRatio, q2XMax, idealRatio);
  ln_q2Comp->SetLineStyle(9);
  ln_q2Comp->SetLineWidth(2);
  ln_q2Comp->SetLineColor(38);
  ln_q2Comp->Draw();

  c_kinComp->cd(6);
  h_w2Ratio->Divide(h_w2MCWgt);
  h_w2Ratio->SetMarkerStyle(22);
  h_w2Ratio->SetMarkerSize(1.25);
  h_w2Ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  h_w2Ratio->GetXaxis()->SetRangeUser(w2XMin, w2XMax);
  h_w2Ratio->GetYaxis()->SetTitle("W^{2}_{data} / W^{2}_{MC}");
  h_w2Ratio->SetTitle("Ratio of W^{2}_{data} to W^{2}_{MC}");
  h_w2Ratio->Draw();

  ln_w2Comp = new TLine(w2XMin, idealRatio, w2XMax, idealRatio);
  ln_w2Comp->SetLineStyle(9);
  ln_w2Comp->SetLineWidth(2);
  ln_w2Comp->SetLineColor(38);
  ln_w2Comp->Draw();

  c_focalComp = new TCanvas("c_focalComp", "Comparison of Focal Plane Quantites", canWidth, canHeight); 
  c_focalComp->Divide(4, 3);
  c_focalComp->cd(1); gPad->SetLogz();
  h2_xVxpFocalMCWgt->Draw("COLZ");

  c_focalComp->cd(2); gPad->SetLogz();
  h2_xVxpFocalData->Draw("COLZ");

  c_focalComp->cd(3); gPad->SetLogz();
  h2_xVyFocalMCWgt->Draw("COLZ");

  c_focalComp->cd(4); gPad->SetLogz();
  h2_xVyFocalData->Draw("COLZ");

  c_focalComp->cd(5); gPad->SetLogz();
  h2_xVypFocalMCWgt->Draw("COLZ");

  c_focalComp->cd(6); gPad->SetLogz();
  h2_xVypFocalData->Draw("COLZ");

  c_focalComp->cd(7); gPad->SetLogz();
  h2_xpVyFocalMCWgt->Draw("COLZ");

  c_focalComp->cd(8); gPad->SetLogz();
  h2_xpVyFocalData->Draw("COLZ");

  c_focalComp->cd(9); gPad->SetLogz();
  h2_xpVypFocalMCWgt->Draw("COLZ");

  c_focalComp->cd(10); gPad->SetLogz();
  h2_xpVypFocalData->Draw("COLZ");

  c_focalComp->cd(11); gPad->SetLogz();
  h2_yVypFocalMCWgt->Draw("COLZ");

  c_focalComp->cd(12); gPad->SetLogz();
  h2_yVypFocalData->Draw("COLZ");

  
  //inFile->Write();
  inFile->cd();
  c_tarComp->Write();
  c_kinComp->Write();
  c_focalComp->Write();
  c_tarComp2->Write();
  /*
  double w2, eprime, xbj, bornBin;
  Double_t dBin, dBinContent, dBinError;
  //Loop Over Range in dp:
  int nDeltaBins;
  nDeltaBins = h_deltaRatio->GetNbinsX();
  TVectorD bornCSErrArr(1000);
  TVectorD bornCSArr(1000);
  TVectorD xbjErrArr(1000);
  TVectorD xbjArr(1000);
  cout << "dBin dBinError xsec xsecErr eprime xbj\n";
  for(int i = 0; i < nDeltaBins; i++) 
  {
    dBin = h_deltaRatio->GetXaxis()->GetBinCenter(i+1);
    dBinContent = h_deltaRatio->GetBinContent(i+1);
    dBinError = h_deltaRatio->GetBinError(i+1);
    eprime = cal_ep(hsec, dBin);
    w2     = cal_w2(protonMass, beamEnergy, thetac, eprime);
    xbj    = cal_xbj2(beamEnergy, eprime, thetac, protonMass);
    bornBin = born->Interpolate(eprime, thetac);
    if(dBinContent != 0 ){
    bornCSArr[i] = dBinContent * bornBin;
    bornCSErrArr[i] = dBinError * bornBin;
    xbjArr[i] = xbj;
    }
    //if(bornBin > 0.0) {
      cout << dBin << " " << dBinError << " " << dBinContent * bornBin << " " << dBinError * bornBin << " " << eprime << " " << xbj << " " << endl;
      //}
    //Grab the error from the delta histo
    //Multiply error of delta histo by born CS
    //Save point-by-point data somehow for CS...
  }
  TGraphErrors *gr = new TGraphErrors(xbjArr, bornCSArr, xbjErrArr, bornCSErrArr);
  TCanvas *c1 = new TCanvas();
  c1->cd();
  gr->Draw();
    //Save results
    */
}
