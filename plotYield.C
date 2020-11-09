TFile *comparisonFile;
TDirectory *dataDir, *mcWgtDir;
TH1D  *h_xbjData, *h_xbjMCWgt;

const char* formatOutFilePath = "./csTables/dataCS_%s.dat";
const char* formatOutFilePath2 = "./csTables/modelCS_%s.dat";
const char* kinFileDir    = "./runList/kinSettingsAll.dat";
const char* csTableDir = "./csTables/kin%d_%s_m%d.out";
const char* mcComparisonFile = "./dataToMC/kin%d_m%dc_%s.root";

double modelCS;
double xbjData, xbjDataError, xbjMCWgt, xbjMCWgtError;
double xbjBinCenter, xbjRatioError, ePrime, model, W2;
int nBins;
double Mp = 0.9382723;
int totalBins = 0;
double central_p;

std::fstream& GetTblLine(std::fstream& file, unsigned int num);

TCanvas *yieldCanvas = new TCanvas("yields","yields",600,600);
TLegend *legend = new TLegend(0.75, 0.70, .875, 0.875);

vector<int> markerStyles = {20,21,22,23,33};
int plotindex = 0;
int plotYield(int model) {

  int palette[5] = {kMagenta+2,kRed-3,kBlue-3,kGreen+2,kCyan-6};
  gStyle->SetPalette(5, palette);


  vector<int> allKins = {0};
  vector<string> allStrs = {"C12","D2","B10","B11","Be9"};

  TGraph2D* xt = new TGraph2D();
  TGraph2D* rt = new TGraph2D();
  TGraph2D* cst = new TGraph2D();
  
//============================================
//   Input Kinematics Settings for MC
//============================================
  double hsec, thetac;
  string kinFile = kinFileDir;
  fstream kinTblFile(kinFile);
  //Kin0 has a 2% offset in hsec.
 
  //Determine momentum based on kinematic.. Switch or if statement for now.
  for(string currentStr : allStrs) {

    comparisonFile = new TFile(Form(mcComparisonFile,allKins[0],model,currentStr.c_str()));
    if(std::strcmp(currentStr.c_str(),"H1")==0 || std::strcmp(currentStr.c_str(),"D2")==0) {
      dataDir = dynamic_cast <TDirectory*> (comparisonFile->Get("subDataDir"));
      //mcWgtDir = dynamic_cast <TDirectory*> (comparisonFile->Get("mcWgtDir"));
      h_xbjData = dynamic_cast <TH1D*> (dataDir->Get("h_xbjSubData"));
      //h_xbjMCWgt = dynamic_cast <TH1D*> (mcWgtDir->Get("h_xbjMCWgt"));
    } else {
      dataDir = dynamic_cast <TDirectory*> (comparisonFile->Get("dataDir"));
      //mcWgtDir = dynamic_cast <TDirectory*> (comparisonFile->Get("mcWgtDir"));
      h_xbjData = dynamic_cast <TH1D*> (dataDir->Get("h_xbjData"));
      //h_xbjMCWgt = dynamic_cast <TH1D*> (mcWgtDir->Get("h_xbjMCWgt"));
    }

    yieldCanvas->cd();
    legend->AddEntry(h_xbjData, currentStr.c_str());
    //h_xbjData->SetTitle("Spring 2018 Yield");
    h_xbjData->SetXTitle("X");
    h_xbjData->SetYTitle("Yield");
    h_xbjData->GetXaxis()->CenterTitle();
    h_xbjData->GetYaxis()->CenterTitle();
    h_xbjData->GetYaxis()->SetTitleOffset(1.2);
    //h_xbjData->GetXaxis()->SetTitleOffset(1.1);
    h_xbjData->SetAxisRange(0.4,2.2,"X");
    h_xbjData->SetAxisRange(1.7E-3,0.6E3,"Y");
    h_xbjData->SetStats(0);
    h_xbjData->SetMarkerSize(0.65);
    h_xbjData->SetMarkerStyle(kFullCircle);
    h_xbjData->Draw("SAME E1 P0 X0 PMC PLC");
    yieldCanvas->SetLogy();
    yieldCanvas->SetTicks(1,1);
    yieldCanvas->SetGrid(1,1);
    yieldCanvas->Update();
    plotindex+=1;
  }
  legend->Draw();

  //TLine *hline = new TLine(1.45,1E-6,1.45,1E0);
  //hline->SetLineColor(kOrange+4);
  //hline->SetLineWidth(2);
  //hline->Draw("SAME");

  //TPaveText *p1 = new TPaveText();
  //p1->SetX1(0.5);
  //p1->SetY1(-.7);
  //p1->AddText("Q^{2}=2.05 GeV/c @ X=1");
  //p1->Draw();


  return 0;
}

std::fstream& GetTblLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}
