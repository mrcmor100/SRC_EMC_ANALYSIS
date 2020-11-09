using namespace ROOT::RDF;
auto snapshotOptions = RSnapshotOptions("UPDATE",ROOT::kZLIB,1,0,99,false);
//const char* kinFileDir = "/volatile/hallc/xem2/cmorean/runlist/kinSettings.dat";
//const char* format2018Path = "/cache/hallc/E12-10-008/cmorean/realpass-3b-data/shms_replay_production_%d_-1.root";
//const char* pathTo2018Replays = "/cache/hallc/E12-10-002/abishek/realpass-3b-shms-data/shms_replay_production_%d_-1.root";
//const char* pathTo2019Replays = "/cache/hallc/E12-10-008/cmorean/realpass-3b-data/shms_replay_production_%d_-1.root";
//const char* dataRunlistDir = "/volatile/hallc/xem2/cmorean/runlist/kin%d_%s_runs.dat";
//const char* runEfficienciesFile = "/volatile/hallc/xem2/cmorean/runlist/emcEfficiencies2.dat";
//const char* mcComparisonFile = "kin%d_m%dc_%s.root";

const char* kinFileDir          = "./runList/kinSettingsAll.dat";
const char* dataRunlistDir      = "./runList/kin%d_%s_runs.dat";
const char* runEfficienciesFile = "./runList/xgt1_2018Efficiencies.dat";
const char* mcComparisonFile    = "./dataToMC/kin%d_m%dc_%s.root";

const char* pathTo2018Replays = "./realpass-3c-shms-data/shms_replay_production_%d_-1.root";
const char* pathTo2019Replays = "/cache/hallc/E12-10-008/cmorean/realpass-3b-data/shms_replay_production_%d_-1.root";

TFile *compFile;
TChain *dataChain;
TDirectory *dataDir;

TFile *ngcerEffFile = TFile::Open("./NGCER_Efficiencies.root"); //Fix rel. dir.
TH2F *h2_ngcerEff = dynamic_cast <TH2F*> (ngcerEffFile->FindObjectAny("histEff"));

int aluminumSubtraction(TFile *compFile, int kin, string targ, int model);

Double_t calc_CalEff(Double_t ePrime);

Double_t getNorm(string fileList, string effTblFile, string targ);

std::fstream& GetTblLine(std::fstream& file, unsigned int num);

TChain *formChain(string list, string format);



int implicitCryo(int kin , string targ, int model, bool snapshot=false) {

  TH1DModel m_xFocalData ("h_xFocalData" , "Data: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm"  , 100,  -40  ,  40  );
  TH1DModel m_xpFocalData("h_xpFocalData", "Data: X'_{fp}; X'_{fp}; Number of Entries / 2 mrad"   , 100, -100.0, 100.0);
  TH1DModel m_yFocalData ("h_yFocalData" , "Data: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm"  , 100,  -40  ,  40  );
  TH1DModel m_ypFocalData("h_ypFocalData", "Data: Y'_{fp}; Y'_{fp}; Number of Entries / 2 mrad"   , 100, -100.0, 100.0);
  TH1DModel m_yTarData   ("h_yTarData"   , "Data: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm", 334,  -10  ,  10  );
  TH1DModel m_xpTarData  ("h_xpTarData"  , "Data: X'_{tar}; X'_{tar}; Number of Entries / 2 mrad" , 100, -100.0, 100.0);
  TH1DModel m_ypTarData  ("h_ypTarData"  , "Data: Y'_{tar}; Y'_{tar}; Number of Entries / 2 mrad" , 100, -100.0, 100.0);
  TH1DModel m_deltaData  ("h_deltaData"  , "Data: #delta; #delta; Number of Entries"              ,  60,  -30  ,  30  );
  TH1DModel m_xbjData    ("h_xbjData"    , "Data: x_{bj}; x_{bj}; Number of Entries "             , 120,    0.0,   3.0);
  TH1DModel m_q2Data     ("h_q2Data"     , "Data: Q^{2}; Q^{2} (GeV^{2}); Number of Entries / 0.025 GeV^{2}", 240,   0.0,  6.0);
  TH1DModel m_w2Data     ("h_w2Data"     , "Data: W^{2}; W^{2} (GeV^{2}); Number of Entries / 0.050 GeV^{2}", 375, -10.0, 20.0);
  TH1DModel m_thetaData  ("h_thetaData"  , "Data: #theta; #theta; Number of Entries / 0.01 deg"             , 200,   8.0, 18.0);

  TH2DModel m_xVxpFocalData ("h2_xVxpFocalData" ,"Data: X_{fp} vs. X'_{fp}; X'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm"   ,100,-100.0,100.0,160, -40  , 40  );
  TH2DModel m_xVyFocalData  ("h2_xVyFocalData"  ,"Data: X_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X_{fp} (cm) / 5 mm"  ,160, -40  , 40  ,160, -40  , 40  );
  TH2DModel m_xVypFocalData ("h2_xVypFocalData" ,"Data: X_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm"   , 60, -60.0, 60.0,160, -40  , 40  );
  TH2DModel m_xpVyFocalData ("h2_xpVyFocalData" ,"Data: X'_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X'_{fp} / 2 mrad"   ,160, -40  , 40  ,100,-100.0,100.0);
  TH2DModel m_xpVypFocalData("h2_xpVypFocalData","Data: X'_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X'_{fp} / 2 mrad"    , 60, -60.0, 60.0,100,-100.0,100.0);
  TH2DModel m_yVypFocalData ("h2_yVypFocalData" ,"Data: Y_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; Y_{fp} (cm) / 5 mm"   , 60, -60.0, 60.0,160, -40  , 40  );
  TH2DModel m_yVxpTarData   ("h2_yVxpTarData"   ,"Data: Y_{tar} vs. X'_{tar}; X'_{tar} / 2 mrad; Y_{tar} / 1 mm"    ,200,-100.0,100.0,100,  -5  ,  5  ); 
  TH2DModel m_yVypTarData   ("h2_yVypTarData"   ,"Data: Y_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; Y_{tar} / 1 mm"    ,200,-100.0,100.0,100,  -5  ,  5  ); 
  TH2DModel m_xpVypTarData  ("h2_xpVypTarData"  ,"Data: X'_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; X'_{tar} / 2 mrad", 60, -60.0, 60.0,100,-100.0,100.0);

  TH1::SetDefaultSumw2(false);
  //============================================
  //   Input Kinematics Settings for MC
  //============================================
  string kinFile = kinFileDir;
  float Ebeam, hsec, thetac;
  fstream kinTblFile(kinFile);
  GetTblLine(kinTblFile, kin+1);  //Changes line pointed to in file.
  kinTblFile >> Ebeam >> hsec >> thetac;
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "                 Beam Info                      \n";
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "Sourcing the kinematics for kin: " << kin << endl
       << kinFile << endl;
  cout << " Beam Energy:         " << Ebeam << 
        "\n Central Momentum:    " << hsec << 
        "\n Central Angle (rad): " << thetac << endl;
  cout << endl << endl;

  //Convert the central angel to radians
  float thetacrad = thetac * TMath::DegToRad();

  string dataRunlistFile = Form(dataRunlistDir, kin, targ.c_str());
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "                 Run List\n";
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "Sourcing RunList from:\n" << dataRunlistFile << endl;
  if(kin==1) {
    dataChain = formChain(dataRunlistFile, pathTo2019Replays);
  }
  else if(kin >= 0) {
    dataChain = formChain(dataRunlistFile, pathTo2018Replays);
  }
  else {
    dataChain = formChain(dataRunlistFile, pathTo2019Replays);
  }
  cout << endl << endl;

  Double_t dataNorm= getNorm(dataRunlistFile, runEfficienciesFile, targ);

  int nEntries;
  nEntries = dataChain->GetEntries();
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "                  Entries \n";
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "The total number of Entries is: " << nEntries << endl;
  cout << endl << endl;

  auto calcCalEff = [=](double deltaData) {
    return 0.9984 - TMath::Exp(-1.98 - 0.2356 * TMath::Power(hsec * (1 + deltaData / 100), 3));
  };
  auto calcCerEff = [&](double yFocalData, double ypFocalData, double xFocalData, double xpFocalData) {
    double xCer = xFocalData - 89.1*xpFocalData;
    double yCer = yFocalData - 89.1*ypFocalData;
    int xBin    = h2_ngcerEff->GetYaxis()->FindBin(xCer);
    int yBin    = h2_ngcerEff->GetXaxis()->FindBin(yCer);
    double cerEff  = h2_ngcerEff->GetBinContent(yBin,xBin);
    return cerEff==0 ? 1 : cerEff;
  };
  auto rad2mrad   = [](double pFocalData) {
    return pFocalData*1000.0;
  };
  auto getTheta   = [=] (double ypTarData) {
    return  (thetac*TMath::DegToRad() + ypTarData)*TMath::RadToDeg();
  };
  auto weightAll  = [=] (double calEff, double cerEff) {
    return (1.0 / calEff / cerEff / dataNorm);
  };

  ROOT::EnableImplicitMT(8);
  //cout << ROOT::GetImplicitMTPoolSize() << endl;
  ROOT::RDataFrame d(*dataChain);
  //Filter all the crap events BEFORE making histograms and new variables.

  auto d2 = d.Filter("P.ngcer.npeSum > 2.0")
    .Filter("P.cal.etottracknorm > 0.7")
    .Filter("P.gtr.dp > -10. && P.gtr.dp < 22.")
    .Filter("P.bcm.bcm4c.AvgCurrent > 5")
    .Filter("abs(P.gtr.ph) < 0.1")
    .Filter("abs(P.gtr.th) < 0.1")
    .Filter("abs(P.gtr.y) < 6.0");
    //.Filter("P.dc.InsideDipoleExit==1")
    //.Filter("T.shms.pEDTM_tdcTimeRaw==0.0");  //Remove because not valuable
 //.Filter("P.gtr.beta > 0.5 && P.gtr.beta < 1.5")
 //.Filter("P.kin.W2 > 6.0 && P.kin.W2 < 13.")
 //.Filter("abs(P.dc.x_fp) <= 40.")
 //.Filter("abs(P.dc.y_fp) <= 40.")
  //Define new variables of interest.
  d2 = d2.Define("P_gtr_xp_mrad" ,rad2mrad, {"P.gtr.th"}  )
         .Define("P_gtr_yp_mrad" ,rad2mrad, {"P.gtr.ph"}  )
         .Define("P_dc_xpfp_mrad",rad2mrad, {"P.dc.xp_fp"})
         .Define("P_dc_ypfp_mrad",rad2mrad, {"P.dc.yp_fp"})
         .Define("P_theta"       ,getTheta, {"P.gtr.ph"}  )
         .Define("caloEff",calcCalEff,{"P.gtr.dp"})
         .Define("cerEff",calcCerEff,{"P.dc.y_fp","P.gtr.ph","P.dc.x_fp","P.gtr.th"})
         .Define("totalWeight", weightAll, {"caloEff","cerEff"});

  auto nEvents = d2.Count();
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "             Good Electron Events\n";
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "Number of good electron events: " << *nEvents << 
    endl << endl << endl;

  if(snapshot) {
    d2.Snapshot("DataTree",Form(mcComparisonFile,kin, model, targ.c_str()),
		{"caloEff", "cerEff", "P.gtr.dp", "P.gtr.ph","P.gtr.th", "P.gtr.y", "P.kin.W2"}
		,snapshotOptions);
  }

  //Make all the relevant histograms
  auto h_xFocalData  = d2.Histo1D(m_xFocalData  ,"P.dc.x_fp","totalWeight");
  auto h_xpFocalData = d2.Histo1D(m_xpFocalData ,"P_dc_xpfp_mrad","totalWeight");
  auto h_yFocalData  = d2.Histo1D(m_yFocalData  ,"P.dc.y_fp","totalWeight");
  auto h_ypFocalData = d2.Histo1D(m_ypFocalData ,"P_dc_ypfp_mrad","totalWeight");
  auto h_yTarData    = d2.Histo1D(m_yTarData , "P.gtr.y","totalWeight");
  auto h_xpTarData   = d2.Histo1D(m_xpTarData, "P_gtr_xp_mrad","totalWeight");
  auto h_ypTarData   = d2.Histo1D(m_ypTarData, "P_gtr_yp_mrad","totalWeight");
  auto h_deltaData   = d2.Histo1D(m_deltaData, "P.gtr.dp","totalWeight");
  auto h_xbjData     = d2.Histo1D(m_xbjData  , "P.kin.x_bj","totalWeight");
  auto h_q2Data      = d2.Histo1D(m_q2Data   , "P.kin.Q2","totalWeight");
  auto h_w2Data      = d2.Histo1D(m_w2Data   , "P.kin.W2","totalWeight");
  auto h_thetaData   = d2.Histo1D(m_thetaData, "P_theta","totalWeight");

  // Data 2D histos
  auto h2_xVxpFocalData  = d2.Histo2D(m_xVxpFocalData , "P_dc_xpfp_mrad","P.dc.x_fp"      ,"totalWeight");
  auto h2_xVyFocalData   = d2.Histo2D(m_xVyFocalData  , "P.dc.y_fp","P.dc.x_fp"           ,"totalWeight");
  auto h2_xVypFocalData  = d2.Histo2D(m_xVypFocalData , "P_dc_ypfp_mrad","P.dc.x_fp"      ,"totalWeight");
  auto h2_xpVyFocalData  = d2.Histo2D(m_xpVyFocalData , "P.dc.y_fp","P_dc_xpfp_mrad"      ,"totalWeight");
  auto h2_xpVypFocalData = d2.Histo2D(m_xpVypFocalData, "P_dc_ypfp_mrad","P_dc_xpfp_mrad" ,"totalWeight");
  auto h2_yVypFocalData  = d2.Histo2D(m_yVypFocalData , "P_dc_ypfp_mrad","P.dc.y_fp"      ,"totalWeight");
  auto h2_yVxpTarData    = d2.Histo2D(m_yVxpTarData   , "P_gtr_xp_mrad","P.gtr.y"         ,"totalWeight");
  auto h2_yVypTarData    = d2.Histo2D(m_yVypTarData   , "P_gtr_yp_mrad","P.gtr.y"         ,"totalWeight");
  auto h2_xpVypTarData   = d2.Histo2D(m_xpVypTarData  , "P_gtr_yp_mrad","P_gtr_xp_mrad"   ,"totalWeight");


  compFile  = new TFile(Form(mcComparisonFile,kin, model, targ.c_str()), "UPDATE");
  if (!compFile || compFile->IsZombie()) { cout <<"The file could not be opened!"; delete compFile; } 
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "                     Output \n";
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "Output File is:\n" << 
    Form(mcComparisonFile,kin, model, targ.c_str()) << 
    endl << endl << endl;

  // Create data directory and descend into it
  dataDir = dynamic_cast <TDirectory*> (compFile->Get("dataDir"));
  if(!dataDir) {dataDir = compFile->mkdir("dataDir"); dataDir->cd();}
  else if(dataDir) {compFile->rmdir("dataDir"); dataDir = compFile->mkdir("dataDir"); dataDir->cd();}

  h_xFocalData->Write();
  h_xpFocalData->Write();
  h_yFocalData->Write();
  h_ypFocalData->Write();
  h_yTarData->Write();
  h_xpTarData->Write();
  h_ypTarData->Write();
  h_deltaData->Write();
  h_xbjData->Write();
  h_q2Data->Write();
  h_w2Data->Write();
  h_thetaData->Write();
  h2_xVxpFocalData->Write();
  h2_xVyFocalData->Write();
  h2_xVypFocalData->Write();
  h2_xpVyFocalData->Write();
  h2_xpVypFocalData->Write();
  h2_yVypFocalData->Write();
  h2_yVxpTarData->Write();
  h2_yVypTarData->Write();
  h2_xpVypTarData->Write();

  aluminumSubtraction(compFile, kin, targ, model);

  compFile->Close();

  return 0;
}

TChain *formChain(string fileList, string format) {

  TChain *chain = new TChain("T");
  ifstream myfile (fileList);
  int currentRun;
  string line;
  while(getline(myfile, line)) {
    stringstream ss(line);
    ss >> currentRun;
    cout << "Adding run: " << currentRun << " to chain\n";
    chain->Add(Form(format.c_str(), currentRun));
  }

  return chain;
}

std::fstream& GetTblLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

Double_t getNorm(string fileList, string effTblFile, string targ) {
  //Vars from eff tbl
  int ps;
  double chg, current, oCurrent, trkEff, clt, elt;
  Double_t sumCharge = 0;
  //Get initial run from runList
  fstream runList(fileList);
  string runEntry;
  getline(runList, runEntry);
  istringstream ss(runEntry);
  int currentRun;
  ss >> currentRun;

  fstream effTbl(effTblFile);
  string runLine;
  int runNo;

  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "            Charge and Efficiencies  \n";
  //cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";

  double slopeVal = 0;
  double LD2_slope = 3.42/10000.; //Casey's slopes !?!?!?
  double LH2_slope = 5.3314/10000.;

  //double LD2_slope = 0.00080029; //Carlos's Slopes
  //double LH2_slope = 0.00063396;

  double boil_corr;
  if(std::strcmp(targ.c_str(),"H1")==0) {slopeVal = LH2_slope;}
  else if(std::strcmp(targ.c_str(),"D2")==0) {slopeVal = LD2_slope;}
  else {slopeVal = 0.0;}
  //cout << "No Boiling Correction for solid Target: " << targ.c_str() << endl;
  //Assume runlist is in order and complete?
  while (getline(effTbl, runLine)) {
    istringstream ss(runLine);
    ss >> runNo;
    if (runNo == currentRun) {
      ss >> ps >> chg >> current >> oCurrent >> trkEff >> clt >> elt;
      cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
      cout << " Run Number: " << runNo << endl <<
	" Cut Charge from Report Output File: " << chg << endl <<
	" Cut Current from Report Output File: " << current << endl <<
	" Tracking Efficiency: " << trkEff << endl <<
	" Computer Live Time: " << clt << endl << 
	" Electronic Live Time: " << elt << endl << 
	" Pre-Scale: " << ps << endl;
      boil_corr = (1. - slopeVal * current);
      if(std::strcmp(targ.c_str(),"H1")==0) {cout << " Boiling Correction for run " << runNo << " is: " << boil_corr << endl;}
      else if(std::strcmp(targ.c_str(),"D2")==0) {cout << " Boiling Correction for run " << runNo << " is: " << boil_corr << endl;}
      sumCharge += (chg*(clt * boil_corr * trkEff * elt/ ps));
      cout << " SumCharge: " << sumCharge << endl;
      //Next run.
      getline(runList, runEntry);
      istringstream ss(runEntry);
      ss >> currentRun;
    }
  }
  cout << endl << endl;
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "               Data Normalization\n";
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << " The total Normalization is: " << sumCharge << endl << endl;
  return sumCharge;
}

Double_t calc_CalEff(Double_t ePrime) {
  return 0.998394241 - TMath::Exp(-1.98 - 0.2356 * TMath::Power(ePrime, 3));
}

int aluminumSubtraction(TFile *compFile, int kin, string targ, int model) {

  // Aluminum Data histos
  TH1D *h_xFocalAlData, *h_xpFocalAlData, *h_yFocalAlData, *h_ypFocalAlData;
  TH1D *h_yTarAlData, *h_xpTarAlData, *h_ypTarAlData, *h_deltaAlData;
  TH1D *h_thetaAlData, *h_q2AlData, *h_w2AlData, *h_xbjAlData;
  
  //Subtracted Data Data histos
  TH1D *h_xFocalSubData, *h_xpFocalSubData, *h_yFocalSubData, *h_ypFocalSubData;
  TH1D *h_yTarSubData, *h_xpTarSubData, *h_ypTarSubData, *h_deltaSubData;
  TH1D *h_thetaSubData, *h_q2SubData, *h_w2SubData, *h_xbjSubData;

  static const Double_t  dFactor_hydrogen  = 1. / 3.789;
  static const Double_t  dFactor_deuterium = 1. / 4.063;
  Float_t dummyFactor;

  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  cout << "            Aluminum Subtraction\n";
  cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  if(std::strcmp(targ.c_str(),"H1")==0) {dummyFactor = dFactor_hydrogen;
  cout << " Subtracting Aluminum from: " << targ.c_str() << endl << endl;
  }
  else if(std::strcmp(targ.c_str(),"D2")==0) {dummyFactor = dFactor_deuterium;
    cout << " Subtracting Aluminum from: " << targ.c_str() << endl << endl;
  }
  else {
    cout << " Not Subtracting Aluminum from: " << targ.c_str() << endl << endl;
    return 0;
  }

  TFile *alumFile;
  TDirectory *alumDir, *unsubDataDir, *alumDataDir;
  TDirectory *subDataDir;
  const char* alumFilePattern;
  const char* alumFilePath;
  alumFilePattern = "./dataToMC/kin%d_m%dc_%s.root";
  alumFilePath = Form(alumFilePattern, kin, model, "Al");
  if(!gSystem->AccessPathName(alumFilePath)) {
    alumFile = new TFile(alumFilePath, "READ");
    alumFile->cd();
    alumDataDir = dynamic_cast <TDirectory*> (alumFile->Get("dataDir"));
    if(!alumDataDir) {
      cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
      cout << "ERROR! No dataDir in Aluminum File\n\n";
      return 0;
    }
  } else {
    cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
    cout << "ERROR: No Aluminum File found for kin:" <<
      kin << "and model:  " << model << endl << endl;
    return 0;
  }

  // Aluminum Data 1D histos
  h_xFocalAlData  = dynamic_cast <TH1D*> (alumDataDir->Get("h_xFocalData")->Clone("h_xFocalAlData"));
  h_xFocalAlData->SetTitle("AlData: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm");
  h_xpFocalAlData = dynamic_cast <TH1D*> (alumDataDir->Get("h_xpFocalData")->Clone("h_xpFocalAlData"));
  h_xpFocalAlData->SetTitle("AlData: X'_{fp}; X'_{fp}; Number of Entries / 2 mrad");
  h_yFocalAlData  = dynamic_cast <TH1D*> (alumDataDir->Get("h_yFocalData")->Clone("h_yFocalAlData"));
  h_yFocalAlData->SetTitle("AlData: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm");
  h_ypFocalAlData = dynamic_cast <TH1D*> (alumDataDir->Get("h_ypFocalData")->Clone("h_ypFocalAlData"));
  h_ypFocalAlData->SetTitle("AlData: Y'_{fp}; Y'_{fp}; Number of Entries / 2 mrad");
  h_yTarAlData    = dynamic_cast <TH1D*> (alumDataDir->Get("h_yTarData")->Clone("h_yTarAlData"));
  h_yTarAlData->SetTitle("AlData: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm");
  h_xpTarAlData   = dynamic_cast <TH1D*> (alumDataDir->Get("h_xpTarData")->Clone("h_xpTarAlData"));
  h_xpTarAlData->SetTitle("AlData: X'_{tar}; X'_{tar}; Number of Entries / 2 mrad");
  h_ypTarAlData   = dynamic_cast <TH1D*> (alumDataDir->Get("h_ypTarData")->Clone("h_ypTarAlData"));
  h_ypTarAlData->SetTitle("AlData: Y'_{tar}; Y'_{tar}; Number of Entries / 2 mrad");
  h_deltaAlData   = dynamic_cast <TH1D*> (alumDataDir->Get("h_deltaData")->Clone("h_deltaAlData"));
  h_deltaAlData->SetTitle("AlData: #delta; #delta; Number of Entries");
  h_thetaAlData   = dynamic_cast <TH1D*> (alumDataDir->Get("h_thetaData")->Clone("h_thetaAlData"));
  h_thetaAlData->SetTitle("AlData: #theta; #theta; Number of Entries / 0.01 deg");
  h_q2AlData      = dynamic_cast <TH1D*> (alumDataDir->Get("h_q2Data")->Clone("h_q2AlData"));
  h_q2AlData->SetTitle("AlData: Q^{2}; Q^{2} (GeV^{2}); Number of Entries / 0.025 GeV^{2}");
  h_w2AlData      = dynamic_cast <TH1D*> (alumDataDir->Get("h_w2Data")->Clone("h_w2AlData"));
  h_w2AlData->SetTitle("AlData: W^{2}; W^{2} (GeV^{2}); Number of Entries / 0.050 GeV^{2}");
  h_xbjAlData     = dynamic_cast <TH1D*> (alumDataDir->Get("h_xbjData")->Clone("h_xbjAlData"));
  h_xbjAlData->SetTitle("AlData: X_{bj}; X_{bj}; Number of Entries;");

  //Normalize based on dummy factor
  h_xFocalAlData->Scale(dummyFactor);
  h_xpFocalAlData->Scale(dummyFactor);
  h_yFocalAlData->Scale(dummyFactor);
  h_ypFocalAlData->Scale(dummyFactor);
  h_yTarAlData->Scale(dummyFactor);
  h_xpTarAlData->Scale(dummyFactor);
  h_ypTarAlData->Scale(dummyFactor);
  h_deltaAlData->Scale(dummyFactor);
  h_thetaAlData->Scale(dummyFactor);
  h_q2AlData->Scale(dummyFactor);
  h_w2AlData->Scale(dummyFactor);
  h_xbjAlData->Scale(dummyFactor);

  compFile->cd();
  alumDir = dynamic_cast <TDirectory*> (compFile->Get("alumDir"));
  if(!alumDir) {
    alumDir = compFile->mkdir("alumDir");
    alumDir->cd();
  } else if(alumDir) {
    compFile->rmdir("alumDir");
    alumDir = compFile->mkdir("alumDir");
    alumDir->cd();
  }

  //Save to the Aluminum Directory
  h_xFocalAlData->Write();
  h_xpFocalAlData->Write();
  h_yFocalAlData->Write();
  h_ypFocalAlData->Write();
  h_yTarAlData->Write();
  h_xpTarAlData->Write();
  h_ypTarAlData->Write();
  h_deltaAlData->Write();
  h_thetaAlData->Write();
  h_q2AlData->Write();
  h_w2AlData->Write();
  h_xbjAlData->Write();

  compFile->cd("../");
  unsubDataDir = dynamic_cast <TDirectory*> (compFile->Get("dataDir"));
  if(!unsubDataDir) {
    cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
    cout << "ERROR!  No data in file.\n";
    cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
    return 0;
  }
  //Soon-to-be subtracted Data 1D histos
  h_xFocalSubData  = dynamic_cast <TH1D*> (unsubDataDir->Get("h_xFocalData")->Clone("h_xFocalSubData"));
  h_xFocalSubData->SetTitle("subData: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm");
  h_xpFocalSubData = dynamic_cast <TH1D*> (unsubDataDir->Get("h_xpFocalData")->Clone("h_xpFocalSubData"));
  h_xpFocalSubData->SetTitle("subData: X'_{fp}; X'_{fp}; Number of Entries / 2 mrad");
  h_yFocalSubData  = dynamic_cast <TH1D*> (unsubDataDir->Get("h_yFocalData")->Clone("h_yFocalSubData"));
  h_yFocalSubData->SetTitle("subData: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm");
  h_ypFocalSubData = dynamic_cast <TH1D*> (unsubDataDir->Get("h_ypFocalData")->Clone("h_ypFocalSubData"));
  h_ypFocalSubData->SetTitle("subData: Y'_{fp}; Y'_{fp}; Number of Entries / 2 mrad");
  h_yTarSubData    = dynamic_cast <TH1D*> (unsubDataDir->Get("h_yTarData")->Clone("h_yTarSubData"));
  h_yTarSubData->SetTitle("subData: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm");
  h_xpTarSubData   = dynamic_cast <TH1D*> (unsubDataDir->Get("h_xpTarData")->Clone("h_xpTarSubData"));
  h_xpTarSubData->SetTitle("subData: X'_{tar}; X'_{tar}; Number of Entries / 2 mrad");
  h_ypTarSubData   = dynamic_cast <TH1D*> (unsubDataDir->Get("h_ypTarData")->Clone("h_ypTarSubData"));
  h_ypTarSubData->SetTitle("subData: Y'_{tar}; Y'_{tar}; Number of Entries / 2 mrad");
  h_deltaSubData   = dynamic_cast <TH1D*> (unsubDataDir->Get("h_deltaData")->Clone("h_deltaSubData"));
  h_deltaSubData->SetTitle("subData: #delta; #delta; Number of Entries");
  h_thetaSubData   = dynamic_cast <TH1D*> (unsubDataDir->Get("h_thetaData")->Clone("h_thetaSubData"));
  h_thetaSubData->SetTitle("subData: #theta; #theta; Number of Entries / 0.01 deg");
  h_q2SubData      = dynamic_cast <TH1D*> (unsubDataDir->Get("h_q2Data")->Clone("h_q2SubData"));
  h_q2SubData->SetTitle("subData: Q^{2}; Q^{2} (GeV^{2}); Number of Entries / 0.025 GeV^{2}");
  h_w2SubData      = dynamic_cast <TH1D*> (unsubDataDir->Get("h_w2Data")->Clone("h_w2SubData"));
  h_w2SubData->SetTitle("subData: W^{2}; W^{2} (GeV^{2}); Number of Entries / 0.050 GeV^{2}");
  h_xbjSubData     = dynamic_cast <TH1D*> (unsubDataDir->Get("h_xbjData")->Clone("h_xbjSubData"));
  h_xbjSubData->SetTitle("subData: X_{bj}; X_{bj}; Number of Entries;");

  //Check, create and descend into subDataDir
  subDataDir = dynamic_cast <TDirectory*> (compFile->Get("subDataDir"));
  if(!subDataDir) {
    subDataDir = compFile->mkdir("subDataDir");
    subDataDir->cd();}
  else if(subDataDir) {
    compFile->rmdir("subDataDir");
    subDataDir = compFile->mkdir("subDataDir");
    subDataDir->cd();
  }

  //Subtract aluminum from these histograms
  h_xFocalSubData->Add( h_xFocalAlData, -1.);
  h_xpFocalSubData->Add( h_xpFocalAlData, -1.);
  h_yFocalSubData->Add( h_yFocalAlData, -1.);
  h_ypFocalSubData->Add( h_ypFocalAlData, -1.);
  h_yTarSubData->Add( h_yTarAlData, -1.);
  h_xpTarSubData->Add( h_xpTarAlData, -1.);
  h_ypTarSubData->Add( h_ypTarAlData, -1.);
  h_deltaSubData->Add( h_deltaAlData, -1.);
  h_thetaSubData->Add( h_thetaAlData, -1.);
  h_q2SubData->Add( h_q2AlData, -1.);
  h_w2SubData->Add( h_w2AlData, -1.);
  h_xbjSubData->Add( h_xbjAlData, -1.);

  //Save the Aluminum subtracted data
  subDataDir->cd();
  h_xFocalSubData->Write();
  h_xpFocalSubData->Write();
  h_yFocalSubData->Write();
  h_ypFocalSubData->Write();
  h_yTarSubData->Write();
  h_xpTarSubData->Write();
  h_ypTarSubData->Write();
  h_deltaSubData->Write();
  h_thetaSubData->Write();
  h_q2SubData->Write();
  h_w2SubData->Write();
  h_xbjSubData->Write();

  alumFile->Close();

  return 0;
}
