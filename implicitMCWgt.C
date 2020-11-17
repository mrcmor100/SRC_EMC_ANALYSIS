using namespace ROOT::RDF;
auto snapshotOptions = RSnapshotOptions("UPDATE",ROOT::kZLIB,1,0,99,false);
//const char* targetInfoDir = "/group/c-xem2/cmorean/E12-10-008/target_info/targetInfo_new.dat";
//const char* mcRawInputDir = "/group/c-xem2/cmorean/E12-10-008/monteCarlos/infiles/kin%d_%s.inp";
//const char* csTableDir = "/group/c-xem2/cmorean/E12-10-008/csTables/kin%d_%s_m%d.out";
//const char* kinFileDir = "/group/c-xem2/cmorean/E12-10-008/runList/kinSettingsAll.dat";
//const char* mcComparisonFile = "kin%d_m%dc_%s.root";  //For testing purposes only.  Need to unify this and implicitCryo.C
//const char* mcRawForm = "/volatile/hallc/xem2/cmorean/mc-single-arm/worksim/kin%d_%s_*.root";  //All files for a given targ and kin (~50M evts)

const char* targetInfoDir = "./targetInfo/target_info_new.dat";
const char* kinFileDir    = "./runList/kinSettingsAll.dat";
const char* mcRawInputDir = "./monteCarlos/infiles/kin%d_%s.inp";
const char* csTableDir    = "./csTables/kin%d_%s_m%d.out";

const char* mcComparisonFile = "./dataToMC/kin%d_m%dc_%s.root";
const char* mcRawForm = "./monteCarlos/worksim/kin%d_%s.root";

TFile *compFile;
TDirectory *mcWgtDir;

double GetInpVar(std::fstream& file, unsigned int num);
std::fstream& GetTblLine(std::fstream& file, unsigned int num);
double phase_space, sigtot, lummc, lumfract, fract;

int implicitMCWgt(int kin, string targ, int model, bool snapshot=false, bool doCSBCorr = false) {

  TH1DModel m_xFocalMCWgt ("h_xFocalMCWgt","Weighted Monte-Carlo: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm",100, -40, 40);
  TH1DModel m_xpFocalMCWgt("h_xpFocalMCWgt","Weighted Monte-Carlo: X'_{fp}; X'_{fp}; Number of Entries / 2 mrad",100, -100.0, 100.0);
  TH1DModel m_yFocalMCWgt ("h_yFocalMCWgt","Weighted Monte-Carlo: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm",100, -40, 40);
  TH1DModel m_ypFocalMCWgt("h_ypFocalMCWgt","Weighted Monte-Carlo: Y'_{fp}; Y'_{fp}; Number of Entries / 2 mrad",100, -100.0, 100.0);
  TH1DModel m_yTarMCWgt   ("h_yTarMCWgt","Weighted Monte-Carlo: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm",334, -10, 10);
  TH1DModel m_xpTarMCWgt  ("h_xpTarMCWgt","Weighted Monte-Carlo: X'_{tar}; X'_{tar}; Number of Entries / 2 mrad",100, -100.0, 100.0);
  TH1DModel m_ypTarMCWgt  ("h_ypTarMCWgt","Weighted Monte-Carlo: Y'_{tar}; Y'_{tar}; Number of Entries / 2 mrad",100, -100.0, 100.0);
  TH1DModel m_deltaMCWgt  ("h_deltaMCWgt","Weighted Monte-Carlo: #delta; #delta; Number of Entries",60, -30, 30);
  TH1DModel m_thetaMCWgt  ("h_thetaMCWgt","Weighted Monte-Carlo: #theta; #theta; Number of Entries / 0.01 deg",400, 5.0, 45.0);
  TH1DModel m_q2MCWgt     ("h_q2MCWgt","Weighted Monte-Carlo: Q^{2}; Q^{2} (GeV^{2}); Number of Entries / 0.025 GeV^{2}",240, 0.0, 6.0);
  TH1DModel m_w2MCWgt     ("h_w2MCWgt","Weighted Monte-Carlo: W^{2}; W^{2} (GeV^{2}); Number of Entries / 0.050 GeV^{2}",375, -10.0, 20.0);
  TH1DModel m_xbjMCWgt    ("h_xbjMCWgt","Weighted Monte-Carlo: X_{bj}; X_{bj}; Number of Entries / 0.050", 120, 0.0, 3.0);

  TH2DModel m_xVxpFocalMCWgt ("h2_xVxpFocalMCWgt","Weighted Monte-Carlo: X_{fp} vs. X'_{fp}; X'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm",
			       100, -100.0, 100.0, 160, -40, 40);
  TH2DModel m_xVypFocalMCWgt ("h2_xVypFocalMCWgt","Weighted Monte-Carlo: X_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm",
			       60, -60.0, 60.0, 160, -40, 40);
  TH2DModel m_xVyFocalMCWgt  ("h2_xVyFocalMCWgt","Weighted Monte-Carlo: X_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X_{fp} (cm) / 5 mm",
			       160, -40, 40, 160, -40, 40);
  TH2DModel m_xpVyFocalMCWgt ("h2_xpVyFocalMCWgt","Weighted Monte-Carlo: X'_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X'_{fp} / 2 mrad",
			       160, -40, 40, 100, -100.0, 100.0);
  TH2DModel m_xpVypFocalMCWgt("h2_xpVypFocalMCWgt","Weighted Monte-Carlo: X'_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X'_{fp} / 2 mrad",
			       60,-60.0,60.0,100,-100.0, 100.0);
  TH2DModel m_yVypFocalMCWgt ("h2_yVypFocalMCWgt","Weighted Monte-Carlo: Y_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; Y_{fp} (cm) / 5 mm",
			       60,-60.0,60.0,160,-40,40);
  TH2DModel m_yVxpTarMCWgt   ("h2_yVxpTarMCWgt","Weighted Monte-Carlo: Y_{tar} vs. X'_{tar}; X'_{tar} / 2 mrad; Y_{tar} / 1 mm",
			       200,-100.0,100.0,100,-5,5);
  TH2DModel m_yVypTarMCWgt   ("h2_yVypTarMCWgt","Weighted Monte-Carlo: Y_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; Y_{tar} / 1 mm",
			       200,-100.0,100.0,100,-5,5);
  TH2DModel m_xpVypTarMCWgt  ("h2_xpVypTarMCWgt","Weighted Monte-Carlo: X'_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; X'_{tar} / 2 mrad",
			       60,-60.0,60.0,100,-100.0, 100.0);

//============================================
//    Defined offsets?? Probably not used
//============================================
  double ypcor = 0.0;

//============================================
//      Input file directories and names
//============================================
  string targetInfoFile  = targetInfoDir;
  string kinFile         = kinFileDir;
  string csTableFile     = Form(csTableDir,kin, targ.c_str(), model);
  TString mcRawInputFile = Form(mcRawInputDir, kin, targ.c_str());
  string mcRawFiles      = Form(mcRawForm, kin, targ.c_str());

//============================================
//        Initialize Physical Constants
//============================================
  float  mp     = 0.9382723;
  float  mp2    = mp*mp;
  double sigave = 0.0;

//============================================
//     input from input file to Raw MC
//============================================
  fstream mcInpFile(mcRawInputFile);
  float deldown = GetInpVar(mcInpFile, 8);
  float delup   = GetInpVar(mcInpFile, 9);
  float dxpup   = GetInpVar(mcInpFile, 11);
  float dxpdown = GetInpVar(mcInpFile, 10);
  float dypup   = GetInpVar(mcInpFile, 13);
  float dypdown = GetInpVar(mcInpFile, 12);
  cout << "=============================================\n";
  cout << "            Phase Space Limits in MC         \n";
  cout << "=============================================\n";
  cout << " Delta: deldown = " << deldown << endl;
  cout << " Delta: delup   = " << delup << endl;
  cout << " xpTar: dxpdown = " << dxpdown << endl;
  cout << " xpTar: dypup   = " << dxpup << endl;
  cout << " ypTar: dypdown = " << dypdown << endl;
  cout << " ypTar: dypup   = " << dypup << endl;
  cout << endl << endl;

//============================================
//   Input Kinematics Settings for MC
//============================================
  float beamEnergy, hsec, thetac;
  fstream kinTblFile(kinFile);
  //Kin0 has a 2% offset in hsec.
  GetTblLine(kinTblFile, kin+1);  //Changes pointer.
  kinTblFile >> beamEnergy >> hsec >> thetac;

  cout << "=============================================\n";
  cout << "           Beam and Kineamtic Info           \n";
  cout << "=============================================\n";
  cout << " Beam Energy:      " << beamEnergy << endl;
  cout << " Central Momentum: " << hsec << endl;
  cout << " Central Angle:    " << thetac << endl;
  cout << endl << endl;

//Convert the central angel to radians
  double thetacrad = thetac * TMath::DegToRad();

  double dep          = (delup - deldown) / 100. * hsec;
  double dxp          = (dxpup - dxpdown);
  double dyp          = (dypup - dypdown);
  phase_space = dxp*dyp*dep / 1000.0;
  cout << "=============================================\n";
  cout << "              Phase Space                    \n";
  cout << "=============================================\n";
  cout << "dE': " << dep << endl;
  cout << "dx': " << dxp << endl;
  cout << "dy': " << dyp << endl;
  cout << "Phase Space: " <<  phase_space << endl;
  cout << endl << endl;

//============================================
//       Read in Target date from file
//============================================
  fstream targInfo(targetInfoFile);
  string targLine, targ2;
  int   tar_id;
  float density, tar_length, at_mass, at_no, aerial_density;

  while (getline(targInfo, targLine)) {
    if (targLine.find(targ.c_str()) != std::string::npos) {
      istringstream ss(targLine);
      ss >> targ2 >> tar_id >> density >> tar_length >>
	at_mass >> at_no >> aerial_density;
    }
  }

  double cCorrFact;
//Insert lumdata
//double lumdata = (aerial_density * 6.022137E-10 / at_mass) * (bcmavecharge / 1.602177E-13);
//  double lumdata = (density * tar_length / at_mass * 3758.72141467515);
  double lumdata = (aerial_density / at_mass * 3758.72141467515);
  //NEED ONLY APPLY TO CRYO TARGETS SILLY!!!!
  if(std::strcmp(targ.c_str(),"H1")==0) {cCorrFact = 0.996;}
  else if(std::strcmp(targ.c_str(),"D2")==0) {cCorrFact = 0.996;}
  else {cout << "No Cryo Correction solid Target: " << targ.c_str() << endl;
   cCorrFact = 1.0;
  }
  lumdata = lumdata * cCorrFact;

  cout << "=============================================\n";
  cout << "              Target Information             \n";
  cout << "=============================================\n";
  cout << " Target:    " << targ2 << endl;
  cout << " Target ID: " << tar_id << endl;
  cout << " Target Density:  " << density << endl;
  cout << " Aerial Density:  " << aerial_density << endl;
  cout << " Target Length:   " << tar_length << endl;
  cout << " Atomic Mass:     " << at_mass << endl;
  cout << " Atomic Number:   " << at_no << endl;
  if(std::strcmp(targ.c_str(),"H1")==0) {
    cout << " Cryo Correction: " << cCorrFact << endl;}
  else if(std::strcmp(targ.c_str(),"D2")==0) {
    cout << " Cryo Correction: " << cCorrFact << endl;}
  cout << " Data Luminosity: " << lumdata << endl;
  cout << endl << endl;


//==============================================
//  Interpolation from a 2D graph using root.
//     and put them in a TGraph2D().
//==============================================
 TGraph2D* xt = new TGraph2D();
 TGraph2D* rt = new TGraph2D();
 TGraph2D* cst = new TGraph2D();

//variables to hold elements in each line
 double ep, theta, x_bj, Q2, amu, radCorrFactor, C_cor;
 double Sig_Born, Sig_Born_In, Sig_Born_QE,
   Sig_Rad, Sig_Rad_EL, Sig_Rad_QE, Sig_Rad_DIS;
 double x_sec, wsqr, radCorrfactor;
 int i = 0;
 string line;
 ifstream crossFile(csTableFile);
//Open file

  int nSkip = 2000;
  cout << "=============================================\n";
  cout << "            Cross Section Table              \n";
  cout << "=============================================\n";
  cout << "Every entry that is a multiple of " << nSkip << ":\n";
  cout << " E'\tTheta\tX\tSigRad\n";
 
  while(getline(crossFile, line)) {
    if(i == 0) {
      i+=1; continue;
    }
    istringstream ss(line);
    i+=1;
    ss >> beamEnergy >> ep >> theta  >> x_bj >> Q2 >>
      Sig_Born >> Sig_Born_In >> Sig_Born_QE >> Sig_Rad >>
      Sig_Rad_EL >> Sig_Rad_QE >> Sig_Rad_DIS >> C_cor;
    if(i%nSkip == 0) {
      cout << ep << "\t" << theta << "\t" << x_bj << "\t" << Sig_Rad  << endl;
    }
    if(x_bj < 1.85) {
      xt->SetPoint(i, ep, theta, x_bj);
      cst->SetPoint(i, ep, theta, Sig_Rad);
    }
  }

  //Open the root tree output file
  compFile  = new TFile(Form(mcComparisonFile,kin, model, targ.c_str()), "UPDATE");     // May need to fix open type and location soon.

  auto rad2mrad   = [](float pFocalData) {return pFocalData*1000.0;};
  auto hseCalc = [=] (float deltai) {return hsec*(1.0 + deltai / 100.);};
  auto thetainiCalc = [=] (float yptari, float xptari) {return TMath::ACos(TMath::Cos(thetacrad+yptari)*TMath::Cos(xptari));};
  auto hsthetaCalc = [=] (float yptar, float xptar) {return  TMath::ACos(TMath::Cos(thetacrad+yptar)*TMath::Cos(xptar));};
  auto sin2Calc = [] (double thetaini) {return TMath::Power(TMath::Sin(thetaini/2.), 2);};
  auto nuCalc = [=] (double hse) {return beamEnergy - hse;};
  auto q2Calc = [=] (double hse, double sin2) {return 4.0*hse*beamEnergy*sin2;};
  auto w2Calc = [=] (double nu, double q2) {return mp2 + 2.*mp*nu - q2;};
  auto xbjCalc = [=] (double q2, double nu) {return q2 / (2*mp*nu);};
  auto phaseSpCorCalc = [] (float psxptar, float psyptar) {return 1./TMath::Power((1+TMath::Power(psxptar,2)+TMath::Power(psyptar,2)),(3/2));};
  auto csbCorCalc = [=] (double thetaini, double nu) {
    if(doCSBCorr) {
      Double_t csb_p0=-2.09 * thetaini*180./TMath::Pi() + 12.47;
      Double_t csb_p1=0.2 * thetaini*180./TMath::Pi() - 0.6338;
      return TMath::Exp(csb_p0)*(TMath::Exp(csb_p1*(nu))-1);
    }
    else return 0.0;
  };

  //Open the unweighted root tree files. chain with 1411
  ROOT::EnableImplicitMT();
  //cout << ROOT::GetImplicitMTPoolSize() << endl;
  ROOT::RDataFrame d("h1411", mcRawFiles); // Interface to TTree and TChain with ALL monte-carlo files '*'
  //Perform cuts on the monte-carlo (fail_id==0, delta cut? AND THATS IT!)
  auto d2 = d.Filter("psdelta < 22 && psdelta > -10")
             .Filter("stop_id == 0");

  //Add releveant branches
  d2 = d2.Define("hse", hseCalc, {"psdelta"})
         .Define("xptari_mrad", rad2mrad, {"psxptar"})
         .Define("yptari_mrad", rad2mrad, {"psyptar"})
         .Define("xptar_mrad", rad2mrad, {"psxptar"})
         .Define("yptar_mrad", rad2mrad, {"psyptar"})
         .Define("xpfp_mrad", rad2mrad, {"psxpfp"})
         .Define("ypfp_mrad", rad2mrad, {"psypfp"})
         .Define("thetaini", thetainiCalc, {"psyptar","psxptar"})
         .Define("hstheta", hsthetaCalc, {"psyptar","psxptar"})
         .Define("sin2", sin2Calc , {"thetaini"})
         .Define("nu", nuCalc , {"hse"})
         .Define("q2", q2Calc, {"hse","sin2"})
         .Define("w2", w2Calc, {"nu","q2"})
         .Define("xbj", xbjCalc, {"q2","nu"})
         .Define("csbCorr",csbCorCalc,{"thetaini","nu"})
         .Define("phasespcor", phaseSpCorCalc, {"psxptar","psyptar"});

  auto radCorCalc = [&] (double w2, double thetaini) {
    double radCorFac = rt->Interpolate(w2, TMath::RadToDeg()*thetaini);
    //cout << "radCorFac failed\n";
    return radCorFac !=0.0 ? radCorFac : radCorFac;
  };
  auto bornEMCCalc = [&] (double w2, double thetaini) {
    double cs = (cst->Interpolate(w2, TMath::RadToDeg()*thetaini));
    //cout << "cs failed\n";
    return cs != 0 ? cs :  cs;
  };

  auto bornCalc = [&] (double hse, double thetaini) {return cst->Interpolate(hse, TMath::RadToDeg()*thetaini);};

  auto d3 = d2;
  d3 = d2.Define("born", bornCalc, {"hse","thetaini"});

  int failCtr = 0;
  auto ngenCutBorn = [] (double ep, double th, double born) {
    //if(born == 0 ) {cout << ep << " " << th*TMath::RadToDeg() << endl;}
    return born > 0.0;
  };
  auto ngenCut = [=] (float yptari, float xptari, float deltai, float ytar) {
    return (xptari < dxpup/1000. && xptari > dxpdown/1000. && yptari < dypup/1000. && 
	    yptari > dypdown / 1000. && deltai > deldown && deltai < delup && abs(ytar) < 10.0);
  };

  cout << "=============================================\n";
  cout << "            Monte Carlo Events               \n";
  cout << "=============================================\n";
  auto ngenPassed = d3.Filter(ngenCut,{"psyptar","psxptar","psdelta","psytar"}).Count();
  cout << *ngenPassed << " events made it to the detector\n  within the phase space.\n";
  auto ngenPassedBorn = d3.Filter(ngenCutBorn,{"hse","thetaini","born"}).Count();
  cout << *ngenPassedBorn << " events many events passed the\n  born cut and are in phase space.\n";
  auto ngen = d.Count();
  cout << *ngen << " events were generated in the\n  gen limits of  the monte-carlo.\n";
  cout << endl << endl;

  fract = lumdata * phase_space / (*ngen + (*ngenPassed - *ngenPassedBorn)) / 1000.00;
  cout << "=============================================\n";
  cout << "           Monte Carlo Weighting             \n";
  cout << "=============================================\n";
  cout << " Fract: " << fract << endl;
  cout << endl << endl;

  auto mcWgtCalc = [=] (double born, double phasespcor, double csb) {return (born + csb)*phasespcor*fract;};

  cout << "=============================================\n";
  cout << "                   Output                    \n";
  cout << "=============================================\n";
  cout << " Writing root file to: \n";
  cout << Form(mcComparisonFile,kin, model, targ.c_str()) << endl;
  cout << endl << endl;

  d3 = d3.Define("mcWgt", mcWgtCalc, {"born","phasespcor","csbCorr"});

  if(snapshot) {
    d3.Snapshot("MCWgtTree",Form(mcComparisonFile,kin, model, targ.c_str()),
		{"mcWgt, born","w2","hse","phasespcor","thetaini","deltai","psyptar","psxptar"}
		,snapshotOptions);
  }
  compFile  = new TFile(Form(mcComparisonFile,kin, model, targ.c_str()), "UPDATE");
  //Define all the histograms
  auto h_xFocalMCWgt  = d3.Histo1D(m_xFocalMCWgt,"psxfp","mcWgt");
  auto h_xpFocalMCWgt = d3.Histo1D(m_xpFocalMCWgt,"xpfp_mrad","mcWgt");
  auto h_yFocalMCWgt  = d3.Histo1D(m_yFocalMCWgt,"psyfp","mcWgt");
  auto h_ypFocalMCWgt = d3.Histo1D(m_ypFocalMCWgt,"ypfp_mrad","mcWgt");
  auto h_yTarMCWgt    = d3.Histo1D(m_yTarMCWgt,"psytar","mcWgt");
  auto h_xpTarMCWgt   = d3.Histo1D(m_xpTarMCWgt,"xptar_mrad","mcWgt");
  auto h_ypTarMCWgt   = d3.Histo1D(m_ypTarMCWgt,"yptar_mrad","mcWgt");
  auto h_deltaMCWgt   = d3.Histo1D(m_deltaMCWgt,"psdelta","mcWgt");
  auto h_thetaMCWgt   = d3.Histo1D(m_thetaMCWgt,"hstheta","mcWgt");
  auto h_q2MCWgt      = d3.Histo1D(m_q2MCWgt,"q2","mcWgt");
  auto h_w2MCWgt      = d3.Histo1D(m_w2MCWgt,"w2","mcWgt");
  auto h_xbjMCWgt     = d3.Histo1D(m_xbjMCWgt,"xbj","mcWgt");
  // Weighted 2D MC histos
  auto h2_xVxpFocalMCWgt = d3.Histo2D(m_xVxpFocalMCWgt,"xpfp_mrad","psxfp","mcWgt");
  auto h2_xVyFocalMCWgt  = d3.Histo2D(m_xVyFocalMCWgt,"psyfp","psxfp","mcWgt");
  auto h2_xVypFocalMCWgt = d3.Histo2D(m_xVypFocalMCWgt,"ypfp_mrad","psxfp","mcWgt");
  auto h2_xpVyFocalMCWgt = d3.Histo2D(m_xpVyFocalMCWgt,"psyfp","xpfp_mrad","mcWgt");
  auto h2_xpVypFocalMCWgt= d3.Histo2D(m_xpVypFocalMCWgt,"ypfp_mrad","xpfp_mrad","mcWgt");
  auto h2_yVypFocalMCWgt = d3.Histo2D(m_yVypFocalMCWgt,"ypfp_mrad","psyfp","mcWgt");
  auto h2_yVxpTarMCWgt   = d3.Histo2D(m_yVxpTarMCWgt,"xptar_mrad","psytar","mcWgt");
  auto h2_yVypTarMCWgt   = d3.Histo2D(m_yVypTarMCWgt,"yptar_mrad","psytar","mcWgt");
  auto h2_xpVypTarMCWgt  = d3.Histo2D(m_xpVypTarMCWgt,"yptar_mrad","xptar_mrad","mcWgt");

  mcWgtDir = dynamic_cast <TDirectory*> (compFile->Get("mcWgtDir"));
  if(!mcWgtDir) {mcWgtDir = compFile->mkdir("mcWgtDir"); mcWgtDir->cd();}
  else if(mcWgtDir) {compFile->rmdir("mcWgtDir"); mcWgtDir = compFile->mkdir("mcWgtDir"); mcWgtDir->cd();}

  //Write all the histograms to the output file.
  //Define all the histograms
  h_xFocalMCWgt->Write();
  h_xpFocalMCWgt->Write();
  h_yFocalMCWgt->Write();
  h_ypFocalMCWgt->Write();
  h_yTarMCWgt->Write();
  h_xpTarMCWgt->Write();
  h_ypTarMCWgt->Write();
  h_deltaMCWgt->Write();
  h_thetaMCWgt->Write();
  h_q2MCWgt->Write();
  h_w2MCWgt->Write();
  h_xbjMCWgt->Write();
  // Weighted 2D MC histos
  h2_xVxpFocalMCWgt->Write();
  h2_xVyFocalMCWgt->Write();
  h2_xVypFocalMCWgt->Write();
  h2_xpVyFocalMCWgt->Write();
  h2_xpVypFocalMCWgt->Write();
  h2_yVypFocalMCWgt->Write();
  h2_yVxpTarMCWgt->Write();
  h2_yVypTarMCWgt->Write();
  h2_xpVypTarMCWgt->Write();

  return 0;
}


double GetInpVar(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    double var;
    file >> var;
    return var;
}

std::fstream& GetTblLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}
