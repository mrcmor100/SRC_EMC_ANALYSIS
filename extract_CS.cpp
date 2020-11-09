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

double cal_ep_delta (double mom_val,double  delta) {
    return (mom_val * (1. + delta/100.));
}

double cal_ep (double mom_val,double  xbj, double eb, double thetai) {
  return Mp*xbj*eb / (Mp*xbj + 2*eb*pow(sin(0.5 * TMath::DegToRad() * thetai),2));
}

double cal_w2 (double eb, double ep, double Mp, double thetai) {
  return pow(Mp,2) + 2.* Mp*(eb - ep) - 4.*eb*ep*pow(sin(0.5*TMath::DegToRad() * thetai),2);
}

double cal_xbj(double ep, double thetai, double eb, double Mp) {
  return 4*eb*ep*pow(sin(TMath::DegToRad()*thetai/2.),2) / (2 * Mp * (eb - ep));
}

double cal_err(double n1, double e1, double n2, double e2) {
  return( (n1/n2) * sqrt( pow((e1/n1),2) + pow((e2/n2),2)));
} 

TVectorD bornCS(1000);
TVectorD bodek(1000);
TVectorD bodek_xval(1000);
TVectorD bodek_xerr(1000);
TVectorD bodek_CS_ratio(1000);
TVectorD bornErr(1000);
TVectorD xErr(1000);
TVectorD x_val(1000);


int extract_CS(int kin, string targ, int model) {

  vector<int> allKins = {0};

  TGraph2D* xt = new TGraph2D();
  TGraph2D* rt = new TGraph2D();
  TGraph2D* cst = new TGraph2D();
  
  string csTableFile     = Form(csTableDir,kin, targ.c_str(), model);
  //variables to hold elements in each line
  double ep, theta, x_bj, Q2, amu, radCorrFactor, C_cor;
  double Sig_Born, Sig_Born_In, Sig_Born_QE,
   Sig_Rad, Sig_Rad_EL, Sig_Rad_QE, Sig_Rad_DIS;
  double x_sec, wsqr, radCorrfactor;
  double beamEnergy;
  int i = 0;
  string line;
  ifstream crossFile(csTableFile);
 //Open file
 while(getline(crossFile, line)) {
   if(i == 0) {
     i+=1; continue;
   }
   istringstream ss(line);
   i+=1;
   ss >> beamEnergy >> ep >> theta  >> x_bj >> Q2 >>
     Sig_Born >> Sig_Born_In >> Sig_Born_QE >> Sig_Rad >>
     Sig_Rad_EL >> Sig_Rad_QE >> Sig_Rad_DIS >> C_cor;
   if(i%10 == 0) {
     cout << ep << " " << theta << " " << x_bj << " " << Sig_Born  << endl;
   }
   if(x_bj < 1.85) {
     xt->SetPoint(i, x_bj, theta, x_bj);
     cst->SetPoint(i, x_bj, theta, Sig_Born);
   }
 }

//============================================
//   Input Kinematics Settings for MC
//============================================
  double hsec, thetac;
  string kinFile = kinFileDir;
  fstream kinTblFile(kinFile);
  //Kin0 has a 2% offset in hsec.
 
  //Determine momentum based on kinematic.. Switch or if statement for now.
  for(int currentKin : allKins) {

    cout << "Doing kinematic: " << currentKin << endl;
    GetTblLine(kinTblFile, currentKin+1);  //Changes pointer.
    kinTblFile >> beamEnergy >> hsec >> thetac;
    cout << "=============================================\n";
    cout << "           Beam and Kineamtic Info           \n";
    cout << "=============================================\n";
    cout << " Beam Energy:      " << beamEnergy << endl;
    cout << " Central Momentum: " << hsec << endl;
    cout << " Central Angle:    " << thetac << endl;
    cout << endl << endl;

    comparisonFile = new TFile(Form(mcComparisonFile,currentKin,model,targ.c_str()));
    if(std::strcmp(targ.c_str(),"H1")==0 || std::strcmp(targ.c_str(),"D2")==0) {
      dataDir = dynamic_cast <TDirectory*> (comparisonFile->Get("subDataDir"));
      mcWgtDir = dynamic_cast <TDirectory*> (comparisonFile->Get("mcWgtDir"));
      h_xbjData = dynamic_cast <TH1D*> (dataDir->Get("h_xbjSubData"));
      h_xbjMCWgt = dynamic_cast <TH1D*> (mcWgtDir->Get("h_xbjMCWgt"));
    } else {
      dataDir = dynamic_cast <TDirectory*> (comparisonFile->Get("dataDir"));
      mcWgtDir = dynamic_cast <TDirectory*> (comparisonFile->Get("mcWgtDir"));
      h_xbjData = dynamic_cast <TH1D*> (dataDir->Get("h_xbjData"));
      h_xbjMCWgt = dynamic_cast <TH1D*> (mcWgtDir->Get("h_xbjMCWgt"));
    }

    //int totalBins = 0;
    nBins = h_xbjData->GetNbinsX(); 
    for(int bin = 0; bin <= nBins; bin++) {

      xbjBinCenter=h_xbjData->GetBinCenter(bin);

      xbjData= (double) h_xbjData->GetBinContent(bin);
      xbjDataError= (double) h_xbjData->GetBinError(bin);

      xbjMCWgt= (double) h_xbjMCWgt->GetBinContent(bin);
      xbjMCWgtError= (double) h_xbjMCWgt->GetBinError(bin);
      if(xbjMCWgt==0) {xbjMCWgt=1;}
      xbjRatioError = cal_err(xbjData, xbjDataError, xbjMCWgt, xbjMCWgtError);

      ePrime=cal_ep(hsec,xbjBinCenter,beamEnergy, thetac);
      W2 = cal_w2(beamEnergy, ePrime, Mp, thetac);
      modelCS=(double) cst->Interpolate(xbjBinCenter,thetac);

      if(xbjBinCenter > 0.3 && xbjBinCenter < 1.85) {
	bornCS[totalBins+bin] = (modelCS*xbjData/xbjMCWgt);
	bodek_CS_ratio[totalBins+bin] = (modelCS*xbjData/xbjMCWgt) / modelCS;
	bornErr[totalBins+bin] = (modelCS*xbjRatioError);
	if(bornCS[totalBins+bin] == 0) {
	  cout << "XbjBin: " << xbjBinCenter << " Model: " << model << 
	    " xbjData: " << xbjData << " xbjMCWgt: " << xbjMCWgt <<
	    " bornCS: " << (modelCS*xbjData/xbjMCWgt) << endl;
	}
	x_val[totalBins+bin] = xbjBinCenter;
	xErr[totalBins+bin] = (0);
      }
    }
    totalBins = totalBins + nBins;
  }
  double xNow; 
  for (int i = 1; i < 1000; i++) {
    xNow = (double) (i * ((1.80 - 0.2) / 1000.) + 0.2);
    bodek_xval[i] = xNow;
    bodek_xerr[i] = 0.;
    bodek[i] =(double) cst->Interpolate(xNow,thetac);
    cout << xNow << " " << bodek[i] << endl;
  }

  
  ofstream outCS;
  outCS.open(Form(formatOutFilePath, targ.c_str()));
  outCS << "xbj\tCS\tCS_stat\n";
  for(int i = 0; i < totalBins; i++) {
    if(bornCS[i] != 0) {
      outCS << x_val[i] << "\t" << bornCS[i] << "\t" << bornErr[i] << "\n";
    }
  }
  outCS.close();

  ofstream outCSmodel;
  outCSmodel.open(Form(formatOutFilePath2, targ.c_str()));
  outCSmodel << "xbj\tmodelCS\n";
  for(int i = 0; i < 1000; i++) {
    if(bodek[i] != 0) {
      outCSmodel << bodek_xval[i] << "\t" << bodek[i] << "\n";
    }
  }
  outCSmodel.close();

  TGraphErrors *gcx= new TGraphErrors(x_val,bornCS,xErr,bornErr);
  TGraphErrors *gbodec= new TGraphErrors(bodek_xval,bodek,bodek_xerr, bodek_xerr);
  TGraphErrors *gbodekRatio= new TGraphErrors(x_val,bodek_CS_ratio,xErr,xErr);
  
  gcx->Draw("ap");
  gbodec->SetMarkerStyle(6);
  gbodec->Draw("L SAME");
  //TCanvas *c2 = new TCanvas();
  //c2->cd();
  //gbodekRatio->SetMarkerStyle(22);
  //gbodekRatio->Draw("ap");

  return 0;
}

std::fstream& GetTblLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}
