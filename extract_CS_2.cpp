
const   double Mp = 0.9382723;

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
  return( (n1/n2) * sqrt( pow(((e1)/n1),2) + pow((e2/n2),2) )  );
} 

void extract_a_ratio(TGraphErrors *num, TGraphErrors *den, string targ) {

  TGraphErrors *ratio = new TGraphErrors();
  int numNum = num->GetN();
  int denNum = den->GetN();

  double numX, numY, denX, denY;
  double ratioX, ratioY, ratioErrX, ratioErrY;
  int largest;
  if(numNum > denNum) {
    largest = numNum;
    cout << largest << endl;
  } else if(numNum < denNum) {
    largest = denNum;
    cout << largest << endl;
  } else {
    largest = numNum;
    cout << largest << endl;
  }
  for(int i = 0; i < largest; i++) {
    num->GetPoint(i,numX, numY);
    den->GetPoint(i,denX, denY);

    if(numX == denX) {
      ratioX = numX;
      if(denY !=0) {
	ratioY = numY / denY;
	cout << ratioX << " " << ratioY << endl;
	ratio->SetPoint(i, 1., 1.*i);
      }
      else {continue;}
      ratioErrX = 0;
      ratioErrY = cal_err(numY,num->GetErrorY(i),denY, den->GetErrorY(i));
      ratio->SetPointError(i,0.1,0.1);
      if(i!=0){
	ratio->GetPoint(i-1,numX,numY);
	cout <<numX << " " << numY << endl;
      }
    }
    else {
      cout << "SOMETHING BAD HAPPENED!\n";
      continue;
    }
  }
  cout << ratio->GetN() << endl;
  ratio->Write(Form("%s_ratio",targ.c_str()));
  return;
}

TGraphErrors *extract_a_CS(int kin, string targ, int model) {
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
  double xbjBinCenter, xbjRatioError, ePrime, W2;
  int nBins;
  int totalBins = 0;
  double central_p;

  vector<double> bornCS = {};
  vector<double> bodek = {};
  vector<double> bodek_xval = {};
  vector<double> bodek_xerr = {};
  vector<double> bodek_CS_ratio = {};
  vector<double> bornErr = {};
  vector<double> xErr = {};
  vector<double> x_val = {};

  int palette[5] = {kMagenta+2,kRed-3,kBlue-3,kGreen+2,kCyan-6};
  gStyle->SetPalette(5, palette);

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
   if(i%1000 == 0) {
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
	bornCS.push_back(modelCS*xbjData/xbjMCWgt);
	bodek_CS_ratio.push_back((modelCS*xbjData/xbjMCWgt) / modelCS);
	bornErr.push_back(2.414*modelCS*(xbjData / xbjMCWgt) * xbjDataError / xbjData);
	//if(bornCS.back() == 0) {
	//cout << "XbjBin: " << xbjBinCenter << " Model: " << model << 
	//" xbjData: " << xbjData << " xbjMCWgt: " << xbjMCWgt <<
	//" bornCS: " << (modelCS*xbjData/xbjMCWgt) << endl;
	//}
	x_val.push_back(xbjBinCenter);
	xErr.push_back(1.E-6);
      }
    }
    totalBins = totalBins + nBins;
  }

  TGraphErrors *gcx= new TGraphErrors();

  for(int i = 0; i < bornCS.size(); i++) {
    gcx->SetPoint(i,x_val[i],bornCS[i]/1000.);
    gcx->SetPointError(i,xErr[i],bornErr[i]/1000.);
  }

  return gcx;
}



TMultiGraph* extract_CS_2(int kin, string targ, int model) {

  int palette[5] = {kMagenta+2,kRed-3,kBlue-3,kGreen+2,kCyan-6};
  gStyle->SetPalette(5, palette);

  vector<string> allTargs = {"C12","D2","B10","B11"};
 
  TGraphErrors *gcx = new TGraphErrors();
  TGraphErrors *deut = new TGraphErrors();
  TCanvas *csCanvas = new TCanvas("csCanvas","csCanvas",600,600);
  TLegend *legend = new TLegend(0.55, 0.70, .875, 0.875);

  TFile *f1 = new TFile("temp.root","RECREATE");

  for(string targ : allTargs) {
    
    gcx = extract_a_CS(kin, targ, model);
    f1->cd();
    gcx->SetTitle(Form("%s Spring 2018 Born Cross-Section",targ.c_str()));
    gcx->GetXaxis()->SetTitle("X");
    gcx->GetYaxis()->SetTitle("#sigma [nb/sr]");
    gcx->GetXaxis()->CenterTitle();
    gcx->GetYaxis()->CenterTitle();
    gcx->GetYaxis()->SetTitleOffset(1.3);
    gcx->SetLineColorAlpha(kBlack,0.0);
    //gcx->GetXaxis()->SetTitleOffset(1.1);
    gcx->GetXaxis()->SetRangeUser(0.85,1.8);
    gcx->GetYaxis()->SetRangeUser(1E-4,1E1);
    gcx->SetMarkerSize(0.65);
    gcx->SetMarkerStyle(kFullCircle);
    gcx->Write(Form("%s_cs",targ.c_str()));
  }

  TMultiGraph *mg = new TMultiGraph();
  TMultiGraph *mg_r = new TMultiGraph();
  TCanvas *c2 = new TCanvas();
  TGraphErrors *num, *den;
  TGraphErrors *ratio;
  int numNum, denNum;
  double numX, numY, denX, denY;
  double ratioX, ratioY, ratioErrX, ratioErrY;
  double rX, rY;
  int f=0, j, offset = 0;


  for(string targ : allTargs) {
    csCanvas->cd();
    gcx = dynamic_cast<TGraphErrors*> (f1->Get(Form("%s_cs",targ.c_str())));
    mg->Add(gcx);
    legend->AddEntry(gcx, Form("%s Cross-Section",targ.c_str()));
  }
  /*
  for(string targ : allTargs) {
    num = dynamic_cast<TGraphErrors*> (f1->Get(Form("%s_cs",targ.c_str())));
    den = dynamic_cast<TGraphErrors*> (f1->Get(Form("%s_cs","D2")));
    if(strcmp(targ.c_str(),"D2") == 0) {
      cout << targ << endl;
    } else {
      numNum = num->GetN();
      denNum = den->GetN();
      cout << "Total N - num: " << numNum << " den: " << denNum << endl;
      for(int i = 0; i < numNum; i++) {
	num->GetPoint(i,numX, numY);
	den->GetPoint(i,denX, denY);
	//cout << offset << endl;
	if(numX < 0.7 || numX > 1.8) {offset+=1; continue;}
	//j = i - offset;
	//cout << j << endl;
	//cout << numX << " " << denX << endl;
	//cout << numY << " " << denY << endl;
	ratioX = (double) numX;
	ratioY = (double) numY / denY + 1E-10;
	cout << f << " " << ratioX << " " << ratioY << endl;
	ratio->SetPoint(f, ratioX, ratioY);
	f+=1;
	//ratioErrX = 0.;
	//ratioErrY = cal_err(numY,num->GetErrorY(i),denY, den->GetErrorY(i));
	//ratio->SetPointError(j,ratioErrX,ratioErrY);
      }
    }
  }
  */

  csCanvas->SetLogy();  
  csCanvas->SetTicks(1,1);
  csCanvas->SetGrid(1,1);

  mg->Draw("APE PLC PMC");
  mg->SetTitle("Spring 2018 Born Cross-Section");
  mg->GetXaxis()->SetTitle("X");
  mg->GetYaxis()->SetTitle("#sigma [nb/MeV/sr]");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.3);
  //gcx->GetXaxis()->SetTitleOffset(1.1);
  mg->GetXaxis()->SetRangeUser(0.85,1.8);
  mg->GetYaxis()->SetRangeUser(1E-4,1E1);
  legend->Draw();  
  csCanvas->Update();

  c2->cd();
  mg_r->GetXaxis()->SetRangeUser(0.85,1.8);
  //mg_r->GetYaxis()->SetRangeUser(1E-8,0.5E0);
  mg_r->Draw("AP PLC PMC");

  //gcx->Draw("ap");
  //gbodec->SetMarkerStyle(6);
  //gbodec->Draw("L SAME");
  //TCanvas *c2 = new TCanvas();
  //c2->cd();
  //gbodekRatio->SetMarkerStyle(22);
  //gbodekRatio->Draw("ap");
  
  return mg_r;
}

std::fstream& GetTblLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}
