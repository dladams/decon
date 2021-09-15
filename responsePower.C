//************************************************************************
// Helpers.
//************************************************************************

using Index = unsigned int;

int fetchTool(string tnam, AdcChannelTool*& ptoo, string prefix) {
  ptoo = ptm->getShared<AdcChannelTool>(tnam);
  if ( ptoo == nullptr ) {
    cout << prefix << "ERROR: Unable to find tool " << tnam << endl;
    return 1;
  }
  cout << prefix << "Found tool " << tnam << endl;
  return 0;
}

void printSamples(const AdcChannelDataMap& acm) {
  Index nsam = acm.begin()->second.samples.size();
  cout << "Printing " << nsam << " samples for " << acm.size() << " channels." << endl;
  Index wsam = 10;
  for ( Index isam=0; isam<nsam; ++isam ) {
    cout << setw(4) << isam << ":";
    for ( auto& kcha : acm ) {
      const AdcChannelData& acd = kcha.second;
      //cout << fixed << setprecision(5) << setw(10) << acd.samples.size();
      if ( isam < acd.samples.size() ) {
        cout << fixed << setprecision(5) << setw(10) << acd.samples[isam];
      } else {
        cout << setw(wsam) << "-";
      }
    }
    cout << endl;
  }
}

template<typename C>
TH2* printPowerByColumn(const C& dat, string sview, string hnam ="", float xmax =1.0, float ymax =1.0) {
  Index n0 = dat.size(0);
  Index n1 = dat.size(1);
  TH2* ph = nullptr;
  string httl = "DFT power " + sview + "; f_{t} [MHz]; f_{wire} [(wire)^{-1}]";
  if ( hnam.size() ) {
    ph = new TH2F(hnam.c_str(), httl.c_str(), n1, 0, xmax, n0, 0, ymax);
    ph->SetDirectory(0);
    ph->SetStats(0);
    ph->SetContour(200);
  }
  cout << "Printing power for " << dat.size(0) << " x " << dat.size(1) << " by column." << endl;
  float pwrsum = 0.0;
  typename C::IndexArray isams;
  Index& i0 = isams[0];
  Index& i1 = isams[1];
  for ( i1=0; i1<n1; ++i1 ) {
    cout << setw(4) << i1 << ":";
    for ( i0=0; i0<n0; ++i0 ) {
      float pwr = std::norm(dat.value(isams));
      pwrsum += pwr;
      if ( ph != nullptr ) ph->SetBinContent(i1+1, i0+1, pwr);
      cout << fixed << setprecision(5) << setw(10) << pwr;
    }
    cout << endl;
  }
  cout << "Sum power: " << fixed << setprecision(5) << pwrsum << endl;
  return ph;
}

//************************************************************************

int responsePower(string sview, Index nsam, Index nnbr) {
  string myname = "responsePower: ";
  // Fetch the tools.
  // Tools that use Root objects in their dtors are created as private and deleted
  // here to ensure those objects are still present.
  string line = "************************************************************";
  cout << myname << line << endl;
  cout << myname << "Fetching tools." << endl;
  cout << myname << "Fetching 2D response convolution tool." << endl;
  AdcChannelTool* pcon = nullptr;
  int nbin = 1;
  string tnam = "conv2dNbin" + to_string(nbin) + sview;
  if ( fetchTool(tnam, pcon, myname) ) return 2;
  AdcChannelTool* ppltSam = nullptr;
  if ( fetchTool("plotSamples", ppltSam, myname) ) return 2;
  cout << myname << "Create channel data." << endl;
  AdcChannelDataMap acm;
  Index irun = 123;
  Index ievt = 1;
  Index icha0 = sview == "x" ? 1800 : sview == "u" ? 200 : sview == "v" ? 1000 : 0;
  if ( icha0 == 0 ) {
    cout << "Invalid view: " << sview << endl;
    return 1;
  }
  Index icha1 = icha0 -nnbr;
  Index icha2 = icha0 + nnbr + 1;
  for ( Index icha=icha1; icha<icha2; ++icha ) {
    AdcChannelData& jacd = acm[icha];
    jacd.setEventInfo(irun, ievt);
    jacd.sampleUnit = "ke/tick";
    jacd.setChannelInfo(icha, 0, 0, 0);;
    jacd.samples.resize(nsam);
    jacd.binSamples.resize(1);
    jacd.binSamples[0].resize(nsam);
  }
  Index isam0 = 30;
  acm[icha0].samples[isam0] = 1.0;
  acm[icha0].binSamples[0][isam0] = 1.0;
  printSamples(acm);
  cout << myname << "Convolute." << endl;
  if ( pcon->updateMap(acm) ) {
    cout << myname << "Convolution failed." << endl;
    return 3;
  }
  printSamples(acm);
  cout << myname << "Plot samples." << endl;
  if ( ppltSam->viewMap(acm) ) {
    cout << myname << "Convolution failed." << endl;
    return 3;
  }

  cout << myname << "Extracting data." << endl;
  Index ncha = icha2 - icha1;
  Fw2dFFT::IndexArray ndxs = {ncha, nsam};
  Fw2dFFT::IndexArray idxs;
  Fw2dFFT::DFT::Norm norm(2,1);
  Fw2dFFT::Data dat(ndxs);
  Fw2dFFT::DFT dft(norm, ndxs);
  Index& dcha = idxs[0];
  Index& isam = idxs[1];
  Index ichk = 0;
  Index* pchk = &ichk;
  set<Index> idxchk;
  for ( dcha=0; dcha<ncha; ++dcha ) {
    Index icha = icha1 + dcha;
    cout << myname << "Setting data for channel " << icha << endl;
    const AdcChannelData& acd = acm[icha];
    if ( acd.samples.size() != nsam ) {
      cout << myname << "ERROR: Channel " << icha << " has a bad sample count: "
           << acd.samples.size() << " != " << nsam << endl;
      return 4;
    }
    Index ichk = 0;
    for ( isam=0; isam<nsam; ++isam ) {
      float val = acd.samples[isam];
      Index idat = dat.setValue(idxs, val, pchk);
      if ( pchk != nullptr && *pchk ) {
        cout << myname << "ERROR: Unable to set data value for channel " << icha
             << " sample " << isam << endl;
        return 5;
      }
      if ( idxchk.count(idat) ) {
        cout << myname << "ERROR: Duplicate set index with channel " << icha
             << " sample " << isam << endl;
        return 6;
      }
      idxchk.insert(idat);
    }
  }
  cout << myname << "Value set count: " << idxchk.size() << endl;
  printPowerByColumn(dat, sview);
  
  cout << myname << "Doing FFT." << endl;
  Index ndft = 2*ncha*(nsam/2+1);
  Index dftopt = 1;
  Fw2dFFT xf(ndft, dftopt);
  xf.fftForward(dat, dft);
  cout << myname << " Data size is " << dat.size(0) << " x " << dat.size(1) << endl;
  cout << myname << "  DFT size is " << dft.size(0) << " x " << dft.size(1) << endl;
  string dftname = "dft-" + sview;
  float xmax = 2.0;  // MHz
  float ymax = 1.0;
  int nx = nsam;
  int ny = ncha;
  float xbin = xmax/float(nx);
  float ybin = ymax/float(ny);
  TH2* phdft = printPowerByColumn(dft, sview, dftname, xmax, ymax);
  TPadManipulator man(1400, 1000);
  man.add(phdft, "colz");
  man.setPalette(2010);
  man.setLogZ();
  man.setLogRangeZ(1.e-10, 1.e-2);
  bool zoom =  true;
  bool symlines = true;
  if ( symlines ) {
    vector<float> xsyms = { 1.0, float((nx+1)/2), float((nx+2)/2) };
    vector<float> ysyms = { 1.0, float((ny+1)/2), float((ny+2)/2) };
    for ( float xsym : xsyms ) {
      float xval = xbin*xsym;
      man.addVerticalLine(xval);
    }
    for ( float ysym : ysyms ) {
      float yval = ybin*ysym;
      man.addHorizontalLine(yval);
    }
  }
  if ( zoom ) {
    int ix = nx/2 + 1;
    float xmax = xbin*ix;
    int iy = ny/2 + 1;
    float ymax = ybin*iy;
    man.setRangeX(0, xmax);
    man.setRangeY(0, ymax);
  }
  man.addAxis();
  man.print(dftname + ".png");
  
  phdft->Print();

  cout << myname << "Done." << endl;
  return 0;
}
